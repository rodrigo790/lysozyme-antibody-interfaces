"""Compare experimental antibody contacts with hen egg-white lysozyme."""
from itertools import combinations

from Bio.PDB import PDBParser
from Bio.SeqUtils import seq1
import numpy as np
import pandas as pd
import py3Dmol
from scipy.spatial import cKDTree

from common import RESULTS, download, plot_style, save_summary

COMPLEXES = {
    "1VFB": ("D1.3", "C", "AB"),
    "1YQV": ("HyHEL-5", "Y", "LH"),
    "1MLC": ("D44.1", "E", "AB"),
}
CUTOFFS = (3.5, 4.0, 4.5, 5.0)


def residue_id(residue):
    return f"{residue.id[1]}{residue.id[2].strip()}"


def heavy_atoms(chain):
    return [
        atom for residue in chain if residue.id[0] == " "
        for atom in residue if atom.element.upper() not in {"H", "D"}
    ]


def residue_contacts(antigen, antibody, radius):
    """The minimum heavy-atom distance for each unique residue pair."""
    neighbors = cKDTree([atom.coord for atom in antibody])
    matches = {}
    for atom, near in zip(
        antigen, neighbors.query_ball_point([atom.coord for atom in antigen], radius)
    ):
        for index in near:
            partner = antibody[index]
            a, b = atom.parent, partner.parent
            key = (residue_id(a), a.resname, b.parent.id, residue_id(b), b.resname)
            distance = float(atom - partner)
            matches[key] = min(distance, matches.get(key, float("inf")))
    return matches


def main():
    RESULTS.mkdir(exist_ok=True)
    parser = PDBParser(QUIET=True)
    columns = ["antigen_residue", "antigen_aa", "antibody_chain", "antibody_residue", "antibody_aa"]
    rows, counts, sets = [], [], {}
    viewer = py3Dmol.view(viewergrid=(1, 3), width=990, height=450, linked=False)
    reference = None
    pml = ["hide everything", "bg_color white"]
    for panel, (pdb, (name, antigen_id, partners)) in enumerate(COMPLEXES.items()):
        path = download(pdb.lower() + ".pdb", f"https://files.rcsb.org/download/{pdb}.pdb")
        model = parser.get_structure(pdb, path)[0]
        antigen = heavy_atoms(model[antigen_id])
        residues = [res for res in model[antigen_id] if res.id[0] == " " and "CA" in res]
        numbering = [(residue_id(res), seq1(res.resname)) for res in residues]
        if reference is None:
            reference = numbering
        if numbering != reference:
            raise ValueError(f"Lysozyme residue mapping differs in {pdb}")
        antibody = [atom for chain in partners for atom in heavy_atoms(model[chain])]
        maximum = residue_contacts(antigen, antibody, max(CUTOFFS))
        for cutoff in CUTOFFS:
            subset = {key for key, distance in maximum.items() if distance <= cutoff}
            sets[pdb, cutoff] = {key[0] for key in subset}
            counts.append({
                "pdb": pdb, "antibody": name, "cutoff_A": cutoff,
                "residue_pairs": len(subset), "antigen_residues": len(sets[pdb, cutoff]),
                "antibody_residues": len({key[2:4] for key in subset}),
            })
        selected = sorted(map(int, sets[pdb, 4.0]))
        rows.extend(
            {"pdb": pdb, "antibody": name, **dict(zip(columns, key)), "distance_A": distance}
            for key, distance in maximum.items()
        )
        keep = set(antigen_id + partners)
        coordinates = "\n".join(
            line for line in path.read_text().splitlines()
            if line.startswith("ATOM  ") and line[21] in keep
        ) + "\nEND\n"
        slot = (0, panel)
        viewer.addModel(coordinates, "pdb", viewer=slot)
        viewer.setStyle({}, {"cartoon": {"color": "#bdc9d4"}}, viewer=slot)
        viewer.setStyle({"chain": antigen_id}, {"cartoon": {"color": "#258579"}}, viewer=slot)
        viewer.addStyle(
            {"chain": antigen_id, "resi": selected}, {"stick": {"color": "#d87737"}}, viewer=slot
        )
        viewer.addLabel(f"{name} · {pdb}", {"position": {"x": 0, "y": 0, "z": 0},
                        "inFront": True, "backgroundColor": "#26384a"}, viewer=slot)
        viewer.zoomTo(viewer=slot)
        pml.extend([
            f"load raw/{pdb.lower()}.pdb, {pdb}",
            f"remove {pdb} and not chain {'+'.join(keep)}",
            f"show cartoon, {pdb}", f"color gray70, {pdb}",
            f"color teal, {pdb} and chain {antigen_id}",
            f"show sticks, {pdb} and chain {antigen_id} and resi {'+'.join(map(str, selected))}",
        ])

    pairs = []
    for cutoff in CUTOFFS:
        for a, b in combinations(COMPLEXES, 2):
            left, right = sets[a, cutoff], sets[b, cutoff]
            pairs.append({
                "a": a, "b": b, "cutoff_A": cutoff,
                "shared_residues": len(left & right),
                "jaccard": len(left & right) / len(left | right),
            })
    pd.DataFrame(counts).to_csv(RESULTS / "counts.csv", index=False)
    pd.DataFrame(rows).sort_values(["pdb", "distance_A"]).to_csv(
        RESULTS / "residue_pairs.csv", index=False, float_format="%.4f"
    )
    pd.DataFrame(pairs).to_csv(RESULTS / "overlap.csv", index=False, float_format="%.4f")
    viewer.write_html(str(RESULTS / "viewer.html"))
    (RESULTS / "view.pml").write_text("\n".join([*pml, "set grid_mode, 1", "zoom"]) + "\n")

    plt = plot_style()
    fig, ax = plt.subplots(figsize=(9, 3.3))
    for level, (pdb, (name, _, _)) in enumerate(COMPLEXES.items()):
        positions = sorted(map(int, sets[pdb, 4.0]))
        ax.scatter(positions, [level] * len(positions), marker="|", s=320, linewidth=2.1)
    ax.set(xlim=(0, 130), xlabel="Lysozyme residue number", yticks=range(len(COMPLEXES)),
           yticklabels=[f"{name} · {pdb}" for pdb, (name, _, _) in COMPLEXES.items()],
           title="Antibody footprints on hen egg-white lysozyme · 4 Å")
    ax.grid(axis="y", visible=False)
    fig.savefig(RESULTS / "footprints.png")
    plt.close(fig)
    save_summary({
        "complexes": len(COMPLEXES), "lysozyme_residues": len(reference),
        "residues_at_4A": {pdb: len(sets[pdb, 4.0]) for pdb in COMPLEXES},
        "jaccard_at_4A": {
            f"{r['a']}–{r['b']}": round(r['jaccard'], 3) for r in pairs if r["cutoff_A"] == 4.0
        }
    })


if __name__ == "__main__":
    main()
