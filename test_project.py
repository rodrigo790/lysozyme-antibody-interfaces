import unittest

import pandas as pd

from analysis import CUTOFFS, COMPLEXES
from common import RESULTS


class InterfaceTests(unittest.TestCase):
    def test_footprints_grow_with_distance_cutoff(self):
        counts = pd.read_csv(RESULTS / "counts.csv")
        self.assertEqual(set(counts.pdb), set(COMPLEXES))
        self.assertEqual(set(counts.cutoff_A), set(CUTOFFS))
        for _, group in counts.sort_values("cutoff_A").groupby("pdb"):
            self.assertTrue(group.antigen_residues.is_monotonic_increasing)
            self.assertTrue(group.residue_pairs.is_monotonic_increasing)
        overlaps = pd.read_csv(RESULTS / "overlap.csv")
        self.assertTrue(overlaps.jaccard.between(0, 1).all())
