"""
 Copyright (C) 2024 Jessie Fielding

 This file is part of simble.

 simble is free software: you can redistribute it and/or modify
 it under the terms of the GNU Affero General Public License as
 published by the Free Software Foundation, either version 3 of the
 License, or (at your option) any later version.

 simble is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU Affero General Public License for more details.

 You should have received a copy of the GNU Affero General Public License
 along with simble.  If not, see <https://www.gnu.org/licenses/>.
 """

import argparse
import json
import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from simble import helper
from simble.cell import Cell, SingleChainCell
from simble.chain import EmptyChain, HeavyChain, LightChain
from simble.location import as_enum
from simble.parsing import (validate_and_process_chains,
                            validate_and_process_target_input)
from simble.settings import Settings, s
from simble.simble import do_simulation
from tests.test_internal_nodes import fasta_records

ENTROPY = 54897022524486695084299880814690718190


def run_clone(chains="paired", uniform=False, neutral=False):
    """Runs one short clone with a fixed seed and returns its results."""
    settings = Settings()
    settings.LOCATIONS[0].sample_times = list(range(0, 11, 5))
    settings.LOCATIONS[1].sample_times = list(range(0, 11, 5))
    settings.LOCATIONS[0].migration_rate = 0.5
    settings.INTERNAL_NODES = True
    settings.CHAINS = chains
    if uniform:
        settings.UNIFORM = True
    if neutral:
        settings.SELECTION = False
    if chains == "heavy":
        settings.TARGET_MUTATIONS_LIGHT = 0
    seed = np.random.SeedSequence(ENTROPY)
    with tempfile.NamedTemporaryFile(mode="w") as tmpf:
        json.dump(settings, tmpf, default=lambda o: o.encode(), indent=4)
        tmpf.flush()
        return do_simulation(1, seed, None, tmpf.name)


def make_args(**kwargs):
    """Returns command line arguments with the defaults of the simble parser."""
    defaults = {
        "single_chain": False,
        "uniform": False,
        "light_shm_per_site": None,
        "target_mutations_heavy": None,
        "target_mutations_light": None,
        "target": None,
        "target_heavy": None,
        "target_light": None,
        "target_join": "cell_id",
        "sequence_length": None,
    }
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


class SettingsSaver(unittest.TestCase):
    """Saves and restores the global settings that these tests change."""
    FIELDS = ["CHAINS", "_UNIFORM", "_SELECTION", "TARGET", "TARGET_MUTATIONS_HEAVY",
              "TARGET_MUTATIONS_LIGHT", "LIGHT_SHM_PER_SITE", "NAIVE_FILE"]

    def setUp(self):
        self.saved = {x: getattr(s, x) for x in self.FIELDS}
        # other tests can leave the global settings changed, so start from the defaults
        defaults = Settings()
        for key in self.FIELDS:
            setattr(s, key, getattr(defaults, key))

    def tearDown(self):
        for key, value in self.saved.items():
            setattr(s, key, value)
        helper.update_helper_tables()


class TestSingleChainSimulation(unittest.TestCase):
    """A single chain simulation uses the heavy chain of the naive pair only."""

    @classmethod
    def setUpClass(cls):
        cls.paired = run_clone("paired")
        cls.single = run_clone("heavy")

    def test_only_heavy_rows(self):
        airr = self.single["airr"]
        self.assertGreater(len(airr), 0)
        self.assertTrue(all(airr["sequence_id"].str.endswith("_heavy")))
        self.assertEqual(set(airr["locus"]), {"IGH"})
        # one row per cell
        self.assertEqual(len(airr), airr["cell_id"].nunique())

    def test_same_founder_as_paired(self):
        """With the same seed, the founder is the same naive heavy chain."""
        paired_heavy = self.paired["airr"]
        paired_heavy = paired_heavy[paired_heavy["sequence_id"].str.endswith("_heavy")]
        self.assertEqual(
            set(self.single["airr"]["germline_alignment"]),
            set(paired_heavy["germline_alignment"])
            )

    def test_heavy_target_has_mutations_light_target_is_empty(self):
        """Selection is on, so the heavy chain has a target. There is no light target."""
        self.assertEqual(self.single["targets"]["light"], "")
        self.assertNotEqual(self.single["targets"]["heavy"], "")
        germline = self.single["airr"]["germline_alignment"].iloc[0]
        naive_aa = helper.translate_to_amino_acid(germline)
        target_aa = self.single["targets"]["heavy"]
        differences = sum(1 for a, b in zip(naive_aa, target_aa) if a != b)
        self.assertEqual(differences, Settings().TARGET_MUTATIONS_HEAVY)

    def test_internal_nodes_one_row_per_node(self):
        internal = self.single["internal_airr"]
        self.assertGreater(len(internal), 0)
        self.assertEqual(len(internal), internal["cell_id"].nunique())
        self.assertTrue(all(internal["sequence_id"].str.endswith("_heavy")))

    def test_fasta_only_heavy(self):
        records = fasta_records(self.single["fasta"])
        self.assertGreater(len(records), 0)
        self.assertTrue(all(x.endswith("_heavy") for x in records))
        records = fasta_records(self.single["internal_fasta"])
        self.assertTrue(all(x.endswith("_heavy") for x in records))


class TestSingleChainNeutral(unittest.TestCase):
    """--single-chain works with --neutral."""

    def tearDown(self):
        # validate_and_process_args changes many settings, so reset them all
        s.update_from_dict(json.loads(
            json.dumps(Settings(), default=lambda o: o.encode()), object_hook=as_enum))
        helper.update_helper_tables()

    def test_parsed_settings_reach_the_workers(self):
        """The command line gives no selection, heavy chain only, and no target
        mutations, and these settings go to the workers in the settings JSON."""
        from simble.parsing import get_parser, validate_and_process_args  # pylint: disable=import-outside-toplevel
        with tempfile.TemporaryDirectory() as folder:
            args = get_parser().parse_args(
                ["-o", os.path.join(folder, "out"), "--single-chain", "--neutral", "-q"])
            warnings = validate_and_process_args(args)
        self.assertEqual(warnings, [])
        encoded = json.loads(json.dumps(s, default=lambda o: o.encode()), object_hook=as_enum)
        worker = Settings()
        worker.update_from_dict(encoded)
        for settings in [s, worker]:
            self.assertFalse(settings.SELECTION)
            self.assertTrue(settings.SINGLE_CHAIN)
            self.assertEqual(settings.TARGET_MUTATIONS_HEAVY, 0)
            self.assertEqual(settings.TARGET_MUTATIONS_LIGHT, 0)

    def test_simulation(self):
        """A neutral single chain clone has only heavy rows, mutates, and its
        target is its own naive heavy chain (no target mutations)."""
        result = run_clone("heavy", neutral=True)
        airr = result["airr"]
        self.assertTrue(all(airr["sequence_id"].str.endswith("_heavy")))
        self.assertTrue((airr["sequence_alignment"] != airr["germline_alignment"]).any())
        germline = airr["germline_alignment"].iloc[0]
        self.assertEqual(result["targets"]["heavy"], helper.translate_to_amino_acid(germline))
        self.assertEqual(result["targets"]["light"], "")


class TestSingleChainCells(SettingsSaver):
    """Cells and chains in a single chain or uniform simulation."""

    def test_empty_sequence_is_an_error(self):
        """An empty chain is made with EmptyChain, not with an empty sequence."""
        with self.assertRaises(ValueError):
            LightChain("", gapped_seq="")
        with self.assertRaises(ValueError):
            HeavyChain("", gapped_seq="")

    def test_empty_chain(self):
        chain = EmptyChain()
        self.assertEqual(chain.nucleotide_seq, "")
        self.assertEqual(chain.similarity, 1)

    def _make_founder(self):
        s._x_RNG = np.random.default_rng(1)  # pylint: disable=protected-access
        try:
            return Cell(None, None, created_at=0)
        finally:
            s._x_RNG = None  # pylint: disable=protected-access

    def test_single_chain_founder(self):
        s.CHAINS = "heavy"
        cell = self._make_founder()
        self.assertIsInstance(cell, SingleChainCell)
        self.assertIsInstance(cell.light_chain, EmptyChain)
        self.assertGreater(len(cell.heavy_chain.nucleotide_seq), 0)

    def test_uniform_founder_light_chain_is_complete(self):
        """The founder's empty light chain is a full EmptyChain (before, it had
        no similarity values)."""
        s.UNIFORM = True
        cell = self._make_founder()
        self.assertIsInstance(cell, SingleChainCell)
        self.assertIsInstance(cell.light_chain, EmptyChain)
        self.assertEqual(cell.light_chain.similarity, 1)

    def test_paired_founder(self):
        cell = self._make_founder()
        self.assertNotIsInstance(cell, SingleChainCell)
        self.assertIsInstance(cell.light_chain, LightChain)

    def test_active_chains(self):
        self.assertEqual(s.ACTIVE_CHAINS, ["heavy", "light"])
        s.CHAINS = "heavy"
        self.assertEqual(s.ACTIVE_CHAINS, ["heavy"])
        # a uniform simulation keeps its own light chain handling
        s.CHAINS = "paired"
        s._UNIFORM = True  # pylint: disable=protected-access
        self.assertFalse(s.SINGLE_CHAIN)
        self.assertEqual(s.ACTIVE_CHAINS, ["heavy", "light"])


class TestValidateChains(SettingsSaver):
    """Tests for the validation of the CHAINS setting."""

    def test_default_is_paired(self):
        warnings = []
        validate_and_process_chains(make_args(), warnings)
        self.assertEqual(s.CHAINS, "paired")
        self.assertEqual(warnings, [])

    def test_single_chain_flag(self):
        s.TARGET_MUTATIONS_LIGHT = 5
        warnings = []
        validate_and_process_chains(make_args(single_chain=True), warnings)
        self.assertEqual(s.CHAINS, "heavy")
        self.assertEqual(s.TARGET_MUTATIONS_LIGHT, 0)
        self.assertEqual(warnings, [])

    def test_light_arguments_are_ignored(self):
        warnings = []
        validate_and_process_chains(
            make_args(single_chain=True, light_shm_per_site=0.1, target_mutations_light=3),
            warnings)
        self.assertEqual(len(warnings), 2)
        self.assertEqual(s.TARGET_MUTATIONS_LIGHT, 0)

    def test_invalid_value(self):
        s.CHAINS = "light"
        with self.assertRaises(ValueError):
            validate_and_process_chains(make_args(), [])

    def test_uniform_is_unchanged(self):
        """In a uniform simulation, the light chain settings are not changed."""
        s.UNIFORM = True
        s.TARGET_MUTATIONS_LIGHT = 4
        warnings = []
        validate_and_process_chains(make_args(uniform=True, single_chain=True), warnings)
        self.assertEqual(s.TARGET_MUTATIONS_LIGHT, 4)
        self.assertEqual(len(warnings), 1)
        # CHAINS goes back to the default, so the uniform target and input
        # handling are the same as without --single-chain
        self.assertEqual(s.CHAINS, "paired")
        self.assertFalse(s.SINGLE_CHAIN)


class TestSingleChainTarget(SettingsSaver):
    """A user-supplied target in a single chain simulation."""

    @classmethod
    def setUpClass(cls):
        cls.folder = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        naive = pd.read_csv(helper.get_data("naive_pairs_filtered.csv"), header=0)
        naive = naive.rename(columns={
            "heavy": "heavy_sequence",
            "light": "light_sequence",
            "heavy_aligned": "heavy_sequence_alignment",
            "light_aligned": "light_sequence_alignment",
            })
        cls.naive_file = os.path.join(cls.folder.name, "naive.csv")
        cls.target_file = os.path.join(cls.folder.name, "target.csv")
        naive.head(300).to_csv(cls.naive_file, index=False)
        naive.head(1).to_csv(cls.target_file, index=False)

    @classmethod
    def tearDownClass(cls):
        cls.folder.cleanup()

    def _pool_size(self, chains):
        s.CHAINS = chains
        s.TARGET = {}
        s.NAIVE_FILE = self.naive_file
        validate_and_process_target_input(make_args(target=self.target_file), [])
        helper.update_helper_tables()
        return helper.NAIVE_ROWS

    def test_target_has_no_light_chain(self):
        self._pool_size("heavy")
        self.assertEqual(s.TARGET["light_aligned"], "")
        self.assertEqual(s.TARGET["light_cdr3_aa_length"], 0)
        self.assertEqual(
            sorted(s.TARGET["signature"].keys()),
            ["heavy_aligned", "heavy_junction", "heavy_locus", "heavy_v"]
            )

    def test_only_heavy_chain_restricts_the_pool(self):
        """The light chain is not used, so it does not restrict the naive pool."""
        paired = self._pool_size("paired")
        single = self._pool_size("heavy")
        self.assertGreaterEqual(paired, 1)
        self.assertGreaterEqual(single, paired)


if __name__ == '__main__':
    unittest.main()
