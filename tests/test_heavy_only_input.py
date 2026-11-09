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

import json
import os
import tempfile
import unittest

import numpy as np
import pandas as pd

from simble import helper
from simble.data.process_naive import create_naive_table
from simble.parsing import (read_target_table, validate_and_process_naive_input,
                            validate_and_process_target_input)
from simble.settings import Settings, s
from simble.simble import do_simulation
from tests.test_single_chain import SettingsSaver, make_args

ENTROPY = 54897022524486695084299880814690718190
N_ROWS = 30


def write_airr_tables(folder):
    """Writes heavy and light AIRR tsv files made from the bundled naive pairs.

    Returns:
        tuple: The heavy file, the light file, and the bundled rows used.
    """
    naive = pd.read_csv(helper.get_data("naive_pairs_filtered.csv"), header=0).head(N_ROWS)
    # the bundled cell ids are not unique, so give each pair its own id
    naive["cell_id"] = [f"cell{i}" for i in range(len(naive))]
    files = []
    for chain in ["heavy", "light"]:
        columns = {f"{chain}_{x}": x for x in [
            c[len(chain) + 1:] for c in naive.columns if c.startswith(f"{chain}_")
            ]}
        columns[chain] = "sequence"
        columns[f"{chain}_aligned"] = "sequence_alignment"
        table = naive[list(columns.keys()) + ["cell_id"]].rename(columns=columns)
        table["productive"] = "T"
        path = os.path.join(folder, f"{chain}.tsv")
        table.to_csv(path, sep="\t", index=False)
        files.append(path)
    return files[0], files[1], naive


class InputFiles(SettingsSaver):
    """Makes paired and heavy chain only input files once for the test class."""

    @classmethod
    def setUpClass(cls):
        cls.folder = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with
        cls.heavy_tsv, cls.light_tsv, cls.bundled = write_airr_tables(cls.folder.name)
        cls.paired = create_naive_table(cls.heavy_tsv, cls.light_tsv, "cell_id")
        cls.heavy_only = create_naive_table(cls.heavy_tsv, None, "cell_id")
        cls.paired_file = os.path.join(cls.folder.name, "paired.csv")
        cls.heavy_file = os.path.join(cls.folder.name, "heavy_only.csv")
        cls.paired.to_csv(cls.paired_file, index=False)
        cls.heavy_only.to_csv(cls.heavy_file, index=False)
        cls.paired_target = os.path.join(cls.folder.name, "paired_target.csv")
        cls.heavy_target = os.path.join(cls.folder.name, "heavy_target.csv")
        cls.paired.head(1).to_csv(cls.paired_target, index=False)
        cls.heavy_only.head(1).to_csv(cls.heavy_target, index=False)
        # one-row AIRR tables for --target-heavy and --target-light
        cls.heavy_target_tsv = os.path.join(cls.folder.name, "heavy_target.tsv")
        cls.light_target_tsv = os.path.join(cls.folder.name, "light_target.tsv")
        pd.read_csv(cls.heavy_tsv, sep="\t").head(1).to_csv(
            cls.heavy_target_tsv, sep="\t", index=False)
        pd.read_csv(cls.light_tsv, sep="\t").head(1).to_csv(
            cls.light_target_tsv, sep="\t", index=False)

    @classmethod
    def tearDownClass(cls):
        cls.folder.cleanup()


class TestProcessNaive(InputFiles):
    """process_naive with and without a light chain table."""

    def test_paired_table(self):
        self.assertEqual(len(self.paired), N_ROWS)
        self.assertTrue(any(x.startswith("light_") for x in self.paired.columns))
        self.assertIn("cell_id", self.paired.columns)

    def test_heavy_only_columns_match_paired(self):
        """The heavy chain only table has the same heavy columns as the paired table."""
        paired_heavy = [x for x in self.paired.columns if not x.startswith("light_")]
        self.assertEqual(sorted(self.heavy_only.columns), sorted(paired_heavy))
        self.assertFalse(any(x.startswith("light_") for x in self.heavy_only.columns))

    def test_heavy_only_values_match_paired(self):
        columns = list(self.heavy_only.columns)
        pd.testing.assert_frame_equal(
            self.heavy_only[columns].reset_index(drop=True),
            self.paired[columns].reset_index(drop=True)
            )

    def test_no_join_column_needed(self):
        table = pd.read_csv(self.heavy_tsv, sep="\t").drop(columns=["cell_id"])
        path = os.path.join(self.folder.name, "heavy_no_id.tsv")
        table.to_csv(path, sep="\t", index=False)
        heavy_only = create_naive_table(path, None, "cell_id")
        self.assertEqual(len(heavy_only), N_ROWS)
        self.assertNotIn("cell_id", heavy_only.columns)

    def test_productive_filter(self):
        table = pd.read_csv(self.heavy_tsv, sep="\t")
        table.loc[0, "productive"] = "F"
        path = os.path.join(self.folder.name, "heavy_one_unproductive.tsv")
        table.to_csv(path, sep="\t", index=False)
        self.assertEqual(len(create_naive_table(path, None, "cell_id")), N_ROWS - 1)

    def test_keep_cols(self):
        table = pd.read_csv(self.heavy_tsv, sep="\t")
        table["cell_label"] = "x"
        path = os.path.join(self.folder.name, "heavy_label.tsv")
        table.to_csv(path, sep="\t", index=False)
        heavy_only = create_naive_table(path, None, "cell_id", ["cell_label"], "heavy")
        self.assertIn("cell_label", heavy_only.columns)
        with self.assertRaises(ValueError):
            create_naive_table(path, None, "cell_id", ["cell_label"], "light")


class TestHeavyOnlyNaiveInput(InputFiles):
    """simble with a heavy chain only naive table."""

    def test_paired_mode_needs_light_columns(self):
        with self.assertRaises(ValueError) as context:
            validate_and_process_naive_input(
                make_args(naive=self.heavy_file, keep_cols=None), [])
        self.assertIn("--single-chain", str(context.exception))

    def test_single_chain_mode_accepts_heavy_only(self):
        s.CHAINS = "heavy"
        validate_and_process_naive_input(make_args(naive=self.heavy_file, keep_cols=None), [])
        self.assertEqual(s.NAIVE_FILE, self.heavy_file)

    def test_paired_input_still_valid(self):
        validate_and_process_naive_input(make_args(naive=self.paired_file, keep_cols=None), [])
        s.CHAINS = "heavy"
        validate_and_process_naive_input(make_args(naive=self.paired_file, keep_cols=None), [])


class TestHeavyOnlySimulation(InputFiles):
    """A heavy chain only naive table gives the same simulation as a paired table."""

    def _run(self, naive_file):
        settings = Settings()
        settings.LOCATIONS[0].sample_times = list(range(0, 9, 4))
        settings.LOCATIONS[1].sample_times = list(range(0, 9, 4))
        settings.LOCATIONS[0].migration_rate = 0.5
        settings.CHAINS = "heavy"
        settings.TARGET_MUTATIONS_LIGHT = 0
        settings.NAIVE_FILE = naive_file
        settings.NAIVE_SELECTION = "ordered"
        seed = np.random.SeedSequence(ENTROPY)
        with tempfile.NamedTemporaryFile(mode="w") as tmpf:
            json.dump(settings, tmpf, default=lambda o: o.encode(), indent=4)
            tmpf.flush()
            return do_simulation(3, seed, None, tmpf.name)

    def test_same_as_paired(self):
        from_paired = self._run(self.paired_file)
        from_heavy = self._run(self.heavy_file)
        drop = ["sequence_id", "cell_id"]
        pd.testing.assert_frame_equal(
            from_paired["airr"].drop(columns=drop).reset_index(drop=True),
            from_heavy["airr"].drop(columns=drop).reset_index(drop=True)
            )
        pd.testing.assert_frame_equal(from_paired["pop_data"], from_heavy["pop_data"])
        self.assertEqual(from_paired["targets"], from_heavy["targets"])


class TestHeavyOnlyTarget(InputFiles):
    """Heavy chain only targets in a single chain simulation."""

    def _target(self, **kwargs):
        s.TARGET = {}
        warnings = []
        validate_and_process_target_input(make_args(**kwargs), warnings)
        return dict(s.TARGET), warnings

    def test_heavy_only_target_file(self):
        s.CHAINS = "heavy"
        from_paired, _ = self._target(target=self.paired_target)
        from_heavy, _ = self._target(target=self.heavy_target)
        self.assertEqual(from_paired, from_heavy)
        self.assertEqual(from_heavy["light_aligned"], "")

    def test_target_heavy_alone(self):
        s.CHAINS = "heavy"
        from_paired, _ = self._target(target=self.paired_target)
        from_airr, warnings = self._target(target_heavy=self.heavy_target_tsv)
        self.assertEqual(from_paired, from_airr)
        self.assertEqual(warnings, [])

    def test_target_light_ignored(self):
        s.CHAINS = "heavy"
        alone, _ = self._target(target_heavy=self.heavy_target_tsv)
        both, warnings = self._target(target_heavy=self.heavy_target_tsv, target_light=self.light_target_tsv)
        self.assertEqual(alone, both)
        self.assertEqual(len(warnings), 1)
        self.assertIn("--target-light", warnings[0])

    def test_target_light_alone_is_an_error(self):
        s.CHAINS = "heavy"
        with self.assertRaises(ValueError):
            read_target_table(make_args(target_light=self.light_target_tsv))

    def test_paired_mode_needs_both(self):
        with self.assertRaises(ValueError) as context:
            read_target_table(make_args(target_heavy=self.heavy_target_tsv))
        self.assertIn("--single-chain", str(context.exception))
        with self.assertRaises(ValueError):
            self._target(target=self.heavy_target)

    def test_paired_mode_unchanged(self):
        """A paired target in a paired simulation keeps both chains."""
        target, _ = self._target(target_heavy=self.heavy_target_tsv, target_light=self.light_target_tsv)
        self.assertNotEqual(target["light_aligned"], "")
        self.assertIn("light_locus", target["signature"])
        from_file, _ = self._target(target=self.paired_target)
        self.assertEqual(target, from_file)


if __name__ == '__main__':
    unittest.main()
