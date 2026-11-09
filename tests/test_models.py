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
import shutil
import tempfile
import unittest

import numpy as np
import pandas as pd

from simble import helper, models
from simble.chain import HeavyChain
from simble.parsing import validate_and_process_model_files
from simble.settings import Settings, s
from simble.simble import do_simulation
from tests.test_single_chain import SettingsSaver

ENTROPY = 54897022524486695084299880814690718190
MODEL_SETTINGS = ["MUTABILITY_HEAVY_FILE", "SUBSTITUTION_HEAVY_FILE",
                  "MUTABILITY_LIGHT_FILE", "SUBSTITUTION_LIGHT_FILE"]


def bundled_tables(chain="heavy"):
    """Returns the bundled mutability and substitution files as raw tables."""
    prefix = "hh" if chain == "heavy" else "hkl"
    return (pd.read_csv(helper.get_data(f"{prefix}_sf5.csv")),
            pd.read_csv(helper.get_data(f"{prefix}_sf5_substitution.csv")))


def apobec_tables():
    """Returns a toy APOBEC-like model: only C in a TC context mutates, always to T."""
    mutability = pd.DataFrame({"Fivemer": models.ACGT_KMERS})
    mutability["Mutability"] = [
        1.0 if x[1] == "T" and x[2] == "C" else 0.0 for x in models.ACGT_KMERS
        ]
    substitution = pd.DataFrame({"Fivemer": models.ACGT_KMERS})
    for nucleotide in models.NUCLEOTIDES:
        substitution[nucleotide] = [
            1.0 if (x[1] == "T" and x[2] == "C" and nucleotide == "T") else 0.0
            for x in models.ACGT_KMERS
            ]
    # leave out the rows of the 5-mers that cannot mutate; they are not needed
    substitution = substitution[mutability["Mutability"] > 0]
    return mutability, substitution


def run_clone(chains="paired", neutral=False, **model_files):
    """Runs one short clone with a fixed seed and returns its results."""
    settings = Settings()
    settings.LOCATIONS[0].sample_times = list(range(0, 9, 4))
    settings.LOCATIONS[1].sample_times = list(range(0, 9, 4))
    settings.LOCATIONS[0].migration_rate = 0.5
    settings.CHAINS = chains
    if chains == "heavy":
        settings.TARGET_MUTATIONS_LIGHT = 0
    if neutral:
        settings.SELECTION = False
    for key, value in model_files.items():
        setattr(settings, key, value)
    seed = np.random.SeedSequence(ENTROPY)
    with tempfile.NamedTemporaryFile(mode="w") as tmpf:
        json.dump(settings, tmpf, default=lambda o: o.encode(), indent=4)
        tmpf.flush()
        return do_simulation(1, seed, None, tmpf.name)


class ModelFiles(SettingsSaver):
    """Writes model files in a temporary folder."""
    FIELDS = SettingsSaver.FIELDS + MODEL_SETTINGS

    @classmethod
    def setUpClass(cls):
        cls.folder = tempfile.TemporaryDirectory()  # pylint: disable=consider-using-with

    @classmethod
    def tearDownClass(cls):
        cls.folder.cleanup()

    def write(self, table, name):
        path = os.path.join(self.folder.name, name)
        table.to_csv(path, index=False)
        return path

    def copy_bundled(self, chain="heavy"):
        prefix = "hh" if chain == "heavy" else "hkl"
        mutability = os.path.join(self.folder.name, f"{chain}_mutability.csv")
        substitution = os.path.join(self.folder.name, f"{chain}_substitution.csv")
        shutil.copy(helper.get_data(f"{prefix}_sf5.csv"), mutability)
        shutil.copy(helper.get_data(f"{prefix}_sf5_substitution.csv"), substitution)
        return mutability, substitution


class TestReadModelFiles(ModelFiles):
    """Reading the bundled files as user files gives the bundled values."""

    def test_bundled_files_are_valid(self):
        for chain in ["heavy", "light"]:
            with self.subTest(chain=chain):
                mutability_file, substitution_file = self.copy_bundled(chain)
                mutability = models.read_mutability(mutability_file)
                substitution = models.read_substitution(substitution_file)
                models.validate_model(mutability, substitution, chain)
                bundled = models.bundled_as_dicts(*helper.BUNDLED_MODEL_TABLES[chain])
                self.assertEqual(mutability, bundled[0])
                self.assertEqual(substitution, bundled[1])

    def test_n_rows_are_made_like_the_bundled_files(self):
        """Without N rows, the N rows are made the same way as in the bundled files."""
        raw_mutability, raw_substitution = bundled_tables()
        acgt = ~raw_mutability["Fivemer"].str.contains("N")
        mutability = models.read_mutability(
            self.write(raw_mutability[acgt], "mutability_no_n.csv"))
        acgt = ~raw_substitution["Fivemer"].str.contains("N")
        substitution = models.read_substitution(
            self.write(raw_substitution[acgt], "substitution_no_n.csv"))
        bundled_mutability, bundled_substitution = models.bundled_as_dicts(
            *helper.BUNDLED_MODEL_TABLES["heavy"])
        self.assertEqual(len(mutability), 3125)
        for kmer in models.N_KMERS:
            self.assertAlmostEqual(mutability[kmer], bundled_mutability[kmer], places=15)
            np.testing.assert_allclose(
                substitution[kmer], bundled_substitution[kmer], rtol=0, atol=1e-12)

    def test_tsv_and_lower_case(self):
        mutability, _ = bundled_tables()
        mutability["Fivemer"] = mutability["Fivemer"].str.lower()
        path = os.path.join(self.folder.name, "mutability.tsv")
        mutability.to_csv(path, sep="\t", index=False)
        self.assertEqual(len(models.read_mutability(path)), 3125)


class TestModelFileErrors(ModelFiles):
    """Each kind of bad model file gives a clear error."""

    def assert_error(self, function, table, text):
        path = self.write(table, "bad.csv")
        with self.assertRaises(models.ModelFileError) as context:
            function(path)
        self.assertIn(text, str(context.exception))

    def test_missing_file(self):
        with self.assertRaises(models.ModelFileError):
            models.read_mutability(os.path.join(self.folder.name, "nothing.csv"))

    def test_missing_column(self):
        mutability, _ = bundled_tables()
        self.assert_error(models.read_mutability,
                          mutability.rename(columns={"Mutability": "rate"}), "no column")

    def test_bad_kmer(self):
        mutability, _ = bundled_tables()
        mutability.loc[0, "Fivemer"] = "AAXAA"
        self.assert_error(models.read_mutability, mutability, "AAXAA")
        mutability.loc[0, "Fivemer"] = "AAAA"
        self.assert_error(models.read_mutability, mutability, "AAAA")

    def test_duplicate_kmer(self):
        mutability, _ = bundled_tables()
        mutability.loc[1, "Fivemer"] = mutability.loc[0, "Fivemer"]
        self.assert_error(models.read_mutability, mutability, "duplicate")

    def test_negative_and_not_a_number(self):
        mutability, _ = bundled_tables()
        mutability["Mutability"] = mutability["Mutability"].astype(object)
        mutability.loc[0, "Mutability"] = -1
        self.assert_error(models.read_mutability, mutability, "negative")
        mutability.loc[0, "Mutability"] = "high"
        self.assert_error(models.read_mutability, mutability, "not numbers")

    def test_missing_acgt_kmer(self):
        mutability, _ = bundled_tables()
        self.assert_error(models.read_mutability, mutability[mutability["Fivemer"] != "ACGTA"],
                          "ACGTA")

    def test_all_zero(self):
        mutability, _ = bundled_tables()
        mutability["Mutability"] = 0
        self.assert_error(models.read_mutability, mutability, "greater than 0")

    def test_self_substitution(self):
        mutability = models.read_mutability(self.copy_bundled()[0])
        _, raw = bundled_tables()
        raw.loc[raw["Fivemer"] == "AAAAA", "A"] = 0.1
        substitution = models.read_substitution(self.write(raw, "self.csv"))
        with self.assertRaises(models.ModelFileError) as context:
            models.validate_model(mutability, substitution, "test")
        self.assertIn("AAAAA", str(context.exception))
        self.assertIn("stay the same", str(context.exception))

    def test_row_sum(self):
        mutability = models.read_mutability(self.copy_bundled()[0])
        _, raw = bundled_tables()
        raw.loc[raw["Fivemer"] == "AAAAC", "C"] = 0.9
        substitution = models.read_substitution(self.write(raw, "sum.csv"))
        with self.assertRaises(models.ModelFileError) as context:
            models.validate_model(mutability, substitution, "test")
        self.assertIn("AAAAC", str(context.exception))

    def test_missing_substitution_row_for_mutable_kmer(self):
        mutability = models.read_mutability(self.copy_bundled()[0])
        _, raw = bundled_tables()
        substitution = models.read_substitution(
            self.write(raw[raw["Fivemer"] != "GGCGG"], "missing.csv"))
        with self.assertRaises(models.ModelFileError) as context:
            models.validate_model(mutability, substitution, "test")
        self.assertIn("GGCGG", str(context.exception))

    def test_rounding_is_renormalized(self):
        """A row a little away from 1 (rounding) is accepted and renormalized."""
        mutability = models.read_mutability(self.copy_bundled()[0])
        _, raw = bundled_tables()
        raw.loc[raw["Fivemer"] == "AAAAC", "C"] += 5e-7
        substitution = models.read_substitution(self.write(raw, "rounding.csv"))
        models.validate_model(mutability, substitution, "test")
        self.assertAlmostEqual(sum(substitution["AAAAC"]), 1, places=14)

    def test_unused_rows_do_not_spoil_n_rows(self):
        """A row of a 5-mer that cannot mutate is not checked, and it is not
        used to make the N rows either."""
        mutability, substitution = apobec_tables()
        # AACAA cannot mutate (its C is not after a T), so its row is not used.
        # Give it a row that sums to 1 but keeps the C.
        bad = pd.DataFrame([{"Fivemer": "AACAA", "A": 0.0, "C": 1.0, "G": 0.0, "T": 0.0}])
        substitution = pd.concat([substitution, bad])
        mutability = models.read_mutability(self.write(mutability, "unused_m.csv"))
        substitution = models.read_substitution(self.write(substitution, "unused_s.csv"))
        # NNCAA can mutate (some of its expansions are TC), so its row is checked;
        # it comes only from the valid rows
        self.assertGreater(mutability["NNCAA"], 0)
        self.assertEqual(substitution["NNCAA"], [0.0, 0.0, 0.0, 1.0])
        models.validate_model(mutability, substitution, "test")

    def test_toy_model_is_valid(self):
        """Rows can be missing for 5-mers that cannot mutate."""
        mutability, substitution = apobec_tables()
        mutability = models.read_mutability(self.write(mutability, "apobec_m.csv"))
        substitution = models.read_substitution(self.write(substitution, "apobec_s.csv"))
        models.validate_model(mutability, substitution, "test")
        # an edge 5-mer with N that can mutate gets a valid substitution row
        self.assertGreater(mutability["NTCAA"], 0)
        self.assertEqual(substitution["NTCAA"], [0.0, 0.0, 0.0, 1.0])


class TestModelSettings(ModelFiles):
    """The command line arguments for model files."""

    def make_args(self, **kwargs):
        from tests.test_single_chain import make_args  # pylint: disable=import-outside-toplevel
        defaults = {x: None for x in [
            "mutability", "substitution", "mutability_heavy", "substitution_heavy",
            "mutability_light", "substitution_light"]}
        defaults.update(kwargs)
        return make_args(**defaults)

    def test_shortcut_sets_both_chains(self):
        mutability, substitution = self.copy_bundled()
        validate_and_process_model_files(
            self.make_args(mutability=mutability, substitution=substitution), [])
        self.assertEqual(s.MUTABILITY_HEAVY_FILE, mutability)
        self.assertEqual(s.MUTABILITY_LIGHT_FILE, mutability)
        self.assertEqual(s.SUBSTITUTION_HEAVY_FILE, substitution)
        self.assertEqual(s.SUBSTITUTION_LIGHT_FILE, substitution)

    def test_specific_overrides_shortcut(self):
        heavy, _ = self.copy_bundled("heavy")
        light, _ = self.copy_bundled("light")
        warnings = []
        validate_and_process_model_files(
            self.make_args(mutability=heavy, mutability_light=light), warnings)
        self.assertEqual(s.MUTABILITY_HEAVY_FILE, heavy)
        self.assertEqual(s.MUTABILITY_LIGHT_FILE, light)
        self.assertEqual(len(warnings), 1)

    def test_mutability_only(self):
        """A mutability file alone is used with the bundled substitution file."""
        mutability, _ = self.copy_bundled()
        validate_and_process_model_files(self.make_args(mutability_heavy=mutability), [])
        self.assertEqual(s.MUTABILITY_HEAVY_FILE, mutability)
        self.assertIsNone(s.SUBSTITUTION_HEAVY_FILE)

    def test_bad_file_stops_the_run(self):
        mutability, _ = bundled_tables()
        path = self.write(mutability[mutability["Fivemer"] != "ACGTA"], "short.csv")
        with self.assertRaises(ValueError):
            validate_and_process_model_files(self.make_args(mutability=path), [])
        with self.assertRaises(FileNotFoundError):
            validate_and_process_model_files(
                self.make_args(mutability=os.path.join(self.folder.name, "nothing.csv")), [])

    def test_uniform_is_an_error(self):
        mutability, _ = self.copy_bundled()
        s.UNIFORM = True
        with self.assertRaises(ValueError):
            validate_and_process_model_files(self.make_args(mutability=mutability), [])

    def test_single_chain_ignores_light_files(self):
        s.CHAINS = "heavy"
        warnings = []
        validate_and_process_model_files(
            self.make_args(mutability_light=os.path.join(self.folder.name, "nothing.csv")),
            warnings)
        self.assertEqual(len(warnings), 1)

    def test_no_files_is_unchanged(self):
        warnings = []
        validate_and_process_model_files(self.make_args(), warnings)
        self.assertEqual(warnings, [])
        for setting in MODEL_SETTINGS:
            self.assertIsNone(getattr(s, setting))


class TestModelSimulation(ModelFiles):
    """Simulations with model files."""

    def test_bundled_copies_give_the_same_simulation(self):
        """Model files with the bundled values give the same seeded results."""
        heavy_mutability, heavy_substitution = self.copy_bundled("heavy")
        light_mutability, light_substitution = self.copy_bundled("light")
        default = run_clone()
        from_files = run_clone(
            MUTABILITY_HEAVY_FILE=heavy_mutability,
            SUBSTITUTION_HEAVY_FILE=heavy_substitution,
            MUTABILITY_LIGHT_FILE=light_mutability,
            SUBSTITUTION_LIGHT_FILE=light_substitution)
        drop = ["sequence_id", "cell_id"]
        pd.testing.assert_frame_equal(
            default["airr"].drop(columns=drop), from_files["airr"].drop(columns=drop))
        pd.testing.assert_frame_equal(default["pop_data"], from_files["pop_data"])

    def test_toy_apobec_model(self):
        """With the toy model, every mutation is C to T."""
        mutability, substitution = apobec_tables()
        result = run_clone(
            chains="heavy", neutral=True,
            MUTABILITY_HEAVY_FILE=self.write(mutability, "apobec_m.csv"),
            SUBSTITUTION_HEAVY_FILE=self.write(substitution, "apobec_s.csv"))
        changes = {}
        for row in result["airr"].itertuples():
            for germline, observed in zip(row.germline_alignment, row.sequence_alignment):
                if germline != observed:
                    changes[(germline, observed)] = changes.get((germline, observed), 0) + 1
        self.assertGreater(sum(changes.values()), 0)
        self.assertEqual(set(changes), {("C", "T")})

    def test_context_comes_from_current_sequence(self):
        """A second mutation in the same division uses the sequence after the first.

        In ATCCA, only the first C (after T) can mutate. When it becomes T, the
        second C is after a T, so it can mutate too. Before this fix, the second
        substitution used the old context (a C before it), which has no
        substitution row in this model, and the simulation failed.
        """
        mutability, substitution = apobec_tables()
        s.MUTABILITY_HEAVY_FILE = self.write(mutability, "apobec_m.csv")
        s.SUBSTITUTION_HEAVY_FILE = self.write(substitution, "apobec_s.csv")
        helper.update_model_tables()
        s._x_RNG = np.random.default_rng(1)  # pylint: disable=protected-access
        try:
            chain = HeavyChain("ATCCA", gapped_seq="ATCCA")
            self.assertEqual(chain.mutate(1, n=2), 2)
            self.assertEqual(chain.nucleotide_seq, "ATTTA")
        finally:
            s._x_RNG = None  # pylint: disable=protected-access

    def test_no_mutable_positions(self):
        """A chain with no mutable 5-mer does not mutate, and does not fail."""
        mutability, substitution = apobec_tables()
        s.MUTABILITY_HEAVY_FILE = self.write(mutability, "apobec_m.csv")
        s.SUBSTITUTION_HEAVY_FILE = self.write(substitution, "apobec_s.csv")
        helper.update_model_tables()
        s._x_RNG = np.random.default_rng(1)  # pylint: disable=protected-access
        try:
            chain = HeavyChain("AAAAGGGGAAAAGGGG", gapped_seq="AAAAGGGGAAAAGGGG")
            self.assertEqual(sum(chain.mutability_map), 0)
            self.assertEqual(chain.mutate(1, n=3), 0)
            self.assertEqual(chain.nucleotide_seq, "AAAAGGGGAAAAGGGG")
        finally:
            s._x_RNG = None  # pylint: disable=protected-access


if __name__ == '__main__':
    unittest.main()
