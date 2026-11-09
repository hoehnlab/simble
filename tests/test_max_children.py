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
import tempfile
import unittest

import numpy as np

from simble.parsing import validate_and_process_max_children
from simble.settings import Settings, s
from simble.simble import do_simulation
from tests.test_internal_nodes import parse_newick

ENTROPY = 54897022524486695084299880814690718190


def run_clone(max_children=None):
    """Runs one short clone with a fixed seed and keeps the full tree."""
    settings = Settings()
    settings.LOCATIONS[0].sample_times = list(range(0, 7, 3))
    settings.LOCATIONS[1].sample_times = list(range(0, 7, 3))
    settings.KEEP_FULL_TREE = True
    if max_children is not None:
        settings.MAX_CHILDREN = max_children
    seed = np.random.SeedSequence(ENTROPY)
    with tempfile.NamedTemporaryFile(mode="w") as tmpf:
        json.dump(settings, tmpf, default=lambda o: o.encode(), indent=4)
        tmpf.flush()
        return do_simulation(1, seed, None, tmpf.name)


def children_columns(pop_data):
    """Returns the number_of_cells_with_N_children columns of the population data."""
    return [x for x in pop_data.columns if x.startswith("number_of_cells_with_")]


class TestMaxChildren(unittest.TestCase):
    """Tests for the maximum number of children of a GC cell."""

    def test_default_columns(self):
        """The default is 100, so there are 100 children columns."""
        result = run_clone()
        columns = children_columns(result["pop_data"])
        self.assertEqual(len(columns), 100)
        self.assertEqual(columns[0], "number_of_cells_with_1_children")
        self.assertEqual(columns[-1], "number_of_cells_with_100_children")

    def test_limit_is_applied(self):
        """No cell has more children than the limit, and some cells reach it."""
        result = run_clone(max_children=3)
        columns = children_columns(result["pop_data"])
        self.assertEqual(
            columns,
            [f"number_of_cells_with_{i}_children" for i in range(1, 4)]
            )
        # there is no migration, so all children come from GC reproduction
        nodes = parse_newick(result["full_tree"])
        most_children = max(len(x["children"]) for x in nodes.values())
        self.assertEqual(most_children, 3)


class TestValidateMaxChildren(unittest.TestCase):
    """Tests for the validation of the maximum number of children."""

    def setUp(self):
        self.saved = s.MAX_CHILDREN

    def tearDown(self):
        s.MAX_CHILDREN = self.saved

    def test_valid(self):
        warnings = []
        validate_and_process_max_children(7, warnings)
        self.assertEqual(s.MAX_CHILDREN, 7)
        self.assertEqual(warnings, [])

    def test_none_keeps_setting(self):
        s.MAX_CHILDREN = 12
        validate_and_process_max_children(None, [])
        self.assertEqual(s.MAX_CHILDREN, 12)

    def test_integer_float_from_config(self):
        """A config file can give 5.0; it becomes the integer 5."""
        s.MAX_CHILDREN = 5.0
        validate_and_process_max_children(None, [])
        self.assertEqual(s.MAX_CHILDREN, 5)
        self.assertIsInstance(s.MAX_CHILDREN, int)

    def test_invalid(self):
        for value in [0, -1, 2.5, "3", True]:
            with self.subTest(value=value):
                with self.assertRaises(ValueError):
                    validate_and_process_max_children(value, [])

    def test_warns_above_max_population(self):
        warnings = []
        validate_and_process_max_children(s.LOCATIONS[0].max_population + 1, warnings)
        self.assertEqual(len(warnings), 1)
        self.assertIn("has no effect", warnings[0])


if __name__ == '__main__':
    unittest.main()
