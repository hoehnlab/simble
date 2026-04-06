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
import pandas as pd

from simble.settings import Settings
from simble.simble import do_simulation


class TestSeeds(unittest.TestCase):
    """Test case for simble simulation with specific seeds."""
    def assert_dataframe_equal(self, a, b, msg):
        """Custom assertion to compare two DataFrames."""
        try:
            pd.testing.assert_frame_equal(a, b)
        except AssertionError as e:
            raise self.failureException(msg) from e

    def setUp(self):
        self.addTypeEqualityFunc(pd.DataFrame, self.assert_dataframe_equal)


    def test_main(self):
        """Test the main simulation with a specific seed."""
        clean_settings = Settings()
        clean_settings.LOCATIONS[0].sample_times = list(range(0, 10, 5))
        clean_settings.LOCATIONS[1].sample_times = list(range(0, 10, 5))
        entropy = 54897022524486695084299880814690718190
        seed = np.random.SeedSequence(entropy)
        with tempfile.NamedTemporaryFile(mode="w") as tmpf:
            json.dump(clean_settings, tmpf, default=lambda o: o.encode(), indent=4)
            tmpf.flush()
            result1 = do_simulation(1, seed, tmpf.name)

        clean_settings = Settings()
        clean_settings.LOCATIONS[0].sample_times = list(range(0, 10, 5))
        clean_settings.LOCATIONS[1].sample_times = list(range(0, 10, 5))
        with tempfile.NamedTemporaryFile(mode="w") as tmpf:
            json.dump(clean_settings, tmpf, default=lambda o: o.encode(), indent=4)
            tmpf.flush()
            result2 = do_simulation(1, seed, tmpf.name)

        drop = ["sequence_id", "cell_id"]
        self.assertEqual(result1["data"], result2["data"])
        self.assertEqual(result1["pop_data"], result2["pop_data"])
        self.assertEqual(result1["clone_id"], result2["clone_id"])
        self.assertEqual(result1["airr"].drop(columns=drop), result2["airr"].drop(columns=drop))
        # cell_ids will not match so we'd need to write custom function to compare trees,
        # but if these dataframes are exactly equal it should be good enough for now

if __name__ == '__main__':
    unittest.main()
