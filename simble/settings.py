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

import logging
from .location import LocationName

logger = logging.getLogger(__package__)


class Encodable():
    """Base class for objects that can be encoded to a dictionary."""
    def encode(self):
        """Encodes the object as a dictionary for serialization."""
        return vars(self)

class LocationSettings(Encodable):
    """Settings for a specific location in the simulation.
    Attributes:
        name (LocationName): The name of the location.
        sample_times (list): Times at which samples are taken.
        mutation_rate (float): The mutation rate for the location.
        max_population (int): The maximum population allowed in the location.
        migration_rate (float): The rate of migration out of this location.
        sample_size (int): The number of cells to sample from the location.
    """
    def __init__(self,
                 name,
                 sample_times=None,
                 mutation_rate=None,
                 max_population=1000,
                 migration_rate=0,
                 sample_size=None):
        self.name = name
        if name == LocationName.GC:
            defaults = {
                "sample_times": list(range(0, 201, 25)),
                "mutation_rate": 1.0,
                "sample_size": 50
            }
        elif name == LocationName.OTHER:
            defaults = {
                "sample_times": [],
                "mutation_rate": 0.0,
                "sample_size": 12
            }
        else:
            raise ValueError(f"Invalid location name: {name}")
        self.sample_times = sample_times if sample_times else defaults["sample_times"]
        self.mutation_rate = mutation_rate if mutation_rate else defaults["mutation_rate"]
        self.max_population = max_population
        self.migration_rate = migration_rate
        self.sample_size = sample_size if sample_size else defaults["sample_size"]

    def __repr__(self):
        return repr(vars(self))


class Settings(Encodable):
    """Global settings for the simulation.
    Attributes:
        LOCATIONS (list): List of LocationSettings for different locations.
        HEAVY_SHM_PER_SITE (float): SHM rate per site of the heavy chain.
        LIGHT_SHM_PER_SITE (float): SHM rate per site of the light chain.
        TARGET_MUTATIONS_HEAVY (int): Number of target mutations for heavy chain.
        TARGET_MUTATIONS_LIGHT (int): Number of target mutations for light chain.
        SELECTION (bool): Whether selection is applied in the simulation.
        UNIFORM (bool): Whether to use a uniform mutation model.
        RESULTS_DIR (str): Directory for saving results.
        MULTIPLIER (float): Multiplier for affinity calculations.
        _x_RNG (random.Random): Random number generator instance.
        DEV (bool): Development mode flag.
        FASTA (bool): Whether to output results in FASTA format.
        VERBOSE (bool): Verbosity level for logging.
        CDR_DIST (str): Distribution type for CDR mutations.
        CDR_VAR (float): Variance for CDR mutations.
        FWR_DIST (str): Distribution type for FWR mutations.
        FWR_VAR (float): Variance for FWR mutations.
        TIME_SWITCH (int): Time at which to switch locations.
        GENERATIONS_PER_DAY (float): Number of generations per day.
        MEMORY_SAVE (bool): Whether to save memory during the simulation.
        KEEP_FULL_TREE (bool): Whether to keep the full tree of cells.
        QUIET (bool): Whether to suppress output.
    """
    # pylint: disable=invalid-name

    def __init__(self):
        # pylint: disable=invalid-name
        self.LOCATIONS = [
            LocationSettings(
                name=LocationName.GC,
                mutation_rate=1.0,
                max_population=1000,
                migration_rate=0.0,
                sample_size=50),
            LocationSettings(
                name=LocationName.OTHER,
                mutation_rate=0.0,
                max_population=1000,
                migration_rate=0.0,
                sample_size=12)
                ]
        # Based on an average heavy chain length of 370, the calculated value is:
        # 0.33/370 = 0.0008908272571108565
        # 0.16/325 (avg light chain length) = 0.0004923076923076923
        self.HEAVY_SHM_PER_SITE = 0.0008908272571108565
        self.LIGHT_SHM_PER_SITE = 0.0004923076923076923
        self.TARGET_MUTATIONS_HEAVY = 5
        self.TARGET_MUTATIONS_LIGHT = 2
        self._UNIFORM = False
        self.RESULTS_DIR = ""
        self.MULTIPLIER = 2
        self._x_RNG = None # pylint: disable=protected-access
        self.DEV = False
        self.FASTA = False
        self.VERBOSE = False
        self.CDR_DIST = "exponential"
        self.CDR_VAR = 0.995
        self.FWR_DIST = "exponential"
        self.FWR_VAR = 0.85
        self.TIME_SWITCH = 50
        self.GENERATIONS_PER_DAY = 0.5
        self.MEMORY_SAVE = False
        self.KEEP_FULL_TREE = False
        self.QUIET = False
        self.SEQUENCE_LENGTH = 370 # default for uniform seqs, avg heavy len
        self._SELECTION = True

    @property
    def RNG(self):
        """Returns the random number generator instance."""
        return self._x_RNG

    @property
    def SELECTION(self):
        """Returns whether selection is applied in the simulation."""
        return self._SELECTION

    @SELECTION.setter
    def SELECTION(self, value):
        """Sets whether selection is applied in the simulation."""
        if self.UNIFORM and value:
            logger.warning(("Uniform mutation/substitution model specified, "
            "ignoring specified selection"))
            self._SELECTION = False
        else:
            self._SELECTION = value

        if not value:
            self.TARGET_MUTATIONS_HEAVY = 0
            self.TARGET_MUTATIONS_LIGHT = 0

    @property
    def UNIFORM(self):
        """Returns whether a uniform mutation model is used."""
        return self._UNIFORM

    @UNIFORM.setter
    def UNIFORM(self, value):
        """Sets whether a uniform mutation model is used."""
        if value and self.SELECTION:
            self.SELECTION = False
        if value and self.LIGHT_SHM_PER_SITE != 0:
            logger.warning(("Uniform mutation and substitution model specified, "
            "ignoring light chain"))
            self.LIGHT_SHM_PER_SITE = 0
        self._UNIFORM = value

    @property
    def END_TIME(self):
        """Calculates the end time of the simulation based on sample times."""
        def _max_sample_time(sample_times):
            if len(sample_times) > 0:
                return max(sample_times)
            else:
                return 0
        return max([_max_sample_time(x.sample_times) for x in self.LOCATIONS]) + 1

    def update_from_dict(self, dictionary):
        """Updates the settings from a dictionary."""
        for key, value in dictionary.items():
            if key == "LOCATIONS":
                self.LOCATIONS = [LocationSettings(**x) for x in value]
                continue
            setattr(self, key, value)


s = Settings()
