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

from .location import LocationName


class Encodable():
    def encode(self):
        return vars(self)

class LocationSettings(Encodable):
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

    def __init__(self):
        self.LOCATIONS = [
            LocationSettings(
                name=LocationName.GC, 
                mutation_rate=1.0, 
                max_population=1000, 
                migration_rate=0,
                sample_size=50), 
            LocationSettings(
                name=LocationName.OTHER, 
                mutation_rate=0.0, 
                max_population=1000, 
                migration_rate=0,
                sample_size=12)
                ]
        self.HEAVY_MUTATE_PROBABILITY = 0.5
        self.LIGHT_MUTATE_PROBABILITY = 0.25
        self.TARGET_MUTATIONS_HEAVY = 5
        self.TARGET_MUTATIONS_LIGHT = 2
        self.SELECTION = True
        self.RESULTS_DIR = ""
        self.MULTIPLIER = 2
        self._RNG = None
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

    @property
    def END_TIME(self):
        def _max_sample_time(sample_times):
            if len(sample_times) > 0:
                return max(sample_times)
            else:
                return 0
        return max([_max_sample_time(x.sample_times) for x in self.LOCATIONS]) + 1
    
    def update_from_dict(self, dict):
        for key, value in dict.items():
            if key == "LOCATIONS":
                self.LOCATIONS = [LocationSettings(**x) for x in value]
            else:
                setattr(self, key, value)


s = Settings()