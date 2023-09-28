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

from enum import Enum


class LocationName(Enum):
    GC = "germinal_center"
    OTHER = "other"

    def encode(self):
        return {"__enum__": str(self)}

class Location():
    def __init__(
            self,
            name,
            settings
    ):
        self.name = name
        self.settings = settings
        # self.mutation_rate = mutation_rate
        # self.MAX_POPULATION = max_population
        # self.migration_rate = migration_rate
        self.current_generation = []
        self.immigrating_population = []
        self.number_of_reproducing_cells = 0
        self.number_of_cells_with_more_than_one_child = 0

    def finish_migration(self):
        def update_cell(node):
            node.cell.location = self.name
            node.cell.mutation_rate = self.settings.mutation_rate

        [update_cell(x) for x in self.immigrating_population]
        self.current_generation.extend(self.immigrating_population)
        self.immigrating_population = []

PUBLIC_ENUMS = {
    'LocationName': LocationName,
}

def as_enum(d):
    if "__enum__" in d:
        name, member = d["__enum__"].split(".")
        return getattr(PUBLIC_ENUMS[name], member)
    else:
        return d