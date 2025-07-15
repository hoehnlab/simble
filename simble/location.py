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
# pylint: disable=expression-not-assigned

from enum import Enum


class LocationName(Enum):
    """Enum representing different locations in the simulation."""
    GC = "germinal_center"
    OTHER = "other"

    def encode(self):
        """Encodes the enum as a dictionary for serialization."""
        return {"__enum__": str(self)}

class Location():
    """Represents a location in the simulation.
    Attributes:
        name (str): The name of the location.
        settings (LocationSettings): The settings for the location.
        current_generation (list): The current population in the location.
        immigrating_population (list): The population that is immigrating to the location.
        number_of_children (list): The number of children produced by the population.
    """
    def __init__(
            self,
            name,
            settings
    ):
        """Initializes a Location instance.
        Args:
            name (str): The name of the location.
            settings (LocationSettings): The settings for the location.
        """
        self.name = name
        self.settings = settings
        self.current_generation = []
        self.immigrating_population = []
        self.number_of_children = []

    def update_cell(self, node):
        """Updates the cell's location and mutation rate."""

        node.cell.location = self.name
        node.cell.mutation_rate = self.settings.mutation_rate

    def finish_migration(self):
        """Finalizes the migration of cells to this location."""

        [self.update_cell(x) for x in self.immigrating_population]
        self.current_generation.extend(self.immigrating_population)
        self.immigrating_population = []

PUBLIC_ENUMS = {
    'LocationName': LocationName,
}

def as_enum(d):
    """Converts a dictionary (e.g. from json) to an enum if it contains an encoded enum."""
    if "__enum__" in d:
        name, member = d["__enum__"].split(".")
        return getattr(PUBLIC_ENUMS[name], member)
    else:
        return d
