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

from .chain import HeavyChain, LightChain, EmptyChain
from .helper import get_start_pair
from .location import LocationName
from .settings import s


class CellType(Enum):
    """Enum representing different cell types in the simulation."""
    DEFAULT = "gc_b_cell"
    PC = "plasma_cell"
    MBC = "memory_b_cell"

class Cell:
    """Represents a cell in the simulation.
    
    Attributes:
        heavy_chain (HeavyChain): The heavy chain of the cell.
        light_chain (LightChain): The light chain of the cell.
        created_at (int): The generation at which the cell was created.
        is_alive (bool): Whether the cell is alive.
        location (LocationName): The location of the cell.
        cell_type (CellType): The type of the cell.
    """
    def __init__(
            self,
            heavy_chain,
            light_chain,
            created_at,
            is_alive=True,
            location=LocationName.GC,
            cell_type=CellType.DEFAULT,
            clone_id=None,
            naive_row_idx=None) -> None:
        self.user_constants = {}
        if heavy_chain is None:
            pair = get_start_pair(clone_id, naive_row_idx=naive_row_idx)
            heavy_chain = HeavyChain(**pair.heavy.chain._asdict())
            heavy_chain.airr_constants = pair.heavy.constants
            if pair.light is None:
                # a single chain or uniform simulation: no light chain
                light_chain = None
            else:
                light_chain = LightChain(**pair.light.chain._asdict())
                light_chain.airr_constants = pair.light.constants
            self.user_constants = pair.user_constants

        if light_chain is None:
            light_chain = EmptyChain()
        if isinstance(light_chain, EmptyChain):
            self.__class__ = SingleChainCell

        self.heavy_chain = heavy_chain
        self.light_chain = light_chain
        self.created_at = created_at
        self.is_alive = is_alive
        self.location = location
        self.cell_type = cell_type

        self.mutation_rate = s.LOCATIONS[0].mutation_rate if self.location == LocationName.GC else 0
        self.affinity = 1
        self.clone_id = clone_id

    def kill_cell(self):
        """Marks the cell as dead."""
        self.is_alive = False


    def mutate_cell(self):
        """Mutates the cell's heavy and light chains."""
        heavy_n = self.heavy_chain.mutate(self.mutation_rate)
        light_n = self.light_chain.mutate(self.mutation_rate)
        return heavy_n, light_n

    def _add_cell_AIRR(self, row): # pylint: disable=invalid-name
        """Adds cell-specific information to the AIRR row."""
        row["sequence_id"]=f'{str(id(self))}_{row["sequence_id"]}'
        row["cell_id"]=str(id(self))
        row["location"]=self.location.value
        row["celltype"]=self.cell_type.value
        row.update(self.user_constants)
        return row

    def as_AIRR(self, generation): # pylint: disable=invalid-name
        """Returns the cell data in AIRR format."""
        heavy = self._add_cell_AIRR(self.heavy_chain.as_AIRR(generation))
        light = self._add_cell_AIRR(self.light_chain.as_AIRR(generation))
        return [heavy, light]

    def as_fasta_helper(self, generation, heavy, name=None, extra_fields=None):
        """Returns the cell and chain data in FASTA format for a single chain.

        Args:
            generation (int): The sample time. If None (for example, for an
                internal node), the header does not include ``sample_time``.
            heavy (bool): Whether to write the heavy chain (else the light chain).
            name (str, optional): The name of the record. Use the node's name
                (``Node.name``), so that the FASTA IDs are the same as the AIRR
                ``sequence_id`` and the Newick names. If None, the cell's id
                is used (without a clone prefix).
            extra_fields (dict, optional): More fields to add to the header.
        """
        chain = self.heavy_chain if heavy else self.light_chain
        row = {"locus": chain.locus}
        if generation is not None:
            row["sample_time"] = generation
        row["location"] = self.location.value
        row["celltype"] = self.cell_type.value
        if extra_fields is not None:
            row.update(extra_fields)
        row = "|".join([f"{x}={y}" for x, y in row.items()])
        name = str(id(self)) if name is None else name
        return f">{name}_{'heavy' if heavy else 'light'}|{row}\n{chain.nucleotide_seq}\n"

    def as_fasta(self, generation, name=None, extra_fields=None):
        """Returns both chains with cell data in FASTA format.

        See ``as_fasta_helper`` for the arguments.
        """
        heavy = self.as_fasta_helper(generation, True, name, extra_fields)
        light = self.as_fasta_helper(generation, False, name, extra_fields)
        return heavy + light

    def remake_self(self):
        """Creates a new Cell instance with the same properties."""
        new = type(self)(
            self.heavy_chain,
            self.light_chain,
            self.created_at,
            self.is_alive,
            self.location,
            self.cell_type,
            self.clone_id
        )
        new.mutation_rate = self.mutation_rate
        new.affinity = self.affinity
        new.user_constants = self.user_constants
        return new

    def calculate_affinity(self, target_pair):
        """Calculates the affinity of the cell's chains to a target pair."""
        self.affinity = (
            self.heavy_chain.calculate_affinity(target_pair)
            * self.light_chain.calculate_affinity(target_pair)
        )

        if not self.heavy_chain.is_functional or not self.light_chain.is_functional:
            self.affinity = 0

    def differentiate(self, cell_type):
        """Changes the cell type to a different type."""
        self.cell_type = cell_type



class SingleChainCell(Cell):
    """Represents a cell with only one chain (heavy)."""

    def as_AIRR(self, generation): # pylint: disable=invalid-name
        """Returns the cell data in AIRR format."""
        heavy = self._add_cell_AIRR(self.heavy_chain.as_AIRR(generation))
        return [heavy]

    def as_fasta(self, generation, name=None, extra_fields=None):
        """Returns the heavy chain with cell data in FASTA format.

        See ``Cell.as_fasta_helper`` for the arguments.
        """
        return self.as_fasta_helper(generation, True, name, extra_fields)

    def calculate_affinity(self, target_pair):
        """Calculates the affinity of the cell's chain to a target pair."""
        self.affinity = self.heavy_chain.calculate_affinity(target_pair)

        if not self.heavy_chain.is_functional:
            self.affinity = 0

    def mutate_cell(self):
        """Mutates the cell's heavy chain."""
        heavy_n = self.heavy_chain.mutate(self.mutation_rate)
        return heavy_n, 0
