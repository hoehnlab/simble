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

from .chain import HeavyChain, LightChain
from .helper import get_random_start_pair
from .location import LocationName
from .settings import s


class CellType(Enum):
    DEFAULT = "default"


class Cell:
    def __init__(
            self, 
            heavy_chain, 
            light_chain, 
            created_at, 
            is_alive=True,
            location=LocationName.GC,
            type=CellType.DEFAULT) -> None:
        if heavy_chain is None and light_chain is None:
            pair = get_random_start_pair()
            heavy_chain = HeavyChain(
                nucleotide_seq=pair.heavy.input.chain, 
                gapped_seq=pair.heavy.input.aligned,
                cdr3_aa_length=pair.heavy.input.cdr3_aa_length,
                junction=pair.heavy.input.junction)
            heavy_chain.airr_constants = pair.heavy.constants
            light_chain = LightChain(
                nucleotide_seq=pair.light.input.chain, 
                gapped_seq=pair.light.input.aligned,
                cdr3_aa_length=pair.light.input.cdr3_aa_length,
                junction=pair.light.input.junction)
            light_chain.airr_constants = pair.light.constants
        self.heavy_chain = heavy_chain
        self.light_chain = light_chain
        self.created_at = created_at
        self.is_alive = is_alive
        self.location = location
        self.type = type
        self.mutation_rate = s.LOCATIONS[0].mutation_rate if self.location == LocationName.GC else 0
        self.affinity = 1

    def kill_cell(self):
        self.is_alive = False


    def mutate_cell(self):
        heavy_n = self.heavy_chain.mutate(self.mutation_rate)
        light_n = self.light_chain.mutate(self.mutation_rate)
        return heavy_n, light_n
    
    def _add_cell_AIRR(self, row):
        row["sequence_id"]=f'{str(id(self))}_{row["sequence_id"]}'
        row["cell_id"]=str(id(self))
        row["location"]=self.location.value
        return row
    
    def as_AIRR(self, generation):
        heavy = self._add_cell_AIRR(self.heavy_chain.as_AIRR(generation))
        light = self._add_cell_AIRR(self.light_chain.as_AIRR(generation))
        return [heavy, light]
    
    def as_fasta_helper(self, generation, heavy):
        row = {
            "locus": "IGH" if heavy else "IGL",
            "sample_time": generation,
            "location": self.location.value,
        }
        row = "|".join([f"{x}={y}" for x, y in row.items()])
        seq = self.heavy_chain.nucleotide_seq if heavy else self.light_chain.nucleotide_seq
        return f">{str(id(self))}_{'heavy' if heavy else 'light'}|{row}\n{seq}\n"
    
    def as_fasta(self, generation):
        heavy = self.as_fasta_helper(generation, heavy=True)
        light = self.as_fasta_helper(generation, heavy=False)
        return heavy + light
    
    def remake_self(self):
        new = Cell(
            self.heavy_chain,
            self.light_chain,
            self.created_at,
            self.is_alive,
            self.location,
            self.type
        )
        new.mutation_rate = self.mutation_rate
        new.affinity = self.affinity
        return new
    
    def calculate_affinity(self, target_pair):
        self.affinity = self.heavy_chain.calculate_affinity(target_pair) * self.light_chain.calculate_affinity(target_pair)

        if not self.heavy_chain.is_functional or not self.light_chain.is_functional:
            self.affinity = 0
