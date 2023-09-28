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

import abc

import numpy as np

from .helper import (get_mutability_of_kmer, get_substitution_probability,
                     translate_to_amino_acid)
from .settings import s


class Chain:
    __metaclass__ = abc.ABCMeta
    def __init__(
            self,
            nucleotide_seq,
            amino_acid_seq=None,
            nucleotide_gaps=None,
            mutability_map=None,
            gapped_seq=None,
            cdr3_aa_length=13,
            junction=None
    ):
        if nucleotide_seq is None:
            nucleotide_seq, gapped_seq = self.create_naive_seq()
        self.nucleotide_seq = nucleotide_seq
        self.nucleotide_gaps = self._get_gaps(gapped_seq, nucleotide_gaps)
        if amino_acid_seq is None:
            amino_acid_seq = translate_to_amino_acid(self.get_gapped_sequence())
        self.amino_acid_seq = amino_acid_seq
        if mutability_map is None:
            mutability_map = self.create_mutability_map()
        self.mutability_map = mutability_map
        self.is_functional = self.get_functionality()
        self.affinity = 1
        self.CDR3_length = cdr3_aa_length
        self.airr_constants = None
        if junction is not None:
            self.junction_start = self.nucleotide_seq.find(junction)
            self.junction_length = len(junction)

    @property
    def junction(self):
        return self.nucleotide_seq[self.junction_start:self.junction_start+self.junction_length]
    
    @property
    def junction_aa(self):
        return translate_to_amino_acid(self.junction)
    

    def copy(self):
        new = type(self)(
            self.nucleotide_seq, 
            amino_acid_seq=self.amino_acid_seq, 
            nucleotide_gaps=self.nucleotide_gaps, 
            mutability_map=self.mutability_map.copy(),
            cdr3_aa_length=self.CDR3_length
            )
        new.airr_constants = self.airr_constants
        new.junction_start = self.junction_start
        new.junction_length = self.junction_length
        return new


    def _get_gaps(self, gapped_seq, nucleotide_gaps):
        if nucleotide_gaps is None and gapped_seq is None:
            return None
        elif nucleotide_gaps is not None:
            return nucleotide_gaps
        else:
            gaps = {}
            curr_gap = None
            for i in range(len(gapped_seq)):
                if gapped_seq[i] == ".":
                    if curr_gap is None:
                        curr_gap = i
                        gaps[curr_gap] = 1
                    else:
                        gaps[curr_gap] += 1
                else:
                    curr_gap = None
            return gaps
                

    def get_functionality(self):
        if "_" in self.amino_acid_seq:
            return False
        return True


    def get_gapped_sequence(self):
        if self.nucleotide_gaps is None:
            return self.nucleotide_seq
        else:
            gapped_seq = self.nucleotide_seq
            for gap_position, gap_length in self.nucleotide_gaps.items():
                gapped_seq = gapped_seq[:gap_position] + "."*gap_length + gapped_seq[gap_position:]
            return gapped_seq


    def create_mutability_map(self):
        mutability_map = []
        padded_seq = "NN" + self.nucleotide_seq + "NN"
        for i in range(len(self.nucleotide_seq)):
            current_weight = get_mutability_of_kmer(
                padded_seq[i:i+5], 
                self.IS_HEAVY
                )
            mutability_map.append(current_weight)
        return mutability_map


    def update_mutability_map(self, mutated_positions):
        padded_seq = "NN" + self.nucleotide_seq + "NN"
        for i in mutated_positions:
            floor = max(0, i-2)
            ceiling = min(len(self.nucleotide_seq), i+3)
            for j in range(floor, ceiling):
                self.mutability_map[j] = get_mutability_of_kmer(
                    padded_seq[j:j+5], 
                    self.IS_HEAVY
                    )

    def mutate(self, cell_mutation_rate, n=None):

        if n is None:
            n = s._RNG.poisson(self.MUTATE_PROBABILITY*cell_mutation_rate)

        if n == 0:
            return n
        
        padded_seq = "NN" + self.nucleotide_seq + "NN"
        for _ in range(n):
            probability_map = np.array(self.mutability_map) / np.sum(self.mutability_map)
            mutated_position = s._RNG.choice(range(len(self.nucleotide_seq)),
                                                size=1,
                                                p=probability_map,
                                                replace=False)[0]
            probabilities = get_substitution_probability(
                padded_seq[mutated_position:mutated_position+5], 
                self.IS_HEAVY
                )
            substitution = s._RNG.choice(["A", "C", "G", "T"], p=probabilities)
            # substitute the nucleotide at that position
            self.nucleotide_seq = self.nucleotide_seq[:mutated_position] + substitution + self.nucleotide_seq[mutated_position+1:]
            self.update_mutability_map([mutated_position])
        

        self.amino_acid_seq = translate_to_amino_acid(self.get_gapped_sequence())
        self.is_functional = self.get_functionality()

        return n
    
    def calculate_affinity(self, target_pair):
        target = self.get_target_from_pair(target_pair)

        similarities = 0
        cdr_similarities = 0
        affinity = 1
        for i in range(len(self.amino_acid_seq)):
            if self.amino_acid_seq[i] == target.amino_acid_seq[i]:
                similarities += 1
                if i in target.CDR_POSITIONS:
                    cdr_similarities += 1
                affinity *= target.all_multipliers[i]
            else:
                affinity *= 1/target.all_multipliers[i]

        self.cdr_similarity = cdr_similarities/len(target.CDR_POSITIONS)
        self.fwr_similarity = (similarities - cdr_similarities)/(len(self.amino_acid_seq) - len(target.CDR_POSITIONS))
        self.similarity = similarities/len(self.amino_acid_seq)
        self.affinity = affinity

        return affinity
    
    def get_observed_mutations(self, germline_gapped, targets):
        cdr = list(range(3*27, 3*39)) + list(range(3*56, 3*66)) + list(range(3*105, 105+3*(self.CDR3_length)))
        observed_mutations = 0
        cdr_mutations = 0
        target_site_mutations = 0
        gapped = self.get_gapped_sequence()
        for i in range(312):
            if gapped[i] != germline_gapped[i]:
                observed_mutations += 1
                if i in targets:
                    target_site_mutations += 1
                if i in cdr:
                    cdr_mutations += 1
        filtered = observed_mutations - target_site_mutations
        fwr_mutations = observed_mutations - cdr_mutations
        return observed_mutations/312, filtered/(312-len(targets)), cdr_mutations/len(cdr), fwr_mutations/(312-len(cdr))

    def as_AIRR(self, generation):
        generated = {
            "sequence_id": "heavy" if self.IS_HEAVY else "light",
            "sequence": self.nucleotide_seq,
            "sequence_alignment": self.get_gapped_sequence(),
            "sample_time": generation,
            "locus": "IGH" if self.IS_HEAVY else "IGL",
            "junction": self.junction,
            "junction_aa": self.junction_aa,
            "junction_length": self.junction_length
        }
        generated.update(self.airr_constants)
        return generated

        

    @abc.abstractclassmethod
    def get_target_from_pair(self, target_pair):
        pass
  


class HeavyChain(Chain):
    MUTATE_PROBABILITY = s.HEAVY_MUTATE_PROBABILITY
    IS_HEAVY = True
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        

    def get_target_from_pair(self, target_pair):
        return target_pair.heavy
    


class LightChain(Chain):
    MUTATE_PROBABILITY = s.LIGHT_MUTATE_PROBABILITY
    IS_HEAVY = False
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
    

    def get_target_from_pair(self, target_pair):
        return target_pair.light
