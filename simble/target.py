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

import numpy as np

from .helper import translate_to_amino_acid
from .settings import s


class TargetAminoPair:
    def __init__(self, heavy_gapped_nucleotide, light_gapped_nucleotide, heavy_cdr3_length, light_cdr3_length):
        self.heavy = TargetAminoAcid(heavy_gapped_nucleotide, CDR3_length=heavy_cdr3_length)
        self.light = TargetAminoAcid(light_gapped_nucleotide, CDR3_length=light_cdr3_length)

    @property
    def max_affinity(self):
        return self.heavy.max_affinity * self.light.max_affinity

    def mutate(self, heavy_n, light_n):
        self.heavy.mutate(heavy_n)
        self.light.mutate(light_n)

        

class TargetAminoAcid:
    def __init__(self, gapped_nucleotide_seq, CDR3_length):
        self.gapped_nucleotide_seq = gapped_nucleotide_seq
        # CDR1 = 27-38 CDR2 = 56-65 CDR3 = 105 to 105+CDR3_length
        self.CDR_POSITIONS = list(range(27, 39)) + list(range(56, 66)) + list(range(105, 105+CDR3_length))
        self.amino_acid_seq = translate_to_amino_acid(self.gapped_nucleotide_seq)
        FWR_POSITIONS = [x for x in range(len(self.amino_acid_seq)) if x not in self.CDR_POSITIONS]
        self.mutation_locations = []
        self.all_multipliers = {}
        if s.CDR_DIST == "exponential":
            exp_mean = -(s.MULTIPLIER-1)/np.log(1-s.CDR_VAR)
            exp_distribution = [1 + x for x in s._RNG.exponential(exp_mean, len(self.CDR_POSITIONS))]
        elif s.CDR_DIST == "constant":
            exp_distribution = [s.CDR_VAR for _ in range(len(self.CDR_POSITIONS))]
        else:
            exp_distribution = [s.MULTIPLIER for _ in range(len(self.CDR_POSITIONS))]
        
        if s.FWR_DIST == "exponential":
            exp_mean = -(s.MULTIPLIER-1)/np.log(1-s.FWR_VAR)
            fwr_distribution = [1 + x for x in s._RNG.exponential(exp_mean, len(FWR_POSITIONS))]
        elif s.FWR_DIST == "constant":
            fwr_distribution = [s.FWR_VAR for _ in range(len(FWR_POSITIONS))]
        elif s.FWR_DIST == "constant-noise":
            fwr_distribution = [s.FWR_VAR + s._RNG.normal(0, 0.1) for _ in range(len(FWR_POSITIONS))]
        else:
            fwr_distribution = [s.MULTIPLIER for _ in range(len(FWR_POSITIONS))]
        
        self.cdr_multipliers = {self.CDR_POSITIONS[i]: exp_distribution[i] for i in range(len(self.CDR_POSITIONS))}
        self.fwr_multipliers = {FWR_POSITIONS[i]: fwr_distribution[i] for i in range(len(FWR_POSITIONS))}
        self.all_multipliers.update(self.cdr_multipliers)
        self.all_multipliers.update(self.fwr_multipliers)

    @property
    def max_affinity(self):
        return np.prod([x for _, x in self.all_multipliers.items()])
        

    def choose_replacement_nucleotide(self, codon, curr_amino_acid):
        # get all one nucleotide mutations
        possible_codon = []
        for i in range(3):
            nucleotide = codon[i]
            for replacement in [x for x in ["A", "C", "G", "T"] if x != nucleotide]:
                new_codon = codon[:i] + replacement + codon[i+1:]
                new_amino_acid = translate_to_amino_acid(new_codon)
                if new_amino_acid == curr_amino_acid or new_amino_acid == "_":
                    continue
                possible_codon.append(codon[:i] + replacement + codon[i+1:])
        
        new_codon = s._RNG.choice(possible_codon)
        new_amino_acid = translate_to_amino_acid(new_codon)
        return new_codon, new_amino_acid
    
    def mutate(self, n):
        CDR_PROB = 1
        OTHER_PROB = 0
        mutate_probability = []
        amino_acid_seq = self.amino_acid_seq
        nucleotide_seq = self.gapped_nucleotide_seq
        for i in range(len(amino_acid_seq)):
            if amino_acid_seq[i] in ["X", "_"]:
                mutate_probability.append(0)
            elif i in self.CDR_POSITIONS:
                mutate_probability.append(CDR_PROB)
            else:
                mutate_probability.append(OTHER_PROB)
        mutate_probability = np.array(mutate_probability) / np.sum(mutate_probability)
        is_nan = np.isnan(mutate_probability)
        if True in is_nan:
            print("NaN in mutate probability!")
            mutate_probability=None
        mutate_positions = s._RNG.choice(
            range(len(amino_acid_seq)), 
            size=n, 
            p=mutate_probability, 
            replace=False)
        for i in mutate_positions:
            new_codon, new_amino_acid = self.choose_replacement_nucleotide(
                nucleotide_seq[i*3:i*3+3], 
                amino_acid_seq[i]
                )
            nucleotide_seq = nucleotide_seq[:i*3] + new_codon + nucleotide_seq[i*3+3:]
            amino_acid_seq = amino_acid_seq[:i] + new_amino_acid + amino_acid_seq[i+1:]

        self.gapped_nucleotide_seq = nucleotide_seq
        self.amino_acid_seq = amino_acid_seq
        self.mutation_locations = mutate_positions
        multipliers = {x: s.MULTIPLIER for x in mutate_positions}
        self.all_multipliers.update(multipliers)