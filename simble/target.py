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

# IMGT conserved sites are 23, 41, 89, 104, 104+cdr3+1, but cdr3 is variable
# python is 0-indexed so we need to subtract 1
CONSERVED_SITES = [23-1, 41-1, 89-1, 104-1]


class TargetAminoPair:
    """Represents a pair of target amino acids for heavy and light chains.
    Attributes:
        heavy (TargetAminoAcid): The target amino acid for the heavy chain.
        light (TargetAminoAcid): The target amino acid for the light chain.
    """
    def __init__(
            self,
            heavy_gapped_nucleotide,
            light_gapped_nucleotide,
            heavy_cdr3_length,
            light_cdr3_length
            ):
        """Initializes a TargetAminoPair instance.
        Args:
            heavy_gapped_nucleotide (str): The gapped nucleotide sequence for the heavy chain.
            light_gapped_nucleotide (str): The gapped nucleotide sequence for the light chain.
            heavy_cdr3_length (int): The length of the CDR3 region for the heavy chain.
            light_cdr3_length (int): The length of the CDR3 region for the light chain.
        """
        self.heavy = TargetAminoAcid(heavy_gapped_nucleotide, cdr3_length=heavy_cdr3_length)
        self.light = TargetAminoAcid(light_gapped_nucleotide, cdr3_length=light_cdr3_length)

    @property
    def max_affinity(self):
        """Calculates the maximum affinity of the target pair."""
        return self.heavy.max_affinity * self.light.max_affinity

    def mutate(self, heavy_n, light_n):
        """Creates target mutations in the target amino acid chains.
        Args:
            heavy_n (int): The number of mutations to apply to the heavy chain.
            light_n (int): The number of mutations to apply to the light chain."""
        self.heavy.mutate(heavy_n)
        self.light.mutate(light_n)



class TargetAminoAcid:
    """Represents a target amino acid sequence.
    Attributes:
        gapped_nucleotide_seq (str): The gapped nucleotide sequence of the target.
        CDR_POSITIONS (list): The positions of the CDR regions in the amino acid 
            sequence.
        amino_acid_seq (str): The amino acid sequence derived from the gapped 
            nucleotide sequence.
        mutation_locations (list): The positions of mutations from germline in 
            the target amino acid sequence.
        all_multipliers (dict): A dictionary of multipliers for each position in 
            the amino acid sequence.
        cdr_multipliers (dict): A dictionary of multipliers specifically for CDR positions.
        fwr_multipliers (dict): A dictionary of multipliers specifically for FWR positions.
    """
    def __init__(
            self,
            gapped_nucleotide_seq,
            cdr3_length
            ):
        self.gapped_nucleotide_seq = gapped_nucleotide_seq
        # IMGT numbering: CDR1 = 27-38 CDR2 = 56-65 CDR3 = 105 to 105+CDR3_length
        # but python is 0-indexed so we need to subtract 1
        self.CDR_POSITIONS = ( # pylint: disable=invalid-name
            list(range(27-1, 39-1))
            + list(range(56-1, 66-1))
            + list(range(105-1, 105-1+cdr3_length))
        )
        self.amino_acid_seq = translate_to_amino_acid(self.gapped_nucleotide_seq)
        FWR_POSITIONS = [ # pylint: disable=invalid-name
            x
            for x in range(len(self.amino_acid_seq))
            if x not in self.CDR_POSITIONS
            ]
        self.mutation_locations = []
        self.all_multipliers = {}
        if s.CDR_DIST == "exponential":
            exp_mean = -(s.MULTIPLIER-1)/np.log(1-s.CDR_VAR)
            exp_distribution = [
                1 + x
                for x in s.RNG.exponential(
                    exp_mean,
                    len(self.CDR_POSITIONS)
                    )
                ]
        elif s.CDR_DIST == "constant":
            exp_distribution = [s.CDR_VAR for _ in range(len(self.CDR_POSITIONS))]
        else:
            exp_distribution = [s.MULTIPLIER for _ in range(len(self.CDR_POSITIONS))]

        if s.FWR_DIST == "exponential":
            exp_mean = -(s.MULTIPLIER-1)/np.log(1-s.FWR_VAR)
            fwr_distribution = [1 + x for x in s.RNG.exponential(exp_mean, len(FWR_POSITIONS))]
        elif s.FWR_DIST == "constant":
            fwr_distribution = [s.FWR_VAR for _ in range(len(FWR_POSITIONS))]
        elif s.FWR_DIST == "constant-noise":
            fwr_distribution = [
                s.FWR_VAR + s.RNG.normal(0, 0.1)
                for _ in range(len(FWR_POSITIONS))
                ]
        else:
            fwr_distribution = [s.MULTIPLIER for _ in range(len(FWR_POSITIONS))]

        self.cdr_multipliers = {
            self.CDR_POSITIONS[i]: exp_distribution[i]
            for i in range(len(self.CDR_POSITIONS))
            }
        self.fwr_multipliers = {
            FWR_POSITIONS[i]: fwr_distribution[i]
            for i in range(len(FWR_POSITIONS))
            }
        self.conserved_sites = CONSERVED_SITES + [CONSERVED_SITES[-1] + cdr3_length + 1]
        conserved_multipliers = {
            x: s.MULTIPLIER * 1.25
            for x in self.conserved_sites
        }
        self.all_multipliers.update(self.cdr_multipliers)
        self.all_multipliers.update(self.fwr_multipliers)
        self.all_multipliers.update(conserved_multipliers)

    @property
    def max_affinity(self):
        """Calculates the maximum affinity of the target amino acid sequence."""
        return np.prod([x for _, x in self.all_multipliers.items()])


    def choose_replacement_nucleotide(self, codon, curr_amino_acid):
        """Chooses a replacement nucleotide for a codon that results in a different amino acid.
        Args:
            codon (str): The codon to mutate.
            curr_amino_acid (str): The current amino acid represented by the codon.
        Returns:
            tuple: A tuple containing the new codon and the new amino acid.
        """
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

        new_codon = s.RNG.choice(possible_codon)
        new_amino_acid = translate_to_amino_acid(new_codon)
        return new_codon, new_amino_acid


    def mutate(self, n):
        """Mutates the target amino acid sequence by replacing nucleotides.
        Args:
            n (int): The number of mutations to apply.
        """
        if self.amino_acid_seq == "" or n == 0:
            self.mutation_locations = []
            return
        CDR_PROB = 1 # pylint: disable=invalid-name
        OTHER_PROB = 0 # pylint: disable=invalid-name
        mutate_probability = []
        amino_acid_seq = self.amino_acid_seq
        nucleotide_seq = self.gapped_nucleotide_seq
        for i, amino_acid in enumerate(amino_acid_seq):
            if amino_acid in ["X", "_"]:
                mutate_probability.append(0)
            elif i in self.conserved_sites:
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
        mutate_positions = s.RNG.choice(
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
