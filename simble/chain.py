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
import logging

import numpy as np

from .helper import (get_mutability_of_kmer, get_substitution_probability,
                     translate_to_amino_acid)
from .settings import s

logger = logging.getLogger(__package__)

class Chain:
    """
    Represents a chain of the BCR.

    Attributes:
        nucleotide_seq (str): The nucleotide sequence of the chain.
        amino_acid_seq (str): The amino acid sequence of the chain.
        nucleotide_gaps (dict): A dictionary mapping gap positions to their lengths.
        mutability_map (list): A list of mutability weights for each nucleotide position.
        CDR3_length (int): The length of the CDR3 region in amino acids.
        junction (str): The junction sequence of the chain.
        cdr_similarity (float): Similarity of the CDR regions to the target.
        fwr_similarity (float): Similarity of the FWR regions to the target.
        similarity (float): Overall similarity of the chain to the target.
    """

    __metaclass__ = abc.ABCMeta

    _mutate_probability = None
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
        """
        Initializes a new Chain object.

        Args:
            nucleotide_seq (str): The nucleotide sequence of the chain.
            amino_acid_seq (str, optional): The amino acid sequence of the chain. If 
                not provided, it will be derived from the nucleotide sequence.
            nucleotide_gaps (dict, optional): A dictionary mapping gap positions 
                to their lengths. If not provided, it will be derived from the
                gapped sequence.
            mutability_map (list, optional): A list of mutability weights for 
                each nucleotide position. If not provided, it will be created 
                based on the nucleotide sequence.
            gapped_seq (str, optional): The gapped nucleotide sequence. 
            cdr3_aa_length (int, optional): The length of the CDR3 region in 
                amino acids. Defaults to 13.
            junction (str, optional): The junction sequence of the chain. If 
                provided, it will be used to set the junction start and 
                length attributes.
        Raises:
            ValueError: If nucleotide_seq is not provided.
        """

        if nucleotide_seq is None:
            raise ValueError("nucleotide_seq must be provided")
        self.nucleotide_seq = nucleotide_seq
        self.nucleotide_gaps = self._get_gaps(gapped_seq, nucleotide_gaps)
        if amino_acid_seq is None:
            amino_acid_seq = translate_to_amino_acid(self.get_gapped_sequence())
        self.amino_acid_seq = amino_acid_seq
        if mutability_map is None:
            mutability_map = self.create_mutability_map()
        self.mutability_map = mutability_map
        self.affinity = 1
        self.cdr3_length = cdr3_aa_length

        self.airr_constants = None

        self.cdr_similarity = None
        self.fwr_similarity = None
        self.similarity = None

        if junction is not None:
            self.junction_start = self.nucleotide_seq.find(junction)
            self.junction_length = len(junction)

    @property
    def junction(self):
        """Returns the junction nucleotide sequence of the chain."""
        return self.nucleotide_seq[self.junction_start:self.junction_start+self.junction_length]

    @property
    def junction_aa(self):
        """Returns the junction amino acid sequence of the chain."""
        return translate_to_amino_acid(self.junction)

    @property
    def is_functional(self):
        """ Checks if the chain is functional based on its amino acid sequence."""
        return self.get_functionality()


    def copy(self):
        """Creates a deep copy of the Chain object"""
        new = type(self)(
            self.nucleotide_seq,
            amino_acid_seq=self.amino_acid_seq,
            nucleotide_gaps=self.nucleotide_gaps,
            mutability_map=self.mutability_map.copy(),
            cdr3_aa_length=self.cdr3_length
            )
        new.airr_constants = self.airr_constants
        new.junction_start = self.junction_start
        new.junction_length = self.junction_length
        new.mutate_probability = self.mutate_probability
        return new


    def _get_gaps(self, gapped_seq, nucleotide_gaps):
        if nucleotide_gaps is None and gapped_seq is None:
            logger.warning("no gap information provided, assuming no gaps")
            return None
        if nucleotide_gaps is not None:
            return nucleotide_gaps
        gaps = {}
        curr_gap = None
        for i, nucleotide in enumerate(gapped_seq):
            if nucleotide == ".":
                if curr_gap is None:
                    curr_gap = i
                    gaps[curr_gap] = 1
                else:
                    gaps[curr_gap] += 1
            else:
                curr_gap = None
        return gaps


    def get_functionality(self):
        """Checks if the chain is functional based on its amino acid sequence"""
        if "_" in self.amino_acid_seq:
            return False
        return True


    def get_gapped_sequence(self):
        """Returns the nucleotide sequence with gaps represented by '.'"""
        if self.nucleotide_gaps is None:
            return self.nucleotide_seq
        else:
            gapped_seq = self.nucleotide_seq
            for gap_position, gap_length in self.nucleotide_gaps.items():
                gapped_seq = gapped_seq[:gap_position] + "."*gap_length + gapped_seq[gap_position:]
            return gapped_seq


    def create_mutability_map(self):
        """Creates a mutability map based on the nucleotide sequence"""
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
        """Updates the mutability map based on mutated positions"""
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
        """
        Mutates the chain based on a mutation rate and returns the number of mutations.

        Args:
            cell_mutation_rate (float): The mutation rate for the cell.
            n (int, optional): The number of mutations to perform. If None, 
                it will be sampled from a Poisson distribution.
        Returns:
            int: The number of mutations performed.
        """

        if n is None:
            n = s._RNG.poisson(self.mutate_probability*cell_mutation_rate) # pylint: disable=protected-access

        if n == 0:
            return n

        padded_seq = "NN" + self.nucleotide_seq + "NN"
        for _ in range(n):
            probability_map = np.array(self.mutability_map) / np.sum(self.mutability_map)
            mutated_position = s._RNG.choice( # pylint: disable=protected-access
                range(len(self.nucleotide_seq)),
                size=1,
                p=probability_map,
                replace=False
                )[0]
            probabilities = get_substitution_probability(
                padded_seq[mutated_position:mutated_position+5],
                self.IS_HEAVY
                )
            substitution = s._RNG.choice( # pylint: disable=protected-access
                ["A", "C", "G", "T"],
                p=probabilities
                )
            # substitute the nucleotide at that position
            self.nucleotide_seq = (
                self.nucleotide_seq[:mutated_position]
                + substitution
                + self.nucleotide_seq[mutated_position+1:]
            )
            self.update_mutability_map([mutated_position])


        self.amino_acid_seq = translate_to_amino_acid(self.get_gapped_sequence())

        return n

    def calculate_affinity(self, target_pair):
        """Calculates the affinity of the chain to a target pair"""
        target = self.get_target_from_pair(target_pair)

        similarities = 0
        cdr_similarities = 0
        affinity = 1
        for i, amino_acid in enumerate(self.amino_acid_seq):
            if amino_acid == target.amino_acid_seq[i]:
                similarities += 1
                if i in target.CDR_POSITIONS:
                    cdr_similarities += 1
                affinity *= target.all_multipliers[i]
            else:
                affinity *= 1/target.all_multipliers[i]

        self.cdr_similarity = cdr_similarities/len(target.CDR_POSITIONS)
        self.fwr_similarity = (
            (similarities - cdr_similarities)
            /(len(self.amino_acid_seq) - len(target.CDR_POSITIONS))
        )
        self.similarity = similarities/len(self.amino_acid_seq)
        self.affinity = affinity

        return affinity

    def get_observed_mutations(self, germline_gapped, targets):
        """ Calculates the observed mutations in the chain compared to a germline sequence.
        Args:
            germline_gapped (str): The gapped germline sequence to compare against.
            targets (list): A list of target positions to exclude from the mutation count.
        Returns:
            tuple: A tuple containing the observed mutations, filtered mutations, 
                CDR mutations, and FWR mutations.
        """
        cdr = (
            list(range(3*27, 3*39))
            + list(range(3*56, 3*66))
            + list(range(3*105, 105+3*(self.cdr3_length)))
        )
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
        return (
            observed_mutations/312,
            filtered/(312-len(targets)),
            cdr_mutations/len(cdr),
            fwr_mutations/(312-len(cdr))
        )

    def as_AIRR(self, generation): # pylint: disable=invalid-name
        """Generates a dictionary representation of the chain in AIRR format"""
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

    @abc.abstractmethod
    def get_target_from_pair(self, target_pair):
        """Returns the appropriate target from a TargetAminoPair based on the chain type"""


    @property
    def mutate_probability(self):
        """Returns the mutation probability for the chain"""
        if self._mutate_probability is None:
            logger.debug("calculating mutate_probability based on settings")
            self._mutate_probability = self.shm_per_site * len(self.nucleotide_seq)
        return self._mutate_probability

    @mutate_probability.setter
    def mutate_probability(self, value):
        """Sets the mutation probability for the chain if it is not already set"""
        if self._mutate_probability is None:
            self._mutate_probability = value

    @property
    @abc.abstractmethod
    def shm_per_site(self):
        """Returns the mutation rate per site per week for the chain"""

    @property
    @abc.abstractmethod
    def IS_HEAVY(self): # pylint: disable=invalid-name
        """Returns whether the chain is a heavy chain"""



class HeavyChain(Chain):
    """
    Represents an IGH (heavy) chain of the BCR.

    Attributes:
        shm_per_site_per_week (float): Mutation rate per site per week for the heavy chain
        IS_HEAVY (bool): Indicates that this is a heavy chain
    """

    @property
    def shm_per_site(self):
        return s.HEAVY_SHM_PER_SITE

    @property
    def IS_HEAVY(self):
        return True

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_target_from_pair(self, target_pair):
        return target_pair.heavy



class LightChain(Chain):
    """
    Represents an IGL or IGK (light) chain of the BCR.

    Attributes:
        shm_per_site_per_week (float): Mutation rate per site per week for the light chain
        IS_HEAVY (bool): Indicates that this is NOT a heavy chain
    """

    @property
    def shm_per_site(self):
        return s.LIGHT_SHM_PER_SITE

    @property
    def IS_HEAVY(self):
        return False

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def get_target_from_pair(self, target_pair):
        return target_pair.light
