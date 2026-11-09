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

 Reads and validates 5-mer mutability and substitution model files.

 The file formats are the same as the bundled S5F files in simble/data:

 - mutability: columns "Fivemer" and "Mutability". One row for each 5-mer.
 - substitution: columns "Fivemer", "A", "C", "G" and "T" (an "N" column is
   allowed and not used). Each row gives the probability that the center
   nucleotide of the 5-mer becomes A, C, G or T.

 A user file must have all 1024 5-mers made of A, C, G and T. Rows with N are
 optional. If they are missing, they are made the same way as in the bundled
 files (and in shazam): a 5-mer with an N at an edge position gets the mean of
 all the 5-mers that replace each N with A, C, G or T. A 5-mer with an N at the
 center has a mutability of 0.

 Only the relative mutability of the 5-mers matters, because simble divides the
 mutability map of a sequence by its sum. The number of mutations comes from
 the SHM rate settings, not from the mutability file.
 """

import itertools
import logging
import os

import numpy as np
import pandas as pd

logger = logging.getLogger(__package__)

NUCLEOTIDES = ["A", "C", "G", "T"]
KMER_LENGTH = 5
CENTER = KMER_LENGTH // 2
ALL_KMERS = ["".join(x) for x in itertools.product("ACGTN", repeat=KMER_LENGTH)]
ACGT_KMERS = [x for x in ALL_KMERS if "N" not in x]
N_KMERS = [x for x in ALL_KMERS if "N" in x]

# A substitution row must sum to 1 within this tolerance. Rows that are within
# the tolerance but further than RENORMALIZE_ABOVE from 1 are renormalized, so
# that the random choice of the new nucleotide does not fail during the
# simulation (numpy accepts a sum within about 1.5e-8 of 1). Rows closer to 1
# than RENORMALIZE_ABOVE are not changed, so a copy of a bundled file (whose rows
# are within about 1e-10 of 1) gives the same results as the bundled file.
ROW_SUM_TOLERANCE = 1e-6
RENORMALIZE_ABOVE = 1e-9

MAX_EXAMPLES = 5


class ModelFileError(ValueError):
    """A 5-mer model file is not valid."""


def _examples(items):
    items = list(items)
    shown = ", ".join(str(x) for x in items[:MAX_EXAMPLES])
    if len(items) > MAX_EXAMPLES:
        shown += f", and {len(items) - MAX_EXAMPLES} more"
    return shown


def _read_table(path, columns, kind):
    """Reads a model file and checks its 5-mers.

    Returns:
        pd.DataFrame: The table, indexed by 5-mer, with only the given columns.
    """
    if not os.path.exists(path):
        raise ModelFileError(f"the {kind} file {path} cannot be found")
    separator = "\t" if str(path).endswith(".tsv") else ","
    try:
        table = pd.read_csv(path, sep=separator, header=0)
    except Exception as e:
        raise ModelFileError(f"cannot read the {kind} file {path}: {e}") from e

    missing = [x for x in ["Fivemer"] + columns if x not in table.columns]
    if missing:
        raise ModelFileError(
            f"the {kind} file {path} has no column {_examples(missing)}. It must have "
            f"the columns: Fivemer, {', '.join(columns)}"
            )
    table = table[["Fivemer"] + columns].copy()
    table["Fivemer"] = table["Fivemer"].astype(str).str.strip().str.upper()

    bad = [x for x in table["Fivemer"] if len(x) != KMER_LENGTH or set(x) - set("ACGTN")]
    if bad:
        raise ModelFileError((
            f"the {kind} file {path} has 5-mers that are not 5 characters of A, C, G, "
            f"T and N: {_examples(bad)}"
            ))
    duplicated = table["Fivemer"][table["Fivemer"].duplicated()].unique()
    if len(duplicated) > 0:
        raise ModelFileError(
            f"the {kind} file {path} has duplicate 5-mers: {_examples(duplicated)}"
            )

    for column in columns:
        values = pd.to_numeric(table[column], errors="coerce")
        not_numbers = table["Fivemer"][values.isna() & table[column].notna()]
        if len(not_numbers) > 0:
            raise ModelFileError((
                f"the {kind} file {path} has values in column {column} that are not "
                f"numbers, for 5-mers: {_examples(not_numbers)}"
                ))
        # an empty value (NA) is 0, the same as for the bundled files
        values = values.fillna(0)
        if not np.all(np.isfinite(values)) or (values < 0).any():
            bad = table["Fivemer"][~np.isfinite(values) | (values < 0)]
            raise ModelFileError((
                f"the {kind} file {path} has values in column {column} that are "
                f"negative or not finite, for 5-mers: {_examples(bad)}"
                ))
        table[column] = values.astype(float)

    table = table.set_index("Fivemer")
    missing_kmers = [x for x in ACGT_KMERS if x not in table.index]
    return table, missing_kmers


def _expansions(kmer):
    """Returns the 5-mers that replace each N in kmer with A, C, G or T."""
    options = [NUCLEOTIDES if x == "N" else [x] for x in kmer]
    return ["".join(x) for x in itertools.product(*options)]


def read_mutability(path):
    """Reads a mutability file.

    Returns:
        dict: 5-mer -> mutability, for all 3125 5-mers of A, C, G, T and N.
    Raises:
        ModelFileError: If the file is not valid.
    """
    table, missing = _read_table(path, ["Mutability"], "mutability")
    if missing:
        raise ModelFileError((
            f"the mutability file {path} must have all 1024 5-mers of A, C, G and T. "
            f"{len(missing)} are missing: {_examples(missing)}"
            ))
    mutability = table["Mutability"].to_dict()
    if not any(mutability[x] > 0 for x in ACGT_KMERS):
        raise ModelFileError(
            f"the mutability file {path} has no 5-mer with a mutability greater than 0"
            )
    for kmer in N_KMERS:
        if kmer in mutability:
            continue
        if kmer[CENTER] == "N":
            mutability[kmer] = 0.0
        else:
            mutability[kmer] = float(np.mean([mutability[x] for x in _expansions(kmer)]))
    return mutability


def read_substitution(path):
    """Reads a substitution file.

    Rows for 5-mers of A, C, G and T can be missing. A missing row is all 0, so
    it is only valid for a 5-mer with a mutability of 0 (see validate_model).

    Returns:
        dict: 5-mer -> [A, C, G, T] probabilities, for all 3125 5-mers.
    Raises:
        ModelFileError: If the file is not valid.
    """
    table, _ = _read_table(path, NUCLEOTIDES, "substitution")
    substitution = {
        kmer: [float(x) for x in row]
        for kmer, row in zip(table.index, table[NUCLEOTIDES].to_numpy())
        }
    for kmer in ACGT_KMERS:
        substitution.setdefault(kmer, [0.0, 0.0, 0.0, 0.0])
    for kmer in N_KMERS:
        if kmer in substitution:
            continue
        if kmer[CENTER] == "N":
            substitution[kmer] = [0.0, 0.0, 0.0, 0.0]
            continue
        # the mean of the valid rows (sum to 1, and the center nucleotide does
        # not stay the same), so that a missing or unused row does not make
        # the mean invalid
        center = NUCLEOTIDES.index(kmer[CENTER])
        rows = [substitution[x] for x in _expansions(kmer)]
        valid = [
            x for x in rows
            if abs(sum(x) - 1) <= ROW_SUM_TOLERANCE and x[center] == 0
            ]
        if valid:
            substitution[kmer] = [float(x) for x in np.mean(valid, axis=0)]
        else:
            substitution[kmer] = [0.0, 0.0, 0.0, 0.0]
    return substitution


def validate_model(mutability, substitution, source):
    """Checks that a mutability and a substitution model can be used together.

    For every 5-mer that can mutate (mutability greater than 0), the
    substitution row must sum to 1, and the probability that the center
    nucleotide stays the same must be 0. Rows within ROW_SUM_TOLERANCE of 1,
    but further than RENORMALIZE_ABOVE from 1, are renormalized. Rows for 5-mers that
    cannot mutate are not used, so they are not checked.

    Args:
        mutability (dict): 5-mer -> mutability.
        substitution (dict): 5-mer -> [A, C, G, T]. Changed in place if rows
            are renormalized.
        source (str): How to refer to the model in an error message.
    Raises:
        ModelFileError: If the models cannot be used together.
    """
    bad_sum = []
    bad_self = []
    for kmer, value in mutability.items():
        if value <= 0:
            continue
        row = substitution[kmer]
        center = kmer[CENTER]
        if center in NUCLEOTIDES and row[NUCLEOTIDES.index(center)] != 0:
            bad_self.append(kmer)
        total = sum(row)
        if abs(total - 1) > ROW_SUM_TOLERANCE:
            bad_sum.append(f"{kmer} ({total:.6g})")
        elif abs(total - 1) > RENORMALIZE_ABOVE:
            substitution[kmer] = [x / total for x in row]
    if bad_self:
        raise ModelFileError((
            f"{source}: the substitution probability for the center nucleotide to "
            f"stay the same must be 0, but it is not 0 for 5-mers: {_examples(bad_self)}"
            ))
    if bad_sum:
        raise ModelFileError((
            f"{source}: the substitution row of each 5-mer that can mutate (mutability "
            "greater than 0) must sum to 1. These rows do not (or are missing): "
            f"{_examples(bad_sum)}"
            ))


def mutability_to_table(mutability):
    """Returns a mutability dict as a table in the bundled file format."""
    values = [mutability[x] for x in ALL_KMERS]
    return pd.DataFrame({"Fivemer": ALL_KMERS, "Mutability": values})


def substitution_to_table(substitution):
    """Returns a substitution dict as a table in the bundled file format."""
    rows = np.array([substitution[x] for x in ALL_KMERS])
    table = pd.DataFrame(rows, columns=NUCLEOTIDES)
    table.insert(0, "Fivemer", ALL_KMERS)
    table["N"] = 0.0
    return table


def bundled_as_dicts(mutability_table, substitution_table):
    """Returns bundled tables (read with helper.read_sf5_table) as dicts."""
    mutability = dict(zip(mutability_table["Fivemer"], mutability_table["Mutability"]))
    substitution = {
        kmer: [float(x) for x in row]
        for kmer, row in zip(
            substitution_table["Fivemer"], substitution_table[NUCLEOTIDES].to_numpy()
            )
        }
    return mutability, substitution
