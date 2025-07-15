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

import logging
import math
import os
from collections import namedtuple

import matplotlib.pyplot as plt
import pandas as pd

from .settings import s

logger = logging.getLogger(__package__)

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    """Returns the absolute path to a data file in the simble package.
    Args:
        path (str): The relative path to the data file.
    Returns:
        str: The absolute path to the data file.
    """
    return os.path.join(_ROOT, 'data', path)

AIRR_REQUIRED_FIELDS = [
    'sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call',
    'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa',
    'v_cigar', 'd_cigar', 'j_cigar', 'np1_length', 'v_germline_start',
    'v_germline_end', 'd_germline_start', 'd_germline_end', 'j_germline_start',
    'j_germline_end', 'locus'
    ]

AIRR_FIELDS_TO_GENERATE = [
    'sequence_id', 'sequence', 'sequence_alignment', 'germline_alignment',
    'location', 'sample_time', 'junction', 'junction_aa', 'junction_length'
    ]

AIRR_FIELDS_TO_KEEP = [x for x in AIRR_REQUIRED_FIELDS if x not in AIRR_FIELDS_TO_GENERATE]

ALL_TREE_NAMES = [
    "full_tree",
    "pruned_tree", "pruned_time_tree",
    "simplified_tree", "simplified_time_tree"
    ]

TREE_NAMES = ["pruned_tree", "pruned_time_tree", "simplified_tree", "simplified_time_tree"]
MEMORY_SAVE_TREE_NAMES = ["simplified_tree", "simplified_time_tree"]

def read_sf5_table(filename):
    """Reads a CSV file containing the SF5 mutability table.
    Args:
        filename (str): The path to the CSV file.
    Returns:
        pd.DataFrame: A DataFrame containing the mutability data.
    """
    data = pd.read_csv(filename, header=0)
    data.fillna(0, inplace=True)
    return data

HEAVY_MUTABILITY_TABLE = read_sf5_table(get_data("hh_sf5.csv"))
LIGHT_MUTABILITY_TABLE = read_sf5_table(get_data("hkl_sf5.csv"))

NAIVE = pd.read_csv(get_data("naive_pairs_filtered.csv"), header=0)

HEAVY_SUBSTITUTION_TABLE = read_sf5_table(get_data("hh_sf5_substitution.csv"))
LIGHT_SUBSTITUTION_TABLE = read_sf5_table(get_data("hkl_sf5_substitution.csv"))


def translate_to_amino_acid(nucleotide_seq):
    """ Translates a nucleotide sequence into an amino acid sequence.
    Args:
        nucleotide_seq (str): The nucleotide sequence to translate.
    Returns:
        str: The translated amino acid sequence.
    """
    amino_acid_seq = ""
    for i in range(0, len(nucleotide_seq), 3):
        if i > len(nucleotide_seq) - 3:
            break
        codon = nucleotide_seq[i:i+3]
        amino_acid = codon_to_amino_acid(codon)
        amino_acid_seq += amino_acid
    return amino_acid_seq


def codon_to_amino_acid(codon):
    """Converts a codon (3-nucleotide sequence) to its corresponding amino acid.
    Args:
        codon (str): A 3-nucleotide sequence representing a codon.
    Returns:
        str: The corresponding amino acid represented by the codon.
    """
    table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',                 
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_',
        'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
        '...': 'X'
    }
    if codon not in table:
        return 'X'
    return table[codon]


def get_substitution_probability(kmer, heavy=True):
    """Gets the substitution probabilities for a given kmer.
    Args:
        kmer (str): The 5-mer sequence for which to get substitution probabilities.
        heavy (bool): Whether to use the heavy chain substitution table.
    Returns:
        list: A list of probabilities for each nucleotide substitution (A, C, G, T).
    """
    table = HEAVY_SUBSTITUTION_TABLE if heavy else LIGHT_SUBSTITUTION_TABLE
    row = table.loc[(table["Fivemer"]==kmer)]
    if row.empty:
        logger.error("%s not found in substitution table", kmer)
        exit(1)
    probabilities = [row['A'].values[0], row['C'].values[0], row['G'].values[0], row['T'].values[0]]
    return probabilities

def get_mutability_of_kmer(kmer, heavy=True):
    """Gets the mutability of a given kmer.
    Args:
        kmer (str): The 5-mer sequence for which to get mutability.
        heavy (bool): Whether to use the heavy chain mutability table.
    Returns:
        float: The mutability value for the kmer.
    """
    table = HEAVY_MUTABILITY_TABLE if heavy else LIGHT_MUTABILITY_TABLE
    mutability =  table.loc[table["Fivemer"]==kmer, "Mutability"].values[0]
    if math.isnan(mutability):
        logger.debug("NaN found")
        return 0

    return mutability


def remove_gaps(aligned):
    """Removes gaps from an aligned sequence.
    Args:
        aligned (str): The aligned sequence with gaps.
    Returns:
        str: The aligned sequence with gaps removed.
    """
    return aligned.replace(".", "")


def get_random_start_pair():
    """Generates a random start pair of heavy and light chains.
    Returns:
        StartPair: A named tuple containing the heavy and light chains.
    """
    row = NAIVE.sample(random_state=s._RNG) # pylint: disable=protected-access
    StartPair = namedtuple("RawStartPair", ["heavy", "light"])
    heavy = _format_random_start_chain(row, "heavy")
    light = _format_random_start_chain(row, "light")
    if len(heavy.input.aligned) < 312 or len(light.input.aligned) < 312:
        logger.warning("aligned sequence length is less than 312")
    return StartPair(heavy, light)


def _format_random_start_chain(row, chain_type):
    """Formats a random start chain from a row of the naive pairs DataFrame.
    Args:
        row (pd.Series): A row from the naive pairs DataFrame.
        chain_type (str): The type of chain ('heavy' or 'light').
    Returns:
        StartInfo: A named tuple containing the input and constants for the chain.
    """
    StartInput = namedtuple("StartInput", ["chain", "aligned", "cdr3_aa_length", "junction"])
    StartInfo = namedtuple("StartInfo", ["input", "constants"])
    get_cdr3_length = lambda x: int(len(x)/3)
    start_input = StartInput(
        remove_gaps(row[f'{chain_type}_aligned'].values[0]),
        row[f'{chain_type}_aligned'].values[0],
        get_cdr3_length(row[f'{chain_type}_cdr3'].values[0]),
        row[f'{chain_type}_junction'].values[0]
        )
    constants = {
        x:row[f'{chain_type}_{x}'].values[0] for x in AIRR_FIELDS_TO_KEEP
    }
    constants["germline_alignment"] = start_input.aligned
    return StartInfo(start_input, constants)


def make_plot(data, times, results_file, ylabel, title, log=False):
    """Creates a plot of the given data and saves it to a file.
    Args:
        data (np.ndarray): The data to plot.
        times (np.ndarray): The time points corresponding to the data.
        results_file (str): The file path to save the plot.
        ylabel (str): The label for the y-axis.
        title (str): The title of the plot.
        log (bool): Whether to use a logarithmic scale for the y-axis.
    """
    fig = plt.figure(title)
    ax = fig.gca()
    ax.plot(times, data)
    ax.set_xlabel("Time (generation)")
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    if log:
        ax.set_yscale("log", base=s.MULTIPLIER)
    fig.savefig(results_file)
    plt.close(fig)

def make_bar_plot(data, results_file, xlabel, title):
    """Creates a bar plot of the given data and saves it to a file.
    Args:
        data (np.ndarray): The data to plot.
        results_file (str): The file path to save the plot.
        xlabel (str): The label for the x-axis.
        title (str): The title of the plot.
    """
    fig = plt.figure(title)
    ax = fig.gca()
    ax.hist(data, bins = 30)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    fig.savefig(results_file)

def snake_case_to_normal(name):
    """Converts a snake_case string to a normal string with spaces.
    Args:
        name (str): The snake_case string to convert.
    Returns:
        str: The converted string with spaces instead of underscores.
    """
    return " ".join([x for x in name.split("_")])

def _axis_label(name):
    """Converts a snake_case string to a more readable axis label.
    Args:
        name (str): The snake_case string to convert.
    Returns:
        str: The converted string with spaces instead of underscores and capitalized.
    """
    return [x for x in name.split("_")][-1]

def make_all_plots(df, result_dir, simulation=False):
    """Creates plots for all columns in the DataFrame and saves them to files.
    Args:
        df (pd.DataFrame): The DataFrame containing the data to plot.
        result_dir (str): The directory to save the plots.
        simulation (bool): Whether the plots are for a simulation (affects title).
    """
    title_suffix = "(across all clones in simulation)" if simulation else ""
    times = df["time"].to_numpy()
    columns = df.columns
    for column in columns:
        if column == "time":
            continue
        if "affinity" in column:
            make_plot(
                df[column].to_numpy(),
                times,
                result_dir + f"/{column}.png", "Average affinity",
                f"Average {snake_case_to_normal(column)} {title_suffix}",
                log=True
                )
        elif "population" in column:
            make_plot(
                df[column].to_numpy(),
                times,
                result_dir + f"/{column}.png",
                "Fraction of population",
                f"{snake_case_to_normal(column).capitalize()} {title_suffix}"
                )
        else:
            make_plot(
                df[column].to_numpy(),
                times,
                result_dir + f"/{column}.png",
                f"Average {_axis_label(column)}",
                f"Average {snake_case_to_normal(column)} {title_suffix}")
