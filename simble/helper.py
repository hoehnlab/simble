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
from .constants import AIRR_REQUIRED_FIELDS, SIMBLE_REQUIRED_FIELDS, AIRR_FIELDS_TO_GENERATE

logger = logging.getLogger(__package__)

_ROOT = os.path.abspath(os.path.dirname(__file__))
def get_data(path):
    return os.path.join(_ROOT, 'data', path)

def get_naive_table():
    if s.NAIVE_FILE:
        naive = pd.read_csv(s.NAIVE_FILE, header=0)
        renamed_columns = {
            "heavy_sequence": "heavy",
            "light_sequence": "light",
            "heavy_sequence_alignment": "heavy_aligned",
            "light_sequence_alignment": "light_aligned"
        }
        naive.rename(columns=renamed_columns, inplace=True)
    else:
        naive = pd.read_csv(get_data("naive_pairs_filtered.csv"), header=0)
    return naive

AIRR_FIELDS_TO_KEEP = [x for x in AIRR_REQUIRED_FIELDS if x not in AIRR_FIELDS_TO_GENERATE]

def read_sf5_table(filename):
    data = pd.read_csv(filename, header=0)
    data.fillna(0, inplace=True)
    return data

HEAVY_MUTABILITY_TABLE = read_sf5_table(get_data("hh_sf5.csv"))
LIGHT_MUTABILITY_TABLE = read_sf5_table(get_data("hkl_sf5.csv"))

NAIVE = get_naive_table()
NAIVE_ROWS = len(NAIVE.index)

HEAVY_SUBSTITUTION_TABLE = read_sf5_table(get_data("hh_sf5_substitution.csv"))
LIGHT_SUBSTITUTION_TABLE = read_sf5_table(get_data("hkl_sf5_substitution.csv"))


def update_helper_tables():
    global NAIVE 
    global NAIVE_ROWS
    NAIVE = get_naive_table()
    NAIVE_ROWS = len(NAIVE.index)

def translate_to_amino_acid(nucleotide_seq):
    amino_acid_seq = ""
    for i in range(0, len(nucleotide_seq), 3):
        if i > len(nucleotide_seq) - 3:
            break
        codon = nucleotide_seq[i:i+3]
        amino_acid = codon_to_amino_acid(codon)
        amino_acid_seq += amino_acid
    return amino_acid_seq


def codon_to_amino_acid(codon):
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
    table = HEAVY_SUBSTITUTION_TABLE if heavy else LIGHT_SUBSTITUTION_TABLE
    row = table.loc[(table["Fivemer"]==kmer)]
    if row.empty:
        logger.error(f"{kmer} not found in substitution table")
        exit(1)
    probabilities = [row['A'].values[0], row['C'].values[0], row['G'].values[0], row['T'].values[0]]
    return probabilities

def get_mutability_of_kmer(kmer, heavy=True):
    table = HEAVY_MUTABILITY_TABLE if heavy else LIGHT_MUTABILITY_TABLE
    mutability =  table.loc[table["Fivemer"]==kmer, "Mutability"].values[0]
    if math.isnan(mutability):
        logger.debug("NaN found")
        return 0

    return mutability

def remove_gaps(aligned):
    return aligned.replace(".", "")


def get_start_pair(i=None):
    if i and not s.NAIVE_RANDOM:
        row = NAIVE.iloc[i % NAIVE_ROWS]
    else:
        row = NAIVE.sample(random_state=s._RNG)
    StartPair = namedtuple("RawStartPair", ["heavy", "light", "user_constants"])
    heavy = _format_start_chain(row, "heavy")
    light = _format_start_chain(row, "light")
    user_constants = {x: row[x] for x in s.USER_FIELDS_TO_KEEP}
    if len(heavy.input.aligned) < 312 or len(light.input.aligned) < 312:
        logger.warning("aligned sequence length is less than 312")
    return StartPair(heavy, light, user_constants)

def _format_start_chain(row, type):
    StartInput = namedtuple("StartInput", ["chain", "aligned", "cdr3_aa_length", "junction"])
    StartInfo = namedtuple("StartInfo", ["input", "constants"])
    get_cdr3_length = lambda x: int(len(x)/3)
    input = StartInput(
        remove_gaps(row[f'{type}_aligned'].values[0]), 
        row[f'{type}_aligned'].values[0], 
        get_cdr3_length(row[f'{type}_cdr3'].values[0]),
        row[f'{type}_junction'].values[0]
        )
    constants = {
        x:row[f'{type}_{x}'].values[0] for x in AIRR_FIELDS_TO_KEEP
    }
    constants["germline_alignment"] = input.aligned
    return StartInfo(input, constants)



def make_plot(data, times, results_file, ylabel, title, log=False):
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
    fig = plt.figure(title)
    ax = fig.gca()
    ax.hist(data, bins = 30)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("Frequency")
    ax.set_title(title)
    fig.savefig(results_file)

def snake_case_to_normal(name):
    return " ".join([x for x in name.split("_")])

def axis_label(name):
    return [x for x in name.split("_")][-1]

def make_all_plots(df, result_dir, simulation=False):
    title_suffix = "(across all clones in simulation)" if simulation else ""
    times = df["time"].to_numpy()
    columns = df.columns
    for column in columns:
        if column == "time":
            continue
        elif "affinity" in column:
            make_plot(df[column].to_numpy(), times, result_dir + f"/{column}.png", "Average affinity", f"Average {snake_case_to_normal(column)} {title_suffix}", log=True)
        elif "population" in column:
            make_plot(df[column].to_numpy(), times, result_dir + f"/{column}.png", "Fraction of population", f"{snake_case_to_normal(column).capitalize()} {title_suffix}")
        else:
            make_plot(
                df[column].to_numpy(), 
                times, 
                result_dir + f"/{column}.png", 
                f"Average {axis_label(column)}", 
                f"Average {snake_case_to_normal(column)} {title_suffix}")
