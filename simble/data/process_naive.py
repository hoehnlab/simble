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

import pandas as pd
import argparse

AIRR_REQUIRED_FIELDS = ['sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'np1_length', 'v_germline_start', 'v_germline_end', 'd_germline_start', 'd_germline_end', 'j_germline_start', 'j_germline_end', 'locus']

AIRR_FIELDS_TO_GENERATE = ['sequence_id', 'sequence', 'sequence_alignment', 'germline_alignment'] 

# CGJ
GERMLINE_LENGTH_FIELDS = ["v_germline_length", "d_germline_length", "j_germline_length"]

# CGJ
ALIGNMENT_FIELDS = GERMLINE_LENGTH_FIELDS + ["np2_length"]

# CGJ
AIRR_FIELDS_TO_KEEP = [x for x in AIRR_REQUIRED_FIELDS + ALIGNMENT_FIELDS if x not in AIRR_FIELDS_TO_GENERATE]

FIELDS_NEEDED_AS_INPUT = ["sequence", "sequence_alignment", "cdr3"]

# CGJ
# AIRR writes booleans as T and F, but some tools write TRUE and FALSE instead
PRODUCTIVE_VALUES = {"T", "TRUE", "True", "true", "1"}

# CGJ
def is_productive(value):
    """Reads an AIRR productive field, whichever way the boolean is spelled.

    Args:
        value (str): The raw value of the productive column.
    Returns:
        bool: True if the row is marked productive.
    """
    return str(value).strip() in PRODUCTIVE_VALUES

# CGJ
def resolve_germline_lengths(filename, use_cols):
    """Works out which germline length columns have to be derived from the file.

    The lengths are optional in AIRR, so a table may carry only the germline
    start and end of each gene. A length that is missing but whose start and end
    are both present is dropped from the columns to read and derived afterwards.
    A length that is missing with no start and end to derive it from is left in
    place, so reading the file fails the way it always has.

    Args:
        filename (str): The path to the AIRR tsv file.
        use_cols (list): The columns to read from the file.
    Returns:
        tuple: The columns to read, and a list of (length, start, end) triples
            naming the lengths to derive once the file has been read.
    """
    available = set(pd.read_csv(filename, sep='\t', header=0, nrows=0).columns)
    derive = []
    for length in [x for x in GERMLINE_LENGTH_FIELDS if x not in available]:
        gene = length.split("_")[0]
        start, end = f"{gene}_germline_start", f"{gene}_germline_end"
        if start in available and end in available:
            derive.append((length, start, end))
    derived_names = [x[0] for x in derive]
    return [x for x in use_cols if x not in derived_names], derive

# CGJ
def add_germline_lengths(table, derive):
    """Adds the germline length columns that were missing from the input.

    Args:
        table (pd.DataFrame): A table read from an AIRR tsv file.
        derive (list): The (length, start, end) triples from
            resolve_germline_lengths.
    """
    for length, start, end in derive:
        # the AIRR coordinates are 1-based and closed, so both ends count
        table[length] = table[end] - table[start] + 1

def create_naive_table(heavy_file, light_file, join_id, keep_cols=None, keep_from=None):
    # CGJ
    use_cols = FIELDS_NEEDED_AS_INPUT + AIRR_FIELDS_TO_KEEP + [join_id]
    use_cols_heavy = use_cols
    use_cols_light = use_cols
    if keep_cols:
        if keep_from == "heavy":
            use_cols_heavy = use_cols + keep_cols
        else:
            use_cols_light = use_cols + keep_cols

    # CGJ
    use_cols_heavy, derive_heavy = resolve_germline_lengths(heavy_file, use_cols_heavy)
    use_cols_light, derive_light = resolve_germline_lengths(light_file, use_cols_light)

    heavy = pd.read_csv(heavy_file, sep = '\t', header=0, usecols=use_cols_heavy, converters={"productive": is_productive})
    light = pd.read_csv(light_file, sep = '\t', header=0, usecols=use_cols_light, converters={"productive": is_productive})
    # CGJ
    add_germline_lengths(heavy, derive_heavy)
    add_germline_lengths(light, derive_light)
    heavy = heavy[heavy['productive'] == True]
    light = light[light['productive'] == True]

    def rename_columns(type):
        columns = {x:f"{type}_{x}" for x in AIRR_FIELDS_TO_KEEP}

        columns["sequence"] = f"{type}_sequence"
        columns["sequence_alignment"] = f"{type}_sequence_alignment"
        columns["cdr3"] = f"{type}_cdr3"

        return columns

    heavy = heavy.rename(columns=rename_columns("heavy"))
    light = light.rename(columns=rename_columns("light"))
    merged = pd.merge(heavy, light, on=join_id)

    return merged

def main():
    parser = argparse.ArgumentParser(
        prog="process_naive_data",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="somewhat improved model of B-cell lineage evolution \n\nsource code available at: www.github.com/hoehnlab/simble",
        epilog="It's that simble!"
        )
    parser.add_argument("-o", "--output", 
                         dest="out", 
                         help="absolute path to output file", 
                         metavar="FILE", 
                         type=str)
    parser.add_argument("--heavy", 
                         dest="heavy", 
                         help="path to heavy chain AIRR tsv file", 
                         metavar="FILE", 
                         type=str)
    parser.add_argument("--light", 
                         dest="light", 
                         help="path to light chain AIRR tsv file", 
                         metavar="FILE", 
                         type=str)
    # CGJ
    parser.add_argument("-j", "--join",
                         dest="id",
                         help="column to join heavy and light tables on (default: cell_id)",
                         metavar="COL",
                         type=str,
                         default="cell_id")
    parser.add_argument("--keep-cols", 
                         dest="keep_cols",
                         help="additional columns from the naive data to keep",
                         metavar="COL",
                         type=str,
                         nargs="+",
                         default=None)
    parser.add_argument("--keep-from",
                         dest="keep_from", 
                         help="keep additional columns from heavy or light table",
                         default="heavy", 
                         choices=["heavy", "light"],
                         type=str)

    args = parser.parse_args()
    naive = create_naive_table(args.heavy, args.light, args.id, args.keep_cols, args.keep_from)
    naive.to_csv(args.out, index=False)


if __name__ == "__main__":
    main()