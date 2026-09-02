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

AIRR_REQUIRED_FIELDS = ['sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'np1_length', 'v_germline_start', 'v_germline_end', 'd_germline_start', 'd_germline_end', 'j_germline_start', 'j_germline_end', 'germline_alignment_d_mask', 'locus']

AIRR_FIELDS_TO_GENERATE = ['sequence_id', 'sequence', 'sequence_alignment', 'germline_alignment'] 

AIRR_FIELDS_TO_KEEP = [x for x in AIRR_REQUIRED_FIELDS if x not in AIRR_FIELDS_TO_GENERATE]

FIELDS_NEEDED_AS_INPUT = ["sequence", "sequence_alignment", "cdr3"]

def create_naive_table(heavy_file, light_file, join_id, keep_cols=None, keep_from=None):
    use_cols = FIELDS_NEEDED_AS_INPUT + AIRR_FIELDS_TO_KEEP
    use_cols_heavy = use_cols
    use_cols_light = use_cols
    if keep_cols:
        if keep_from == "heavy":
            use_cols_heavy = use_cols + keep_cols
        else:
            use_cols_light = use_cols + keep_cols

    heavy = pd.read_csv(heavy_file, sep = '\t', header=0, usecols=use_cols_heavy, converters={"productive": lambda x: x == "TRUE"})
    light = pd.read_csv(light_file, sep = '\t', header=0, usecols=use_cols_light, converters={"productive": lambda x: x == "TRUE"})
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
    parser.add_argument("-j", "--join", 
                         dest="id", 
                         help="column to join heavy and light tables on", 
                         metavar="COL",
                         type=str)
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