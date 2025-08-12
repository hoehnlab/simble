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

AIRR_REQUIRED_FIELDS = ['sequence_id', 'sequence', 'rev_comp', 'productive',
                        'v_call', 'd_call', 'j_call', 'sequence_alignment',
                        'germline_alignment', 'junction', 'junction_aa', 'v_cigar',
                        'd_cigar', 'j_cigar', 'np1_length', 'v_germline_start',
                        'v_germline_end', 'd_germline_start', 'd_germline_end',
                        'j_germline_start', 'j_germline_end',
                        'germline_alignment_d_mask', 'locus']

ALIGNMENT_FIELDS = ["v_germline_length", "d_germline_length", "j_germline_length",
                    'np2_length']

AIRR_FIELDS_TO_GENERATE = ['sequence_id', 'sequence', 'sequence_alignment', 'germline_alignment']

AIRR_FIELDS_TO_KEEP = [x for x in AIRR_REQUIRED_FIELDS + ALIGNMENT_FIELDS if x not in AIRR_FIELDS_TO_GENERATE]

FIELDS_NEEDED_AS_INPUT = ["new_cell_id", "sequence", "sequence_alignment", "cdr3"]
def create_naive_table(heavy_file, light_file):
    use_cols = FIELDS_NEEDED_AS_INPUT + AIRR_FIELDS_TO_KEEP
    heavy = pd.read_csv(heavy_file, sep = '\t', header=0, usecols=use_cols)
    light = pd.read_csv(light_file, sep = '\t', header=0, usecols=use_cols, converters={"productive": lambda x: x == "TRUE"})
    heavy = heavy[heavy['productive'] == True]
    light = light[light['productive'] == True]

    def rename_columns(type):
        columns = {x:f"{type}_{x}" for x in AIRR_FIELDS_TO_KEEP}

        columns["new_cell_id"] = "cell_id"
        columns["sequence"] = f"{type}"
        columns["sequence_alignment"] = f"{type}_aligned"
        columns["cdr3"] = f"{type}_cdr3"

        return columns

    heavy = heavy.rename(columns=rename_columns("heavy"))
    light = light.rename(columns=rename_columns("light"))
    merged = pd.merge(heavy, light, on="cell_id")

    return merged



NAIVE = create_naive_table("naive_heavy_chain.tsv", "naive_light_chain.tsv")

NAIVE.to_csv("naive_pairs_filtered.csv", index=False)