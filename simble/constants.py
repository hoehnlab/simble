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

AIRR_REQUIRED_FIELDS = ['sequence_id', 'sequence', 'rev_comp', 'productive', 'v_call', 'd_call', 'j_call', 'sequence_alignment', 'germline_alignment', 'junction', 'junction_aa', 'v_cigar', 'd_cigar', 'j_cigar', 'np1_length', 'v_germline_start', 'v_germline_end', 'd_germline_start', 'd_germline_end', 'j_germline_start', 'j_germline_end', 'germline_alignment_d_mask', 'locus']

SIMBLE_REQUIRED_FIELDS = ['cdr3'] + AIRR_REQUIRED_FIELDS
SIMBLE_REQUIRED_FIELDS.remove('sequence_id')
SIMBLE_REQUIRED_FIELDS.remove('germline_alignment')

AIRR_FIELDS_TO_GENERATE = ['sequence_id', 'sequence', 'sequence_alignment', 'germline_alignment', 'location', 'sample_time', 'junction', 'junction_aa', 'junction_length']