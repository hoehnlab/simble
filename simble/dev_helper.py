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

from .settings import s

def get_data_points(current_generation, time, heavy_germline_gapped, light_germline_gapped, heavy_targets, light_targets):
    similarity_heavy = 0
    similarity_light = 0
    cdr_similarity_heavy = 0
    cdr_similarity_light = 0
    fwr_similarity_heavy = 0
    fwr_similarity_light = 0
    affinity = 0
    heavy_chain_affinity = 0
    population = 0
    population_heavy = 0
    population_light = 0
    heavy_shm = 0
    heavy_shm_filtered = 0
    light_shm = 0
    light_shm_filtered = 0
    heavy_cdr_shm = 0
    heavy_fwr_shm = 0
    light_cdr_shm = 0
    light_fwr_shm = 0
    n = len(current_generation)
    for node in current_generation:
        similarity_heavy += node.cell.heavy_chain.similarity
        similarity_light += node.cell.light_chain.similarity
        cdr_similarity_heavy += node.cell.heavy_chain.cdr_similarity
        cdr_similarity_light += node.cell.light_chain.cdr_similarity
        fwr_similarity_heavy += node.cell.heavy_chain.fwr_similarity
        fwr_similarity_light += node.cell.light_chain.fwr_similarity
        affinity += node.cell.affinity
        heavy_chain_affinity += node.cell.heavy_chain.affinity
        population += 1 if node.cell.affinity == s.MULTIPLIER**(s.TARGET_MUTATIONS_HEAVY+s.TARGET_MUTATIONS_LIGHT) else 0
        population_heavy += 1 if node.cell.heavy_chain.affinity == s.MULTIPLIER**s.TARGET_MUTATIONS_HEAVY else 0
        population_light += 1 if node.cell.light_chain.affinity == s.MULTIPLIER**s.TARGET_MUTATIONS_LIGHT else 0
        a, b, c, d = node.cell.heavy_chain.get_observed_mutations(heavy_germline_gapped, heavy_targets)
        heavy_shm += a
        heavy_shm_filtered += b
        heavy_cdr_shm += c
        heavy_fwr_shm += d
        a, b, c, d = node.cell.light_chain.get_observed_mutations(light_germline_gapped, light_targets)
        light_shm += a
        light_shm_filtered += b
        light_cdr_shm += c
        light_fwr_shm += d
        

    return {
        "time": time,
        "affinity": affinity/n,
        "heavy_chain_affinity": heavy_chain_affinity/n,
        "similarity": (similarity_heavy + similarity_light)/(2*n),
        "heavy_similarity": similarity_heavy/n,
        "light_similarity": similarity_light/n,
        "heavy_cdr_similarity": cdr_similarity_heavy/n,
        "light_cdr_similarity": cdr_similarity_light/n,
        "heavy_fwr_similarity": fwr_similarity_heavy/n,
        "light_fwr_similarity": fwr_similarity_light/n,
        "population_with_matching_sequence": population/n,
        "population_with_matching_heavy_sequence": population_heavy/n,
        "population_with_matching_light_sequence": population_light/n,
        "heavy_shm": heavy_shm/n,
        "light_shm": light_shm/n,
        "heavy_(non-target)_shm": heavy_shm_filtered/n,
        "light_(non-target)_shm": light_shm_filtered/n,
        "heavy_cdr_shm": heavy_cdr_shm/n,
        "heavy_fwr_shm": heavy_fwr_shm/n,
        "light_cdr_shm": light_cdr_shm/n,
        "light_fwr_shm": light_fwr_shm/n,
    }