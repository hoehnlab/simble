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
from collections import Counter

import numpy as np
import pandas as pd
from tqdm import tqdm

from .cell import Cell, CellType
from .dev_helper import get_data_points
from .helper import make_all_plots, make_bar_plot
from .location import Location, LocationName
from .settings import s
from .target import TargetAminoPair
from .tree import Node, simplify_tree

logger = logging.getLogger(__package__)

def get_population_data(location, time):
    population = len(location.current_generation)

    children_counter = Counter(location.number_of_children)
    children_dict = {
        f"number_of_cells_with_{i}_children": children_counter[i] 
        for i in range(1, 10)
        }
    pop_data = {
        "time": time,
        "location": location.name.value,
        "population": population,
        "number_of_reproducing_cells": population - children_counter[0],
        "average_affinity": np.mean([x.cell.affinity for x in location.current_generation]) if population > 0 else 0,
    }
    pop_data.update(children_dict)
    return pop_data

def do_differentiation(location, time):
    current_generation = location.current_generation
    to_migrate = []
    migrate_size = min(int(s._RNG.poisson(location.settings.migration_rate)), len(current_generation)//2)
    def get_mbc_pc_size(migrate_size, time):
        current_day = s.GENERATIONS_PER_DAY*time
        if current_day < 8:
            mbc_size = int(migrate_size * 0.99)
            pc_size = migrate_size - mbc_size
        elif current_day < 16:
            days = current_day - 8
            percentage_mbc = 0.99 - (days * 0.98/8)
            mbc_size = int(migrate_size * percentage_mbc)
            pc_size = migrate_size - mbc_size
        elif current_day < 41:
            pc_size = int(migrate_size * 0.99)
            mbc_size = migrate_size - pc_size
        else:
            pc_size = migrate_size
            mbc_size = 0
        return mbc_size, pc_size
        

    mbc_size, pc_size = get_mbc_pc_size(migrate_size, time)
    if mbc_size > 0:
        mbcs = s._RNG.choice(
            current_generation, 
            size=mbc_size, 
            replace=False)
        current_generation = [x for x in current_generation if x not in mbcs]
        [mbc.cell.differentiate(CellType.MBC) for mbc in mbcs]
        to_migrate.extend(mbcs)
    if pc_size > 0:
        affinities = [x.cell.affinity for x in current_generation]
        p = np.array(affinities) / np.sum(affinities) if s.SELECTION else None
        pcs = s._RNG.choice(
            current_generation, 
            size=pc_size,
            p=p,
            replace=False)
        current_generation = [x for x in current_generation if x not in pcs]
        [pc.cell.differentiate(CellType.PC) for pc in pcs]
        to_migrate.extend(pcs)
    return to_migrate

def non_gc_population_control(current_generation):
    # right now the other location will have no reproduction or death
    # TODO: add a max population
    # TODO: potentially tweak reproduction 
    new_generation = []
    for node in current_generation:
        child_node = Node(node.cell.remake_self(), parent=node, generation=node.generation+1)
        node.add_child(child_node)
        new_generation.append(child_node)
    return new_generation


def simulate(clone_id, TARGET_PAIR, gc_start_generation, root, time=0):
    dev_data_rows = []
    pop_data_rows = []
    naive = root.cell
    locations = [Location(x.name, x) for x in s.LOCATIONS]
    GC = [x for x in locations if x.name == LocationName.GC][0]
    OTHER = [x for x in locations if x.name == LocationName.OTHER][0]
    GC.current_generation = gc_start_generation
    # fasta_string = naive.as_fasta(time)
    airr = []
    sampled_ids = []
    sampled = []
    # TARGET_PAIR.mutate(s.TARGET_MUTATIONS_HEAVY, s.TARGET_MUTATIONS_LIGHT)

    def make_new_child(node):
        child_cell = Cell(
            node.cell.heavy_chain.copy(), 
            node.cell.light_chain.copy(), 
            location=node.cell.location,
            created_at=time)
        heavy_n, light_n = child_cell.mutate_cell()
        child_node = Node(child_cell, parent=node, heavy_mutations=heavy_n, light_mutations=light_n, generation=node.generation+1)
        child_cell.calculate_affinity(TARGET_PAIR)
        node.add_child(child_node)
        return child_node
        

    def make_new_generation(location):
        new_generation = [] 

        current_generation = location.current_generation
        if len(current_generation) == 0:
            return []

        if location.name == LocationName.GC:
            to_migrate = do_differentiation(location, time)
        else:
            # TODO: potentially allow other locations to migrate in the future
            to_migrate = []

        current_generation = [x for x in current_generation if x not in to_migrate]

        for node in to_migrate:
            child_node = Node(node.cell.remake_self(), parent=node, generation=node.generation+1)
            node.add_child(child_node)
            OTHER.immigrating_population.append(child_node)
        
        if location.name == LocationName.OTHER:
            new_generation = non_gc_population_control(current_generation)
            return new_generation
            
        available_antigen = location.settings.max_population

        if s.SELECTION and location.name == LocationName.GC:
            affinities = [x.cell.affinity for x in current_generation]
            p = np.array(affinities) / np.sum(affinities)
    
        else:
            p = None
            
        for _ in range(available_antigen):
            current_node = s._RNG.choice(
                current_generation,
                p=p
                )
            current_node.antigen += 1

        location.number_of_children = [min(x.antigen, 10) for x in current_generation]
        for node in current_generation:
            node.cell.kill_cell()
            children = [make_new_child(node) for _ in range(min(node.antigen, 10))]
            live_children = [x for x in children if x.cell.is_alive]
            new_generation.extend(live_children)
            if node.antigen == 0 and (s.MEMORY_SAVE or not s.KEEP_FULL_TREE):
                node.prune_up_tree()
        

        return new_generation
    
    bar = tqdm(total=s.END_TIME-1, initial=0, desc=f"Clone {clone_id}", position=clone_id, leave=True)
    while time<s.END_TIME:
        for location in locations:
            location.current_generation = make_new_generation(location)

        targets = lambda x: (3*x, 3*x+1, 3*x+2)
        
        row = get_data_points(
            GC.current_generation, 
            time, 
            naive.heavy_chain.get_gapped_sequence(), 
            naive.light_chain.get_gapped_sequence(),
            [i for x in TARGET_PAIR.heavy.mutation_locations for i in targets(x)],
            [i for x in TARGET_PAIR.light.mutation_locations for i in targets(x)])

        dev_data_rows.append(row)

        if time % 25 ==0:
            logger.debug(f"Time: {time}, population: {len(GC.current_generation)}")
        
        [x.finish_migration() for x in locations]

        for location in locations:
            pop_data_rows.append(get_population_data(location, time))
            if time in location.settings.sample_times:
                if time == 0:
                    # make sure we don't remove the naive cell
                    continue
                if time == s.END_TIME-1:
                    sample_size = min(len(location.current_generation), location.settings.sample_size)
                else:
                    sample_size = min(len(location.current_generation)//2, location.settings.sample_size)
                current_sample = s._RNG.choice(location.current_generation, size=sample_size, replace=False)
                location.current_generation = [x for x in location.current_generation if x not in current_sample]
                for node in current_sample:
                    sampled_ids.append(id(node.cell))
                    node.sampled_time = time
                    sampled.append(node)
                    airr.extend(node.cell.as_AIRR(time))

        time += 1
        if time<s.END_TIME:
            bar.update()


    [x.finish_migration() for x in locations]
    for location in locations:
        [node.prune_up_tree() for node in location.current_generation]

    df = pd.DataFrame(dev_data_rows)
    pop_data = pd.DataFrame(pop_data_rows)
    pop_data["clone_id"] = clone_id

    bar.bar_format = "{desc}: |{bar}| {n}/{total} in {elapsed}"
    bar.refresh()
    bar.close()
    
    return sampled, pop_data, df

    
def run_simulation(i, result_dir):
    time = 0
    clone_id = i+1
    naive = Cell(None, None, created_at=time)
    root = Node(naive, clone_id=clone_id)
    airr = []
    TARGET_PAIR = TargetAminoPair(
        naive.heavy_chain.get_gapped_sequence(), 
        naive.light_chain.get_gapped_sequence(), 
        naive.heavy_chain.CDR3_length, 
        naive.light_chain.CDR3_length)
    TARGET_PAIR.mutate(s.TARGET_MUTATIONS_HEAVY, s.TARGET_MUTATIONS_LIGHT)

    sampled, pop_data, dev_df = simulate(clone_id, TARGET_PAIR, [root], root)

    sampled_ids = [id(x.cell) for x in sampled]
    fasta_string = "".join([x.cell.as_fasta(x.sampled_time) for x in sampled])

    airr = [x for node in sampled for x in node.cell.as_AIRR(node.sampled_time)]
    airr = pd.DataFrame(airr)
    airr["sequence_id"] = airr["sequence_id"].apply(lambda x: f"{clone_id}_{x}")
    airr["cell_id"] = airr["cell_id"].apply(lambda x: f"{clone_id}_{x}")
    airr["clone_id"] = clone_id


    if s.MEMORY_SAVE:
        # in memory saving mode we don't keep the full tree
        newick = ""
        pruned_newick = ""
        pruned_time_tree = ""
        pruned = root
    elif s.KEEP_FULL_TREE:
        newick = f'{root.write_newick()};'
        # TODO (jf): add an interative function that can prune the tree
        pruned = root.prune_subtree(sampled_ids)
        pruned_newick = f'{pruned.write_newick()};'
        pruned_time_tree = f'{pruned.write_newick(time_tree=True)};'
    else:
        newick = ""
        pruned = root
        pruned_newick = f'{root.write_newick()};'
        pruned_time_tree = f'{root.write_newick(time_tree=True)};'

    simplified_tree = simplify_tree(pruned)
    simplified_tree_newick = f'{simplified_tree.write_newick()};'
    simplified_time_tree_newick = f'{simplified_tree.write_newick(time_tree=True)};'

    # TODO (jf): clean up dev code and logging with new tree options
    if s.DEV:
        with open(result_dir + "/all_samples.fasta", "w") as f:
            f.write(fasta_string)

        # logger.info("writing newick tree")
        # with open(result_dir + "/true_tree.tree", "w") as f:
        #     f.write(newick)

        logger.info("writing pruned newick tree")
        with open(result_dir + "/pruned_tree.tree", "w") as f:
            f.write(pruned_newick)

        logger.info("writing simplified newick tree")
        with open(result_dir + "/simplified_time_tree.tree", "w") as f:
            f.write(simplified_tree_newick)
        
        logger.info("writing simplified newick time tree")
        with open(result_dir + "/simplified_time_tree.tree", "w") as f:
            f.write(simplified_time_tree_newick)
    
    if s.DEV:
        logger.info(f"max affinity was: {TARGET_PAIR.max_affinity}")
        logger.info("making plots")
        # make plots
        make_all_plots(dev_df, result_dir)
        make_bar_plot(list(TARGET_PAIR.heavy.cdr_multipliers.values()), result_dir + "/cdr_multiplier.png", "CDR multiplier value", "CDR multiplier distribution")
        make_bar_plot(list(TARGET_PAIR.heavy.fwr_multipliers.values()), result_dir + "/fwr_multiplier.png", "FWR multiplier value", "FWR multiplier distribution")

    return {
        "airr": airr, 
        "fasta": fasta_string, 
        "full_tree": newick, 
        "pruned_tree": pruned_newick,
        "pruned_time_tree": pruned_time_tree, 
        "simplified_tree": simplified_tree_newick,
        "simplified_time_tree": simplified_time_tree_newick,
        "data": dev_df, 
        "clone_id": clone_id, 
        "pop_data": pop_data,
        "targets": {"clone_id": clone_id, "heavy": TARGET_PAIR.heavy.amino_acid_seq, "light": TARGET_PAIR.light.amino_acid_seq}
        }


if __name__ == "__main__":
    run_simulation(1, f's.RESULTS_DIR/results1/')