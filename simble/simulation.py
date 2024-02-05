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

from collections import Counter
import logging

import numpy as np
import pandas as pd

from .cell import Cell
from .dev_helper import get_data_points
from .helper import make_all_plots, make_bar_plot
from .location import Location, LocationName
from .settings import s
from .target import TargetAminoPair
from .tree import Node

logger = logging.getLogger(__package__)

def get_population_data(location, time):
    population = len(location.current_generation)
    children_counter = Counter(location.number_of_children)
    get_more_than = lambda x, y: sum([v for k, v in x.items() if k > y])
    children_dict = {
        f"number_of_cells_with_more_than_{i}_children": get_more_than(children_counter, i) 
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

    
def run_simulation(i, result_dir):
    dev_data_rows = []
    pop_data_rows = []
    time = 0
    clone_id = i+1
    naive = Cell(None, None, created_at=time)
    root = Node(naive, clone_id=clone_id)
    locations = [Location(x.name, x) for x in s.LOCATIONS]
    GC = [x for x in locations if x.name == LocationName.GC][0]
    OTHER = [x for x in locations if x.name == LocationName.OTHER][0]
    GC.current_generation = [root]
    fasta_string = naive.as_fasta(time)
    airr = []
    sampled_ids = []
    TARGET_PAIR = TargetAminoPair(
        naive.heavy_chain.get_gapped_sequence(), 
        naive.light_chain.get_gapped_sequence(), 
        naive.heavy_chain.CDR3_length, 
        naive.light_chain.CDR3_length)
    TARGET_PAIR.mutate(s.TARGET_MUTATIONS_HEAVY, s.TARGET_MUTATIONS_LIGHT)

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
        
        if time >= 25:
            migrate_size = int(s._RNG.poisson(location.settings.migration_rate))
            to_migrate = s._RNG.choice(
                current_generation, 
                size=migrate_size, 
                replace=False)
            current_generation = [x for x in current_generation if x not in to_migrate]


            for node in to_migrate:
                child_node = Node(node.cell.remake_self(), parent=node, generation=node.generation+1)
                node.add_child(child_node)
                OTHER.immigrating_population.append(child_node)
        
        if location.name == LocationName.OTHER:
            # other location reproduces very slowly
            # we only want these cells to have 2 children, unlike the GC
            # get number of cells to reproduce
            cells_to_reproduce = int(s._RNG.poisson(0.25))

            #randomly select cells to reproduce
            to_divide = s._RNG.choice(
                current_generation, 
                size=cells_to_reproduce, 
                replace=False)
            
            not_dividing = [x for x in current_generation if x not in to_divide]
            
            cells_to_live = min(location.settings.max_population - cells_to_reproduce*2, len(not_dividing))
            to_live = s._RNG.choice(
                not_dividing,
                size=cells_to_live,
                replace=False
                )

            def update_antigen(x, i):
                x.antigen = i
            [update_antigen(x, 1) for x in to_live]
            [update_antigen(x, 2) for x in to_divide]
            
        else:
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

        location.number_of_children = [x.antigen for x in current_generation]
        for node in current_generation:
            node.cell.kill_cell()
            children = [make_new_child(node) for _ in range(min(node.antigen, 10))]
            live_children = [x for x in children if x.cell.is_alive]
            new_generation.extend(live_children)
        

        return new_generation
    
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
                sample_size = min(len(location.current_generation)//2, location.settings.sample_size)
                current_sample = s._RNG.choice(location.current_generation, size=sample_size, replace=False)
                location.current_generation = [x for x in location.current_generation if x not in current_sample]
                for node in current_sample:
                    sampled_ids.append(id(node.cell))
                    fasta_string += node.cell.as_fasta(time)
                    airr.extend(node.cell.as_AIRR(time))

        time += 1

    newick = f'({root.write_newick()});'
    pruned = root.prune_subtree(sampled_ids)
    pruned_newick = f'({pruned.write_newick()});'
    pruned_time_tree = f'({pruned.write_newick(time_tree=True)});'
    if s.DEV:
        with open(result_dir + "/all_samples.fasta", "w") as f:
            f.write(fasta_string)

        logger.info("writing newick tree")
        with open(result_dir + "/true_tree.tree", "w") as f:
            f.write(newick)

        logger.info("writing pruned newick tree")
        with open(result_dir + "/pruned_tree.tree", "w") as f:
            f.write(pruned_newick)


    df = pd.DataFrame(dev_data_rows)
    pop_data = pd.DataFrame(pop_data_rows)
    airr = pd.DataFrame(airr)
    airr["sequence_id"] = airr["sequence_id"].apply(lambda x: f"{clone_id}_{x}")
    airr["cell_id"] = airr["cell_id"].apply(lambda x: f"{clone_id}_{x}")
    airr["clone_id"] = clone_id
    pop_data["clone_id"] = clone_id
    
    if s.DEV:
        logger.info(f"max affinity was: {TARGET_PAIR.max_affinity}")
        logger.info("making plots")
        # make plots
        make_all_plots(df, result_dir)
        make_bar_plot(list(TARGET_PAIR.heavy.cdr_multipliers.values()), result_dir + "/cdr_multiplier.png", "CDR multiplier value", "CDR multiplier distribution")
        make_bar_plot(list(TARGET_PAIR.heavy.fwr_multipliers.values()), result_dir + "/fwr_multiplier.png", "FWR multiplier value", "FWR multiplier distribution")

    return {
        "airr": airr, 
        "fasta": fasta_string, 
        "true_tree": newick, 
        "pruned_tree": pruned_newick,
        "pruned_time_tree": pruned_time_tree, 
        "data": df, 
        "clone_id": clone_id, 
        "pop_data": pop_data
        }


if __name__ == "__main__":
    run_simulation(1, f's.RESULTS_DIR/results1/')