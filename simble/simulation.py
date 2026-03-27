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
# pylint: disable=expression-not-assigned

import logging
from collections import Counter

import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import os
import warnings

from matplotlib.animation import FuncAnimation

from .cell import Cell, CellType
from .dev_helper import get_data_points
from .helper import make_all_plots, make_bar_plot
from .location import Location, LocationName
from .settings import s
from .target import TargetAminoPair
from .tree import Node, simplify_tree
from .animation import SimbleEvent, make_animation_html

logger = logging.getLogger(__package__)

MAX_CHILDREN = 100

def get_population_data(location, time):
    """Calculates population data for a given location at a specific time.

    Args:
        location (Location): The location for which to calculate population data.
        time (int): The current time in the simulation.
    Returns:
        dict: A dictionary containing population data, including the number of cells with children.
    """
    population = len(location.current_generation)

    children_counter = Counter(location.number_of_children)
    children_dict = {
        f"number_of_cells_with_{i}_children": children_counter[i]
        for i in range(1, MAX_CHILDREN + 1)
        }
    pop_data = {
        "time": time,
        "location": location.name.value,
        "population": population,
        "number_of_reproducing_cells": sum(children_counter.values()) - children_counter[0],
        "average_affinity": np.mean(
            [
                x.cell.affinity
                for x in location.current_generation
                ]
            ) if population > 0 else 0,
    }
    pop_data.update(children_dict)
    return pop_data

def do_differentiation(location, time):
    """Handles the differentiation of cells as they leave a location.
    Currently, this is only implemented for the germinal center (GC) location.

    Args:
        location (Location): The germinal center location.
        time (int): The current time in the simulation.
    Returns:
        list: A list of nodes that are migrating out of the germinal center.
    """
    current_generation = location.current_generation
    if s.SELECTION:
        current_generation = [x for x in current_generation if x.cell.affinity != 0]
    to_migrate = []
    migrate_size = min(
        int(
            s.RNG.poisson(location.settings.migration_rate)
            ),
            len(current_generation)//2
            )
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
        mbcs = s.RNG.choice(
            current_generation,
            size=mbc_size,
            replace=False)
        current_generation = [x for x in current_generation if x not in mbcs]
        [mbc.cell.differentiate(CellType.MBC) for mbc in mbcs]
        to_migrate.extend(mbcs)
    if pc_size > 0:
        affinities = [x.cell.affinity for x in current_generation]
        p = np.array(affinities) / np.sum(affinities) if s.SELECTION else None
        pcs = s.RNG.choice(
            current_generation,
            size=pc_size,
            p=p,
            replace=False)
        current_generation = [x for x in current_generation if x not in pcs]
        [pc.cell.differentiate(CellType.PC) for pc in pcs]
        to_migrate.extend(pcs)

    for node in to_migrate:
        node.last_migration = time

    return to_migrate


def do_reentry(location, time):
    """Handles the reentry of cells into the GC location.
    Currently, this is only implemented for the other (OTHER) location.

    Args:
        location (Location): The germinal center location.
        time (int): The current time in the simulation.
    Returns:
        list: A list of nodes that are migrating back into the germinal center.
    """
    if location.name != LocationName.OTHER:
        return []
    if s.PLANE is None:
        return []
    
    # cells_to_reenter = []
    reasonable_gens = 15
    # TODO (jf): make this a setting and allow users to specify
    reentry_size = int(s.RNG.poisson(0.02))
    potential_reentry = []
    for i in range(1, reasonable_gens, 1):
        potential_reentry.extend([x for x in location.current_generation if x.last_migration == time - i])
        if i >= 10 and len(potential_reentry) >= reentry_size:
            break

    # migration_times = [time - x.last_migration for x in potential_reentry]
    if len(potential_reentry) == 0:
        return []
    
    cells_to_reenter = s.RNG.choice(
                potential_reentry,
                size=min(reentry_size, len(potential_reentry)),
                replace=False
            )

    # print(f"len(cells_to_reenter) = {len(cells_to_reenter)} and reentry_size = {reentry_size}")
    [x.cell.differentiate(CellType.DEFAULT) for x in cells_to_reenter]
    return cells_to_reenter

def non_gc_population_control(current_generation):
    """Handles population control for non-GC locations.

    Args:
        current_generation (list): The current population in the non-GC location.
    Returns:
        list: A new generation of nodes, where each node is a child of the original nodes.
    """
    # TODO (jf): implement ability of other location to reproduce
    new_generation = []
    for node in current_generation:
        child_node = Node(node.cell.remake_self(), parent=node, generation=node.generation+1)
        node.add_child(child_node)
        new_generation.append(child_node)
    return new_generation


def simulate(clone_id, TARGET_PAIR, gc_start_generation, root, time=0): # pylint: disable=invalid-name
    """Runs the simulation for a single clone.

    Args:
        clone_id (int): The ID of the clone.
        TARGET_PAIR (TargetAminoPair): The target amino acid pair for the simulation.
        gc_start_generation (list): The initial population in the germinal center.
        root (Node): The root node of the simulation tree.
        time (int): The current time in the simulation.
    Returns:
        tuple: A tuple containing the sampled nodes, population data, and development data.
    """
    dev_data_rows = []
    pop_data_rows = []
    naive = root.cell
    locations = [Location(x.name, x) for x in s.LOCATIONS]
    GC = [x for x in locations if x.name == LocationName.GC][0] # pylint: disable=invalid-name
    OTHER = [x for x in locations if x.name == LocationName.OTHER][0] # pylint: disable=invalid-name
    GC.current_generation = gc_start_generation
    # fasta_string = naive.as_fasta(time)
    airr = []
    sampled_ids = []
    sampled = []

    event_handler = SimbleEvent(min_interval=200)

    if s.PLANE:
        # if not s.DEV:
        #     plt.ioff()
        if not s.DEV:
            logging.getLogger('matplotlib').setLevel(level=logging.CRITICAL)
            plt.switch_backend("Agg")
        simulation_video_location = f"{s.RESULTS_DIR}/{clone_id}/video"
        simulation_frame_location = f"{simulation_video_location}/animation_frames/"
        if not os.path.exists(simulation_frame_location):
            os.makedirs(simulation_frame_location)
        with open(f"{simulation_video_location}/animation.html", "w", encoding="utf-8") as f:
            f.write(make_animation_html(clone_id, s.END_TIME+1))
        
        fig, (ax, ax2) = plt.subplots(2, 1, height_ratios=[0.9, 0.1], num=f"Clone {clone_id}")
        for patch in s.PLANE.get_patches():
            ax.add_patch(patch)

        curr_positions = [[nd.cell.get_current_position() for nd in loc.current_generation] for loc in locations]
        x_gc = [pos.x for pos in curr_positions[0]]
        y_gc = [pos.y for pos in curr_positions[0]]
        x_other = [pos.x for pos in curr_positions[1]]
        y_other = [pos.y for pos in curr_positions[1]]

        data_positions_gc = ax.scatter(x_gc, y_gc, c="r", s=10, edgecolor="black")
        data_positions_other = ax.scatter(x_other, y_other, c="y", s=10, edgecolor="black")
        ax.set_aspect('equal', adjustable='box')

        ax2.set_ylim(-1,1)
        ax2.yaxis.set_visible(False)
        ax2.spines[["left", "top", "right", "bottom"]].set_visible(False)
        ax2.spines[["bottom"]].set_position(("axes", 0.5))
        ax2.set_xlim(-0.5, s.END_TIME-0.5)
        ax2.set_xticks(range(0, s.END_TIME, 10))
        ax2.set_xlabel("Generations")
        ax2.plot([0, s.END_TIME], [0, 0], "-k")
        gen = ax2.plot([time], [0], "-o", color="black", markerfacecolor="white")[0]

        title = ax.set_title(f'Generation {time}')
        def update(i):
            # if time % step != 0:
            #     return data_positions, title, gen,
            # curr_positions = [nd.cell.get_current_position() for loc in locations for nd in loc.current_generation]
            # x = [pos.x for pos in curr_positions]
            # y = [pos.y for pos in curr_positions]
            curr_positions = [[nd.cell.get_current_position() for nd in loc.current_generation] for loc in locations]
            x_gc = [pos.x for pos in curr_positions[0]]
            y_gc = [pos.y for pos in curr_positions[0]]
            x_other = [pos.x for pos in curr_positions[1]]
            y_other = [pos.y for pos in curr_positions[1]]
            data_positions_gc.set_offsets(np.stack([x_gc, y_gc]).T)
            data_positions_other.set_offsets(np.stack([x_other, y_other]).T)

            title.set_text(f'Generation {time}')
            gen.set_data([time], [0])
            fig.canvas.flush_events()
            plt.savefig(f"{simulation_frame_location}/frame{time}.png")
            return data_positions_gc, data_positions_other, title, gen,

        ani = FuncAnimation(fig, update, frames=s.END_TIME, blit=False, event_source=event_handler)
        # plt.ion()
        # plt.legend()
        if not s.DEV:
            fig.canvas.draw_idle()
        # plt.pause(1)
        if s.DEV:
            plt.pause(0.1)
            plt.ion()
            plt.show(block=False)
            plt.ion()

    # TARGET_PAIR.mutate(s.TARGET_MUTATIONS_HEAVY, s.TARGET_MUTATIONS_LIGHT)

    def make_new_child(node):
        child_cell = Cell(
            node.cell.heavy_chain.copy(),
            node.cell.light_chain.copy(),
            location=node.cell.location,
            created_at=time,
            start_point=node.cell.get_current_position(),
            current_polygon=node.cell.current_polygon,
            last_polygon=node.cell.last_polygon
        )
        heavy_n, light_n = child_cell.mutate_cell()
        child_node = Node(
            child_cell,
            parent=node,
            heavy_mutations=heavy_n,
            light_mutations=light_n,
            generation=node.generation+1
            )
        child_cell.calculate_affinity(TARGET_PAIR)
        node.add_child(child_node)
        return child_node


    def make_new_generation(location):
        new_generation = []

        current_generation = location.current_generation
        if len(current_generation) == 0:
            return []
        
        # if s.SELECTION:
        #     non_functional = [x for x in current_generation if x.cell.affinity == 0]
        # if len(non_functional) < len(current_generation):
        #     current_generation = [x for x in current_generation if x not in non_functional]

        if location.name == LocationName.GC:
            to_migrate = do_differentiation(location, time)
        else:
            # potentially allow other locations to migrate in future versions of simble
            to_migrate = do_reentry(location, time)

        current_generation = [x for x in current_generation if x not in to_migrate]

        for node in to_migrate:
            child_node = Node(node.cell.remake_self(), parent=node, generation=node.generation+1)
            node.add_child(child_node)
            if location.name == LocationName.GC:
                OTHER.immigrating_population.append(child_node)
            if location.name == LocationName.OTHER:
                GC.immigrating_population.append(child_node)

        emigration_location = None
        if location.name == LocationName.GC:
            emigration_location = OTHER
        if location.name == LocationName.OTHER:
            emigration_location = GC
        [x.cell.do_random_walk(transition=True) for x in emigration_location.immigrating_population] if emigration_location is not None else None

        if location.name == LocationName.OTHER:
            new_generation = non_gc_population_control(current_generation)
            [x.cell.do_random_walk(transition=False) for x in new_generation]
            return new_generation

        available_antigen = location.settings.max_population

        if s.SELECTION and location.name == LocationName.GC and s.PLANE is None:
            affinities = [x.cell.affinity for x in current_generation]
            p = np.array(affinities) / np.sum(affinities)

        else:
            p = None

        if s.PLANE is None:
            for _ in range(available_antigen):
                current_node = s.RNG.choice(
                    current_generation,
                    p=p
                    )
                current_node.antigen += 1
        else:
            # if time % 10 == 0:
            #     temp_poly = [x.cell.current_polygon.name for x in current_generation if x.cell.current_polygon is not None]
            #     print(Counter(temp_poly), flush=True)
            for poly in s.PLANE.polygons:
                poly_nodes = [x for x in current_generation if x.cell.current_polygon is not None and x.cell.current_polygon.name == poly.name]
                if len(poly_nodes) == 0:
                    # print(f"{poly.name} has no cells")
                    continue
                poly_available_antigen = s.PLANE.antigen_allocation.get(poly.name)
                # print(f"{poly.name} has {len(poly_nodes)} current gen at time {time}, has {poly_available_antigen}", flush=True)
                poly_affinities = [x.cell.affinity for x in poly_nodes]
                p_poly = np.array(poly_affinities)/np.sum(poly_affinities)
                if np.any(np.isnan(p_poly)):
                    logger.warning(f"probabilities include nan values: {p_poly}")
                    p_poly = None
                for _ in range(int(poly_available_antigen)):
                    current_node = s.RNG.choice(
                        poly_nodes,
                        p=p_poly
                        )
                    current_node.antigen += 1


        location.number_of_children = [min(x.antigen, MAX_CHILDREN) for x in current_generation]
        for node in current_generation:
            node.cell.kill_cell()
            children = [make_new_child(node) for _ in range(min(node.antigen, MAX_CHILDREN))]
            live_children = [x for x in children if x.cell.is_alive]
            new_generation.extend(live_children)
            if node.antigen == 0 and (s.MEMORY_SAVE or not s.KEEP_FULL_TREE):
                node.prune_up_tree()

        [x.cell.do_random_walk(transition=False) for x in new_generation]
        return new_generation

    progress_bar = tqdm(
        total=s.END_TIME-1,
        initial=0,
        desc=f"Clone {clone_id}",
        position=clone_id,
        leave=True,
        disable=s.QUIET
        )
    while time<s.END_TIME:
        for location in locations:
            location.current_generation = make_new_generation(location)

        targets = lambda x: (3*x, 3*x+1, 3*x+2) #pylint: disable=unnecessary-lambda-assignment

        row = get_data_points(
            GC.current_generation,
            time,
            naive.heavy_chain.get_gapped_sequence(),
            naive.light_chain.get_gapped_sequence(),
            [i for x in TARGET_PAIR.heavy.mutation_locations for i in targets(x)],
            [i for x in TARGET_PAIR.light.mutation_locations for i in targets(x)])

        dev_data_rows.append(row)

        if time % 25 ==0:
            logger.debug("Time: %d, population: %d", time, len(GC.current_generation))

        [x.finish_migration() for x in locations]

        for location in locations:
            pop_data_rows.append(get_population_data(location, time))
            if time in location.settings.sample_times:
                if time == 0:
                    # make sure we don't remove the naive cell
                    continue
                if time == s.END_TIME-1:
                    sample_size = min(
                        len(location.current_generation),
                        location.settings.sample_size
                        )
                else:
                    sample_size = min(
                        len(location.current_generation)//2,
                        location.settings.sample_size
                        )
                current_sample = s.RNG.choice(
                    location.current_generation,
                    size=sample_size,
                    replace=False
                    )
                location.current_generation = [
                    x
                    for x in location.current_generation
                    if x not in current_sample
                    ]
                for node in current_sample:
                    sampled_ids.append(id(node.cell))
                    node.sampled_time = time
                    sampled.append(node)
                    airr.extend(node.cell.as_AIRR(time))

        time += 1
        if time<s.END_TIME:
            event_handler.update_callbacks()
            progress_bar.update()


    [x.finish_migration() for x in locations]
    for location in locations:
        [node.prune_up_tree() for node in location.current_generation]

    df = pd.DataFrame(dev_data_rows)
    pop_data = pd.DataFrame(pop_data_rows)
    pop_data["clone_id"] = clone_id

    progress_bar.bar_format = "{desc}: |{bar}| {n}/{total} in {elapsed}"
    progress_bar.refresh()
    progress_bar.close()
    event_handler.update_callbacks()

    return sampled, pop_data, df


def run_simulation(i, result_dir):
    """Runs the simulation for a single iteration.

    Args:
        i (int): The iteration number of the simulation.
        result_dir (str): The directory where results will be saved.
    Returns:
        dict: A dictionary containing the results of the simulation, 
            including AIRR data, FASTA sequences, trees, and population data.
    """
    time = 0
    clone_id = i+1
    naive = Cell(None, None, created_at=time)
    root = Node(naive, clone_id=clone_id)
    airr = []
    TARGET_PAIR = TargetAminoPair( # pylint: disable=invalid-name
        naive.heavy_chain.get_gapped_sequence(),
        naive.light_chain.get_gapped_sequence(),
        naive.heavy_chain.cdr3_length,
        naive.light_chain.cdr3_length)
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
        with open(result_dir + "/all_samples.fasta", "w", encoding="utf-8") as f:
            f.write(fasta_string)

        # logger.info("writing newick tree")
        # with open(result_dir + "/true_tree.tree", "w") as f:
        #     f.write(newick)

        logger.info("writing pruned newick tree")
        with open(result_dir + "/pruned_tree.tree", "w", encoding="utf-8") as f:
            f.write(pruned_newick)

        logger.info("writing simplified newick tree")
        with open(result_dir + "/simplified_time_tree.tree", "w", encoding="utf-8") as f:
            f.write(simplified_tree_newick)

        logger.info("writing simplified newick time tree")
        with open(result_dir + "/simplified_time_tree.tree", "w", encoding="utf-8") as f:
            f.write(simplified_time_tree_newick)

    if s.DEV:
        logger.info("max affinity was: %s", TARGET_PAIR.max_affinity)
        logger.info("making plots")
        # make plots
        make_all_plots(dev_df, result_dir)
        make_bar_plot(
            list(TARGET_PAIR.heavy.cdr_multipliers.values()),
            result_dir + "/cdr_multiplier.png",
            "CDR multiplier value",
            "CDR multiplier distribution"
            )
        make_bar_plot(
            list(TARGET_PAIR.heavy.fwr_multipliers.values()),
            result_dir + "/fwr_multiplier.png",
            "FWR multiplier value",
            "FWR multiplier distribution"
            )

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
        "targets": {
            "clone_id": clone_id,
            "heavy": TARGET_PAIR.heavy.amino_acid_seq,
            "light": TARGET_PAIR.light.amino_acid_seq
            }
        }


if __name__ == "__main__":
    run_simulation(1, f'{s.RESULTS_DIR}/results1/')
