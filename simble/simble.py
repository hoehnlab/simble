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

import json
import logging
import os
# CGJ 8/4/26
import sys
import tempfile
import time
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd
from tqdm import tqdm

# CGJ 
from .helper import (ALL_TREE_NAMES, MEMORY_SAVE_TREE_NAMES, TREE_NAMES,
                     get_unique_founder_indices, make_all_plots, update_helper_tables)
from .location import as_enum
from .parsing import get_parser, validate_and_process_args
from .settings import s
from .simulation import run_simulation
# CGJ
from .target import TargetAminoPair

logger = logging.getLogger(__package__)
# CGJ 8/4/26 
RECURSION_LIMIT_MULTIPLIER = 1.5

class TqdmLoggingHandler(logging.Handler):
    """Custom logging handler to write logs to tqdm output."""
    def __init__(self, level=logging.NOTSET):
        super().__init__(level)

    def emit(self, record):
        try:
            msg = self.format(record)
            tqdm.write(f"\033[K{msg}")
            self.flush()
        except Exception:
            self.handleError(record)

def set_logger():
    """Sets up the logger for the simulation."""
    if s.DEV:
        logger.setLevel(logging.DEBUG)
        log_format = '%(asctime)s %(process)d \t%(levelname)s: %(message)s'
    elif s.VERBOSE:
        logger.setLevel(logging.INFO)
        log_format = '%(asctime)s %(levelname)s: %(message)s'
    else:
        logger.setLevel(logging.WARNING)
        log_format = '%(levelname)s: %(message)s'

    handler = TqdmLoggingHandler()
    formatter=logging.Formatter(log_format, datefmt='%H:%M:%S')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.propagate = False


# CGJ
def build_target_pair(seed):
    """Builds the target that every clone in the run is scored against.

    Args:
        seed (np.random.SeedSequence): The seed for the multiplier distributions.
    Returns:
        TargetAminoPair: The target built from the supplied target sequence.
    """
    previous_rng = s._x_RNG
    s._x_RNG = np.random.default_rng(seed)
    try:
        return TargetAminoPair(
            s.TARGET["heavy_aligned"],
            s.TARGET["light_aligned"],
            s.TARGET["heavy_cdr3_aa_length"],
            s.TARGET["light_cdr3_aa_length"])
    finally:
        # the settings are serialized to the workers, and the generator is not
        # serializable, so leave it as it was found
        s._x_RNG = previous_rng


# CGJ
def do_simulation(clone_id, seed, founder_idx, filename, target=None):
    """Runs a single simulation with the given seed and settings."""
    with open(filename, "r", encoding="utf-8") as f:
        settings = json.load(f, object_hook=as_enum)
    s.update_from_dict(settings)
    # CGJ
    # workers do not inherit the tables built when the package was imported
    update_helper_tables()
    # CGJ 8/4/26
    # Each multiprocessing worker needs the recursion limit to be raised
    sys.setrecursionlimit(max(1000, int(s.END_TIME * RECURSION_LIMIT_MULTIPLIER)))
    s._x_RNG = np.random.default_rng(seed) # pylint: disable=protected-access
    set_logger()
    logger.info(f"Starting simulation for clone {clone_id}")
    folder = s.RESULTS_DIR
    curr_results = f'{folder}/results{clone_id}/'
    if s.DEV and not os.path.exists(curr_results):
        os.mkdir(curr_results)
    start = time.time()
    # CGJ
    data = run_simulation(clone_id, curr_results, founder_idx, target)
    end = time.time()

    logger.debug("Time taken: %s", end - start)
    return data


def process_results(results):
    """Processes the results of the simulations and saves them to files."""
    all_results = {}
    for key in results[0].keys():
        all_results[key] = [x[key] for x in results]

    if s.FASTA:
        fasta_string = "\n".join(all_results["fasta"])
        with open(s.RESULTS_DIR + "/all_samples.fasta", "w", encoding="utf-8") as f:
            f.write(fasta_string)
    airr = pd.concat(all_results["airr"])
    if "d_germline_start" in airr.columns:
        airr["d_germline_start"] = airr["d_germline_start"].astype(pd.Int64Dtype())
    if "d_germline_end" in airr.columns:
        airr["d_germline_end"] = airr["d_germline_end"].astype(pd.Int64Dtype())
    airr.to_csv(s.RESULTS_DIR + "/all_samples_airr.tsv", sep="\t", index=False)
    pop_data = pd.concat(all_results["pop_data"])
    pop_data.to_csv(s.RESULTS_DIR + "/population_data.csv", index=False)
    if s.DEV:
        df = pd.concat(all_results["data"])
        df.to_csv(s.RESULTS_DIR + "/all_data.csv")
        grouped = df.groupby(['time']).mean().reset_index()
        make_all_plots(grouped, s.RESULTS_DIR, True)

    if s.MEMORY_SAVE:
        tree_names = MEMORY_SAVE_TREE_NAMES
    elif s.KEEP_FULL_TREE:
        tree_names = ALL_TREE_NAMES
    else:
        tree_names = TREE_NAMES
    nexus = ["#NEXUS\n" + "BEGIN TREES;\n" for _ in tree_names]

    for clone in results:
        for i, tree_name in enumerate(tree_names):
            nexus[i] += f'\tTree {clone["clone_id"]} = {clone[tree_name]}\n'

    for i, tree_name in enumerate(tree_names):
        nexus[i] += "END;\n"
        with open(s.RESULTS_DIR + f"/all_{tree_name}s.nex", "w", encoding="utf-8") as f:
            f.write(nexus[i])

    targets = pd.DataFrame(all_results["targets"])
    targets.to_csv(s.RESULTS_DIR + "/all_targets.csv", index=False)


def main():
    """ Main function to run the simulation. """
    parser = get_parser()

    args = parser.parse_args()
    try:
        warnings = validate_and_process_args(args)
    except Exception as e:
        raise SystemExit(e)
    update_helper_tables()

    set_logger()
    for warning in warnings:
        logger.warning(warning)

    if args.seed is not None:
        seed = args.seed
        ss = np.random.SeedSequence(seed)
    else:
        ss = np.random.SeedSequence()
    seeds = ss.spawn(args.n)
    print(f"Seed: {ss.entropy}")

    # CGJ
    if args.unique_founders and not s.UNIFORM:
        founder_rng = np.random.default_rng(ss.spawn(1)[0])
        founder_indices = get_unique_founder_indices(args.n, founder_rng)
    else:
        founder_indices = [None] * args.n

    clone = args.clone
    # CGJ
    target = build_target_pair(ss.spawn(1)[0]) if s.TARGET else None

    with tempfile.NamedTemporaryFile(mode="w") as tmpf:
        json.dump(s, tmpf, default=lambda o: o.encode(), indent=4)
        tmpf.flush()
        start = time.time()
        logger.info("Starting simulation")
        if args.processes > 1:
            with Pool(processes=args.processes) as pool:
                result = pool.starmap(
                    # CGJ
                    partial(do_simulation, filename=tmpf.name, target=target),
                    zip(range(clone, clone + args.n), seeds, founder_indices)
                    )
        else:
            result = []
            for i in range(args.n):
                result.append(
                    # CGJ
                    do_simulation(clone + i, seeds[i], founder_indices[i], tmpf.name, target)
                    )

    process_results(result)

    end = time.time()
    logger.debug("Program finished! Total time taken: %s", end-start)
