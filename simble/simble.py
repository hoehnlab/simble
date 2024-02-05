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
import tempfile
import time
from functools import partial
from multiprocessing import Pool

import numpy as np
import pandas as pd

from .helper import make_all_plots
from .location import as_enum
from .parsing import get_parser, validate_and_process_args
from .settings import s
from .simulation import run_simulation

logger = logging.getLogger(__package__)

def set_logger():
    if s.DEV:
        logger.setLevel(logging.DEBUG)
        log_format = '%(asctime)s %(process)d \t%(levelname)s: %(message)s'
    elif s.VERBOSE:
        logger.setLevel(logging.INFO)
        log_format = '%(asctime)s %(levelname)s: %(message)s'
    else:
        logger.setLevel(logging.WARNING)
        log_format = '%(levelname)s: %(message)s'

    logging.basicConfig(format=log_format, datefmt='%H:%M:%S')


def do_simulation(i, seed, filename):
    with open(filename, "r") as f:
        settings = json.load(f, object_hook=as_enum)
    s.update_from_dict(settings)
    s._RNG = np.random.default_rng(seed)
    set_logger()
    logger.info(f"Starting simulation {i}")
    folder = s.RESULTS_DIR
    curr_results = f'{folder}/results{i}/'
    if s.DEV and not os.path.exists(curr_results):
        os.mkdir(curr_results)
    start = time.time()
    data = run_simulation(i, curr_results)
    end = time.time()

    logger.debug(f"Time taken: {end-start}")
    return data


def process_results(results):
    all_results = {}
    for key in results[0].keys():
        all_results[key] = [x[key] for x in results]

    if s.FASTA:
        fasta_string = "\n".join(all_results["fasta"])
        with open(s.RESULTS_DIR + "/all_samples.fasta", "w") as f:
            f.write(fasta_string)
    airr = pd.concat(all_results["airr"])
    airr.to_csv(s.RESULTS_DIR + "/all_samples_airr.tsv", sep="\t", index=False)
    pop_data = pd.concat(all_results["pop_data"])
    pop_data.to_csv(s.RESULTS_DIR + "/population_data.csv", index=False)
    if s.DEV:
        df = pd.concat(all_results["data"])
        df.to_csv(s.RESULTS_DIR + "/all_data.csv")
        grouped = df.groupby(['time']).mean().reset_index()
        make_all_plots(grouped, s.RESULTS_DIR, True)

    nexus = "#NEXUS\n" + "BEGIN TREES;\n"
    for clone in results:
        nexus += f'\tTree true_tree_{clone["clone_id"]} = {clone["true_tree"]}\n'
        nexus += f'\tTree pruned_tree_{clone["clone_id"]} = {clone["pruned_tree"]}\n'
        nexus += f'\tTree pruned_time_tree_{clone["clone_id"]} = {clone["pruned_time_tree"]}\n'
    nexus += "END;\n"
    with open(s.RESULTS_DIR + "/all_trees.nex", "w") as f:
        f.write(nexus)


def main():
    parser = get_parser()

    args = parser.parse_args()
    warnings = validate_and_process_args(args)

    set_logger()

    for warning in warnings:
        logger.warning(warning)

    if args.seeds is not None:
        seeds = [np.random.SeedSequence(seed) for seed in args.seeds]
    else:
        ss = np.random.SeedSequence()
        seeds = ss.spawn(args.n)
    print(f"Seeds: {[ss.entropy for ss in seeds]}")

    with tempfile.NamedTemporaryFile(mode="w") as tmpf:
        json.dump(s, tmpf, default=lambda o: o.encode(), indent=4)
        tmpf.flush()
        start = time.time()
        logger.info(f"Starting simulation")
        if args.processes > 1:
            with Pool(processes=args.processes) as pool:
                result = pool.starmap(partial(do_simulation, filename=tmpf.name), zip(range(args.n), seeds))

        else:
            result = []
            for i in range(args.n):
                result.append(do_simulation(i, seeds[i], tmpf.name))

    process_results(result)

    end = time.time()
    logger.debug(f"Program finished! Total time taken: {end-start}")