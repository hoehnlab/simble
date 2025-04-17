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
from tqdm import tqdm

from .helper import (ALL_TREE_NAMES, MEMORY_SAVE_TREE_NAMES, TREE_NAMES,
                     make_all_plots)
from .location import as_enum
from .parsing import get_parser, validate_and_process_args
from .settings import s
from .simulation import run_simulation

logger = logging.getLogger(__package__)

class TqdmLoggingHandler(logging.Handler):
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
    airr["d_germline_start"] = airr["d_germline_start"].astype(pd.Int64Dtype())
    airr["d_germline_end"] = airr["d_germline_end"].astype(pd.Int64Dtype())
    airr.to_csv(s.RESULTS_DIR + "/all_samples_airr.tsv", sep="\t", index=False)
    pop_data = pd.concat(all_results["pop_data"])
    pop_data.to_csv(s.RESULTS_DIR + "/population_data.csv", index=False)
    if s.DEV:
        df = pd.concat(all_results["data"])
        df.to_csv(s.RESULTS_DIR + "/all_data.csv")
        grouped = df.groupby(['time']).mean().reset_index()
        make_all_plots(grouped, s.RESULTS_DIR, True)

    tree_names = MEMORY_SAVE_TREE_NAMES if s.MEMORY_SAVE else ALL_TREE_NAMES if s.KEEP_FULL_TREE else TREE_NAMES
    nexus = ["#NEXUS\n" + "BEGIN TREES;\n" for _ in tree_names]

    for clone in results:
        for i, tree_name in enumerate(tree_names):
            nexus[i] += f'\tTree {clone["clone_id"]} = {clone[tree_name]}\n'

    for i, tree_name in enumerate(tree_names):
        nexus[i] += "END;\n"
        with open(s.RESULTS_DIR + f"/all_{tree_name}s.nex", "w") as f:
            f.write(nexus[i])
    
    targets = pd.DataFrame(all_results["targets"])
    targets.to_csv(s.RESULTS_DIR + "/all_targets.csv", index=False)


def main():
    parser = get_parser()

    args = parser.parse_args()
    warnings = validate_and_process_args(args)

    set_logger()
    for warning in warnings:
        logger.warning(warning)

    if args.seed is not None:
        # TODO: fix input for seeds
        seed = args.seed
        ss = np.random.SeedSequence(seed)
    else:
        ss = np.random.SeedSequence()
    seeds = ss.spawn(args.n)
    print(f"Seed: {ss.entropy}")

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