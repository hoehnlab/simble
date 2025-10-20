Using SimBLE
=============

.. contents::
    :local:
    :depth: 1

.. toctree::
    :hidden:

    all_arguments


Using SimBLE from the command line
---------------------------

SimBLE is designed to be run out-of-the-box from the command line.

.. include:: basic_command_line_usage.rst

For more information on command line arguments, please see :ref:`all arguments<all_arguments>`.

Using SimBLE as a Python library
---------------------------

SimBLE can also be used as a Python library for extending functionality. Below is a simple example of how to run a neutral and a selection simulation for the same naive cell from Python:

.. code-block:: python

    import os

    import numpy as np
    import pandas as pd

    from simble.cell import Cell
    from simble.parsing import get_parser, validate_and_process_args
    from simble.settings import s
    from simble.simble import logger, process_results, set_logger
    from simble.simulation import run_simulation, simulate
    from simble.target import TargetAminoPair
    from simble.tree import Node, simplify_tree


    def do_selection_simulation(i, naive):
        # set the results directory to be a subfolder for selection results
        main_dir = s.RESULTS_DIR
        s.RESULTS_DIR = s.RESULTS_DIR + "/selection/"
        if not os.path.exists(s.RESULTS_DIR):
            os.mkdir(s.RESULTS_DIR)

        # do the simulation
        results = do_simulation(i, naive)
        process_results([results])

        # reset results directory
        s.RESULTS_DIR = main_dir

    def do_neutral_simulation(i, naive):
        # set the results directory to be a subfolder for neutral results
        main_dir = s.RESULTS_DIR
        s.RESULTS_DIR = s.RESULTS_DIR + "/neutral/"
        # set selection to false so this is neutral
        s.SELECTION = False
        if not os.path.exists(s.RESULTS_DIR):
            os.mkdir(s.RESULTS_DIR)

        # do the simulation
        results = do_simulation(i, naive)
        process_results([results])

        # reset results directory
        s.RESULTS_DIR = main_dir


    def do_simulation(i, naive):
        clone_id = i+1 # 1-indexed clone id rather than 0-indexed
        # set up the root node and target sequences
        root = Node(naive, clone_id=clone_id)
        TARGET_PAIR = TargetAminoPair(
            naive.heavy_chain.get_gapped_sequence(),
            naive.light_chain.get_gapped_sequence(),
            naive.heavy_chain.cdr3_length,
            naive.light_chain.cdr3_length)
        TARGET_PAIR.mutate(s.TARGET_MUTATIONS_HEAVY, s.TARGET_MUTATIONS_LIGHT)

        # run the simulation with our starting root
        sampled, pop_data, df = simulate(clone_id, TARGET_PAIR, [root], root)

        # process the results
        sampled_ids = [id(x.cell) for x in sampled]
        fasta_string = "".join([x.cell.as_fasta(x.sampled_time) for x in sampled])

        airr = [x for node in sampled for x in node.cell.as_AIRR(node.sampled_time)]
        airr = pd.DataFrame(airr)
        airr["sequence_id"] = airr["sequence_id"].apply(lambda x: f"{clone_id}_{x}")
        airr["cell_id"] = airr["cell_id"].apply(lambda x: f"{clone_id}_{x}")
        airr["clone_id"] = clone_id

        newick = ""
        pruned = root.prune_subtree(sampled_ids)
        pruned_newick = f'{pruned.write_newick()};'
        pruned_time_tree = f'{pruned.write_newick(time_tree=True)};'

        simplified_tree = simplify_tree(pruned)
        simplified_tree_newick = f'{simplified_tree.write_newick()};'
        simplified_time_tree_newick = f'{simplified_tree.write_newick(time_tree=True)};'

        return {
            "airr": airr, 
            "fasta": fasta_string, 
            "full_tree": newick, 
            "pruned_tree": pruned_newick,
            "pruned_time_tree": pruned_time_tree, 
            "simplified_tree": simplified_tree_newick,
            "simplified_time_tree": simplified_time_tree_newick,
            "data": df, 
            "clone_id": clone_id, 
            "pop_data": pop_data,
            "targets": {
                "clone_id": clone_id,
                "heavy": TARGET_PAIR.heavy.amino_acid_seq,
                "light": TARGET_PAIR.light.amino_acid_seq
                }
            }


    def main():

        # use simble's argument parsing and validation for ease of setting arguments
        parser = get_parser()
        args = parser.parse_args()
        warnings = validate_and_process_args(args)

        set_logger()

        for warning in warnings:
            logger.warning(warning)

        if args.seed is not None:
            seed = args.seed
            ss = np.random.SeedSequence(seed)
        else:
            ss = np.random.SeedSequence()
        print(f"Seed: {ss.entropy}")

        # enforce that we start with selection simulation
        s.UNIFORM = False
        s.SELECTION = True

        s._x_RNG = np.random.default_rng(seed) # set RNG from seed

        if not os.path.exists(s.RESULTS_DIR):
            os.mkdir(s.RESULTS_DIR)

        # make a random naive cell to use for both simulations
        naive = Cell(None, None, created_at=0)

        logger.info("Starting simulations")
        do_selection_simulation(0, naive)
        do_neutral_simulation(0, naive)


    if __name__ == "__main__":
        main()


Some more advanced examples will be available soon.

..  in the `examples folder <https://github.com/hoehnlab/simble/example-extensions>`_. These include adoptive transfer experiments and simulations designed to test cross-reactivity.

