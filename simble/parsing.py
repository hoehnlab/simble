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

import argparse
import json
import os

from .location import as_enum
from .settings import s


def get_parser():
    """Creates and returns an argument parser for the simble program.
    Returns:
        argparse.ArgumentParser: The argument parser for the simble program.
    """
    parser = argparse.ArgumentParser(
        prog="simble",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=(
            "somewhat improved model of B-cell lineage evolution \n\n"
            "source code available at: www.github.com/hoehnlab/simble"
            ),
        epilog="It's that simble!"
        )
    program = parser.add_argument_group(
        "program-level settings", 
        description="e.g. output directory or verbosity"
        )

    model = parser.add_argument_group(
        "simulation model settings", 
        description="e.g. mutation probabilities or whether to run neutral simulation"
        )

    sampling = parser.add_argument_group("sampling settings")

    program.add_argument("-v", "--verbose",
                         dest="verbose",
                         help="verbose output",
                         action="store_true")
    program.add_argument("-o", "--output",
                         dest="results",
                         help="absolute path to output directory",
                         metavar="DIR",
                         type=str)
    program.add_argument("-n", "--num",
                         dest="n",
                         help="number of simulations to run",
                         type=int,
                         default=1)
    program.add_argument("-p", "--processes",
                         dest="processes",
                         metavar="P",
                         help="number of processes to run",
                         type=int,
                         default=1)
    program.add_argument("--dev",
                         dest="dev",
                         help="development mode",
                         action="store_true")
    program.add_argument("--fasta",
                         dest="fasta",
                         help="output as fasta",
                         action="store_true")
    program.add_argument("--config",
                         dest="config",
                         help=(
                             "untested! path to config file, will be overwritten "
                             "by command line arguments"
                             ),
                         metavar="FILE",
                         type=str,
                         default=None)
    program.add_argument("--seed",
                         dest="seed",
                         help="seed to seed the random number generators",
                         type=int,
                         default=None)
    program.add_argument("--memory-save",
                         dest="memory_save",
                         help="save memory by only using simplified trees",
                         action="store_true")
    program.add_argument("--full-tree",
                         dest="keep_full_tree",
                         help="keep the full trees, including non-sampled tips and their ancestors",
                         action="store_true")
    program.add_argument("-q", "--quiet",
                         dest="quiet",
                         help="don't display progress bar",
                         action="store_true")

    sampling.add_argument("-s", "--samples",
                          dest="sample_info",
                          metavar=("start", "stop", "step"),
                          help="specify sample times other than the default",
                          nargs=3,
                          default=None,
                          type=int)
    sampling.add_argument("--other-samples",
                          dest="other_sample_info",
                          metavar=("start", "stop", "step"),
                          help="specify sample times for only the 'Other' location",
                          nargs=3,
                          default=None)
    sampling.add_argument("--sample-size",
                          dest="sample_size",
                          metavar="N",
                          help="specify sample size for 'GC' location",
                          default=None,
                          type=int)
    sampling.add_argument("--sample-size-other",
                         dest="sample_size_other",
                         metavar="N",
                         help="specify sample size for the 'Other' location",
                         default=None,
                         type=int)

    model.add_argument("--neutral",
                       dest="neutral",
                       help="neutral simulation (no selection in germinal center)",
                       action="store_true")
    model.add_argument("--uniform",
                       dest="uniform",
                       help="use a uniform mutation and substitution model",
                       action="store_true")
    model.add_argument("-a", "--antigen",
                       dest="antigen",
                       help="amount of antigen",
                       metavar="A",
                       default=None,
                       type=int)
    model.add_argument("-m", "--multiplier",
                       dest="multiplier",
                       help="selection multiplier",
                       metavar="M",
                       default=None,
                       type=int)
    model.add_argument("--migration-rate",
                       dest="migration_rate",
                       help=(
                           "migration rate from the GC (expected value of cells "
                           "that migrate per generation)"
                           ),
                       metavar="R",
                       default=None,
                       type=float)
    model.add_argument("--heavy-mutate-probability",
                       dest="heavy_mutate_probability",
                       help="probability of heavy chain mutation",
                       metavar="P",
                       default=None,
                       type=float)
    model.add_argument("--light-mutate-probability",
                       dest="light_mutate_probability",
                       help="probability of light chain mutation",
                       metavar="P",
                       default=None,
                       type=float)
    model.add_argument("--target-mutations-heavy",
                       dest="target_mutations_heavy",
                       help="number of mutations in heavy chain target",
                       metavar="N",
                       default=None,
                       type=int)
    model.add_argument("--target-mutations-light",
                       dest="target_mutations_light",
                       help="number of mutations in light chain target",
                       metavar="N",
                       default=None,
                       type=int)
    model.add_argument("--cdr-dist",
                       dest="cdr_dist",
                       help="cdr distribution",
                       default=None,
                       choices=["constant", "exponential"],
                       type=str)
    model.add_argument("--cdr-var",
                       dest="cdr_var",
                       help="cdr variable",
                       metavar="V",
                       default=None,
                       type=float)
    model.add_argument("--fwr-dist",
                       dest="fwr_dist",
                       help="fwr distribution",
                       default=None,
                       choices=["constant", "exponential"],
                       type=str)
    model.add_argument("--fwr-var",
                       dest="fwr_var",
                       help="fwr variable",
                       metavar="V",
                       default=None,
                       type=float)

    return parser

def _update_setting(name, value):
    """Updates a setting in the global settings object, leaving it unchanged 
    if the value is None."""
    if value is not None:
        setattr(s, name, value)

def validate_samples(sample_info):
    """Validates the sampling settings.
    Args:
        sample_info (list): A list containing start, stop, and step values.
    Raises:
        ValueError: If the sampling setting is invalid.
    """
    if sample_info[0] > sample_info[1]:
        raise ValueError("sample start must be less than or equal to stop")
    if sample_info[2] <= 0:
        raise ValueError("sample step must be greater than 0")
    if sample_info[2] > sample_info[1] - sample_info[0]:
        raise ValueError("sample step must be less than or equal to stop - start")

def validate_and_process_args(args):
    """Validates and processes command line arguments and updates the simulation settings.
    Args:
        args (argparse.Namespace): The parsed command line arguments.
    Returns:
        list: A list of warnings, if any.
    """
    warnings = []
    if args.config is not None:
        # TODO (jf): validate that file exists
        warnings.append("Using a config file is not fully tested. Use at your own risk!")
        data = read_from_json(args.config)
        s.update_from_dict(data)

    if args.results is not None:
        s.RESULTS_DIR = args.results
    elif s.RESULTS_DIR == "":
        s.RESULTS_DIR = os.getcwd() + "/results"

    if not os.path.exists(s.RESULTS_DIR):
        os.mkdir(s.RESULTS_DIR)

    if args.neutral:
        if args.multiplier is not None:
            warnings.append("Neutral simulation specified, ignoring selection multiplier")
        s.SELECTION = False

    # if uniform selection is specified, ignore selection
    if args.uniform and s.SELECTION:
        warnings.append("Uniform mutation and substitution model specified, ignoring selection")
        s.SELECTION = False
        _update_setting("UNIFORM", args.uniform)

    if args.sample_info:
        validate_samples(args.sample_info)
        start = args.sample_info[0]
        stop = args.sample_info[1]
        step = args.sample_info[2]
        for location in s.LOCATIONS:
            location.sample_times = list(range(start, stop+1, step))
    if args.other_sample_info is not None:
        validate_samples(args.other_sample_info)
        start = args.other_sample_info[0]
        stop = args.other_sample_info[1]
        step = args.other_sample_info[2]
        s.LOCATIONS[1].sample_times = list(range(start, stop+1, step))

    if args.migration_rate:
        s.LOCATIONS[0].migration_rate = args.migration_rate

    if args.memory_save and args.keep_full_tree:
        warnings.append("Memory save and full tree options are incompatible. Using memory save.")
        args.keep_full_tree = False

    _update_setting("MULTIPLIER", args.multiplier)
    _update_setting("HEAVY_MUTATE_PROBABILITY", args.heavy_mutate_probability)
    _update_setting("LIGHT_MUTATE_PROBABILITY", args.light_mutate_probability)
    _update_setting("TARGET_MUTATIONS_HEAVY", args.target_mutations_heavy)
    _update_setting("TARGET_MUTATIONS_LIGHT", args.target_mutations_light)
    _update_setting("DEV", args.dev)
    _update_setting("FASTA", args.fasta)
    _update_setting("CDR_DIST", args.cdr_dist)
    _update_setting("CDR_VAR", args.cdr_var)
    _update_setting("FWR_DIST", args.fwr_dist)
    _update_setting("FWR_VAR", args.fwr_var)
    _update_setting("MAX_POPULATION", args.antigen)
    _update_setting("MEMORY_SAVE", args.memory_save)
    _update_setting("KEEP_FULL_TREE", args.keep_full_tree)
    _update_setting("QUIET", args.quiet)

    if args.sample_size:
        s.LOCATIONS[0].sample_size = args.sample_size
    if args.sample_size_other:
        s.LOCATIONS[1].sample_size = args.sample_size_other

    if s.LOCATIONS[1].sample_times is None:
        # if no sample times are specified for the "Other" location, use the same as the GC
        s.LOCATIONS[1].sample_times = s.LOCATIONS[0].sample_times

    return warnings


def validate_location(location):
    """Validates a location dictionary.
    Args:
        location (dict): A dictionary representing a location.
    Raises:
        ValueError: If the location dictionary contains invalid fields or types.
    """
    valid_fields = {x: type(y) for x, y in vars(s.LOCATIONS[0]).items() if not x.startswith("_")}
    for key, value in location.items():
        if key not in valid_fields:
            raise ValueError(f"invalid LOCATION field: {key}")
        if not isinstance(value, valid_fields[key]):
            raise ValueError(f"invalid type for LOCATION field {key}: {type(value)}")


def validate_json(json_input):
    """Validates the JSON input against the global settings object.
    Args:
        json_input (dict): The JSON input to validate.
    Raises:
        ValueError: If the JSON input contains invalid fields or types.
    """
    valid_fields = {x: type(y) for x, y in vars(s).items() if not x.startswith("_")}
    for key, value in json_input.items():
        if key == "LOCATIONS":
            for location in value:
                validate_location(location)
        if key not in valid_fields:
            raise ValueError(f"invalid field: {key}")
        valid_type = valid_fields[key]
        if valid_type == int:
            # float is fine if it's actually an integer
            if isinstance(value, float) and value.is_integer():
                continue
            if not isinstance(value, valid_type):
                raise ValueError(f"invalid type for field {key}: {type(value)}")
        elif valid_type == float:
            # integer is fine for floats
            if not isinstance(value, valid_type) and not isinstance(value, int):
                raise ValueError(f"invalid type for field {key}: {type(value)}")
        elif not isinstance(value, valid_type):
            raise ValueError(f"invalid type for field {key}: {type(value)}")
        if key == "CDR_DIST" or key == "FWR_DIST":
            if value not in ["constant", "exponential"]:
                raise ValueError(f"invalid value for field {key}: {value}")


def read_from_json(filename):
    """Reads a JSON file and returns its contents as a dictionary.
    Args:
        filename (str): The path to the JSON file.
    Returns:
        dict: The contents of the JSON file as a dictionary.
    """
    with open(filename, "r", encoding="utf-8") as f:
        data = json.load(f, object_hook=as_enum)
    validate_json(data)
    return data
