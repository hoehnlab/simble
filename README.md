# simble

simble is a BCR evolution simulator. It starts with a naive heavy and light chain pair,
and models evolution and selection in the germinal center as well as migration out of 
germinal center. simble also has customizable sampling.

:globe_with_meridians: [simble website](https://simble.readthedocs.io)

To cite the SimBLE package in publications, please use:

Fielding J, Wu S, Melton H, Du Plessis L, Fisk N, Hoehn K (2025). "TyCHE enables time-resolved lineage tracing of heterogeneously-evolving populations." *bioRxiv 2025.10.21.683591* doi:10.1101/2025.10.21.683591 https://doi.org/10.1101/2025.10.21.683591

# Documentation

For full, up-to-date documentation, information, and tutorials, please visit [simble.readthedocs.io](https://simble.readthedocs.io).


## Installing simble

<!-- The easiest way to install simble is with pip:

```sh
pip install simble
``` -->

Clone this repo and install from local files:

```sh
# clone the repo however you like, eg
git clone https://github.com/hoehnlab/simble

# in the repo directory, install simble from local files
pip install .

# run
python3 -m simble -flags
# or
simble -flags
```


## Quick Start

>[!NOTE]
>The default simulation runs with selection, with no migration, sampling every 25 generations for 200 generations. 

To specify an output folder for default simulation:
```sh
simble -o <path-to-folder>
```

To run a neutral-selection BCR simulation (using naive BCRs, heavy and light chains, and S5F mutation/substitution model):
```sh
simble --neutral [other args]
```

To run a single chain simulation (only the heavy chain of each naive pair, with the usual S5F mutation/substitution model and selection):
```sh
simble --single-chain [other args]
```
In a single chain simulation, the output has only heavy chain rows, and
`--light-shm` and `--target-mutations-light` are ignored. If you give a
`--target`, only its heavy chain is used, and only the heavy chain has to match
for a naive pair to be compatible with the target. The naive and target inputs
can be paired tables or heavy chain only tables (see
[Using your own naive sequences](#using-your-own-naive-sequences)).

To run a uniformly neutral simulation (no selection, randomly generated starting nucleotide sequence, and uniform mutations/substitutions):
```sh
simble --uniform [other args]
```
To run a uniformly neutral simulation with a specified sequence length of 100:
```sh
simble --uniform --sequence-length 100 [other args]
```

To run with expected migration of one cell every 25 generations:
```sh
simble --migration-rate 0.04 [other args]
```

To run 5 clones in parallel across 2 processes, with expected migration of one cell every 10 generations
with selection, and sampling every 10 generations for 100 generations:
```sh
simble -o ./current-results -n 5 --processes 2 --migration-rate 0.1 --samples 0 100 10
```

which is equivalent to
```sh
simble -o ./current-results -n 5 -p 2 --migration-rate 0.1 -s 0 100 10
```

>[!TIP]
>Flags can be provided any order.

Frequently used arguments:

| argument | abbr | default | description |
| -------- | ----------   | ------- |   -------   |
| --output | -o           | cwd/results| folder for results |
| --number | -n | 1 | number of clones to simulate |
| --processes| -p | 1 | number of processes (multiprocessing) |
| --neutral | | | if provided, runs a neutral simulation|
| --uniform | | | if provided, runs a uniform neutral simulation |
| --migration-rate| | 0 | expected number of cells that leave the germinal center each generation|
| --samples | -s | [0 200 25] | start, stop, step, to specify sample times other than the default|
| --quiet | -q | | don't display progress bar


## Development

Clone the repo and install necessary packages, which can be found in pyproject.toml:

```sh
# clone the repo however you like, eg
git clone https://github.com/hoehnlab/simble

# install requirements, up to date requirements can be found in pyproject.toml
pip install <requirement>

# run
python3 -m simble -flags
```

## Using your own naive sequences

To use your own naive data, you can use the `--naive` argument:

```sh
simble --naive <naive_file.csv> 
```

Additionally, you may wish to keep columns from your naive dataset:

```sh
simble --naive <naive_file.csv> --keep-cols cell_barcode cell_label
```

These columns may not have the same name as any required fields, and additionally
may not be called any of the following:
-sequence_id
-sequence
-sequence_alignment
-germline_alignment
-location
-sample_time
-junction
-junction_aa
-junction_length

The `--naive-sampling` argument controls how each clone's starting pair is
chosen from the naive data:

- `ordered` (the default whenever `--naive` is given): each clone uses the
  naive pair at row clone_id modulo (rows in naive dataset). For example, if
  the naive dataset has 10 entries, clone 6 will use the pair of sequences in
  the 6th row, and clone 13 will use the 3rd row.
- `random`: each clone draws an independent, random starting pair from the
  table. The same pair can be reused by more than one clone. This is always
  the mode used when `--naive` is not given.
- `unique`: one starting pair is drawn per clone, without replacement, so no
  two clones in the run share a naive pair. This requires at least as many
  naive pairs as clones (`--number`/`-n`).

```sh
simble --naive <naive_file.csv> --naive-sampling random
simble --naive <naive_file.csv> --naive-sampling unique
```

If you would like to specify a clone_id to start at:
```sh
simble --naive <naive_file.csv> --clone-id 4
```


Your naive data should have the following AIRR columns for both heavy and 
light chains (prefixed with 'heavy_' and 'light_' respectively):

-sequence
-rev_comp
-productive
-v_call
-d_call
-j_call
-sequence_alignment
-junction
-junction_aa
-v_cigar
-d_cigar
-j_cigar
-np1_length
-v_germline_start
-v_germline_end
-d_germline_start
-d_germline_end
-j_germline_start
-j_germline_end
-germline_alignment_d_mask
-locus
-cdr3

For a single chain simulation (`--single-chain`), the naive data can have only
the heavy chain columns (prefixed with 'heavy_'). A paired table also works;
its light chain columns are not used.

For convenience, once `simble` is installed, you can also run the command 
`process_naive` to compile the correct csv from a heavy AIRR tsv and a light AIRR tsv:

```sh
process_naive -o <outfile> --heavy <heavy_airr.tsv> --light <light_airr.tsv> --join <column_to_join_on>
```

For example, to join on "cell_id":
```sh
process_naive -o my_naive_data.csv --heavy IGH_naive_airr.tsv --light IGL_IGK_naive_airr.tsv -j cell_id
```

To keep columns "cell_barcode" and "cell_label" from the heavy table:
```sh
process_naive -o my_naive_data.csv --heavy IGH_naive_airr.tsv --light IGL_IGK_naive_airr.tsv -j cell_id --keep-cols cell_barcode cell_label
```

Or to keep columns "cell_barcode" and "cell_label" from the light table:
```sh
process_naive -o my_naive_data.csv --heavy IGH_naive_airr.tsv --light IGL_IGK_naive_airr.tsv -j cell_id --keep-cols cell_barcode cell_label --keep-from light
```

To make a heavy chain only table for a single chain simulation, do not give
`--light`. The join column is kept if it is in the heavy table, but it is not
required:
```sh
process_naive -o my_heavy_naive_data.csv --heavy IGH_naive_airr.tsv
simble --single-chain --naive my_heavy_naive_data.csv [other args]
```

In a single chain simulation, a target can also be given as only a heavy chain
AIRR tsv, with `--target-heavy` and no `--target-light`:
```sh
simble --single-chain --naive my_heavy_naive_data.csv --target-heavy IGH_target_airr.tsv
```

## Using your own 5-mer mutation model

By default, simble uses the S5F mutability and substitution models (see
`simble/data/README.md`). To use your own 5-mer model for both chains:

```sh
simble --mutability my_mutability.csv --substitution my_substitution.csv [other args]
```

To give a model for one chain only, use `--mutability-heavy`,
`--substitution-heavy`, `--mutability-light` or `--substitution-light`. A file
for one chain is used instead of `--mutability` or `--substitution` for that
chain. A chain (or a file) that you do not give uses the bundled S5F file, so
you can, for example, give only a mutability file. In a single chain
simulation, the light chain files are not used. Model files cannot be used with
`--uniform`.

The files have the same format as the bundled files (comma separated, or tab
separated if the name ends in `.tsv`):

- mutability file: columns `Fivemer` and `Mutability`. It must have all 1024
  5-mers of A, C, G and T. Only the relative values matter: the number of
  mutations comes from `--heavy-shm` and `--light-shm`.
- substitution file: columns `Fivemer`, `A`, `C`, `G` and `T`. Each row gives
  the probability that the center nucleotide of the 5-mer becomes each
  nucleotide. For each 5-mer with a mutability greater than 0, the row must sum
  to 1, and the probability for the center nucleotide to stay the same must be
  0 (an empty value is 0). Rows for 5-mers with a mutability of 0 are not used,
  so they can be left out.

Rows for 5-mers with N (used at the ends of a sequence) are optional. If they
are left out, simble makes them the same way as in the bundled files: a 5-mer
with N at an edge position gets the mean of the 5-mers that replace each N with
A, C, G or T. simble checks the files before the simulation starts, and stops
with an error that names the bad 5-mers.

To export a model from [shazam](https://shazam.readthedocs.io/en/stable/) in R
(`model` is a `TargetingModel`, for example from `createTargetingModel`):

```r
write.csv(data.frame(Fivemer = names(model@mutability), Mutability = model@mutability),
          "my_mutability.csv", row.names = FALSE)
substitution <- t(model@substitution)
write.csv(data.frame(Fivemer = rownames(substitution), substitution),
          "my_substitution.csv", row.names = FALSE)
```

## Internal node (ancestral) sequences

simble records the true sequence of every cell in the lineage. To also write the
true sequences of the internal nodes of the simplified trees, use the
`--internal-nodes` argument:

```sh
simble --internal-nodes [other args]
```

This writes `all_internal_nodes_airr.tsv`. With `--fasta`, it also writes
`all_internal_nodes.fasta`. These files are not written by default.

These sequences are not reconstructions. They are the known (true) sequences
from the simulation, so you can use them to test ancestral sequence
reconstruction (ASR) methods.

The internal nodes are the root (the naive cell) and each node where the
lineage splits. The simplified tree and the simplified time tree have the same
nodes, so one file agrees with both trees. The `cell_id` of each row is the same
as the node's name in the Newick trees. Each chain has its own row, as in
`all_samples_airr.tsv`. The file has these columns in addition to the usual
AIRR columns:

| column | description |
| ------ | ----------- |
| `node_type` | `root` or `internal` |
| `parent_id` | `cell_id` of the parent node (empty for the root) |
| `generation` | generation of the node (its position in the time tree) |
| `branch_generations` | branch length to the parent in the time tree |
| `branch_mutations` | mutations in this chain on the branch to the parent (the genetic distance tree's branch length is the sum for both chains) |
| `n_descendant_tips` | number of sampled cells below this node |
| `descendant_tips` | `cell_id`s of the sampled cells below this node, separated by `;` |

`sample_time` is empty, because internal nodes are not sampled.

In `all_samples.fasta` and `all_internal_nodes.fasta`, each record ID is the same
as the `sequence_id` in the AIRR output (`<cell_id>_heavy` or `<cell_id>_light`).
The FASTA headers do not include `descendant_tips`.

A reconstructed tree usually has a different topology from the true tree. To
compare a reconstructed node with a true node, match them by the set of sampled
cells below each node (`descendant_tips`).

## All arguments


Available arguments:

<table>
    <tr>
        <th>argument</th>
        <th>abbr</th>
        <th>default</th>
        <th>description</th>
    </tr>
    <tr>
    <td colspan=4> <b><i>Frequently used</i></b> </td>
    </tr>
    <tr>
        <td>--output</td>
        <td>-o</td>
        <td>cwd/results</td>
        <td>folder for results</td>
    </tr>
    <tr>
        <td>--number</td>
        <td>-n</td>
        <td>1</td>
        <td>number of clones to simulate</td>
    </tr>
    <tr>
        <td>--processes</td>
        <td>-p</td>
        <td>1</td>
        <td>number of processes (multiprocessing)</td>
    </tr>
    <tr>
        <td>--neutral</td>
        <td></td>
        <td></td>
        <td>if provided, runs a neutral simulation</td>
    </tr>
        <tr>
        <td>--uniform</td>
        <td></td>
        <td></td>
        <td>if provided, runs a uniform neutral simulation</td>
    </tr>
    <tr>
        <td>--migration-rate</td>
        <td></td>
        <td>0</td>
        <td>expected number of cells that leave the germinal center each generation</td>
    </tr>
    <tr>
        <td>--samples</td>
        <td>-s</td>
        <td>[0 200 25]</td>
        <td>start, stop, step, to specify sample times for germinal center</td>
    </tr>
    <tr>
    <td colspan=4> <b><i>Sampling</i></b> </td>
    </tr>
    <tr>
        <td>--other-samples</td>
        <td></td>
        <td>GC sample times</td>
        <td>start, stop, step, to specify &quot;Other&quot; location sample times</td>
    </tr>
    <tr>
        <td>--sample-size</td>
        <td></td>
        <td>50</td>
        <td>specify sample size for 'GC' location</td>
    </tr>
    <tr>
        <td>--sample-size-other</td>
        <td></td>
        <td>12</td>
        <td>specify sample size for the 'Other' location</td>
    </tr>
    <tr>
    <td colspan=4> <b><i>Model parameters</i></b> </td>
    </tr>
    <tr>
        <td>--sequence-length</td>
        <td></td>
        <td>370</td>
        <td>length of the sequence to simulate if uniform</td>
    </tr>
    <tr>
        <td>--antigen</td>
        <td>-a</td>
        <td>1000</td>
        <td>amount of antigen</td>
    </tr>
    <tr>
        <td>--single-chain</td>
        <td></td>
        <td></td>
        <td>if provided, simulates only the heavy chain of each naive pair (with the usual mutation model and selection)</td>
    </tr>
    <tr>
        <td>--mutability</td>
        <td></td>
        <td>bundled S5F</td>
        <td>path to a 5-mer mutability file for both chains</td>
    </tr>
    <tr>
        <td>--substitution</td>
        <td></td>
        <td>bundled S5F</td>
        <td>path to a 5-mer substitution file for both chains</td>
    </tr>
    <tr>
        <td>--mutability-heavy, --mutability-light</td>
        <td></td>
        <td>bundled S5F</td>
        <td>path to a 5-mer mutability file for one chain (used instead of --mutability for that chain)</td>
    </tr>
    <tr>
        <td>--substitution-heavy, --substitution-light</td>
        <td></td>
        <td>bundled S5F</td>
        <td>path to a 5-mer substitution file for one chain (used instead of --substitution for that chain)</td>
    </tr>
    <tr>
        <td>--max-children</td>
        <td></td>
        <td>100</td>
        <td>maximum number of children that one GC cell can have in one generation. Antigen above this limit is not used. This value also sets the <code>number_of_cells_with_N_children</code> columns in <code>population_data.csv</code></td>
    </tr>
    <tr>
        <td>--heavy-shm</td>
        <td></td>
        <td>0.0008908272571108565</td>
        <td>expected number of heavy chain mutations each division per site</td>
    </tr>
    <tr>
        <td>--light-shm</td>
        <td></td>
        <td>0.0004923076923076923</td>
        <td>expected number of light chain mutations each division per site</td>
    </tr>
    <tr>
        <td>--target-mutations-heavy</td>
        <td></td>
        <td>5</td>
        <td>number of amino acid mutations the target heavy chain should have</td>
    </tr>
    <tr>
        <td>--target-mutations-light</td>
        <td></td>
        <td>2</td>
        <td>number of amino acid mutations the target light chain should have</td>
    </tr>
        <tr>
        <td>--cdr-dist</td>
        <td></td>
        <td></td>
        <td>cdr distribution (constant or exponential)</td>
    </tr>
    <tr>
        <td>--cdr-var</td>
        <td></td>
        <td></td>
        <td>cdr variable</td>
    </tr>
    <tr>
        <td>--fwr-dist</td>
        <td></td>
        <td></td>
        <td>fwr distribution (constant or exponential)</td>
    </tr>
    <tr>
        <td>--fwr-var</td>
        <td></td>
        <td></td>
        <td>fwr variable</td>
    </tr>
    <tr>
        <td>--multiplier</td>
        <td>-m</td>
        <td>2</td>
        <td>selection multiplier</td>
    </tr>
    <tr>
    <td colspan=4> <b><i>Program settings</i></b> </td>
    </tr>
    <tr>
        <td>--quiet</td>
        <td>-q</td>
        <td></td>
        <td>if present, progress bar suppressed</td>
    </tr>
    <tr>
        <td>--verbose</td>
        <td>-v</td>
        <td></td>
        <td>if present, verbose information provided</td>
    </tr>
    <tr>
        <td>--fasta</td>
        <td></td>
        <td></td>
        <td>if present, also write a fasta file</td>
    </tr>
    <tr>
        <td>--internal-nodes</td>
        <td></td>
        <td></td>
        <td>if present, also write the true sequences of the simplified tree's internal nodes</td>
    </tr>
    <tr>
        <td>--config</td>
        <td></td>
        <td></td>
        <td>path to a config file (still in development)</td>
    </tr>
    <tr>
        <td>--dev</td>
        <td></td>
        <td></td>
        <td>if present, run in dev mode (not recommended)</td>
    </tr>
    <tr>
        <td>--seed</td>
        <td></td>
        <td></td>
        <td>an RNG seed to reproduce specific simulations</td>
    </tr>
    <tr>
        <td>--clone_id</td>
        <td></td>
        <td>1</td>
        <td>specify a starting clone id (1-indexed)</td>
    </tr>
    <tr>
        <td>--naive</td>
        <td></td>
        <td></td>
        <td>path to input naive sequences</td>
    </tr>
    <tr>
        <td>--keep-cols</td>
        <td></td>
        <td></td>
        <td>additional columns from the naive data to keep</td>
    </tr>
    <tr>
        <td>--naive-sampling</td>
        <td></td>
        <td>ordered if --naive is given, else random</td>
        <td>how to pick each clone's naive pair: random, ordered, or unique</td>
    </tr>
</table>

Thank you for using simble!
