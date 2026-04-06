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

When providing a naive dataset, the default behavior is that each clone will
use the naive pair of sequences at row clone_id modulo (rows in naive dataset).
For example, if the naive dataset has 10 entries, clone 6 will use the pair of 
sequences in the 6th row, and clone 13 will use the 3rd row. 

If you would like each clone to use a random starting pair from the table:

```sh
simble --naive <naive_file.csv> --naive-random
```

or if you would like to specify a clone_id to start at:
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
        <td>--naive-random</td>
        <td></td>
        <td></td>
        <td>randomly sample from naive input rather than by clone id, always true if naive file is not specified</td>
    </tr>
</table>

Thank you for using simble!
