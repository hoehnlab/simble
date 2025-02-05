# simble

simble is a BCR evolution simulator. It starts with a naive heavy and light chain pair,
and models evolution and selection in the germinal center as well as migration out of 
germinal center. simble also has customizable sampling.

:construction: This readme is still under construction. :construction:


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

To run a neutral simulation:
```sh
simble --neutral [other args]
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
| --migration-rate| | 0 | expected number of cells that leave the germinal center each generation|
| --samples | -s | [0 200 25] | start, stop, step, to specify sample times other than the default|


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
    <td colspan=4> <b><i>Model parameters</i></b> </td>
    </tr>
    <tr>
        <td>--multiplier</td>
        <td>-m</td>
        <td>2</td>
        <td>selection multiplier</td>
    </tr>
    <tr>
        <td>--heavy-mutate-probability</td>
        <td></td>
        <td>0.5</td>
        <td>expected number of heavy chain mutations each division</td>
    </tr>
    <tr>
        <td>--light-mutate-probability</td>
        <td></td>
        <td>0.25</td>
        <td>expected number of light chain mutations each division</td>
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
    <td colspan=4> <b><i>Program settings</i></b> </td>
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
        <td>--seeds</td>
        <td></td>
        <td></td>
        <td>a list of RNG seeds to reproduce specific simulations</td>
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
