
Run a default simulation with a specified output folder:

.. code:: sh

   simble -o <path-to-folder>

.. note:: The default simulation runs with selection, with no
   migration, sampling every 25 generations for 200 generations.

At any time, you can run simble with the ``-h`` or ``--help`` flag to see all available 
arguments (also see :ref:`all_arguments`):

.. code:: sh

   simble -h

To run a neutral-selection BCR simulation (using naive BCRs, heavy and
light chains, and S5F mutation/substitution model):

.. code:: sh

   simble --neutral [other args]

To run a uniformly neutral simulation (no selection, randomly generated
starting nucleotide sequence, and uniform mutations/substitutions):

.. code:: sh

   simble --uniform [other args]

To run a uniformly neutral simulation with a specified sequence length
of 100:

.. code:: sh

   simble --uniform --sequence-length 100 [other args]

To run with expected migration of one cell every 25 generations:

.. code:: sh

   simble --migration-rate 0.04 [other args]

To run 5 clones in parallel across 2 processes, with expected migration
of one cell every 10 generations with selection, and sampling every 10
generations for 100 generations:

.. code:: sh

   simble -o ./current-results -n 5 --processes 2 --migration-rate 0.1 --samples 0 100 10

which is equivalent to

.. code:: sh

   simble -o ./current-results -n 5 -p 2 --migration-rate 0.1 -s 0 100 10


.. tip:: Flags can be provided any order.

Frequently used arguments:

+-----------------+---------------------+--------------+-------------------+
| argument        | abbr                | default      | description       |
+=================+=====================+==============+===================+
| –output         | -o                  | cwd/results  | folder for        |
|                 |                     |              | results           |
+-----------------+---------------------+--------------+-------------------+
| –number         | -n                  | 1            | number of clones  |
|                 |                     |              | to simulate       |
+-----------------+---------------------+--------------+-------------------+
| –processes      | -p                  | 1            | number of         |
|                 |                     |              | processes         |
|                 |                     |              | (multiprocessing) |
+-----------------+---------------------+--------------+-------------------+
| –neutral        |                     |              | if provided, runs |
|                 |                     |              | a neutral         |
|                 |                     |              | simulation        |
+-----------------+---------------------+--------------+-------------------+
| –uniform        |                     |              | if provided, runs |
|                 |                     |              | a uniform neutral |
|                 |                     |              | simulation        |
+-----------------+---------------------+--------------+-------------------+
| –migration-rate |                     | 0            | expected number   |
|                 |                     |              | of cells that     |
|                 |                     |              | leave the         |
|                 |                     |              | germinal center   |
|                 |                     |              | each generation   |
+-----------------+---------------------+--------------+-------------------+
| –samples        | -s                  | [0 200 25]   | start, stop,      |
|                 |                     |              | step, to specify  |
|                 |                     |              | sample times      |
|                 |                     |              | other than the    |
|                 |                     |              | default           |
+-----------------+---------------------+--------------+-------------------+
| –quiet          | -q                  |              | don’t display     |
|                 |                     |              | progress bar      |
+-----------------+---------------------+--------------+-------------------+