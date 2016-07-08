CAMSA
==

``CAMSA``  - is a tool for **C**omparative **A**nalysis and **M**erging of **S**caffold **A**ssemblies.

``CAMSA`` is developed using Python 3.5+ programming language and requires Python 3.5+ interpreter to work properly.

To install ``CAMSA``: download the repository snapshot from the master branch, install all the python packages required for ``CAMSA`` work (listed in requirements.txt).
We suggest using a ``virtual environment`` of python, to isolate the installation.

    python3.5 -m venv camsa_env
    source camsa_env/bin/activate
    pip install -r requirements.txt

``CAMSA`` usage instructions can be found by running:

    python3.5 camsa.py --help

Input
--

``CAMSA`` expect as an input a set of different assemblies on the same set of scaffolds.
Assemblies are ought to be provided in CSV files, with the following format:

    origin, ctg1,   ctg1_or,  ctg2,    ctg2_or, gap_size,   cw
    A1,     s1,     +,        s2,       -,      ?,          ?
    A1,     s2,     -,        s3,       +,      ?,          ?
    ...


Fields description:

* ``origin``:   ``id`` of the assembly, that provides an assembly point
* ``ctg1``:     ``id`` of the first fragment, that participates in the assembly point
* ``ctg1_or``:  relative orientation, of the first fragment in the assembly point. (+/-/?)
* ``ctg2``:     ``id`` of the second fragment, that participates in the assembly point
* ``ctg2_or``:  relative orientation, of the second fragment in the assembly point. (+/-/?)
* ``gap_size``: integer value, determining a gap size between two assembled fragments (>=0 / ?)
* ``cw``:       confidence weight, of the reported assembly point ([0, 1]/?)

Output
--
``CAMSA`` stores all the performed analysis in the created folder.
By default this folder has name of the following form: ``camsa_report_%b_%d_%Y__%H_%M``.

The comparative analysis ``report`` is stored at the top level of the analysis folder.
The folder also contains all the required libraries for the report to properly work.
We drew our inspiration for the structure of the output from the [QUAST software](http://bioinf.spbau.ru/en/quast) workflow ideology.

Contributing / Bug reports
--
Please submit any information about identified bugs to the corresponding [GitHub bug tracker system](https://github.com/aganezov/CAMSA/issues).

Thank you!



