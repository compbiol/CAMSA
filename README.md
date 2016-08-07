CAMSA
==
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/hyperium/hyper/master/LICENSE)
[![Build Status](https://travis-ci.org/aganezov/CAMSA.svg?branch=dev)](https://travis-ci.org/aganezov/CAMSA)


``CAMSA``  - is a tool for **C**omparative **A**nalysis and **M**erging of **S**caffold **A**ssemblies.

``CAMSA`` is distributed under the ***MIT license***. Refer to [LICENSE.txt](https://github.com/aganezov/CAMSA/blob/master/LICENSE.txt) file for more details.

``CAMSA`` is developed using Python programming language and is compatible with both python2 and python3 interpreter to work properly.

To install ``CAMSA`` download the repository snapshot from the master branch, install all the python packages required for ``CAMSA`` work (listed in requirements.txt).
We suggest using a ``virtual environments`` of python, to isolate the installation.

    virtualenv camsa-env
    source camsa-env/bin/activate
    pip install -r requirements.txt

``CAMSA`` usage instructions can be found by running:

    run_camsa.py --help

With any questions about installation or usage, please, contact **Sergey Aganezov** *[aganezov(at)gwu.edu]*

Input
--

``CAMSA`` expect as an input a set of different assemblies on the same set of scaffolds.
Assemblies are ought to be provided in CSV files, with the following format:

    origin    seq1    seq1_or    seq2    seq2_or    gap_size    cw
    A1        s1      +          s2      -          ?           ?
    A1        s2      -          s3      +          ?           ?
    ...


Fields description:

* ``origin``:   ``id`` of the assembly, that provides an assembly point
* ``seq1``:     ``id`` of the first fragment, that participates in the assembly point
* ``seq1_or``:  relative orientation, of the first fragment in the assembly point. (+/-/?)
* ``seq2``:     ``id`` of the second fragment, that participates in the assembly point
* ``seq2_or``:  relative orientation, of the second fragment in the assembly point. (+/-/?)
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



