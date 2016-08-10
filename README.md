CAMSA
==
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/hyperium/hyper/master/LICENSE)
[![Build Status](https://travis-ci.org/aganezov/CAMSA.svg?branch=dev)](https://travis-ci.org/aganezov/CAMSA)


CAMSA  - is a tool for **C**omparative **A**nalysis and **M**erging of **S**caffold **A**ssemblies, distributed as a standalone software under the [MIT license]((https://github.com/aganezov/CAMSA/blob/master/LICENSE.txt)).

CAMSA is developed using Python programming language and is compatible with both Python2 (2.7+) and Python3 (3.3+). This means CAMSA can be ran on any modern operating system properly.
>Note: only Linux and MacOS were tested with respect to CAMSA installation and usage. If you want to use CAMSA on windows, please contact us.

For detailed instructions and documentation please refer to [CAMSA wiki]().

Installation:
---

The simplest way to install CAMSA is by using `pip`. Open your terminal and execute the following command: 

    pip install camsa

For more details please refer to the [CAMSA installation wiki page]().

Input
--

``CAMSA`` expect as an input a set of different assemblies on the same set of scaffolds.
Eash assembly must be represented as a set of assemblies points in TSF files, using the following format:

    origin    seq1    seq1_or    seq2    seq2_or    gap_size    cw
    A1        s1      +          s2      -          ?           ?
    A1        s2      -          s3      +          ?           ?
    ...


Fields description:

* ``origin``:   ``id`` of the assembly, that produced a corresponding assembly point.
* ``seq1``:     ``id`` of the first fragment, that participates in the assembly point.
* ``seq1_or``:  relative orientation `(+/-/?)`, of the first fragment in the assembly point. 
* ``seq2``:     ``id`` of the second fragment, that participates in the assembly point.
* ``seq2_or``:  relative orientation `(+/-/?)`, of the second fragment in the assembly point. 
* ``gap_size``: integer value `(>=0/?)`, determining a gap size between two assembled fragments.
* ``cw``:       confidence weight, of the reported assembly point `([0, 1]/?)`


For more details please refer to the [CAMSA input wiki page]().

Output
--
``CAMSA`` comparative analysis results are stored in the automatically generated comprehensive report, powered by HTML, CSS and JavaScript.
This report format allows one to process the report on any operating system, easily share results with colleges and collaborators, as well as ensure reproducibility of each experiment. The enclosing folder also contains all the required libraries for the report to properly work, so internet connection is required.

For more details please refer to the [CAMSA output wiki page]().

Contributing / Bug reports
--
Please submit any information about identified bugs to the corresponding [GitHub bug tracker system](https://github.com/aganezov/CAMSA/issues).

Thank you!



