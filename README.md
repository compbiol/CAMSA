CAMSA
==
[![MIT licensed](https://img.shields.io/badge/license-MIT-blue.svg)](https://raw.githubusercontent.com/hyperium/hyper/master/LICENSE)
[![Build Status](https://travis-ci.org/compbiol/CAMSA.svg?branch=master)](https://travis-ci.org/aganezov/CAMSA)
[![Build status](https://ci.appveyor.com/api/projects/status/u9f7fwsx86k7y8r4/branch/master?svg=true)](https://ci.appveyor.com/project/aganezov/camsa/branch/master)


CAMSA  - is a tool for **C**omparative **A**nalysis and **M**erging of **S**caffold **A**ssemblies, distributed both as a standalone software package and as Python library under the [MIT license]((https://github.com/aganezov/CAMSA/blob/master/LICENSE.txt)).

Main CAMSA features: 

1. works with any number of scaffold assemblies in de-novo non-progressive fashion
2. allows to simultaniously work with scaffold assemblies obtained from any *in silico* and *in vitro* techniques, supporting multiple existing formats via built-in converters
3. creates an extensive report with several comparative quality metrics (both on assembly level and on the level of individual assembly points)
4. constructs a merged combined scaffold assembly
5. provides an interactive framework for a visual comparative analysis of the given assemblies

CAMSA is developed using Python programming language and is compatible with both Python2 (2.7+) and Python3 (3.3+). This means CAMSA can be run on any modern operating system properly.
>Note: only Linux and MacOS were tested extensively with respect to CAMSA installation and usage. If you want to use CAMSA on windows, please try to do so on your own (we have continuous integration setup on Windows platform as well, but no as detailed as on MacOS/Linux), and if you run into any problems - contact us.


Installation:
---

The simplest way to install CAMSA is by using `pip`. Open your terminal and execute the following command: 

    pip install camsa

For more details and other ways to install CAMSA please refer to the [installation wiki page](https://github.com/compbiol/CAMSA/wiki/Installation).

Input
--

``CAMSA`` expect as an input a set of different assemblies on the same set of scaffolds.
Eash assembly must be represented as a set of assemblies points in TSF files, using the following format:

    origin    seq1    seq1_or    seq2    seq2_or    gap_size    cw
    A1        s1      +          s2      -          ?           ?
    A1        s2      -          s3      +          ?           ?
    ...

CAMSA also provides a set of built-in scripts, that allow one to translate other scaffold assembly formats (i.e., [FASTA](https://en.wikipedia.org/wiki/FASTA_format), [AGPv2.0](https://www.ncbi.nlm.nih.gov/assembly/agp/AGP_Specification/), [GRIMM](http://grimm.ucsd.edu/GRIMM/), etc) into CAMSA format, and vice-versa (when possible). For more details please refer to the [input wiki page](https://github.com/compbiol/CAMSA/wiki/Input).

Usage
--
CAMSA is very straightforward to use: installation process adds several executable scripts to your (python)path, which either execute CAMSA itself, or some of its utilities. Basic CAMSA usage is

    run_camsa.py f1.camsa.points f2.camsa.points ... -o output_dir

This command would run CAMSA in both comparative and merging modes, producing extensive comparative assembly reports as well as a merged assembly. For more details on running CAMSA please refer to the [usage wiki page](https://github.com/aganezov/CAMSA/wiki/Usage).

Output
--
``CAMSA`` comparative analysis results are stored in the automatically generated a set of text-based reports as well as a single comprehensive *interactive* report, powered by HTML, CSS, and JavaScript.
This report format allows one to process the output on any operating system, easily share results with colleges and collaborators, as well as ensure reproducibility of each experiment. The output folder also contains all the required libraries for the interactive report to properly work, so internet connection is required.

For more details regarding both text-based and interactive reports please refer to the [output wiki page](https://github.com/aganezov/CAMSA/wiki/Output).

Contributing / Bug reports
--
Please submit any information about identified bugs to the corresponding [GitHub bug tracker system](https://github.com/compbiol/CAMSA/issues).

Thank you!