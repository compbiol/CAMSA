#!/usr/bin/env bash

if [ "$TRAVIS_OS_NAME" = 'osx' ]; then
    if [ "$PYTHON" = "2.7" ]; then
        wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
    else
        wget https://repo.continuum.io/miniconda/Miniconda3-latest-MacOSX-x86_64.sh -O ~/miniconda.sh
    fi
    bash ~/miniconda.sh -b -p $HOME/miniconda
    export PATH="$HOME/miniconda/bin:$PATH"
    conda install -c bioconda mummer
    nucmer -h
    brew update
#    brew tap homebrew/science
    brew install pyenv
#    brew install mercurial
#    brew install mummer
#    brew install gcc
    case "$PYTHON" in
    "2.7")
     pyenv install 2.7.14
     pyenv global 2.7.14
     ;;
     "3.5")
     pyenv install 3.5.4
     pyenv global 3.5.4
     ;;
     "3.6")
     pyenv install 3.6.3
     pyenv global 3.6.3
     ;;
     "pypy.2")
     pyenv install pypy-5.3
     pyenv global pypy-5.3
     ;;
     "pypy.3")
     pyenv install pypy3.3-5.2-alpha1
     pyenv global pypy3.3-5.2-alpha1
     ;;
    esac
    `pyenv which pip` install virtualenv
    `pyenv which virtualenv` camsa-env
fi

if [ "$TRAVIS_OS_NAME" = 'linux' ]; then
    sudo apt-get update
    sudo apt-get install mummer
    pip install virtualenv
    virtualenv camsa-env
fi