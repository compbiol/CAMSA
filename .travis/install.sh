#!/usr/bin/env bash

if [ "$TRAVIS_OS_NAME" = 'osx' ]; then
    brew update
    brew tap homebrew/science
    brew install pyenv
    brew install mercurial
    brew install mummer
    brew install gcc
    case "$PYTHON" in
    "2.7")
     pyenv install 2.7.11
     pyenv global 2.7.11
     ;;
     "3.4")
     pyenv install 3.4.3
     pyenv global 3.4.3
     ;;
     "3.5")
     pyenv install 3.5.1
     pyenv global 3.5.1
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