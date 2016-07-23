#!/usr/bin/env bash

if [ "$TRAVIS_OS_NAME" = 'osx' ]; then
    brew update
    brew tap homebrew/science
    brew install pyenv
    case "$PYTHON" in
     "3.4")
     pyenv install 3.4.3
     pyenv global 3.4.3
     ;;
     "3.5")
     pyenv install 3.5.1
     pyenv global 3.5.1
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