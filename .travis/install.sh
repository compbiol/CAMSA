#!/usr/bin/env bash

if [ "$TRAVIS_OS_NAME" = 'osx' ]; then
    brew update
    brew tap homebrew/science
    brew install pyenv
    case "$PYTHON" in
     "3.4")
     pyenv install 3.4 ;;
     "3.5")
     pyenv install 3.5 ;;
    esac
fi

if [ "$TRAVIS_OS_NAME" = 'linux' ]; then
    apt-get update
    apt-get install mummer
fi