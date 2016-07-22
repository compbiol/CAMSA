#!/usr/bin/env bash

if [ "$TRAVIS_OS_NAME" = 'osx' ]; then
    brew update
    brew tap homebrew/science
    brew install mummer
    case "$PYTHON" in
     "3.4")
     brew install python34 ;;
     "3.5")
     brew install python35 ;;
    esac
fi

if [ "$TRAVIS_OS_NAME" = 'linux' ]; then
    apt-get update
    apt-get install mummer
fi