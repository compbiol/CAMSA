#!/usr/bin/env bash

if [ "$TRAVIS_OS_NAME" = 'osx' ]; then
    brew update
    brew install mummer
    case "$PYTHON" in
     "3.4")
     brew install python3.4 ;;
     "3.5")
     brew install python3.5 ;;
    esac
fi

if [ "$TRAVIS_OS_NAME" = 'linux' ]; then
    apt-get update
    apt-get install mummer
fi