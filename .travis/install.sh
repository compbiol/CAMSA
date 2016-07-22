#!/usr/bin/env bash

if [[ $TRAVIS_OS_NAME == "osx" ]]; then
    brew update
    brew install mummer
    case "${PYVER}" in
        py34)
            brew install python3.4
        ;;
        py35)
            brew install python3.5
        ;;
    esac
else
    apt-get update
    apt-get install mummer
fi