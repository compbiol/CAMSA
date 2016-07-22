#!/usr/bin/env bash

echo "$TRAVIS_OS_NAME"

if [ "$TRAVIS_OS_NAME" = 'osx' ]; then
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
fi

if [ "$TRAVIS_OS_NAME" = 'linux' ]; then
    apt-get update
    apt-get install mummer
fi