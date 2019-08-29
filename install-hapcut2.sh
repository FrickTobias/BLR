#!/bin/bash

set pipefail -eox

# Based on https://github.com/vibansal/HapCUT2/releases.

git clone --recursive https://github.com/vibansal/HapCUT2.git

cd HapCUT2

# For versions 1.1
git checkout 3cb79c1

make

# Test install

./build/extractHAIRS

./build/HAPCUT2