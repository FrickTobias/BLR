#!/bin/bash

set pipefail -eox

# Based on https://github.com/vibansal/HapCUT2/releases.

git clone https://github.com/vibansal/HapCUT2.git

cd HapCUT2

make

# Test install

./build/extractHAIRS

./build/HAPCUT2