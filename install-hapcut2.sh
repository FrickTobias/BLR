#!/bin/bash

set pipefail -eox

# Based on https://github.com/vibansal/HapCUT2/releases.

git clone https://github.com/vibansal/HapCUT2.git

cd HapCUT2

# For versions 1.1
git checkout 3cb79c1

make