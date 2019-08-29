#!/bin/bash

set pipefail -eox

# Based on https://github.com/vibansal/HapCUT2/releases.

git clone https://github.com/vibansal/HapCUT2.git

cd HapCUT2

make

if [[ "$OSTYPE" == "linux-gnu" ]]; then
    sudo make install-hairs

    sudo make install-hapcut2
fi

# Test install

./build/extractHAIRS

./build/HAPCUT2