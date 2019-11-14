#!/bin/bash

# This script creates environment.(osx|linux).lock.yml from environment.yml.
# It needs to be run manually whenever environment.yml is changed. See also
# the developer documentation.

set -euo pipefail

if [[ $# -gt 0 ]]; then
    case $1 in
        osx|linux)
            targets=$1
            ;;
        *)
            echo "Target must be osx or linux"
            exit 1
            ;;
    esac
else
    targets="linux osx"
fi

if [[ -e ~/.condarc ]]; then
    mv ~/.condarc ~/.condarc.condalock.bak
    trap "mv ~/.condarc.condalock.bak ~/.condarc" EXIT
else
    trap "rm ~/.condarc" EXIT
fi

for os in ${targets}; do
    env=blrtmp-$RANDOM-$os
    env_yml=environment.$os.lock.yml
    printf "channels:\n  - conda-forge\n  - bioconda\n  - defaults\n" > ~/.condarc
    printf "subdir: %s-64\nsubdirs:\n  - %s-64\n  - noarch\n" $os $os >> ~/.condarc
    conda env create -n $env -f environment.yml
    conda env export -n $env | grep -Ev '^(name|prefix):' > ${env_yml}
    conda env remove -n $env
    echo "Created ${env_yml}"
done

# Original .condarc will be restored by the exit trap
