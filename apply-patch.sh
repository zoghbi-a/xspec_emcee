#!/usr/bin/env bash

set -e
set -o pipefail

heasoft_src=$1
if [ -z "$heasoft_src" ]; then
    echo "Usage: $0 <path_to_heasoft_src>"
    exit 1
fi

get_path() {
    local fname=$1
    if [ "$fname" == "xsChain.cxx" ]; then
        echo "Xspec/src/XSUser/Handler"
    else
        echo "Xspec/src/XSFit/MCMC"
    fi
}

if [ $(ls *.*.patch | wc -l | awk '{print $1}') -ne 5 ]; then
    echo "ERROR: Patches do not exist. Please run the updater script first."
    exit 1
fi

# Apply patches
for file in `find . -name '*.patch'`; do
    fname=$(echo $file | sed 's|\.patch$||' | sed 's|./||')
    path=$(get_path $fname)
    echo $path $fname
    patch -d "$heasoft_src/$path" < "$file"
done