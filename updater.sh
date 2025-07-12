#!/usr/bin/env bash

# globale variables.
CURRENT_VERSION=6.34
TARGET_VERSION=6.35

HEA_URL="https://heasarc.gsfc.nasa.gov/FTP/software/lheasoft"
FILES=("xsChain.cxx" "Chain.h" "Chain.cxx" "ChainManager.h" "ChainManager.cxx")


get_path() {
    local fname=$1
    if [ "$fname" == "xsChain.cxx" ]; then
        echo "Xspec/src/XSUser/Handler"
    else
        echo "Xspec/src/XSFit/MCMC"
    fi
}

generate_patches() {
    # generate patches for from modified files. Before 6.35, we had the modified
    # files in the repository, this generates patches from them.
    echo "****\nGenerating patches for the current version...\n***\n"
    local update_from=$1
    for file in `find . -maxdepth 1 -name '*.cxx' -o -name '*.h'`; do
        fname=$(basename $file)
        path=$(get_path $fname)
        echo "Found $fname"
        
        curl "$hea_url/lheasoft${update_from}/src/$path/$fname" -o $update_from/$fname
        diff -u $update_from/$fname $fname > $fname.patch
        # replace top two lines with name only
        echo "--- $fname" > $fname.patch.tmp
        echo "+++ $fname" >> $fname.patch.tmp
        tail -n +3 $fname.patch >> $fname.patch.tmp
        mv $fname.patch.tmp $fname.patch
    done
}

download_version() {
    # Given a file name in FILES, return the path in Xspec source
    local version=$1
    mkdir -p $version
    echo "Downloading files for version $version..."
    for fname in ${FILES[@]}; do
        path=$(get_path $fname)
        if [ ! -f "$version/$fname" ]; then
            curl -s "$HEA_URL/lheasoft${version}/src/$path/$fname" -o $version/$fname
        fi
    done
}

apply_patches() {
    # apply patch to a new version
    local update_from=$1
    local update_to=$2
    for fname in ${FILES[@]}; do
        path=$(get_path $fname)
        diff -u $update_from/$fname $update_to/$fname > $fname.diff

        if [ $(wc -l $fname.diff | awk '{print $1}') -gt 0 ]; then
            printf "%.20s | CHANGED $update_from to $update_to\n" $fname
            unpatched="$unpatched $fname"
        else
            printf "%.20s | no change between $update_from to $update_to\n" $fname
            cat $fname.patch | sed "s|$fname|$update_to/$fname|g" > tmp.patch
            echo "    apply patch for $fname"
            patch -b < tmp.patch
            rm tmp.patch $fname.diff
        fi
    done
    echo "+++++++++++++++++++"
    echo "Unpatched files: $unpatched"
}

# Start of the script
download_version $CURRENT_VERSION
download_version $TARGET_VERSION
apply_patches $CURRENT_VERSION $TARGET_VERSION
