#!/usr/bin/env bash
set -e
set -o pipefail

# download the files
update_from=6.34
update_to=6.35

hea_url="https://heasarc.gsfc.nasa.gov/FTP/software/lheasoft"

# download the files for the current version
mkdir -p $update_from $update_to

get_path() {
    local fname=$1
    if [ "$fname" == "xsChain.cxx" ]; then
        echo "Xspec/src/XSUser/Handler"
    else
        echo "Xspec/src/XSFit/MCMC"
    fi
}

# Generate patches for the current version if they don't exist
rm *tmp* > /dev/null 2>&1 || true
if [ $(ls *.*.patch | wc -l | awk '{print $1}') -ne 5 ]; then
    echo "Generating patches for the current version..."
    for file in `find . -maxdepth 1 -name '*.cxx' -o -name '*.h'`; do
        fname=$(basename $file)
        path=$(get_path $fname)
        
        if [ ! -f $update_from/$fname ]; then
            curl "$hea_url/lheasoft${update_from}/src/$path/$fname" -o $update_from/$fname
        fi
        diff -u $update_from/$fname $fname > $fname.patch
        # replace top two lines with name only
        echo "--- $fname" > $fname.patch.tmp
        echo "+++ $fname" >> $fname.patch.tmp
        tail -n +3 $fname.patch >> $fname.patch.tmp
        mv $fname.patch.tmp $fname.patch
    done
else
    echo "Patches already exist, skipping patch generation."
    unpatched=""
    for file in `find . -name '*.patch'`; do
        fname=$(echo $file | sed 's|\.patch$||' | sed 's|./||')
        path=$(get_path $fname)
        curl -s "$hea_url/lheasoft${update_from}/src/$path/$fname" -o $update_from/$fname
        curl -s "$hea_url/lheasoft${update_to}/src/$path/$fname" -o $update_to/$fname
        diff -u $update_from/$fname $update_to/$fname > $fname.diff
        
        if [ $(wc -l $fname.diff | awk '{print $1}') -gt 0 ]; then
            echo "+++++ Changes found in $fname"
            unpatched="$unpatched $fname"
        else
            echo "No changes in $fname between $update_from and $update_to"
            cat $fname.patch | sed "s|$fname|$update_to/$fname|g" > tmp.patch
            patch -b < tmp.patch
            rm tmp.patch $fname.diff
            # now get the new patch file
            # replace top two lines with name only
            echo "--- $fname" > $fname.patch
            echo "+++ $fname" >> $fname.patch
            diff -u $update_to/$fname.orig $update_to/$fname | tail -n +3 >> $fname.patch
        fi
    done
    echo "+++++++++++++++++++"
    echo "Unpatched files: $unpatched"
fi