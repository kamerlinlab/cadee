#!/usr/bin/env bash

# This script is will compress a simpack lossy, by deleting all dcd files and removing information on us runs.

# Author: Beat Amrein, beat.amrein@gmail.com
# This script is part of CADEE.



echo "WARNING: This will *lossy* compress tarballs in $PWD"
echo 'press ctrl+c to abort within 5 secs'

sleep 5


set -e

wd=$PWD
for fil in $(ls *.tar)
do
    echo unpack $fil
    cd $wd
    mkdir -p /dev/shm/tmp/$$
    cd /dev/shm/tmp/$$
    tar xf $wd/$fil
    echo compress
    rm -f *dcd
    rm -f *{1..9}_fep.re
    rm -f *{1..9}_fep.en.gz
    echo hash
    md5sum * > hashes.md5
    echo repack
    tar cf $fil *
    set +e
    mv $fil $wd/$fil || true
    set -e
    echo clean
    rm -r /dev/shm/tmp/$$
done


