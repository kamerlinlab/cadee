#!/usr/bin/env bash

# This script is supposed to repair a run that went out of storage space (hard quoata limit).

# Author: Beat Amrein, beat.amrein@gmail.com
# This script is part of CADEE.

function usage(){
echo "Usage:
    $0 /path/to/tararchive.tar [ --force ]

    1) Unpack tarball
    2) Check for logfile that without restartfile
    3) delete logfiles without restartfile
    4) repack

    --force will continue even if the tar-utility has a non-0 exitcode, and always repack
"
    exit 1
}


set -e

if [ -z $1 ]
then
    echo "need a simpack (.tar) to fix (got noting)!"
    usage
fi

if [ ! -f $1 ]
then
    echo "need a simpack (.tar) to fix (got $1)!"
    usage
fi


mkdir /tmp/$$
cd /tmp/$$

function finish {
    /bin/rm -r /tmp/$$/
  }

trap finish EXIT

pwd

set +e
tar xf $1
tarexit=$?
set -e

found=0

if [ $tarexit -ne 0 ]
then
echo "tar exit code: $tarexit"
if [ "$2" == '--force' ]
then
    echo Specified --force flag.
    echo Will now repack the tar archive.
    found=1
else
    echo There was an error unpacking.
    echo use --force as second argument to force unpacking
    exit 2
fi
fi

for file in $(ls *.log || true)
do
    if [ -f "${file}.gz" ]
    then
        echo "both compressed and uncompressed logfile found. removing uncompressed $file "
        /bin/rm $file
        let found+=1
    fi
done

for file in $(ls *_dyn*.log.gz ||true) $(ls *_eq.log.gz||true) $(ls *_fep.log.gz||true)
do
    if [ ! -f "${file%.log.gz}.re" ]
    then
        echo "logfile exists, but no restartfile. removing outpuffiles of $file ..."
        rm -v $file
        rm -f -v "${file%.log.gz}.dcd"
        rm -f -v "${file%.log.gz}.en.gz"
        rm -f -v "${file%.log.gz}.log"
        let found+=1
    fi
done

echo -n "Simpack: $1: "
if [ $found -gt 0 ]
then
    echo "Repacking Simpack..."
    tar cf $1 *
else
    echo "No problem with archived files found!"
fi


