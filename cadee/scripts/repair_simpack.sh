#!/usr/bin/env bash

# This script is supposed to repair a run that went out of storage space (hard quoata limit).
# It can also help, if an instance was killed while writing the simpack to disk.

# Author: Beat Amrein, beat.amrein@gmail.com
# This script is part of CADEE.

# Version 0.8.1


if [ -d /scratch ]
then
    echo "Using /scratch for temporary files."
    TEMP=/scratch
else
    echo "Using /tmp for temporary files."
    TEMP=/tmp
fi


function usage(){
echo "Usage:
    $0 /path/to/tararchive.tar [ [ --force ] || [ --checkene ] ]

    This script will ...

    1. Search duplicate logfiles.
    2. Search duplicate energy files.
    3. Search missing restartfiles.
    4. Search damaged logfiles.
    5. Search for logfiles lacking 'terminated normally'.
    6. Search for gzipped logfiles lacking 'terminated normally'.
    
    ( --checkene ) 
    => 7. Search for energy files with invalid sizes and delete them including log,re and dcd.

    => 8. Repack the tarball. (If problems have been found).

    Options (exclusive):
    --force will continue even if the tar-utility has a non-0 exitcode, and always repack
    --checkene will unpack, check the energyfile - sizes and repack.
"
    exit 1
}


set -e

found=0

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

TMPFOLDER="/$TEMP/$$"

if [ -d $TMPFOLDER ]
then
    echo "TEMPORARY FOLDER EXISTS. $TMPFOLDER"
    echo "STOP. (You may try to remove it)."
    exit
fi

simpack="$(cd "$(dirname "$1")"; pwd)/$(basename "$1")"

mkdir -p $TMPFOLDER
cd $TMPFOLDER

function finish {
    /bin/rm -rf $TMPFOLDER
  }

trap finish INT TERM

pwd

set +e
tar xf $simpack
tarexit=$?
set -e


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
    


echo "1. Searching duplicate logfiles:"
for file in $(/bin/ls *.log || true)
do
    if [ -f "${file}.gz" ] 
    then
        let found+=1
        set +e
        zcmp $file "${file}.gz"
        if [ $? -eq 0 ]
        then
            echo 'log is equal to log.gz. delete log'
            /bin/rm $file
        else
            echo 'log is not equal to log.gz!'
            zcat $file.gz | tail | grep -q "terminated normally"
            if [ $? -eq 0 ]
            then
                echo "$file.gz seems correct"
                echo "==> keeping $file.gz"
                /bin/rm $file
            else
                echo "$file.gz seems bad. keeping $file."
                echo "keeping $file"
                /bin/rm $file.gz
            fi
        fi
        set -e
    fi
done

echo "2. Searching duplicate energy files:"
for file in $(/bin/ls *.en || true)
do
    if [ -f "${file}.gz" ]
    then
        let found+=1
        size=$(cat $file | wc -c)
        set +e
        gize=$(zcat ${file}.gz | wc -c)
        [ $? -ne 0 ] && gize=0
        echo "Duplicate energy files: $file/.gz"
        if [ $size -gt $gize ]
        then
            echo "Uncompressed is bigger than compressed one! Removing compressed one!"
            /bin/rm "${file}.gz"
            gzip "${file}"
        else
            echo "Removing uncompressed file."
            /bin/rm -v "${file}"
        fi
    fi
    set -e
done


echo "3. Searching missing restartfiles:"
for file in $(/bin/ls *_dyn*.log.gz ||true) $(/bin/ls *_eq.log.gz||true) $(/bin/ls ????_fep.log.gz||true)
do
    if [ ! -f "${file%.log.gz}.re" ]
    then
        echo "The logfile exists, but there is no restartfile. Removing outpuffiles of $file ..."
        /bin/rm -v $file
        /bin/rm -f -v "${file%.log.gz}.dcd"
        /bin/rm -f -v "${file%.log.gz}.en.gz"
        /bin/rm -f -v "${file%.log.gz}.log"
        let found+=1
    fi
done

echo "4. Searching damaged logfiles:"
for file in $(/bin/ls *_dyn*.log.gz ||true) $(/bin/ls *_eq.log.gz||true) $(/bin/ls ????_fep.log.gz||true)
do
    set +e
    zcat $file > /dev/null
    if [ "$?" -eq "1" ]
    then
        echo zcat error
        /bin/rm $file
        fn="${file%.*}"
        fn="${fn%.*}"
        /bin/rm $fn.re
        /bin/rm $fn.dcd
        /bin/rm $fn.en
        let found+=1
    fi
done
set -e

echo "5. Search for logfiles lacking 'terminated normally':"
for file in $(/bin/ls *_dyn*.log ||true) $(/bin/ls *_eq.log||true) $(/bin/ls ????_fep.log||true)
do
    ok=1
    set +e
    (tail $file | grep -q "terminated normally") || ok=0
    set -e
    if [ $ok -eq 0 ]
    then
        echo 'Bad Logfile:' $file 'deleting ...'
        /bin/rm $file
        let found+=1
    fi
done

echo "6. Search for gzipped logfiles lacking 'terminated normally':"
for file in $(/bin/ls *_dyn*.log.gz ||true) $(/bin/ls *_eq.log.gz||true) $(/bin/ls ????_fep.log.gz||true)
do
    ok=1
    set +e
    (zcat $file | tail | grep -q "terminated normally") || ok=0
    set -e
    if [ $ok -eq 0 ]
    then
        echo 'Bad Logfile:' $file 'deleting ...'
        /bin/rm $file
        let found+=1
    fi
done

if [ ! -z $2 ] && [ "$2" == '--checkene' ]
then
    echo "7. Search for energy files with invalid sizes:"
    echo "    i) Unzipping energy files"
    gunzip -f *.en.gz || true
    echo "   ii) collecting files sizes ..."
    /bin/ls -l *.en | awk '{ print $5}' | sort | uniq -c | sort -g 
    
    numsizes=$(/bin/ls -l *.en | awk '{print $5}' | sort | uniq -c | sort -g | wc -l)
    
    if [ $numsizes -gt 2 ]
    then
        echo "Epected 2 files sizes, found $numsizes ."
        echo "There is a problem in this simpack ..."
        size1=$(/bin/ls -l *.en | awk '{print $5}' | sort | uniq -c | sort -g | awk '{ print $2}' | tail -n 1 )
        size2=$(/bin/ls -l *.en | awk '{print $5}' | sort | uniq -c | sort -g | awk '{ print $2}' | tail -n 2 | head -n 1)

        echo "  iii) Comparing to assigned sizes; Size1=$size1 and Size2=$size2."
        for fil in $(/bin/ls *.en)
        do  
            size=$(cat $fil | wc -c)
            if [[ $size -ne $size1 ]] && [[ $size -ne $size2 ]]
            then    
                let found+=1
        	noext="${fil/.en/}"
    	        echo "Energy File Size for $fil is invalid ($size). Removing $noext/en/re/log/dcd ..."
        	echo $noext
                /bin/rm -rfv "${noext}.dcd"
                /bin/rm -rfv "${noext}.log"
                /bin/rm -rfv "${noext}.log.gz"
                /bin/rm -rfv "${noext}.re"
        	/bin/rm -rfv "$fil"
            fi  
        done
    else
        echo "File sizes OK."
    fi
    echo "    iv) ReZipping EnergyFiles ..."    
    gzip *en
fi

echo -n "$simpack: "
if [ $found -gt 0 ]
then
    echo "Found $found issues! Overwriting simpack ..."
    if [ $found -gt 100 ]
    then
	echo ""
	echo "FOUND MANY PROBLEMS. SEE ABOVE."
	echo ""
	echo "THIS WILL OVERWRITE YOUR SIMPACK."
	echo
	echo "DO YOU HAVE A BACKUP?"
	echo
	select yn in "Yes" "No"; do
    	case $yn in
        	Yes ) 
			break;;
        	No ) 
			echo "Stop."; 
			exit;;
    	esac
	done
    fi
    tar cf $simpack *
    echo "Repacking Done. Ready for resubmission."
else
    echo "No Problems with this simpack. Awesome!"
fi

finish
