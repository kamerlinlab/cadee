#!/bin/bash
# Compute one Simpack with qdyn5p.
#
# Author: Beat Amrein
# Email: beat.amrein@gmail.com
# Date: 05.Mar 2017
# Version: 0.2
# 
# Description: Run a simulation with qdyn5p (parallel).
#              This script can be launched multiple times,
#              to allow filling a compute node.
#              The script extracts the simpack and then iterates
#              the input files (for simpack in *inp).
#
# Installation: Put this script and the parent-script (pcadee.sh) into your $PATH.
#               This script needs a SNIC enviroment to work correctly (see SETTINGS).
#
#
# This script is part of CADEE.
# If you use this script, please cite: 
#     Amrein et al. (2017), CADEE: Computer-Aided Directed Evolution of Enzymes, JUCrJ, p50-64
#     https://doi.org/10.1107/S2052252516018017
#
# Usage:
#    srunq.sh $SIMPACK
#
# TODO: 
#    on-the-fly mapping.
#
# This script is adjusted for usage on SNIC resources.
#    (see also section SETTINGS and CLUSTER CUSTOMIZATION)
#
#

############
# SETTINGS #
############


export MYNAME=$(basename "$1")
function write {
    echo -e "$MYNAME @$SECONDS > $*"
}

# Prepare Tempary Directory
if [ -z $SCRATCH_FOLDER ]
then
    # WARNING:
    #   If you adjust this for non SNIC systems, make sure
    #   you give the absolute path, or adjust the else clause
    echo "FATAL: CONFIGURATION ERROR \$SCRATCH_FOLDER not defined in environment($0)."
    exit 3
else
    TEMPDIR="$SCRATCH_FOLDER/$$"
    mkdir -p "$TEMPDIR"
    function cleanup {
        /bin/rm -rf "$TEMPDIR"
	echo "Cleanup Done."
    }
    trap cleanup EXIT
fi


##################
# INITIALIZATION #
##################

write "Start: $(date)"

# Evaluation of Input-Files
# no argument. stop.
if [ $# -eq 0 ]; then
    write "Fatal: No parameter suplied. Need a SIMPACK. Stop."
    write "Usage:"
    write "       $0 /path/to/simpack.tar"
    exit 4

# one argument. is it file or a file mask?
elif [ $# -eq 1 ]; then
    SIMPACK=$(readlink -f "$1")
    if [ ${SIMPACK##*.} == "tar" ]; then
        write "You supplied a SIMPACK: $SIMPACK"
    else
        write "Fatal: Simpacks have to have a .tar extension. Stop."
        exit 5
    fi 
# more than one argument. not ok.
else
    write "Fatal: Currently only 1 SIMPACK can be supplied. Stop."
    write "Usage:"
    write "       $0 /path/to/simpack.tar"
    exit 6
fi

# is SIMPACK variable assigned?
if [ -z $SIMPACK ]
then
    write "Fatal. Can not have happened error. Stop."
    exit 7
fi



if [ -z "$EXE" ]
then
    write "Fatal. \$EXE not assigned. Stop."
    exit 1
fi


##################
# UNPACK SIMPACK #
##################

write "Unpacking Simpack ($SIMPACK) to tmpdir ($TEMPDIR)."
cd "$TEMPDIR"
tar xf "$SIMPACK"
tarexit=$?
if [ $tarexit -ne 0 ]
then
    write "Fatal: There was a problem unpacking the simpack! (exitcode $tarexit). Stop."
    exit $tarexit
fi

#############
# Functions #
#############


# Search for changed files and append them to the simpack.
function backup {
    # if an argument is given, always backup.
    if [ -z $1 ]
    then
        BACKUP=0
    else
        BACKUP=1
    fi

    if [ $BACKUP -eq 0 ]
    then

        if [ ! -f timestamp ]
        then
            touch timestamp
            BACKUP=1
        fi

        TIME_PASSED=$(( $(date +%s) - $(date +%s -r timestamp) ))

        if [ $TIME_PASSED -gt $BACKUPINTERVAL ]
        then
            BACKUP=1
        fi
    fi

    if [ $BACKUP -ne 0 ]
    then
        tbs=$(date +%s)
        find . -newer timestamp | xargs tar --no-recursion --file=$SIMPACK --append
        tarexit=$?
        if [ $tarexit -ne 0 ]
        then
            write "FATAL: (backup) tar exited with non-zero. Stop."
            exit 1
        fi
        touch timestamp
        write "Backup Complete, Duration: $((`date +%s`-$tbs)), [ $(date) ]"
    else
        write "Backup Skipped."
    fi
}


function isnan {
    # TEST FOR NaN
    NAN=$( grep "SUM.*NaN.*NaN.*NaN"  $1 | wc -l )
    if [ $NAN -ne 0 ]
    then
        write "FATAL: We've got a NaN! Stop."
        ERROR=1
    fi
    return 0
}


################
# Print Config #
################

write ""
write ""
write "###########"
write "# CONFIG: #"
write "###########"
write " bkp int: $BACKUPINTERVAL"
write " simpack: $SIMPACK"
write " cores:   $CORES"
write " exe:     $EXE"
write "  md5sum: $(md5sum $(echo $EXE | rev | cut -d' '  -f1 | rev))"
write " workdir: $PWD"
write ""
write ""
write ""


#############
# Work Loop #
#############

ERROR=0

touch timestamp

for inp in *.inp
do
    file=${inp%.inp}
    write "Working Directory;  `hostname`:`pwd`"
    write "Preparing $file ... "

    # check for logfile ...
    if [ -f "$file.log" ]
    then
        write "LogFile $file.log exists! "
        if [ $(tail -n 10 $file.log | grep "terminated normally" | wc -l ) -eq 1 ]
        then
            write "Last run ended successfully. Skipping this run."
            continue
        else
            write "Will rerun, logfile is lacking 'terminated normally'."
        fi
    fi
    
    # check for gzipped log ...
    if [ -f "$file.log.gz" ]
    then
        write "LogFile $file.log.gz exists! "
        if [ $(zcat ${file}.log.gz | tail -n 10 | grep "terminated normally" | wc -l ) -eq 1 ]
        then
            write "Last run ended successfully. Skipping this run."
            continue
        else
            write "Will rerun, logfile is lacking 'terminated normally'."
        fi
    fi

    #RUN IT
    write "Running MD Simulation on $file.inp ..."
    
    string=$(grep "restart.*.rest.re" ${file}.inp)

    unset restraintname
    if [ ! -z $string ]
    then
        stringarray=($string)
        restraintname=${stringarray[-1]}
        restartname="${restraintname%.*}"
        if [ ! -f $restraintname ]
        then
            if [ -f ${restartname}.re ]
            then
                cp -v "${restartname}.re" "${restraintname}"
            else
                write "ERROR: Restartfile not found: ${file}.re"
                write "ERROR: Don't know what to do!"
                ERROR=1
                break
            fi
        fi
    fi

    if $EXE ${file}.inp > ${file}.log
    then
        write "Finished: ${file}.log"
        # check if NAN?
        isnan ${file}.log
        # TODO: on-the-fly mapping
        write "Zipping."
        gzip ${file}.log
        [ -f  "${file}.en" ] && gzip ${file}.en
        [[ ! -z $restraintname ]] && [[ -f "$restraintname" ]] && /bin/rm "$restraintname"
        backup
    else
        ERROR=1
    fi

    if [ $ERROR -ne 0 ]
    then
        write "Critical: An error occured with ${file}.inp"
        write "          Will now stop."
        break
    fi
done

backup 1

echo ""
echo ""
echo ""
echo ""
if [ $ERROR -eq 1 ]
then
    write "Abnormal Termination"
elif [ $ERROR -eq 0 ]
then
    write "All OK."
fi

write "End: $(date)"
write "Duration: $SECONDS s"

exit $ERROR
