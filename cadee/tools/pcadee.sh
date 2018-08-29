#!/bin/bash
#SBATCH -A snic2016-34-27
#SBATCH -p node
#SBATCH --nodes 1
#SBATCH -n 16
#SBATCH -t 24:00:00
#
# Compute a folder with Simpacks with parallel Qdyn6p.
#
# Author: Beat Amrein
# Email: beat.amrein@gmail.com
# Date: 22.Feb 2017
# Version: 0.1
# 
# Description: Iterate trough a folder with Simpacks, using Qdyn6p.
#              Make sure you first test and adjust the child-script (srunq.sh).
#              Also note, that the child-script is expected to be placed in the
#              same folder like this script (see DIR-var need to adjust this).
#
# Installation: Put this script and the child-script (srunq.sh) into your $PATH.
#               This script needs a SLURM enviroment to work correctly.
#
#
# This script is part of CADEE.
# If you use this script, please cite: 
#     Amrein et al. (2017), CADEE: Computer-Aided Directed Evolution of Enzymes, JUCrJ, p50-64
#     https://doi.org/10.1107/S2052252516018017
#
# Usage:
#    pcadee.sh /folder/with/simpacks
#
# TODO: If too few simpacks, auto-adjust
#       (say CORES=4, but AVAIL=16 and only 2 Simpacks, increase CORES to 8)
# TODO: Before exiting, check if new simpacks available and run those, too.
#
#
# This script is adjusted for usage on SLURM resources.
#    (see also child-script srunq.sh)
#################
# USER SETTINGS #
#################

# If you are not running SLURM, you may adjust following lines yourself:
# SLURM_NTASKS      # number of cores to use
# SLURM_NNODES      # number of nodes to use (only 1 supported)

SLURM_NTASKS=4      # number of cores to use
SLURM_NNODES=1      # number of nodes to use (only 1 supported)

export CORES=4
export MACHINE_NAME="$(hostname)"
export SCRATCH_FOLDER=/tmp
export BACKUPINTERVAL=540 # DEFAULT: 540, 9 minutes


# folder where srunq.sh is located:
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

#########################
# CLUSTER CUSTOMIZATION #
#########################

echo "This is $MACHINE_NAME. Loading Modules:"
case "$MACHINE_NAME" in
"rackham")
    ml intel/17.1 intelmpi/17.1 python/2.7.11
    export QPATH="/home/fabst747/qsource/bin"
    export EXE="mpiexec -n $CORES -bind-to none $QPATH/Qdyn6p"
    ;;


"abisko")
    if [ $CORES -ne 6 ]
    then
        echo "WARNING: You run a computation with \$CORES=$CORES on abisko."
        echo "This is discouraged. You should run it with 6 CORES."
        exit
    fi
    module -v load pgi/14.3-0
    module -v load openmpi/pgi/1.8.1
    export EXE="srun -n $CORES Qdyn6p"
    ;;


"")
    echo "This job is not running in SNIC environment."
    ;;


*)
    echo "THIS CLUSTER IS UNKNOWN!"
    echo "I will not add modules"

    export EXE="mpiexec -n $CORES $(which Qdyn6p)"
    ;;


esac

echo "Adjusted Q PATH! Executables:"
echo $(/bin/ls $QPATH)

if [ -z "$EXE" ]
then
    echo "FATAL:"
    echo "      You must configure the pcadee script properly."
    echo "      \$EXE is not defined."
    exit 1
fi




##############
# INITIALIZE #
##############

# slurm?
if [ -z $SLURM_NTASKS ]
then
    echo "Fatal: Need SLURM environment. Stop."
    exit 1
fi

# more than 2 nodes assignned?
if [ $SLURM_NNODES -ne 1 ]
then
    echo "FATAL: User Error"
    echo "       This script can distribute jobs to up to 1 nodes, you asked for $SLURM_NNODES ."
    echo "       The scrpit stops now, so you do not waste compute time. Bye!."
    exit 1
fi

# simpack-folder existing?
if [ -z $1 ]
then
    echo "Fatal: Missing Argument: Folder with Simpacks. Stop."
    echo "Usage:"
    echo "       $0 /path/to/folder/with/simpacks/"
    exit 1
fi


export MAXTASK=$(($SLURM_NTASKS/$CORES))

SIMPACK_FOLDER=$(readlink -f "$1")

################
# PRINT CONFIG #
################

echo "Simpack Folder $SIMPACK_FOLDER"

echo "Will use $CORES per simpack."
echo "Will run at most $MAXTASK simpacks at one time."
echo "This will use $(($CORES*$MAXTASK)) cores from $SLURM_NTASKS"

echo ""
echo ""
echo "   Will Distribute Jobs and Start Work in 1 Second"
echo "   ==============================================="
echo ""

sleep 1


find $SIMPACK_FOLDER -name "*.tar" | xargs -i --max-procs=$MAXTASK bash -c "echo {}; $DIR/srunq.sh {}; echo {}; exit"

echo ""
echo ""
echo ""
echo "No Simpacks left. Terminating after $SECONDS."

