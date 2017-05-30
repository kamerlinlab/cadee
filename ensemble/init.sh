#!/usr/bin/env bash

# Load and check libraries for ensemble run.

# Author: Beat Amrein, beat.amrein@gmail.com
# This script is part of CADEE.


# tempoary file shall be stored in ...
export CADEE_TMP=$SNIC_TMP

# set to 'ignore', to ignore this.
CADEE_SITE=$SNIC_RESOURCE
CADEE_SITE=ignore

echo -n "Initializing ... "
ERR=0
if [ -z $CADEE_SITE ]
then
    echo "Environment variable CADEE_SITE is not set!"
    echo "UNKNOWN MACHINE! Please adjust settings in $0"
    ERR=1
elif [ "$CADEE_SITE" == "triolith" ]
then
    echo  "for triolith ..."
    module add impi/5.1.3
    module add intel/16.0.2
    module add python/2.7.6
elif [ "$CADEE_SITE" == "abisko" ]
then
    echo -n "for abisko ..."
    module add impi/5.1.3
    module add intel/16.0.2
elif [ "$CADEE_SITE" == "tintin" ] 
then
    echo -n "for tintin ..."
    module add intelmpi/5.1.3
    module add intel/16.2
    module add python/2.7.6
elif [ "$CADEE_SITE" == "ignore" ] 
then
    echo "Ignoring site-settings $(hostname) etc. ..."
else
    echo "NO SETTINGS FOR $CADEE_SITE!"
    ERR=1
fi


scriptdir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

qdynexe="$scriptdir/executables/q/qdyn5"

if ! [ -x "$qdynexe" ]
then    
    ERR=2
    echo
    echo 'CANT FIND executable qdyn5!'
    echo ''
    echo "You must copy it into place and flag it executable:"
    echo "$qdynexe"
else
    echo "Getting MD5 of qdyn5: $qdynexe) ... "
    md5sum $qdynexe || md5 $qdynexe
fi

echo " checking if mpi4py is available ... "
python -c "import mpi4py"
if [ $? -ne 0 ]
then
    ERR=3
    echo ''
    echo 'ERROR:'
    echo 'mpi4py is not available!'
    echo
    echo 'You have to install it!'
    echo '  curl -O https://pypi.python.org/packages/26/b4/1a9678ec113b5c654477354169131c88be3f65e665d7de7c5ef306f2f2a5/mpi4py-1.3.1.tar.gz'
    echo '  tar xf mpi4py-1.3.1.tar.gz'
    echo '  cd mpi4py-1.3.1'
    echo '  python setup.py install --user --prefix=""'
    echo '  cd -'

    echo 'I will try to do that for you in 5s (press ctrl+c to interrupt):'

    sleep 5

    curl -O https://pypi.python.org/packages/26/b4/1a9678ec113b5c654477354169131c88be3f65e665d7de7c5ef306f2f2a5/mpi4py-1.3.1.tar.gz
    tar xf mpi4py-1.3.1.tar.gz
    cd mpi4py-1.3.1
    python setup.py install --user --prefix=""
    cd -

    echo 'Done! Checking if it works:'
    python -c "import mpi4py"
    if [ $? -eq 0 ]
    then
        echo 'Great, mpi4py works now!'
        ERR=0
    fi
fi

# checking if qscripts are there...
if [ ! -d qscripts ]
then
    ERR=4
    echo ''
    echo 'ERROR:'
    echo 'qscripts not found'
else
    if [ ! -f 'qscripts/qscripts.cfg' ]
    then
        ERR=5
        echo ''
        echo 'ERROR:'
        echo 'qscripts not configured (missing qscripts.cfg)'
    fi
fi

python -c 'import sys; print(sys.version_info[:])'  | grep -q "(2, 7,"
if [ $? -ne 0 ]
then
    echo 'ERROR: Need python 2.7'
    ERR=6 
fi

if [ $ERR -eq 0 ]
then
    echo " ... OK! you are set!"
else
    if [ $ERR -eq 1 ]
    then
        echo "There is an error with your configuration $scriptdir/init.sh"
    fi

    echo "Please fix and rerun!"
    
    if [ -z "$PS1" ]
    then
        echo "This script runs with login-shell"
        exit $ERR
    else
        [[ $_ != $0 ]] && echo "Script is being sourced" || exit $ERR
    fi
fi
