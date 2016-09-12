#!/usr/bin/env bash
# Check configuration for running the ensemble trajectories

echo -n "This script will check your configuration on this machine: "
hostname

sleep 1

CADEE_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "checking CADEE Ensemblesimulation setup ..."

echo "check python version"

python -c 'import sys; print(sys.version_info[:])'  | grep -q "(2, 7,"
if [ $? -ne 0 ]
then
    echo 'ERROR: Need python 2.7'
fi

echo "check if qscripts are in place..."
echo "================================="
cd $CADEE_DIR

cd qscripts

function fix_sq(){
    echo 'Fatal: qscripts is not configured properly. please do so first:'
    echo "run: qscripts_config.py in the qscript folder"
    exit 1
}

if ! [ -e "qscripts.cfg" ]
then
    fix_sq
fi

qfep_sq=$(grep "^qfep = "  qscripts.cfg | wc -l )

if [ $qfep_sq -eq 0 ]
then
    echo "qscripts has no qfep-executable defined. please fix this."
    fix_sq
elif [ $qfep_sq -ne 1 ]
then
    echo "qscripts has two or more qfep-executables defined. i dont know which to check. please fix this."
    fix_sq
fi


echo "check if qscripts / qfep5 is in place..."

function fix_qfep(){
    echo 'Fatal: qscripts is not configured properly. (qfep5)'
    exit 1
}


qfep_exe=$(grep "^qfep = " qscripts.cfg | sed "s/.*qfep = //g")

if [ ! -f $qfep_exe ]
then
    echo "qfep exe does not exist ($qfep_exe)"
    fix_qfep
fi

if [ ! -x $qfep_exe ]
then
    echo "qfep exe is not executable, please fix this ($qfep_exe)"
    fix_qfep
fi


echo "check if openbabel is in path"
echo "============================="

babel=$(which babel)
if [ $? -eq 0 ]
then
    echo 'Found a babel executable.'
else
    echo 'Fatal: could not find babel in the environment!'
    echo 'Please install openbabel and make sure it is in $PATH'
    exit 1
fi


echo "check if executables are in place..."
echo "===================================="

cd $CADEE_DIR


function fix_q(){
    echo 'Fatal: Please copy Q-executables to ./mutate/executables/q/'
    exit 1
}

if [ ! -f "./mutate/executables/q/qfep5" ]
then
    fix_q
fi

if [ ! "./mutate/executables/q/qfep5" -ef "$qfep_exe" ]
then
    echo 'qscripts.cfg missconfigured:'
    echo 'I found that you are not using the qfep-executable for both qscripts and cadee.'
    echo ''
    echo " Please use $CADEE_DIR/mutate/executables/q/qfep5 as qfep executable for qscripts! "
    echo "(non fatal configuruation problem)"
fi

if [ ! -x "./mutate/executables/q/qdyn5" ]
then
    echo 'Not an executable:'
    echo " $CADEE_DIR/mutate/executables/q/qdyn5"
    fix_q
fi

echo "will now source the init.sh script in ./ensemble:"
echo "================================================="

cd $CADEE_DIR/ensemble

source init.sh || (echo 'There was a problem in script $CADEE_DIR/ensemble/init.sh. You must fix this.'; exit 1)


