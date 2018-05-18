#!/bin/bash

# This script is very quickly testing, if the executables seem to work.
# It will find out, when the tools have wrong names, for example.

echo "" | ./Qcalc6 2>&1 | grep -q 'qcalc'
if [ $? -ne 0 ]
then
    echo "Qcalc6 is invalid"
    exit 1
fi

echo "" | ./Qfep6 2>&1 | grep -q '# Qfep'
if [ $? -ne 0 ]
then
    echo "Qfep6 is invalid"
    exit 1
fi

echo "" | ./Qprep6 2>&1 | grep -q 'Qprep>'
if [ $? -ne 0 ]
then
    echo "Qprep6 is invalid"
    exit 1
fi

echo "" | ./Qdyn6 2>&1 | grep -q "ABNORMAL TERMINATION of Qdyn5"
if [ $? -ne 0 ]
then
    echo "Qdyn6 is invalid"
    exit 1
fi


echo "" | ./Qdyn6p 2>&1 | grep -q "MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD"
if [ $? -ne 0 ]
then
    echo "Qdyn6p is invalid"
    echo Qcalc6, Qfep6, Qprep6 and Qdyn6 are valid.
else
    echo Qcalc6, Qfep6, Qprep6, Qdyn6 and Qdyn6p are valid.
fi
