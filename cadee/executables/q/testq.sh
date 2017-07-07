echo "" | ./qcalc5 2>&1 | grep -q 'qcalc'
if [ $? -ne 0 ]
then
    echo "qcalc5 is invalid"
    exit 1
fi

echo "" | ./qfep5 2>&1 | grep -q '# Qfep'
if [ $? -ne 0 ]
then
    echo "qfep5 is invalid"
    exit 1
fi

echo "" | ./qprep5 2>&1 | grep -q 'Qprep>'
if [ $? -ne 0 ]
then
    echo "qprep5 is invalid"
    exit 1
fi

echo "" | ./qdyn5 2>&1 | grep -q "ABNORMAL TERMINATION of Qdyn5"
if [ $? -ne 0 ]
then
    echo "qdyn5 is invalid"
    exit 1
fi


echo "" | ./qdyn5p 2>&1 | grep -q "MPI_ABORT was invoked on rank 0 in communicator MPI_COMM_WORLD"
if [ $? -ne 0 ]
then
    echo "qdyn5p is invalid"
    echo qcalc5, qfep5, qprep5 and qdyn5 are valid.
else
    echo qcalc5, qfep5, qprep5, qdyn5 and qdyn5p are valid.
fi
