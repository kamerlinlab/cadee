#!/usr/bin/env bash
find . -name "*.pyc" -exec rm {} \;
rm -f rawmut.pdb
rm -f mutant.???
rm -f qprep.???
rm -f proper.???
rm -f proper.fasta

if [ -e 'libmut' ]
then
    echo -n "Delete libmut folder (y/N)? >"
    read ans
    [ $ans = 'y' ] && rm -r libmut && echo "deleted!"
fi

if [ -e 'ala_scan' ]
then
    echo -n "Delete ala_scan folder (y/N)? >"
    read ans
    [ $ans = 'y' ] && rm -r ala_scan && echo "deleted!"
fi
