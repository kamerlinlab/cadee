#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# MIT License
# 
# Copyright (c) 2016  Miha Purg <miha.purg@gmail.com>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
#
# changes the placeholders inside template files ("$593.CA$") to pdb indexes  ("1123")
# takes in three arguments: PDB (after qprep), file containing the placeholders (q_makefep.py generated FEP file, input templates), output filename
# extra keyword that can be used instead of the placeholder is 'LAST.ID' (no explanation needed)

import re
from lib.common import backup_file

PLACEHOLDER_RE = re.compile("\$\S+\.\S+\$")

def convertToIndexes(inputstring, pdbfile, ignore_comments=True):
    atoms = {}

    with open(pdbfile, 'r') as f:
        for line in f.readlines():
            if line.startswith("ATOM") or line.startswith("HETATM"):
                index = line[6:12].strip()
                pdb_id = line[22:26].strip() + "." + line[12:17].strip() 
                atoms[pdb_id] = index

    last_pdbid = pdb_id

    outputstring=""
    for line in inputstring.split("\n"):
        comment=""
        if ignore_comments and "#" in line:
            i = line.index("#")
            line,comment = line[:i], line[i:]

        c = findPlaceholders(line)
        for pid in c:
            pid = pid.strip("$")
            pid2 = pid.replace("LAST.ID", last_pdbid)  
            try:
                padding = (len(pid2)+2  - len(atoms[pid2])) * " "
            except KeyError:
                raise KeyError("Atom '$%s$' does not exist in the pdb structure." % pid2)
            line = re.sub("\$" + pid + "\$", atoms[pid2] + padding, line) 

        outputstring += line + comment + "\n"
    return outputstring
    
def findPlaceholders(inputstring):
    return PLACEHOLDER_RE.findall(inputstring)


if __name__ == "__main__":

    import sys
    import os
    try:
        import argparse
    except ImportError:
        import lib.argparse as argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('pdb', help = 'pdb structure file (created with qprep)')
    parser.add_argument('inp',  help = 'input/fep file containing the placeholders')
    parser.add_argument('out',  help = 'output filename')
    
    if len(sys.argv) < 3:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    
    for k,v in vars(args).iteritems():
        if k in ['inp', 'pdb'] and not os.path.lexists(v):
            print "FATAL! File %s doesn't exist." % v
            sys.exit(1)
    
    inpstr = open(args.inp, "r").read()

    try:
        outstring = convertToIndexes(inpstr,args.pdb)
    except KeyError as e:
        print "FATAL! Exception raised: %s" % str(e)
        sys.exit(1)

    backup = backup_file(args.out)
    if backup:
        print "Backed up '%s' to '%s'" % (args.out, backup)
    open(args.out, 'w').write(outstring)
    print "Created file '%s'..." % args.out

