#!/usr/bin/env python

"""Merge PDB-Files

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function

import logging
import sys
import time

import tools as tools
from config import NLC, ERR_USAGE, NATURAL_AA

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logger = logging.getLogger('prep.mergepdb')


# not used, using babel instead
def get_sequence(pdbfile):
    """Return Aminoacid Sequence as list"""
    seq = []
    oresnum = -1
    for line in open(pdbfile, 'r'):
        if line[:4] != 'ATOM':
            continue
        resnum = int(line[22:26])
        if resnum != oresnum:
            seq.append(line[17:20])
        oresnum = resnum
    return seq


def mergepdb(oldq, scwrlpdb, newfile, pos):
    """
    INPUT:
        oldq: pdbfile
        scwrlpdb: pdbfile
        pos: list with positions that were mutated

    OUTPUT: merged oldq/scwrlpdb

    LIMITATIONS:
        ONLY WORKS WHEN NOT_RESIDUES ARE AT THE END OF THE PDBFILES.
    """

    ores = {}
    postfix = []
    for line in tools.read_pdbatoms(oldq):

        # if line[11:14].strip()=='H':
        #     print('skip proton')
        #     continue

        resi = line[22:27]
        # resi=int(resi)

        if int(resi)-1 in pos:    # mutated residue
            resn = None
        else:
            resn = line[17:20]
        if resn not in NATURAL_AA:
            postfix.append(line[:-1])

        if not line[12:15] == ' N ':
            continue
        coords = line[31:54]
        ores[coords] = resn

    resnumnam = None
    newfile = open(newfile, 'w')
    for line in tools.read_pdbatoms(scwrlpdb):

        # kill hydrogen:
        if line.split()[-1] == "H":
            continue

        if line[12:15] == ' N ':
            coords = line[31:54]
            resnumnam = line[17:27]
            resn = ores[coords]
        if resnumnam == line[17:27] and (resn is not None):
            line = line[:17] + resn + line[20:-1]
        else:
            pass
        print (line.rstrip(), file=newfile, end=NLC)

    for fix in postfix:
        print (fix.rstrip(), file=newfile, end=NLC)

    newfile.close()

    # os.system('vimdiff {0} {1} {2}'.format(oldq, scwrlpdb, newfile))

if __name__ == "__main__":
    # Parse Command Line
    def usage():
        """Print Usage and Exit"""
        print('')
        print('Usage:')
        print('      ' + sys.argv[0] + ' oldq.pdb nowq.pdb output.pdb rs1 [ res2 ... ]')  # NOPEP8
        print('')
        sys.exit(ERR_USAGE)

    def get_residues():
        """Return Arg4+ as list of integers"""
        positions = []
        for pos in sys.argv[4:]:
            try:
                positions.append(int(pos))
            except ValueError:
                usage()
            except TypeError:
                usage()

    if len(sys.argv) < 5:
        usage()

    START = time.time()

    if len(sys.argv) == 3:
        mergepdb(sys.argv[1], sys.argv[2], None, get_residues())
    else:
        mergepdb(sys.argv[1], sys.argv[2],
                 open(sys.argv[3], 'w'), get_residues())

    print('time', time.time()-START)
