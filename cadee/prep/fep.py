#!/usr/bin/env python

"""Module providing methods to create a new FEPfile

Author: {0} ({1})

This module is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function

import logging
import sys

from tools import get_pdb_atom, check_qprep_pdb, read_pdbatoms

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logger = logging.getLogger('prep.fep')

# ENUMS (easy reading)
ATOMS = 1
BONDS = 2
ANGLES = 3
DIHEDRALS = 4
IMPROPERS = 5


def get_section(line):
    """Assign section by reading FEP-File line with '['"""
    if "[atoms]" in line.lower():
        return ATOMS
    elif "[change_bonds]" in line.lower():
        return BONDS
    elif "[change_angles]" in line.lower():
        return ANGLES
    elif "[change_torsions]" in line.lower():
        return DIHEDRALS
    elif "[change_impropers]" in line.lower():
        return IMPROPERS
    else:
        return None


def _rewrite_fep_section(wtlist, mutdict, line, output, section):
    """Rewrite FEPFILE line
    @param wtlist:  list of atoms in wtpdb
    @param mutdict: dict of atoms in mutantpdb
    @param line:    line of wtfep
    @param output:  file or None, used to print output to
    @param section: enum (integer)
    """
    # strip off comment in line
    if "!" in line[1:]:
        (line, comment) = line.split("!", 1)
    else:
        comment = ""

    parts = line.split()
    if section == ATOMS:
        try:
            parts[1] = mutdict[wtlist[int(parts[1])]]
        except KeyError:
            logger.error('FEPFile/PDBFile missmatch:')
            logger.error('Could not find %s in wildtypepdb!',
                         wtlist[int(parts[1])])
            raise
        except IndexError:
            logger.error('FEPFile/PDBFile missmatch:')
            logger.error('Could not find atomnr %s in mutated pdbfile.',
                         int(parts[1]))
            raise
        except ValueError:
            logger.error('FEPFile/PDBFile missmatch:')
            logger.error('An error happened, (%s is not a qatom-number), while reading FEPfile-line: %s', parts[1], line)  # NOPEP8
            raise
        except Exception:
            logger.error('FEPFile/PDBFile missmatch:')
            logger.error('Could not look up FEP atoms: %s, %s, %s, %s',
                         parts,
                         parts[1],
                         int(parts[1]),
                         wtlist[int(parts[1])])
            raise

    elif section == BONDS:
        parts[0] = mutdict[wtlist[int(parts[0])]]
        parts[1] = mutdict[wtlist[int(parts[1])]]

    elif section == ANGLES:
        parts[0] = mutdict[wtlist[int(parts[0])]]
        parts[1] = mutdict[wtlist[int(parts[1])]]
        parts[2] = mutdict[wtlist[int(parts[2])]]

    elif section == DIHEDRALS or section == IMPROPERS:
        parts[0] = mutdict[wtlist[int(parts[0])]]
        parts[1] = mutdict[wtlist[int(parts[1])]]
        parts[2] = mutdict[wtlist[int(parts[2])]]
        parts[3] = mutdict[wtlist[int(parts[3])]]

    for part in parts:
        print(str(part), end="\t", file=output)

    if comment != "":
        print("!", end="", file=output)

    print(comment, file=output)


def create_fep(wtpdb, wtfep, mutpdb, outfep=None):
    """Create FEP file with wildtype pdb, wildtype fep, and mutant pdb.
    NOTE: Both wtpdb and mutpdb must be output of Qprep5.

    If outfep is None, the fepfile is printed to stdout.
    """

    wtlist = [0]
    mutdict = {}

    # open output file
    if outfep is None:
        output = None
    else:
        output = open(outfep, "w")

    # parse wt-pdb (the starting pdb)
    i = 0
    for line in read_pdbatoms(wtpdb):
        check_qprep_pdb(line)
        anum, code = get_pdb_atom(line)
        i += 1
        if i != anum:
            raise Exception("Out of sync: wt-pdb atom numbers are out of order.")  # NOPEP8
        wtlist.append(code)

    # parse mutant pdb
    for line in read_pdbatoms(mutpdb):
        check_qprep_pdb(line)
        anum, code = get_pdb_atom(line)
        mutdict.update({code: anum})

    # Parse and Rewrite the fep file.
    section = None
    for line in open(wtfep):

        line = line[:-1]    # strip off \n
        line = line.replace('#', '!')   # unify coment str

        if "[" in line:
            section = get_section(line)

        if len(line.split()) < 2:
            print(line, file=output)
            continue

        if line[0] == "!":
            print(line, file=output)
            continue

        if section is None:
            print(line, file=output)
        else:
            try:
                _rewrite_fep_section(wtlist, mutdict, line, output, section)
            except KeyError as err:
                logger.error('Error creating the (mutated) fepfile for %s', mutpdb)        # NOPEP8
                logger.error('CADEE was unable to find the atom named "%s"', err.message)  # NOPEP8
                logger.error('              Please fix your fepfile and rerun.')           # NOPEP8
                raise


if __name__ == "__main__":
    # Parse Command Line

    def usage():
        """Print Usage and Exit"""
        print('')
        print('Usage:')
        print('      '+sys.argv[0]+' qprep-wt.pdb wt.fep qprep-mutant.pdb [mutant.fep]')  # NOPEP8
        print('      Optional:')
        print('               [mutant.fep] is optional. if neglected, print to stdout.')  # NOPEP8
        print('')
        sys.exit()

    if len(sys.argv) < 4 or len(sys.argv) > 5:
        usage()

    else:
        if len(sys.argv) == 5:
            create_fep(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
        else:
            create_fep(sys.argv[1], sys.argv[2], sys.argv[3])
