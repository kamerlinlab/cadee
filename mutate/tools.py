#!/usr/bin/env python

"""Module providing tools for pdb, fep and general inputfile creation

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function
import os
import math
from executables import exe

import logging

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger('mutate.tools')

# ENUMS (easy reading)
ATOMS = 1
BONDS = 2
ANGLES = 3
DIHEDRALS = 4
IMPROPERS = 5

BACKBONE = ['H', 'N', 'O', 'C', 'CA']

NLC = '\n'


def check_version():
    from sys import version_info
    v = int("{0}{1}{2}".format(*version_info[:3]))
    if v < 270:
        logger.warning("Found Python Version %s.%s.%s", *version_info[:3])
        logger.warning('You should use python 2.7 or newer!')
        return False
    else:
        return True


def isint(txt):
    """Return True if @param txt, is integer"""
    try:
        int(txt)
        return True
    except TypeError:
        return False
    except ValueError:
        return False


def isnum(txt):
    """Return True if @param txt, is float"""
    try:
        float(txt)
        if math.isinf(float(txt)) or math.isnan(float(txt)):
            return False
        return True
    except TypeError:
        return False
    except ValueError:
        return False


def bool_log_exists(log):
    """Return true if logfile exists.
       If logfile does not exist: return false and log error message."""
    if not os.path.exists(log):
        logger.info('logfile does not exist: %s', log)
        return False
    return True


def get_executable(program):
    """Locate executable with name program"""
    program = exe.which(program)
    if program is None:
        raise Exception("Executable not found")
    return program


def is_pdb_water_only(pdbfile):
    """bool: only water molecules (HOH, WAT, H2O) in pdbfile"""
    for line in read_pdbatoms(pdbfile):
        if 'HOH' in line:
            continue
        elif 'WAT' in line:
            continue
        elif 'H2O' in line:
            continue
        return False
    return True


def strip_comments(line):
    """strip comments ('! blubb', '# blabb') off a line"""
    line = line.replace('#', '!')
    return line.split('!')[0]


def get_pdb_atom_info(line):
    """Split and read pdb-file line (QPREP-OUT).

     Extract and return:
        atomnumber, atomname, resname, resid)
    """
    # UNIT: check if the correct numbers are extracted
    atomnum = int(line[6:11])
    atmname = line[13:17].strip()
    resname = line[17:20].strip()
    resinum = int(line[22:27])

    return (atomnum, atmname, resname, resinum)


def coords_of_atomnr(pdbfile, atomnr):
    """Input: pdbfile, atomnr
       Output: [x,y,z] (list of floats)"""
    # UNIT: check if the correct numbers are extracted
    for line in open(pdbfile):
        line = line.upper()
        if line[:4] == 'ATOM' or line[:6] == 'HETATM':
            num = int(line[6:11])
            if num == atomnr:
                return get_coords(line)


def get_pdb_atom(line):
    """Split and read pdb-file line (QPREP-OUT).

     Extract and return:
       atomnumber, code
    with:
       code = atomname - residue name - residue number"""
    # UNIT: check if the correct numbers are extracted
    try:
        atomnum = int(line[6:11])
        atmname = line[13:17].strip()
        resname = line[17:20].strip()
        resinum = int(line[22:27])
    except:
        logger.fatal('Error with line: %s', line)
        raise

    code = atmname+"-"+resname+"-"+str(resinum)
    return (atomnum, code)


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


def get_coords(line):
    """Input: pdbline
       Output: [x,y,z] (list of floats)"""
    x_val = line[30:38]
    y_val = line[38:46]
    z_val = line[46:54]
    return [float(x_val), float(y_val), float(z_val)]


def is_backbone(line):
    """Is ATOM-Name a backbone atom?"""
    return line[13:17].strip() in BACKBONE


def is_hydrogen(line):
    """Check if (pdb)-line is Hydrogen-Record.

       If it's a proper-pdb line (77+ characters, extract element)
       If it's a Qprep-pdb line (<54 chars, take from atomname)
    """
    if len(line) >= 78:
        return bool(line[77] == 'H')

    line = line.strip()
    if len(line) < 54:
        aname = line[13:17].strip()
        return bool(aname[0] == 'H')


def read_pdbatoms(pdbfile):
    """Read pdbfile and return list of ATOM entries"""
    pdblines = []
    for line in open(pdbfile, 'r'):
        if line[:4] == "ATOM":
            pdblines.append(line)
    return pdblines


def euklid_dist(coord1, coord2):
    """
    Return euklidian distance between atoms coord1 and coord2.
    """
    xd2 = (coord1[0]-coord2[0]) ** 2
    yd2 = (coord1[1]-coord2[1]) ** 2
    zd2 = (coord1[2]-coord2[2]) ** 2
    return (xd2 + yd2 + zd2)**0.5


def get_atomnumber(pdbfile, g_atomname, g_resname, g_resid):
    """Return AtomNumber of Atom
    @params:
        pdbfile:    string,  path to file
        g_atomname: string,  atom-name
        g_resname:  string,  residue-name;
                      None,  ignore resname, find resid:atomname
        g_resid:    string,  residue-number
                        -1,  ignore resid, find atomname/resname

    """

    """
    PDB SPECS: http://deposit.rcsb.org/adit/docs/pdb_atom_format.html
     1 -  6        Record name     "ATOM  "
     7 - 11        Integer         Atom serial number.
    13 - 16        Atom            Atom name.
    17             Character       Alternate location indicator.
    18 - 20        Residue name    Residue name.
    22             Character       Chain identifier.
    23 - 26        Integer         Residue sequence number.
    27             AChar           Code for insertion of residues.
    31 - 38        Real(8.3)       Orthogonal coordinates for X in Angstroms.
    39 - 46        Real(8.3)       Orthogonal coordinates for Y in Angstroms.
    47 - 54        Real(8.3)       Orthogonal coordinates for Z in Angstroms.
    55 - 60        Real(6.2)       Occupancy.
    61 - 66        Real(6.2)       Temperature factor (Default = 0.0).
    73 - 76        LString(4)      Segment identifier, left-justified.
    77 - 78        LString(2)      Element symbol, right-justified.
    79 - 80        LString(2)      Charge on the atom.
    """
    debug = False
    g_atomname = g_atomname.upper()
    if g_resname is not None:
        g_resname = g_resname.upper()
    g_resid = int(g_resid)
    # counters, for debugging
    resn = 0
    resi = 0
    for line in open(pdbfile):
        line = line.upper()
        record = line[0:6].strip()

        if record == "ATOM" or record == "HETATM":
            if debug:
                logger.debug("!have an atom entry")

            atomnumber = int(line[6:11])
            atomname = line[12:17].strip().upper()  # nonstandard mod
            resname = line[17:20].strip().upper()
            # chain = line[21].upper()
            resid = int(line[22:26])

            if g_resname is None or resname == g_resname:
                logger.debug("!we have the residue name")
                resn += 1
                if resid == g_resid or g_resid == -1:  # changed, 2015-Mar-14
                    if debug:
                        logger.debug("!we have the residuenumber")
                        logger.debug("!search for atom %s", g_atomname)
                        logger.debug("!current atom %s", atomname)
                    resi += 1
                    if atomname == g_atomname:
                        if debug:
                            logger.debug("!we found the atom!!!")
                        return atomnumber

    # atom was not found in pdb:(
    error = "Atom {0}, {1}, {2} not found in {3}".format(g_atomname, g_resname,
                                                         g_resid, pdbfile)
    logger.exception(error)
    raise Exception(error)


def rename_pdb_res(residue, new_resname):
    """Rename a pdb-residue
    @param:         residue, fullresidue
    @new_resname:   3-character-string with name of new residue

    @return: residue with new_resname"""

    if len(new_resname) != 3:
        raise Exception('PDB-Residue must be exactly 3 characters.')

    for i in range(0, len(residue)):
        prefix = residue[i][0:17]
        suffix = residue[i][20:]
        residue[i] = prefix + new_resname + suffix

    return residue


def check_qprep_pdb(line):
    """Check if line is ATOM or HETATM entry.
    @return True if so
    @return False if TER or GAP entry
    @raise  Exception, if other entry"""
    if len(line) > 55:
        raise Exception("Not a Qprep5-ed pdbfile (line too long)", line)
    if line[:6] == 'ATOM  ' or line[:6] == 'HETATM':
        return True
    if line[:3] == 'TER' or line[:3] == 'GAP':
        return False

    raise Exception("Not Qprep5-ed pdbfile. Neither ATOM, HETATM nor TER, GAP entry.")  # NOPEP8


def get_ranges(indexes):
    """Form start-end list, made form indexes:
    eg. indexes=[1,2,3,4,10,11]:
        return [ [1,4], [10,11] ]
    """

    ## UNITTEST: in [2540, 2541, 2542, 2543, 2544, 2545, 2546, 2547, 2548, 2549, 7543, 7544, 7545, 7546, 7547, 7548, 7549, 7550, 7551, 7552, 7553, 7554, 7555, 7556, 7557]  # NOPEP8
    ##           out: [[2540, 2549], [7543, 7557]]

    ranges = []
    previdx = None
    start = None
    for idx in sorted(indexes):
        idx = int(idx)
        if previdx is None or previdx != idx-1:
            if previdx is not None:
                ranges.append([start, previdx])
            start = idx
        previdx = idx

    if start is not None:
        ranges.append([start, idx])

    return ranges


def get_non_backbone_atoms(wtpdb, resids):
    """Return sorted list of atom-numbers that are part of resids and not backbone
    
    WARNING: BACKBONE IS NOT INCLUDING 'CA'

    :param wtpdb: wild-type pdbfile
    :param resids: residue numbers in pdbfile to search
    :type wtpdb: str, path to file
    :type resids: list of integers
    """

    indexes = []
    for line in open(wtpdb):
        if line[:6] == 'ATOM  ' or line[:6] == 'HETATM':
            resid = int(line[22:26].strip())
            if resid in resids:
                if not line[13:17].strip() in ['H', 'N', 'O', 'C']:
                    serial = int(line[6:11])
                    indexes.append(serial)
    return indexes



def get_fep_atom_pdbindexes(wtfep):
    """Return list of FEP-atoms

    Parse fepfile and pdbfile to find residue numbers of FEP-atoms.
    
    :param wtpdb: wild-type pdbfile
    :param wtfep: wild-type fepfile
    :type wtpdb: str
    :type wtfep: str
    :return: list of atoms that are in fepfile
    """
    s_is_at = False
    indexes = []   
    for line in open(wtfep):
        line = line.replace('#', '!').split('!')[0].strip()
        if line == '':
            continue
        if line.lower() == '[atoms]':
            s_is_at = True
            continue
        elif '[' in line:
            s_is_at = False
        elif s_is_at:
            qnum, pdbnum = line.split()
            indexes.append(int(pdbnum))

    return indexes        


def get_fep_resids(wtpdb, wtfep):
    """Return list of residue-numbers containing FEP-atoms

    Parse fepfile and pdbfile to find residue numbers of FEP-atoms.

    :param wtpdb: wild-type pdbfile
    :param wtfep: wild-type fepfile
    :type wtpdb: str
    :type wtfep: str
    :return: list of residues that are in fepfile
    """

    indexes = get_fep_atom_pdbindexes(wtfep)

    fep_resids = []
    for line in open(wtpdb):
        if line[:6] == 'ATOM  ' or line[:6] == 'HETATM':
            serial = int(line[6:11])
            if serial in indexes:
                resid = int(line[22:26].strip())
                fep_resids.append(resid)

    return list(set(fep_resids))


def get_last_pdb_atom_number(pdb):
    """Find in last COORDINATE line and return its atomnumber
    """
    for line in open(pdb):
        if line[:6] == 'ATOM  ' or line[:6] == 'HETATM':
            lastline = line
    if lastline == '':
        raise(Exception, 'could not find any atom line in pdbfile')

    if isint(lastline[6:11]):
        index = int(lastline[6:11])
    else:
        errmsg = 'could not get integer for last ATOM-line in pdbfile'
        logger.error(errmsg)
        raise(Exception, errmsg)
    return index


def get_solute_and_solvent_ranges(pdb):
    """Find in pdb solute and solvent ranges
    return [ [1, end_solute], [startwat, lastatom] ]
    """
    for line in open(pdb):
        if line[:6] == 'ATOM  ' or line[:6] == 'HETATM':
            if line[17:20] == 'HOH':
                startwat = int(get_pdb_atom_info(line)[0])
                endsolute = startwat-1
                break
    solute = [1, endsolute]
    solvent = [startwat, get_last_pdb_atom_number(pdb)]
    return [solute, solvent]
