#!/usr/bin/env python

"""Module to locate Clashes in a pdbfile.

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function

import logging
import sys

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logger = logging.getLogger('mutate.clash')

BACKBONE = ['H', 'N', 'O', 'C', 'CA']
IGNORE_BB = True
IGNORE_H = True

RADIUS = 2.5   # Default: 2.5 A.


def clash_score(distance):
    """
    Return clash_score offset equlibrium distance^6 (distance^6).
    """
    if distance >= RADIUS:
        return 0
    ooe = (RADIUS-distance)
    return ooe**6


def clashscore_and_residues(resnums, pdblines):
    """
    Input:          resids: List of residues
    Input:          pdblines: List of ATOM-records OR filename

    Output: List of residues clashing with those indicated in resids
    """

    logfile = open('clashscore', 'a')

    if not isinstance(pdblines, list):
        import tools
        pdblines = tools.read_pdbatoms(pdblines)

    if isinstance(resnums, (int, str)):
        resnums = [int(resnums)]

    score = 0.
    clashes = []
    for resnum in resnums:
        incr_score, new_clashes = clashes_of_resnum(pdblines,
                                                    resnum, logfile=logfile)
        score += incr_score
        clashes.extend(new_clashes)

    clashes = list(set(clashes))

    logfile.write(str(clashes))
    logfile.write(str(score))

    return score, clashes


def get_resnum_coords(resnum, pdblines):
    """
    input: resnum: residue number
           pdblines: list of ATOM entries
    return: list of [x,y,z] coordinate
    """
    import mutate.tools as tools
    coords_a = []
    anames_a = []
    for line in pdblines:
        try:
            rnum = line[22:28]
            rnum = int(rnum)
            if IGNORE_BB and tools.is_backbone(line):
                continue
            if IGNORE_H and tools.is_hydrogen(line):
                continue
            if rnum == resnum:
                aname = line[13:17].strip()
                coords_a.append(tools.get_coords(line))
                anames_a.append(aname)
        except Exception as e:
            logger.fatal('Fatal: Exception in line %s', line)
            for line in pdblines:
                logger.fatal(line.strip())
            logger.fatal('Exception was: %s, %s, %s', e, e.args, e.message)
            raise e
    return coords_a, anames_a


def clashes_of_resnum(pdblines, resnum, logfile=None):
    """
    @param pdblines: list of pdbfile ATOM-records
    @param resnum:   integer

    @return: score, residues
    """
    import mutate.tools as tools

    coords_a, anames_a = get_resnum_coords(resnum, pdblines)
    clashing_resids = []

    score = 0.
    for line in pdblines:
        if IGNORE_H and tools.is_hydrogen(line):
            continue
        coord_b = tools.get_coords(line)
        for coord in coords_a:
            dist = tools.euklid_dist(coord, coord_b)
            if tools.euklid_dist(coord, coord_b) < RADIUS:
                if int(line[22:28]) != resnum:
                    score += clash_score(dist)
                    logger.debug('found clash: %s with resid: %s atom: %s dist: %s', line[:28], resnum, anames_a[coords_a.index(coord)], round(dist, 3))  # NOPEP8
                clashing_resids.append(int(line[23:26]))
    return score, list(set(clashing_resids))


def main(pdbfile, resnums):
    """Return clash-scores of resnums in pdbfile
    @param pdbfile: path to pdbfile
    @param resnums: int or list of ints
    """
    import mutate.tools as tools
    print(clashscore_and_residues(resnums, tools.read_pdbatoms(pdbfile)))

if __name__ == "__main__":
    def usage():
        """ Print usage info and exit"""
        print("Usage:")
        print(sys.argv[0], ' pdbfile residue-number ')
        sys.exit()

    if len(sys.argv) != 3:
        usage()
    main(sys.argv[1], sys.argv[2])
