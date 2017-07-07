#!/usr/bin/env python

"""
Generate Sequence from a pdbfile and to modify the squences.

Author: {0} ({1})

This module is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function

import logging
import os
import sys
import time

import config

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logger = logging.getLogger('mutate.genseqs')

# ERROR/EXIT CODES
ERR_USAGE = 1
ERR_OUTPUTFOLDER_EXISTS = 2
ERR_TOPO_GENERATION_WT = 3
ERR_QPREP5_INEXISTENT = 4
ERR_MKTOP_INEXISTENT = 5
ERR_NO_BABEL = 6

# CONSTANTS
NLC = '\n'


def genseq2(wtseq, mutations, keepdupes=False):
    """ generate a sequences library based of wtseq
    @param: list of tupel, [ (resid, library), (resid, library), ...]

    @returns: list of sequences
    """

    def estimator(mutations):
        est = 1
        for mut in mutations:
            lib = mut[1]
            est *= (len(lib)+1)
        return est

    logger.info('will mutate wtseq %s and create about %s mutations',
                wtseq, estimator(mutations))

    seqo = list(wtseq)
    sequences = [seqo]
    while len(mutations) > 0:
        newseqs = sequences[:]
        res, lib = mutations.pop()
        for seqo in sequences:
            res = int(res)
            if res < 1:
                raise ValueError('Impossible: resid < 1!', res)
            pos = res - 1
            for aa in lib:
                if len(aa) != 1:
                    raise ValueError('Impossible 1-letter aminoacid',
                                     aa, 'in lib', lib)
                seqn = seqo[:]
                seqn[pos] = aa
                if keepdupes or seqn not in newseqs:
                    newseqs.append(seqn)
        sequences = newseqs

    return sequences


def combine(lib, pos):
    """generate combinations of up to 7.
    @param lib: library
    @param pos: positions to mutate
    # TODO: implement in readable (recursively)
    """
    numseqs = 1
    for each in lib:
        numseqs *= len(each)
    logger.info('Generating %s %s', numseqs, 'sequeces. Please wait.')
    seqlib = []

    logger.info('Library %s, Positions %s', lib, pos)

    for every in lib[0]:
        if len(pos) > 1:
            for every2, in lib[1]:
                if len(pos) > 2:
                    for every3, in lib[2]:
                        if len(pos) > 3:
                            for every4, in lib[3]:
                                if len(pos) > 4:
                                    for every5, in lib[4]:
                                        if len(pos) > 5:
                                            for every6, in lib[5]:
                                                if len(pos) > 6:
                                                    for every7 in lib[6]:
                                                        seqlib.append([every,
                                                                       every2,
                                                                       every3,
                                                                       every4,
                                                                       every5,
                                                                       every6,
                                                                       every7])
                                                else:
                                                    seqlib.append([every,
                                                                   every2,
                                                                   every3,
                                                                   every4,
                                                                   every5,
                                                                   every6])
                                        else:
                                            seqlib.append([every,
                                                           every2,
                                                           every3,
                                                           every4,
                                                           every5])
                                else:
                                    seqlib.append([every, every2, every3,
                                                   every4, every4])
                        else:
                            seqlib.append([every, every2, every3])
                else:
                    seqlib.append([every, every2])
        else:
            seqlib.append([every])

    return seqlib


def gen_seqlib(sequence, pos, lib):
    """
        Generates sequences, mutating at pos[x] to all as in lib[x]
        Generates sequences, mutating at pos[x] if len(lib)==1,
        the same lib will be used for all
        Return sequences
    """
    # is lib a string?
    if isinstance(lib, str):
        lib = [lib]

    # when only 1 library is given, reuse it
    if len(lib) == 1:
        while range(1, len(pos)):
            lib.append(lib[0])

    if len(pos) != len(lib):
        msg = 'Bad Input: Dimensions of pos and lib must be equal: '
        msg += 'found: #pos: {0}, #lib {1}'.format(len(pos), len(lib))
        raise (Exception, msg)

    seqlib = combine(lib, pos)

    # insert combinations into sequence
    sequences_1d = {}
    for i in range(0, len(seqlib)):
        nfa = list(sequence)
        for j, posj in pos:
            if nfa[posj].upper() != seqlib[i][j].upper():
                nfa[posj] = seqlib[i][j]
        modseq = ''.join(nfa)
        sequences_1d[modseq] = 1

    return sequences_1d


def get_fasta(wtpdb):
    """Return fasta code of wtpdb"""

    # preparations

    from pyscwrl import babel_pdb_for_scwrl

    babel_pdb_for_scwrl(wtpdb)

    # read fasta
    fasta = ''
    for line in open('proper.fasta'):
        line = line[:-1]
        if line[0] == '>':
            # fasta-comment, ignore line
            continue
        for char in line:
            fasta += char.lower()

    return fasta


def get_sequences(wtpdb, resids, library):
    """Return list of sequences for resids, created with library"""
    print(wtpdb, resids)
    # Get the fasta sequence from pdbfile
    fasta = get_fasta(wtpdb)

    posids = []
    # position - ids start from 0 (not 1), so we have to convert
    for resid in resids:
        posids.append(int(resid)-1)

    # generate sequences:
    sequences = gen_seqlib(fasta, posids, [library])

    return sequences

if __name__ == "__main__":
    # Parse Command Line
    LIB = config.SatLibs.ALL

    def usage():
        """Print Usage and exit"""
        print('')
        print('Usage:')
        print('      ' + sys.argv[0] + ' qprep-wt.pdb  res1 [ res2 ...] ]')
        print('')
        sys.exit(ERR_USAGE)

    def get_resnumbers(args):
        """Return residue-numbers as list-of-integers"""
        resids = []
        for resid in args:
            try:
                resids.append(int(resid))
            except ValueError:
                print('ValueError with ', resid, ' expected: Integer')
                usage()
        if len(resids) > 7:
            print('FATAL:')
            print('You ask me to mutate more than 7 residues at one time.')
            print('This is NOT IMPLEMENTED...    ...probably a BAD IDEA :')
            print('This is a bad idea, because we grow with LIBRARY^{#RES}!')
            print('In your case ', len(LIB), '^', len(LIB), '=',
                  len(LIB)**len(resids), '!')
            usage()
        return resids

    START = time.time()

    if len(sys.argv) < 3:
        usage()

    if len(get_resnumbers) > 7:
        usage()

    get_sequences(os.path.abspath(sys.argv[1]),
                  get_resnumbers(sys.argv[2:]), LIB)

    print('time', round(time.time()-START, 2), 's')
