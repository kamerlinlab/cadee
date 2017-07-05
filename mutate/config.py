#!/usr/bin/env python

""" Configuration Module
Author: {0} ({1})

This module is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function

import logging

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logger = logging.getLogger('mutate.config')

# ERROR/EXIT CODES
ERR_USAGE = 1
ERR_OUTPUTFOLDER_EXISTS = 2
ERR_TOPO_GENERATION_WT = 3
ERR_QPREP5_INEXISTENT = 4
ERR_MKTOP_INEXISTENT = 5
ERR_NO_BABEL = 6

# CONSTANTS
NLC = '\n'


class SatLibs(object):
    """Saturation Libraries"""
    ALL = ['A', 'R', 'N', 'D', 'C', 'E', 'Q', 'G', 'H', 'I', 'L', 'K',
           'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V']
    NDT = ['F', 'L', 'I', 'V', 'Y', 'H', 'N', 'D', 'C', 'R', 'S', 'G']
    SPECIAL = ['C', 'G', 'P']
    HYDROPHOBIC = ['A', 'V', 'I', 'L', 'M', 'F', 'Y', 'W']
    MINUS = ['D', 'E']
    PLUS = ['R', 'H', 'K']
    CHARGED = MINUS[:]
    CHARGED.extend(PLUS)
    NEUTRAL = ['S', 'T', 'N', 'Q']
    POLAR = CHARGED[:]
    POLAR.extend(CHARGED)

    @staticmethod
    def get_lib(name):
        """ check if library with name exists, alternatively check
        if the provided library is fine

        @raises Exception if an AA code is used 2x,
        @raises Exception if an len(name) is <1,
        @raises Exception if invalid AA code is provided.

        """
        # TODO: MAKE SURE THAT THE RESIDUES ARE ACTUALLY CHARGED CORRECTLY,
        #       I.E. ASP => ASH etc... (check with propka)
        name = name.upper().strip()
        if name == 'ALL' or name == 'SATURATE':
            return SatLibs.ALL
        elif name == 'NDT':
            return SatLibs.NDT
        elif name == 'SPECIAL':
            return SatLibs.SPECIAL
        elif name == 'HYDROPHOBIC' or name == 'APOLAR':
            return SatLibs.HYDROPHOBIC
        elif name == 'CHARGED':
            return SatLibs.CHARGED
        elif name == 'CHARGED+' or name == 'POSITIVE':
            return SatLibs.PLUS
        elif name == 'CHARGED-' or name == 'NEGATIVE':
            return SatLibs.MINUS
        elif name == 'POLAR' or name == 'HYDROPHILE':
            return SatLibs.POLAR
        elif name == 'NEUTRAL' or name == 'UNCHARGED':
            return SatLibs.NEUTRAL
        else:
            customlib = []
            for char in name:
                if char in SatLibs.ALL:
                    if char in customlib:
                        raise Exception('1-Letter Code used twice:', char)
                    else:
                        customlib.append(char)
                else:
                    raise Exception('Invalid 1-Letter Code:', char)
            if len(customlib) < 1:
                raise Exception('Invalid: Empty Library', char)
            return customlib


RESRENAME = [
    ["HID", "HIS"],
    ["HIE", "HIS"],
    ["HIP", "HIS"],
    ["ASH", "ASP"],
    ["GLH", "GLU"],
    ["ARN", "ARG"],
    ["LYN", "LYS"],
    ["DNI", "ASP"],
    ["HNI", "HIS"]
]

# residues compatible with scwrl4, only these residue will be scwrl-ed
NATURAL_AA = ['ARG', 'HIS', 'LYS', 'ASP', 'GLU', 'SER', 'THR', 'ASN', 'GLN',
              'CYS', 'GLY', 'PRO', 'ALA', 'VAL', 'ILE', 'LEU', 'MET', 'PHE',
              'TYR', 'TRP']

# USER-DEFINED CONSTANTS
SOLVENT_MOLECULES = ['WAT', 'HOH', 'H2O', 'T3P']
