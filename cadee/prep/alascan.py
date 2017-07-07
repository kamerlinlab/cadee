#!/usr/bin/env python

""" Create Alanine Scan Inputs

Usage: python alascan.py qprep-wt.pdb wt.fep qprep.inp

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function

import logging
import os
import sys

import qprep5 as qprep5
import tools as tools

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logger = logging.getLogger('mutate.alascan')

# CONSTANTS (easy reading)
ALA = ['N', 'H', 'CA', 'HA', 'CB', 'HB1', 'HB2', 'HB3', 'C', 'O']
MINALA = ['N', 'CA', 'CB', 'C', 'O']


def get_number_mutatable(wtpdb, wtfep, qprep5inp, mutate_radius=None,
                         center_xyz=None, immutable_resids=[]):
    """Prepare AlaScan Inputs.
    @param wtpdb: pdbfile of reference run
    """

    return len(get_mutatable_residues(wtpdb, wtfep, qprep5inp, mutate_radius,
               center_xyz, immutable_resids)[0])


def get_mutatable_residues(wtpdb, wtfep, qprep5inp, mutate_radius=None,
                           center_xyz=None, immutable_resids=[]):
    """parse wt-pdb (the starting pdb) and return mutatable residues
    mutatable residues are within  mutate_radius of center_xyz.
    """

    def round_xyz(xyz):
        """round xyz coordinates for printing"""
        ret = str(round(xyz[0],3)) + " "
        ret += str(round(xyz[1],3)) + " "
        ret += str(round(xyz[2],3)) + " "
        return ret

    msgs = {}

    if (mutate_radius is None) and (center_xyz is None):
        mutate_radius, center_xyz = qprep5.get_sphere_size_center(qprep5inp, wtpdb)  # NOPEP8
        mutate_radius *= 0.85
        msg = 'Looked up radius (scaled by 85%) and center_xyz in qprep5inp -> {0}, {1}'.format(mutate_radius, round_xyz(center_xyz))
    elif mutate_radius is None:
        mutate_radius = qprep5.get_sphere_size_center(qprep5inp, wtpdb)[0]  # NOPEP8
        mutate_radius *= 0.85
        msg = 'Look up radius (scaled by 85%) in qprep5inp -> {0}'.format(mutate_radius)
    elif center_xyz is None:
        center_xyz = qprep5.get_sphere_size_center(qprep5inp, wtpdb)[1]  # NOPEP8
        msg = 'Look up center_xyz in qprep5inp -> {0}'.format(round_xyz(center_xyz))
    else:
        raise 'WTF'

    msgs[msg] = 1

    immutable_fepresids = tools.get_fep_resids(wtpdb, wtfep)

    i = 0
    resids = []

    # TODO: Check if the backbone of the residue is within the cutoff region
    for line in open(wtpdb):
        if not tools.check_qprep_pdb(line):
            continue

        anum, aname, resname, resnum = tools.get_pdb_atom_info(line)
        i += 1
        if i != anum:
            raise Exception("Out of Sync")
        
        if resname.upper() == 'ALA':
            continue

        if resname.upper() == 'GLY':
            continue

        msg = ''
        if int(resnum) in immutable_resids:
            msg = "Won't mutate residue, immutable: {0}".format(resnum)
            msgs[msg] = 1
        elif int(resnum) in immutable_fepresids:
            msg = "Won't mutate residue, contains FEPatoms: {0}".format(resnum)
            msgs[msg] = 1
        if msg:
            logger.debug("%s: %s|%s|%s|%s", msg, anum, aname, resname, resnum)
            continue

        if aname.upper().strip() == 'CA':
            logger.debug('Found CA %s', resnum)
            coords = tools.get_coords(line)
            distance_from_center = tools.euklid_dist(coords, center_xyz)  # NOPEP8
            if distance_from_center < mutate_radius:
                resids.append((resname, str(resnum)))

    resids = sorted(list(set(resids)))
    return resids, msgs


def main(wtpdb, wtfep, qprep5inp, outputfolder,
         mutate_radius=None, center_xyz=None, immutable_resids=[]):
    """Prepare AlaScan Inputs.
    :param wtpdb: pdbfile of reference run
    :param wtfep: fepfile of reference run
    :param qprep5inp: qprep5-input of reference run
    :param outputfolder: to store mutants in
    :param mutate_radius: radius_to_mutate
    :param center_xyz:    center of simulation sphere
    """

    def write_raw_ala_mutant(raw_mut_name, resname2ala, resnum2ala, wtpdb):
        """write output pdb"""

        output = open(raw_mut_name, 'w')

        for line in open(wtpdb):
            line = line[:-1]
            if not tools.check_qprep_pdb(line):
                print(line, file=output)
                continue

            aname, resname, resnum = tools.get_pdb_atom_info(line)[1:4]

            if resname2ala == resname and resnum2ala == str(resnum):
                line = line[:17] + 'ALA' + line[20:]
                if aname.strip() not in ALA:
                    continue  # remove this atom

            print(line, file=output)
        output.close()

    if not os.path.exists(qprep5inp):
        logger.warning('WARNING: inexistent qprep5-input: %s', qprep5inp)
        raise Exception('Qprep-input file does not exist')

    # do actual mutation

    mutatable, msgs = get_mutatable_residues(wtpdb, wtfep, qprep5inp,
                                             mutate_radius,
                                             center_xyz, immutable_resids=[])
    for each in msgs:
        logger.info(each)

    for resid in mutatable:
        resname2ala, resnum2ala = resid

        if resname2ala.upper() == 'ALA':
            logger.warning('%s is ALA, skip', resid)
            continue
        elif resname2ala.upper() == 'GLY':
            logger.warning('%s is GLU, skip', resid)
            continue
        else:
            # previous 2 cases should not happen
            logger.info('Mutate %s to ALA', resid)

        # prepare output folder
        raw_mut_folder = outputfolder + '/' + resname2ala + resnum2ala + 'ALA'
        raw_mut_folder = os.path.abspath(raw_mut_folder)
        raw_mut_name = os.path.abspath(raw_mut_folder + '/raw_mutant.pdb')
        try:
            os.makedirs(raw_mut_folder)
        except OSError:
            logging.info("Skipping %s %s", raw_mut_folder, "exists!")
            continue

        write_raw_ala_mutant(raw_mut_name, resname2ala, resnum2ala, wtpdb)

        try:
            qprep5.create_top_and_fep(qprep5inp, raw_mut_folder,
                                      in_pdb=raw_mut_name, wtpdb=wtpdb,
                                      wtfep=wtfep)
        except KeyError:
            import shutil
            failedfolder = raw_mut_folder+'_fep_failed'
            while os.path.isdir(failedfolder):
                failedfolder += "1"
            shutil.move(raw_mut_folder, failedfolder)
            logger.warning('FEP-generation failed, moved to %s', failedfolder)
            continue


if __name__ == "__main__":
    # Parse Command Line

    def usage():
        """Print usage and exit"""
        print('')
        print('Usage:')
        print('      ' + sys.argv[0] + ' qprep-wt.pdb wt.fep qprep.inp')
        print('')
        sys.exit(1)

    if len(sys.argv) != 4:
        print(sys.argv)
        usage()
    else:
        # Parse Command Line
        main(os.path.abspath(sys.argv[1]),
             os.path.abspath(sys.argv[2]),
             os.path.abspath(sys.argv[3]),
             os.getcwd()+"/ala_scan")
