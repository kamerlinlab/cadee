#!/usr/bin/env python

"""Generate Inputfiles


Author: {0} ({1})

This module is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function

from random import randint
from fractions import gcd

import logging
import os
import sys

from tools import get_atomnumber

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logger = logging.getLogger('mutate.inputgen')

# TO CONVERT A ODS INTO CSV, USE:
#  soffice --headless --convert-to csv test.ods

NLC = '\n'


def lcm(numbers):
    '''kleinstes gemeinsames vielfaches'''
    return reduce(lambda x, y: (x*y)/gcd(x, y), numbers, 1)


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
        return True
    except TypeError:
        return False
    except ValueError:
        return False


def parse_pdbref(pdbref, pdbfilename):
    """
    Parse a pdbref - string:
    END or 'HB2:HIS:299'

    Returns the atom-number of the specified atom.
    """
    if pdbref == '':
        raise (Exception, 'empty line not allowed')

    # seperate / strip off comments
    parts = pdbref.split('!')[0].split('#')[0].split()

    new_keyout = ''

    # split into parts
    for part in parts:
        if part.upper() == 'END':
            # extract number of last atom
            lastline = ''
            for line in open(pdbfilename):
                if line[:4] == 'ATOM':
                    lastline = line
            if lastline == '':
                raise(Exception, 'could not find any atom line in pdbfile')

            if isint(lastline[6:11]):
                part = int(lastline[6:11])
            else:
                errmsg = 'could not get integer for last ATOM-line in pdbfile'
                logger.error(errmsg)
                raise(Exception, errmsg)

        elif part[0] == "(" or not isnum(part):
            part = part.replace(',', ':')
            part = part.replace('.', ':')
            part = part.replace(';', ':')
            part = part.replace('(', '')
            part = part.replace(')', '')
            if len(part.split(':')) != 3:
                print(part.split(':'), part)
                logger.error('Could not split %s to 3 (atomname:residuename:residuenumber)', part)  # NOPEP8
                raise (Exception, 'could not part into 3')
            aname, resn, resi = part.split(':')
            part = get_atomnumber(pdbfilename, aname, resn, resi)
        else:
            pass

        new_keyout = new_keyout + " " + str(part)

    return new_keyout


def walk(pdbfilename, csvfile, folder):
    """ Walk trough subfolders and create inputs"""
    wd = os.getcwd()
    os.chdir(folder)
    for fol in os.listdir('.'):
        os.chdir(folder)
        try:
            os.chdir(fol)
        except:
            # not a folder
            continue
        try:
            if os.path.exists('4670_fep.inp'):
                logger.info('Skipping %s', fol)
                continue
            main(pdbfilename, csvfile)
        except Exception as e:
            print('exception ', e, 'happened in', fol)
            pass
    os.chdir(wd)


def main(pdbfilename, csvfile):
    """Genereate Inputfiles from csvfile"""
    # TODO: Rewrite in Clean
    def add_to_section(section, key, fn, value=None):
        """Add a key/value pair to section"""
        if key is None or key == '':
            if value != '#lambdas':
                print('warning, skip empty key, with value', value)
            return section

        if value is None or value.strip().lower() == 'none' or value == '':
            value = None
            if key[0] != '[':
                logger.debug('# skip empty key {KEY}'.format(KEY=key))
                return section
        elif '$' in value:
            if value == '$seed' or value == '$random':
                value = randint(1, 999999)
            elif '$prev' in value:
                value = value.replace('$prev', prev_name)
            elif '$name' in value:
                value = value.replace('$name', name)

        if '[' not in key[0] or '#' in key[0] or '!' in key[0]:

            # trying to parse atom-descriptor
            if len(key.split('!')[0].split('#')[0].split()) > 1:
                key = parse_pdbref(key, pdbfilename)

            section.append(' {KEY:<20} {VALUE}'.format(KEY=key, VALUE=value))

        else:
            if '[' in key[0]:
                # This is a new section. We have to write all the stuff in
                # 'section' to disk.
                if len(section) > 1:
                    logger.debug("%s, %s", len(section), section[0])
                    for line in section:
                        fn.write(line + NLC)
                    fn.write(NLC)
                    fn.write(NLC)
                elif len(section) == 1:
                    logger.debug('wont write to disk: %s', section)

                section = []
                section.append('{KEY:<20}'.format(KEY=key))
                if value is not None:
                    section.append(' {VALUE}'.format(VALUE=value))
            else:
                # this key is a comment. we will print it as debug msg
                logger.debug('Comment: {KEY:<20} {VALUE}'.format(
                    KEY=key, VALUE=value))

        return section

    print(os.getcwd())

    fep_md_stepsizes = []
    fep_energy_inter = []
    fep_energy_fname = []
    hij = None
    gps = None

    section = []

    name = ''
    prev_name = ''
    for line in open(csvfile, 'r'):
        cells = line.split(',')
        if cells[0].lower() == 'description':
            # create qinputobj
            dsc = 'Description' + NLC + '=========' + NLC
            for cell in cells[1:]:
                if cell != '':
                    dsc += cell + NLC
            # ToDo: create obj

        if cells[0].lower() == 'header':
            # read header
            header = []
            header.append(cells[0])
            for head in cells[1:]:
                if head == "None":
                    head = None
                header.append(head)

        if cells[0].lower() == 'mapping_hij':
            hij = float(cells[1])

        if cells[0].lower() == 'mapping_gps':
            gps = float(cells[1])

        if cells[0].lower() == 'file' or cells[0].lower() == 'fep':
            if name != '':
                prev_name = name
            name = "{0}_{1}".format(cells[1], cells[2])
            fn = open(name + '.inp', 'w')
            section = []
            ready = False
            for i, head in enumerate(header):
                if head == 'steps' and cells[0].lower() == 'fep':
                    fep_md_stepsizes.append(int(cells[i]))
                if head == 'energy' and cells[0].lower() == 'fep':
                    try:
                        fep_energy_inter.append(int(cells[i]))
                    except ValueError:
                        temp = cells[i]
                        # next 5 lines are douplicates (from above)
                        if '$' in temp:
                            if '$prev' in temp:
                                temp = temp.replace('$prev', prev_name)
                            if '$name' in temp:
                                temp = temp.replace('$name', name)
                        fep_energy_fname.append(temp)

                if len(head) == 0:
                    continue
                if head[0] == '[' and ready:
                    section = add_to_section(section, NLC, fn)

                if head == '[MD]':
                    ready = True

                if head == "[Lambdas]":
                    # the column containing 'Lambdas' holds them, so
                    # we write our [Lambdas]
                    section = add_to_section(section, head, fn)
                    # and right after, the value of the cell
                    section = add_to_section(section, cells[i], fn, '#lambdas')
                    continue

                if ready:
                    if cells[i] == 'None':
                        cell = None
                    if cells[i] == '':
                        cell = None
                    else:
                        cell = cells[i]
                    section = add_to_section(section, head, fn, cell)

            # "Flush"; cheap unclean way
            # section = add_to_section('[]')
            # section = []

    if gps is not None and hij is not None:
        generate_mapping_input(gps, hij, fep_energy_inter,
                               fep_md_stepsizes, fep_energy_fname)


def generate_mapping_input(gps, hij, fep_energy_inter, fep_md_stepsizes, fep_energy_fname):  # NOPEP8
    """Generate Mapping Inputs, given gps&hij are not None"""
    # We have to write a mapping input!

    if len(fep_energy_inter) != len(fep_md_stepsizes):
        print('WARNING: Number of energyfiles and fepfiles is inequal.')
        print('         I am unable to write a les anhe same.')

    max_steps = max(fep_md_stepsizes)
    min_steps = min(fep_md_stepsizes)

    if len(set(fep_energy_inter)) != 1:
        print('WARNING: Unequal energy-saving intervals not implemented.')
    else:
        ene_interval = list(set(fep_energy_inter))[0]

    # drop the 1st 10% of the smallest energy file, but at least 1000 steps
    drop_10pz = int(round(min_steps / ene_interval / 10, 1))
    if drop_10pz < 1000:
        drop_10pz = 1000

    logger.debug('10 perz: %s', drop_10pz)

    if max_steps != min_steps:
        # the user chose inequal sampling for different windows.
        # the map-file has to be balanced for this.
        num_ene = []
        for stepsize in fep_md_stepsizes:
            num_ene.append((stepsize/ene_interval)-drop_10pz)

        kgv = lcm(num_ene)
        logger.debug('kgv: %s', kgv)

        enefiles = ''
        num_ene_files = 0
        for i, fil in enumerate(fep_energy_fname):
            reps = kgv / num_ene[i]
            enefiles += (fil + NLC) * reps
            num_ene_files += reps

        # write the mapping-input file
        mappingfile = """{NUM_ENE_FILES}
2 0
.596163 {POINTS_TO_SKIP}
{NUM_BINS}
{MIN_POINTS_PER_BIN}
{GPS}
1
1 2 {HIJ} 0 0 0
1 -1
{ENEFILES}
stop

        """.format(NUM_ENE_FILES=num_ene_files,
                   POINTS_TO_SKIP=drop_10pz,
                   NUM_BINS=2 * len(fep_energy_fname),
                   MIN_POINTS_PER_BIN=int(kgv / min(num_ene) * 50),
                   GPS=gps,
                   HIJ=hij,
                   ENEFILES=enefiles)

        # print(mappingfile)
        open('mapping.map', 'w').write(mappingfile)

    else:
        enefiles = ''
        num_ene_files = len(fep_energy_fname)
        for fil in fep_energy_fname:
            enefiles += fil + NLC

        mappingfile = """{NUM_ENE_FILES}
2 0
.596163 {POINTS_TO_SKIP}
{NUM_BINS}
{MIN_POINTS_PER_BIN}
{GPS}
1
1 2 {HIJ} 0 0 0
1 -1
{ENEFILES}
stop

        """.format(NUM_ENE_FILES=num_ene_files,
                   POINTS_TO_SKIP=drop_10pz,
                   NUM_BINS=2 * len(fep_energy_fname),
                   MIN_POINTS_PER_BIN=int(50),
                   GPS=gps,
                   HIJ=hij,
                   ENEFILES=enefiles)

        # print(mappingfile)
        open('mapping.map', 'w').write(mappingfile)


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage:")
        print(sys.argv[0], 'pdbfilepdb. csv-file.csv')
        sys.exit(1)
    else:
        main(sys.argv[1], sys.argv[2])
