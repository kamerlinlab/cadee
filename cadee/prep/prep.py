#!/usr/bin/env python
"""
CADEE - User - Interface for Mutagenesis

usage: python cadee.py wt.pdb fep.pdb qp.qpinp /libaries --alascan
       for examples, use python cadee.py --help

Author: {0} ({1})

This program is part of CADEE, a framework for
Computer-Aided Directed Evolution of Enzymes

"""


from __future__ import print_function

from random import randint
import argparse
import os
import shutil
import tarfile
import tempfile
import time
import qprep5 as qprep5
import tools as tools
import genseqs as genseqs
import pyscwrl as scwrl
import config as config
import alascan as alascan

import logging

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

VERSION = '0.7.1'

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

NLC = '\n'

START = time.time()

if not tools.check_version():
    time.sleep(5)

# compatibility of raw_input/input for python2/3
try: input = raw_input
except NameError: pass

class cd:
    """Context manager for changing the current working directory
    """
    # http://stackoverflow.com/questions/431684/how-do-i-cd-in-python/13197763#1319776

    # TODO: MOVE TO TOOLS
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def pack_tarballs(parentdir, seeds=1):
    """Pack and create tarball of cadee run"""
    def _find_seed_file():
        """Locate random_seed key in file"""
        for inp in os.listdir('.'):
            if not os.path.isfile(inp):
                continue
            for line in open(inp):
                if 'random_seed' in line:
                    return inp
        return None

    def _pack(tarfil):
        """Add a files to the tarball"""
        tar = tarfile.open(name=tarfil, mode='w')
        for fil in os.listdir('.'):
            tar.add(fil)
        tar.close()

    os.chdir(parentdir)
    for mutant in os.listdir('.'):
        os.chdir(parentdir)
        try:
            os.chdir(mutant)
        except OSError:
            continue
        logger.info('Packing %s', mutant)
        seedfile = _find_seed_file()
        if seedfile is not None:
            tempsf = seedfile + '.tmp'
            for ctr in range(seeds):
                of = open(tempsf, 'w')
                for line in open(seedfile):
                    if 'random_seed' in line:
                        seed = str(randint(1, 999999))
                        line = line.replace(line.split()[1], seed)
                    of.write(line)
                of.close()
                os.remove(seedfile)
                shutil.move(tempsf, seedfile)
                odir = os.path.dirname(os.path.abspath('.'))
                ofil = odir + '/' + str(mutant) + '_' + str(ctr) + '.tar'
                ofil = str(ofil)
                if os.path.exists(ofil):
                    logger.info('File exists! Do NOT overwrite %s', ofil)
                    continue
                logger.info('Pack # %s, Seed: %s', ctr, seed)
                _pack(ofil)
        else:
            if seeds <= 1:
                _pack('../' + mutant + '.tar')
            else:
                logger.warning('Couldnt find "random_seed" in inputs of %s', mutant)
                pass
        
        # Now removing the temporary files from the harddrive ...
        for fil in os.listdir('.'):
            os.remove(fil)
        os.chdir(parentdir)
        os.removedirs(mutant)


def check_int_or_float(value):
    """Check if a value is castable to float or rise ArgumentTypeError"""
    try:
        return float(value)
    except:
        raise argparse.ArgumentTypeError("%s is neither integer nor float" % value)  # NOPEP8


def check_int_oneplus(value):
    """Check if a value is castable to |N (!=0) or raise Argument Type Error"""
    try:
        if int(value) >= 1:
            return int(value)
        else:
            raise argparse.ArgumentTypeError("%s is not >=1" % value)
    except:
        raise argparse.ArgumentTypeError("%s is not integer" % value)


def outputfolder_exists(shortdir, fulldir):
    # type: (str, str) -> None
    logger.warning('The {0} folder exists! Full path:'.format(shortdir))
    logger.warning(fulldir)
    logger.warning('Warning: Do NOT continue, if a simulation is running (i.e. cadee dyn runs in this folder). Else, your simulation may get corrupted.')
    ans = input('\n\n          Do you want to continue (y/N)?').lower()
    ans = ans.lower()
    if ans != 'y':
        logger.fatal('Fatal: {0} folder exists: Please (re)move {1}.'.format(shortdir, outfolder))
        import sys
        sys.exit(1)


def main():
    """ The main - function in this Module """

    def has_ext(fname, expected_ext):
        """has fname extension expected_ext?
        :return True: if extension is equal
        :return False: else
        """
        base, ext = os.path.splitext(fname)
        if ext != expected_ext:
            return False
        return True

    # TODO: load defaults from somewhere
    parser = argparse.ArgumentParser('CADEE: Simpack Preparation.')

    # Minimum Inputfiles needed
    parser.add_argument('wtpdb', action='store', type=argparse.FileType('r'),
                        help='a reference or wildtype pdbfile')
    parser.add_argument('wtfep', action='store', type=argparse.FileType('r'),
                        help='a reference or wildtype fepfile')
    parser.add_argument('qpinp', action='store', type=argparse.FileType('r'),
                        help='the qprep5 inputfile to use')
    parser.add_argument('qplib', action='store',
                        help='path to folder with the libraries for qpinp')

    # Saturate Options
    parser.add_argument('--libmut', action='append', nargs='+', default=None,
                        help="""
                        --libmut 137:SATURATE (20AA)
                        --libmut 137:'CGP' (3AA)
                        --libmut 137:ALL 138:'AG' (20AAx2=40AA)
                        """)

    # Alanine Scan Options
    parser.add_argument('--alascan', action='store_true', default=False,
                        help='perform alanine scan.')
    parser.add_argument('--radius', action='store', type=check_int_or_float,
                        default=None,
                        help='radius to perform alanine scan within, from central atom')
    parser.add_argument('--nummuts', action='store', default=None, type=int,
                        help='make sure that exactly N mutants are created')

    # TODO: check that this works:
    parser.add_argument('--no_wildtype', action='store_false', default=False,
                        help='dont run the wildtype alanine')

    # Input Generator Options
    parser.add_argument('--template', action='store', default=None, type=argparse.FileType('r'),
                        help='Create simpacks based on the provided template.')
    parser.add_argument('--trajcsv', action='store', default=False,
                        help='A traj-csv file to use for inputfile generation. Experimental.')
    parser.add_argument('--numseeds', action='store', default=4,
                        type=check_int_oneplus, help='number of seeds')

    args = parser.parse_args()

    wtpdb = os.path.abspath(args.wtpdb.name)
    if not has_ext(wtpdb, '.pdb'):
        raise argparse.ArgumentTypeError('File must have ".pdb" extension: {0}'.format(wtpdb))
    if os.path.basename(wtpdb) == 'rawmut.pdb':
        raise argparse.ArgumentTypeError("wtpdb must not be named rawmut.pdb (reserved name)")
    wtfep = os.path.abspath(args.wtfep.name)
    if not has_ext(wtfep, '.fep'):
        raise argparse.ArgumentTypeError('File must have ".fep" extension: {0}'.format(wtfep))
    qpinp = os.path.abspath(args.qpinp.name)
    if not has_ext(qpinp, '.qpinp'):
        raise argparse.ArgumentTypeError('File must have name ".qpinp" extension: {0}'.format(qpinp))
    qplib = os.path.abspath(args.qplib)

    if args.alascan and args.radius is None:
        if args.nummuts is None:
            print('You do a alascan without specifying --radius')
        else:
            print('Determining Radius needed to accommodate', args.nummuts, 'mutants. Please wait...', end='')
            immutable_resids = tools.get_fep_resids(wtpdb, wtfep)
            def search(start, divisor):
                for i in range(0, 100):
                    radius = start + (1. * i / divisor)
                    nummut = alascan.get_number_mutatable(wtpdb, wtfep, qpinp, radius, None, immutable_resids)
                    if int(nummut) == int(args.nummuts):
                        return radius
                    elif nummut > int(args.nummuts):
                        return search(radius - 1./divisor, divisor*10)

                raise Exception('could not determine radius. use --radius instead of --nummuts')

            radius = search(0., 1.)

            print ('Done! Radius is', radius)
            args.radius = radius


    if not os.path.isdir(qplib):
        raise Exception('Not a folder', qplib)

    if args.trajcsv:
        args.trajcsv = os.path.abspath(args.trajcsv)

    tempdir = tempfile.mkdtemp()

    odir = os.getcwd()

    with cd(tempdir):
        try:
            qprep5.check_input(qpinp, qplib)
        except qprep5.QprepInpError:
            logger.info('Will attempt to fix qprep5-input file: ')
            qpinp_parts = qpinp.split('.')
            qpinp_parts.insert(-1, 'new')
            qpinp_new = '.'.join(qpinp_parts)
            logger.info(' The new qprep5 input file will be named {0}.'.format(qpinp_new))
            if os.path.exists(qpinp_new):
                logger.critical('Cannot write new qpinp file {0}, exist! Remove it!'.format(qpinp_new))
                logger.fatal('User interaction required.')
                import sys
                sys.exit(1)

            qprep5.check_input(qpinp, qplib, qpinp_new)
            try:
                qprep5.check_input(qpinp_new, qplib)
            except qprep5.QprepInpError:
                logger.critical(' The newly written qpinp file {0} is invalid!'.format(qpinp_new))
                logger.fatal('User interaction required.')
                import sys
                sys.exit(1)

            logger.warn('Adjusted the user-provided qprep5-input file {0}, using new version {1} instead.'.format(qpinp, qpinp_new))
            qpinp=qpinp_new

        qprep5.create_top(qpinp, tempdir, wtpdb=wtpdb)
        out_pdb, out_top, out_fep = qprep5.create_top_and_fep(qpinp,
                                                              tempdir,
                                                              wtpdb=wtpdb,
                                                              wtfep=wtfep)
        # check if reproduced the pdbfile:
        import filecmp
        if not filecmp.cmp(out_pdb, wtpdb, shallow=False):
            raise Exception('Unable to reproduce the provided pdbfile (wtpdb). Please make sure that the pdbfile your provided (%s) was generated with qprep!', wtpdb)

        def is_same_fepfile(file1, file2):
            """check if file1 and file2 are the same fepfiles?

            :param file1: fepfile1
            :param file2: fepfile2
            :type file1: str
            :type file2: str
            :return: bool

            ::info::
                remove comments (text after # and !)
                and ignores whitespace
            """

            with open(file1) as fil1:
                with open(file2) as fil2:
                    for lin1 in fil1:
                        lin1 = lin1.strip()
                        if lin1 != '':
                            lin1 = lin1.replace('#', '!').split('!')[0]
                        try:
                            lin2 = fil2.readline().strip()
                        except ValueError:
                            logger.debug('Different FEPFile-Lenght!')
                            return False
                        if lin2 != '':
                            lin2 = lin2.replace('#', '!').split('!')[0]
                        lin1 = lin1.replace('\t', ' ')
                        lin2 = lin2.replace('\t', ' ')
                        parts1 = lin1.split()
                        parts2 = lin2.split()
                        if parts1 != parts2:
                            logger.debug('FEPFile Lines are not equal: \n%s  \n%s', lin1, lin2)
                            return False
            return True

        # check if we reproduced the fepfile:
        if not filecmp.cmp(out_fep, wtfep, shallow=False):
            if not is_same_fepfile(out_fep, wtfep):
                debugfep = os.path.join(odir, out_fep)
                shutil.move(out_fep, debugfep)
                logger.warning('Unable to reproduce the provided fepfile! provided: %s re-created: %s', wtfep, os.path.join(os.getcwd(), debugfep))

        # TODO: catch all kind of errors and handle them here, tell the user

        shutil.rmtree(tempdir)
    
    outfolder = None

    if args.alascan:
        logger.info('Preparing alascan.')
        radius = args.radius
        outfolder = os.getcwd() + '/' + 'ala_scan'

        if os.path.exists(outfolder):
            outputfolder_exists('ala_scan', outfolder)
        else:
            os.makedirs(outfolder)


        alascan.main(wtpdb, wtfep, qpinp, outfolder, radius)

    elif args.libmut is not None and len(args.libmut) > 0:
        logger.info('Preparing libmut.')
        mutants = []

        immutable = tools.get_fep_resids(wtpdb, wtfep)

        for seperatemut in range(len(args.libmut)):
            for mut in args.libmut[seperatemut]:
                print(mut)
                parts = mut.split(':')
                resid = int(parts[0])
                if resid in immutable:
                    print("Fatal: Can not mutate %s, atoms of residue are in fepfile!", resid)
                    import sys
                    sys.exit(1)
                lname = parts[1].upper()
                aalib = config.SatLibs.get_lib(lname)
                logger.debug('Will mutate residue %s to %s', resid, aalib)
                mutants.append((resid, aalib))

            fasta = genseqs.get_fasta(wtpdb)

            logger.info(fasta)
            logger.info(mutants)

            sequences = genseqs.genseq2(fasta, mutants)

            outfolder = os.getcwd()+'/libmut'

            if os.path.exists(outfolder):
                outputfolder_exists('libmut', outfolder)
            else:
                os.mkdir(outfolder)

            def get_scwrl_seq(wt, new):
                name = ''

                digits = 0
                lenwt = len(wt)
                while lenwt > 0:
                    digits += 1
                    lenwt = int(lenwt/10)

                template = '{WT}{NUM:0'+str(digits)+'d}{MUT}-'
                nseq = list(wt)
                for i in range(len(new)):
                    if wt[i] == new[i]:
                        continue
                    nseq[i] = new[i].upper()
                    name += template.format(WT=wt[i].upper(), NUM=i+1, MUT=new[i].upper())
                if name == '':
                    name = 'wt'
                else:
                    name = name[:-1]
                return name, nseq

            for seq in sequences:

                name, nseq = get_scwrl_seq(fasta, seq)

                dirname = os.path.join(outfolder, name)

                if name == 'wt' and os.path.isdir(dirname):
                    logger.info('skipping wt, exists: %s', dirname)
                    continue

                try:
                    os.mkdir(dirname)
                except OSError:
                    # TODO: proper error handling
                    logger.info('Folder %s exists. Skipping ...', dirname)
                    continue

                logger.info('Working on %s', dirname)

                seq = ''.join(nseq)

                shutil.copy(wtpdb, os.path.join(dirname, os.path.basename(wtpdb)))
                shutil.copy(wtfep, os.path.join(dirname, os.path.basename(wtfep)))
                shutil.copy(qpinp, os.path.join(dirname, os.path.basename(qpinp)))

                with cd(dirname):
                    scwrl.scwrl2(os.path.basename(wtpdb), seq,
                                 os.path.basename(wtfep), os.path.basename(qpinp))
    else:
        logger.info('No parameters provided. Will prepare simpacks from input; "wt".')
        outfolder = os.getcwd() + '/wt'

        if os.path.exists(outfolder):
            logger.fatal('Cannot continue: Folder {0} exists. Please (re)move it.'.format(outfolder))
            import sys
            sys.exit(1)
        else:
            os.makedirs(outfolder+'/wt')
            with cd(outfolder+'/wt'):
                shutil.copy(wtfep, 'mutant.fep')
                qprep5.create_top_and_fep(qpinp, str(os.getcwd()), in_pdb=wtpdb, out_pdb='mutant.pdb')

    if args.trajcsv:
        if outfolder is None:
            print('You provided a trajcsv, but it is not clear which outputfolder to use.')
            print('Please enter the folder here:')
            print('(If you dont the program will terminate.)')
            outfolder = input('Folder:').rstrip()
            if outfolder.strip() == '':
                outfolder = None
                import sys
                sys.exit(1)
            if not os.path.isdir(outfolder):
                print('Directory does not exist:', outfolder)
                outfolder = None
        if os.path.isdir(outfolder):
            os.chdir(outfolder)
            import mutate.inputgen as inputgen
            inputgen.walk('mutant.pdb', args.trajcsv, outfolder)
            pack_tarballs(outfolder, seeds=args.numseeds)
            print('Success! You find your simpacks in', outfolder)
    else:
        if outfolder is None:
            print('If you want to create inputs, I dont know where they are!')
            print('Please enter the folder:')
            print('(Else the program will terminate.)')
            outfolder = input('Folder:').rstrip()
            if outfolder.strip() == '':
                outfolder = None
                import sys
                sys.exit(1)
            if not os.path.isdir(outfolder):
                print('Directory does not exist:', outfolder)
                outfolder = None

        if os.path.isdir(outfolder):
            os.chdir(outfolder)
            import cadee.prep.create_template_based_simpack as inputgen
            inputgen.main(outfolder, args.template.name)
            pack_tarballs(outfolder, seeds=args.numseeds)
            print('Success! You find your simpacks in', outfolder)


if __name__ == "__main__":
    main()
