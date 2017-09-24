#!/usr/bin/env python
"""
This module checks qprep5 input for errors and problems.

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function

import logging
import os
import shutil

import tools as tools
from fep import create_fep

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logger = logging.getLogger('prep.qprep5')


NLC = '\n'


def check_input(qprep5_inp, qprep5_lib, outfile=False):
    """
    Check if qprep5 input file is ready to be used.

    @param qprep5_inp:    qprep - input file to be tested.
    @param qprep5_lib:    qprep - library - folder
    @param outfile:       file to write suggested fix

    @return

    @raise Exception: if libraries cant be found.
                      if rewrite needed, but outfile is False.
    """
    def normalize(path):
        """Kill artefacts from a path"""
        return path.replace('//', '/').replace('/./', '/')

    def clean_line(line):
        """Clean a line; split into line and comment"""
        line = line.strip()
        if line == "":
            return '', ''
        line = line.replace('#', '!')
        line = line.split('#', 1)
        if len(line) == 1:
            comment = ''
            line = line[0]
        else:
            line, comment = line
        return line.strip(), comment.strip()

    needs_mod = ''
    quit_use = False
    output = ''

    qprep5_out_top = None
    qprep5_out_pdb = None
    qprep5_maketop = None
    qprep5_in_pdb = 'rawmut.pdb'

    lines = open(qprep5_inp).readlines()

    while len(lines) > 0:
        line, comment = clean_line(lines.pop(0))
        if line == '':
            if comment != '':
                output += '! ' + comment + NLC
            continue

        key = line.split()[0].strip().lower()

        if key == 'readlib' or key == 'rl':
            lib = line.split(None, 1)[1]
            if not os.path.isfile(lib) or qprep5_lib != os.path.abspath(lib)[0:len(qprep5_lib)]:
                needs_mod += "readlib: path adjust"+NLC
                lib = lib.split('/')[-1]
                combo = qprep5_lib + '/' + lib
                combo = normalize(combo)
                if not os.path.isfile(combo):
                    logger.error('Library Not Found: %s', lib)
                    logger.error('Library Not Found: %s', combo)
                    raise Exception
                elif len(combo) > 80:
                    logger.error('Full Library path must not exceed 80 characters. Found: %s', len(combo))
                    logger.error('Path: %s (%s)', len(qprep5_lib), qprep5_lib)
                    logger.error('Filename of prmlib: %s', len(lib))
                    raise Exception
                else:
                    line = 'readlib ' + combo

        elif key == 'readprm' or key == 'rprm':
            prm = line.split(None, 1)[1]
            if not os.path.isfile(prm) or qprep5_lib != os.path.abspath(prm)[0:len(qprep5_lib)]:
                needs_mod += "readprm: path adjust"+NLC
                prm = prm.split('/')[-1]
                combo = qprep5_lib + '/' + prm
                combo = normalize(combo)
                if not os.path.isfile(qprep5_lib + '/' + prm):
                    logger.error('Parameter Not Found: %s', prm)
                    logger.error('Parameter Not Found: %s', combo)
                    raise Exception
                elif len(combo) > 80:
                    logger.error('Full path to prmlib must not exceed 80 characters. Found: %s', len(combo))
                    logger.error('Path: %s (%s)', len(qprep5_lib), qprep5_lib)
                    logger.error('Filename of prmlib: %s', len(prm))
                    raise Exception
                else:
                    line = 'readprm ' + combo

        elif key == 'readpdb' or key == 'rp':
            if qprep5_in_pdb is not None:
                if qprep5_in_pdb != line.replace(key, '').strip():
                    needs_mod += "readpdb: change name"+NLC
                    line = 'readpdb ' + qprep5_in_pdb
            else:
                qprep5_in_pdb = line.split(None, 1)[1]

        elif key == 'writetop' or key == 'wt':
            if line.strip().lower() == key:
                line, comment = clean_line(lines.pop(0))
            else:
                line = line.replace(key, '')
            parts = line.strip().split()
            qprep5_out_top = parts[0]
            line = line.replace(qprep5_out_top, '')
            parts = line.strip().split()
            line = key + '    mutant.top'

        elif key == 'writepdb' or key == 'wp':
            if line.strip().lower() == key:
                line, comment = clean_line(lines.pop(0))
            else:
                line = line.replace(key, '')
            parts = line.strip().split()
            qprep5_out_pdb = parts[0]
            line = line.replace(qprep5_out_pdb, '')
            parts = line.strip().split()
            if len(parts) == 0:
                line, comment = clean_line(lines.pop(0))
            line = key + '     mutant.pdb    ' + line

        elif key == 'maketop' or key == 'mt':
            if line.strip().lower() == key:
                line, comment = clean_line(lines.pop(0))
            else:
                line = line.replace(key, '')
            parts = line.strip().split()
            qprep5_maketop = parts[0]
            line = line.replace(qprep5_maketop, '')
            line = key + '     mutant.top    ' + line

        elif key == 'quit':
            quit_use = True

        else:
            pass

        output += line
        if comment:
            output += '        !' + comment
        output += NLC

    if not quit_use:
        needs_mod += "quit keyword is missing"+NLC
        output += 'quit' + NLC

    if qprep5_in_pdb is None:
        logger.error("no pdbfile is read (readpdb)")
        raise Exception

    if qprep5_out_pdb == qprep5_in_pdb:
        logger.error("writepdb overwrites readpdb. Strongly Discouraged!")
        raise Exception

    if qprep5_out_pdb is None:
        logger.error("pdbfile is not written to disk (writepdb)")
        raise Exception

    if qprep5_out_top is None:
        logger.error("topology is not, uncorrectly written to disk (writetop)")
        raise Exception

    if qprep5_maketop is None:
        logger.error("topology is not created (maketop)")
        raise Exception

    if qprep5_out_pdb != 'mutant.pdb':
        needs_mod += 'writepdb must write to mutant.pdb' + NLC

    if qprep5_out_top != 'mutant.top':
        needs_mod += 'writetop must write to mutant.top' + NLC

    if outfile:
        open(outfile, 'w').write(output)
        qprep5_inp = os.path.abspath(outfile)
    else:
        if needs_mod:
            logger.warn(NLC + 'ERROR')
            logger.info('Qprep5 - The input file (*) must be rewritten like this:')
            logger.info('                  (*)' + qprep5_inp)
            logger.info(output)
            raise QprepInpError
            
    qprep5 = {}
    qprep5['input'] = qprep5_inp
    qprep5['inpdb'] = qprep5_in_pdb
    qprep5['outpdb'] = qprep5_out_pdb
    qprep5['outtop'] = qprep5_out_top
    return qprep5


class QprepInpError(Exception):
    """ QprepInpError; an Error that requires QprepInp-file to be re-written """
    pass


def run_qprep5(qprep5inp, log, qprep5_exe=None):
    """run qprep5 with qprep5inp, output to log.

    @qprep5_exe: if None, try to locate.
    @raise Exception if ERROR found in logfile."""

    if qprep5_exe is None:
        qprep5_exe = tools.get_executable('qprep5')

    if not os.path.isfile(qprep5_exe):
        logger.error('qprep5 does not exist: %s', qprep5_exe)
        raise Exception('qprep5 does not exist:', qprep5_exe)

    if not os.path.isfile(qprep5inp):
        logger.error('qprep5 input file does not exist: %s', qprep5inp)
        raise Exception('qprep5inp (mktop) file does not exist', qprep5inp)

    cmd = qprep5_exe + ' ' + qprep5inp + ' > ' + log + ' 2>&1'

    os.system(cmd)

    # check if logfile was created:
    if not tools.bool_log_exists(log):
        raise Exception('Qprep5 - logfile was not created. Is executable ok?')

    # check if qprep logfile has error messages:
    err = []
    for line in open(log):
        line = line.strip()
        if 'ERROR' in line:
            err.append(line)

    if len(err) > 0:
        logger.error('Errors occured: %s', err)
        raise Exception('error message found in log')


def topogen(qprep5inp, wtpdb, qprep5exe=None):
    """generate topology with qprep5inp and to reproduce wtpdb.
    @param qprep5inp: qprep5 input file
    @param wtpdb: qprep5-ed pdbfile
    @param in_pdb: inputpdb
    @param out_pdb: outputpdb
    @param qprep5exe: qprep5-executable, if None: lookup

    @raise Exception if not ok

    """

    in_pdb = get_qprep5_inputpdb(qprep5inp)
    out_pdb = get_qprep5_outputpdb(qprep5inp)

    try:
        if out_pdb == wtpdb:
            raise Exception('wtpdb would be overwritten by qprep5')
        if in_pdb == wtpdb:
            raise Exception('wtpdb equals the pdb read by qprep5. bad idea.')

        log = os.path.basename(qprep5inp) + '.log'

        run_qprep5(qprep5inp, log, qprep5exe)

        # compare wtpdb with out_pdb:
        wtpdb_lines = open(wtpdb).readlines()
        output_lines = open(out_pdb).readlines()

        if len(wtpdb_lines) != len(output_lines):
            logger.error('wtpdb, regenerated pdb have not same length: %s, %s',
                         len(wtpdb_lines), len(output_lines))

            difference = list(set(wtpdb_lines), set(output_lines))

            if difference != 0:
                logger.error('last 5 differing lines:', difference[-5:])

            raise Exception('wtpdb and regenerated pdb to not match in length')

    except Exception as err:
        print("FATAL: Error while testing re-creation topology of wt." + NLC +
              "       Make sure, that " + qprep5inp + " is fine." + NLC +
              "       *AND* that you pdbfile was made with qprep" + NLC +
              "       *AND* that your libraries exist." + NLC + NLC +
              "Additional Information" + NLC +
              "----------------------" + NLC +
              err)
        import sys
        sys.exit(1)


def get_qprep5_inputpdb(qprep5inp):
    """extract qprep5 inputpdb(readpdb)"""
    parts = split_qprep5_inputfile(qprep5inp)
    for i in range(0, len(parts)):
        if parts[i].lower() == 'readpdb':
            return parts[i+1].strip()
    for i in range(0, len(parts)):
        if parts[i].lower() == 'rp':
            return parts[i+1].strip()


def get_qprep5_outputtop(qprep5inp):
    """extract qprep5 output-top (writetop)"""
    parts = split_qprep5_inputfile(qprep5inp)
    for i in range(0, len(parts)):
        if parts[i].lower() == 'writetop':
            return parts[i+1].strip()
    for i in range(0, len(parts)):
        if parts[i].lower() == 'wt':
            return parts[i+1].strip()


def get_qprep5_outputpdb(qprep5inp):
    """extract qprep5 outputpdb(writpdb)"""
    parts = split_qprep5_inputfile(qprep5inp)
    for i in range(0, len(parts)):
        if parts[i].lower() == 'writepdb':
            return parts[i+1].strip()
    for i in range(0, len(parts)):
        if parts[i].lower() == 'wp':
            return parts[i+1].strip()


def split_qprep5_inputfile(qprep5inp):
    """return comment-freed input file"""
    bline = ''
    with open(qprep5inp) as fil:
        for line in fil:
            line = tools.strip_comments(line)
            line = line.strip()
            bline += " " + line

    return bline.split()


def get_sphere_size_center(qprep5inp, pdbfile=None):
    """return sphere-radius, sphere-center [x,y,z]
    raise Exceptions if not spherical boundaries,
                        spherecenter not x,y,z coordinates
    """

    parts = split_qprep5_inputfile(qprep5inp)

    for i in range(0, len(parts)):
        if parts[i].lower() in ['bc', 'boundary']:
            # found boundarysection
            if parts[i+1] != '1':
                # not sperical boundaries!
                raise Exception('ERROR in qprep5inp: CADEE supports sperical boundaries only.')  # NOPEP8
            if ':' in parts[i+2]:
                # residue:atomname
                if pdbfile is None:
                    raise Exception('ERROR in qprep5inp: No numerical spherecenter, but residue:atomname definition ==> Need pdbfile (not provided)!')  # NOPEP8
                g_resid, g_atomname = parts[i+2].split(':')
                g_resname = None
                anum = tools.get_atomnumber(pdbfile, g_atomname, g_resname,
                                            g_resid)
                center = tools.coords_of_atomnr(pdbfile, anum)
                radius = float(parts[i+3])
            elif parts[i+2] == 'mass':
                raise Exception('ERROR in qprep5inp: spherecenter defined by mass. this is NOT supported by CADEE')  # NOPEP8
            else:
                x_val = float(parts[i+2])
                y_val = float(parts[i+3])
                z_val = float(parts[i+4])
                center = [x_val, y_val, z_val]
                radius = float(parts[i+5])
    return radius, center


def create_top(qprep5inp, outfolder, wtpdb=None):
    """Create topology."""
    out_pdb, out_top, out_fep = create_top_and_fep(qprep5inp, outfolder, wtpdb)
    return out_pdb, out_top


def create_top_and_fep(qprep5inp, outfolder, in_pdb=None, out_pdb=None,
                       out_top=None, wtpdb=None, wtfep=None):
    """Create topology and fep.
    As of now i/o to cwd and then move files to outfolder, if all OK."""

    if in_pdb is None:
        in_pdb = get_qprep5_inputpdb(qprep5inp)
    elif in_pdb != get_qprep5_inputpdb(qprep5inp):
        shutil.copy(in_pdb, get_qprep5_inputpdb(qprep5inp))

    if out_pdb is None:
        out_pdb = get_qprep5_outputpdb(qprep5inp)
    out_top = get_qprep5_outputtop(qprep5inp)
    out_log = os.path.splitext(
                os.path.basename(qprep5inp)
                )[0] + ".log"

    if out_log == qprep5inp:
        raise 'Fatal: Will not overwrite the qprep5inp'

    if wtfep is not None:
        out_fep = out_pdb.replace('.pdb', '.fep')
    else:
        out_fep = None

    if wtpdb is not None:
        shutil.copy(wtpdb, in_pdb)

    try:
        run_qprep5(qprep5inp, out_log)
    except Exception as err:
        logger.error(err)
        # raise Exception('failed to generate topo & pdb')
        raise

    if not os.path.samefile(outfolder, os.getcwd()):
        shutil.move(out_log, outfolder)
        shutil.move(out_top, outfolder)
        shutil.move(out_pdb, outfolder)

    if wtfep is not None:
        cwd = os.getcwd()
        wtpdb = os.path.abspath(wtpdb)
        wtfep = os.path.abspath(wtfep)
        try:
            os.chdir(outfolder)
            create_fep(wtpdb, wtfep, out_pdb, out_fep)
        except:
            logger.error('failed to create fep file')
            raise
        finally:
            os.chdir(cwd)

    if out_pdb != get_qprep5_outputpdb(qprep5inp):
        shutil.move(get_qprep5_outputpdb(qprep5inp), out_pdb)

    return out_pdb, out_top, out_fep
