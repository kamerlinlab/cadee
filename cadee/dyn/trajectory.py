#!/usr/bin/env python

"""
Module calling Q, and analysing logfile / Simulation

Author: {0} ({1})

This module is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function
from filecmp import cmp as comparefiles
from platform import node as hostname
import gzip
import os
import subprocess
import shutil
import time
import analysis
import tools

from tools import File

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logger = tools.getLogger('dyn.traj')

# This class is running a simulation
# Its essentially my python-implementation of submit.sh

# RELEVANT FOR PYTHON - MPI Interaction:
# http://stackoverflow.com/questions/10590064/interact-with-a-mpi-binary-via-a-non-mpi-python-script

# TODO  1: Custom Exceptions / Handling
#       2: if hot-atoms in first 4 steps, restart with smaller stepsize
#       3: automatic mapping
#       4: compress/uncompress pdbfile
#       5: compress/uncompress topology
#       6: compress/uncompress input files
#       7: fastforward md.QdynPackage, when loaded
#       8: mutation starting from just a sequence
#       9: Add flag: has_failed use instead of rising exception
#


# User-Defined Constants
DONT_LOG_LINES_WITH = ["--- Nonbonded pair list generation -----",
                       "---------- Timing ----------",
                       "Seconds per step (wall-clock):"]
NLC = '\n'

# CONSTANTS
NORMALTERM = "terminated normally"
NAN_INDICATOR = 'SUM                    NaN       NaN       NaN'
SHAKE_TERM = "Terminating due to shake failure"
WARNING_HOT_ATOM = ">>> WARNING: hot atom, i ="

ERR_LOG_TOO_SHORT = 1
ERR_ABNORMAL_TERM = 2
ERR_LOG_INEXISTENT = 4
ERR_SHAKE = 8
ERR_NAN = 16


class WorkUnit(object):
    """ container for 1 qdyn-simuluation """
    # class WorkUnitException(Exception):
    #     pass

    # class NaNException(WorkUnitException):
    #     pass

    # class ShakeFailure(WorkUnitException):
    #     pass

    # class AbnormalTermination(WorkUnitException):
    #     pass

    # class NoLogFile(WorkUnitException):
    #     pass

    # class TooShortLogFile(WorkUnitException):
    #     pass
    DEBUG = False

    def __init__(self, unitnumber, inputfile, topology,
                 pdbfile=None, fepfile=None, restraintfile=None,
                 restartfile=None):

        if isinstance(inputfile, str):
            inputfile = [inputfile, '']
        elif isinstance(inputfile, list):
            pass
        else:
            raise 'WTF'

        self.unitnumber = unitnumber

        # error Msg
        self.errMsg = ''

        # shared inputs
        self.topology = topology
        self.pdbfile = pdbfile  # opt
        self.fepfile = fepfile  # opt

        # inputs
        self.inputfile = inputfile
        logger.debug('Input file: %s', self.inputfile[0])
        self.restartfile = restartfile      # opt
        self.restraintfile = restraintfile  # opt

        # outputs
        self.logfile = None
        self.velocityfile = None            # opt
        self.dcdfile = None                 # opt
        self.energyfile = None              # opt

        # stats
        self.time = 0
        self.status = None
        self.q_exitcode = None

        self._parse_inputfile()

        if self.logfile is None or self.logfile[1] == '':
            log = os.path.splitext(self.inputfile[0])[0]+".log"
            loggz = os.path.splitext(self.inputfile[0])[0]+".log.gz"
            if os.path.exists(log) and os.path.exists(loggz):
                # this should not happen, so we give a Msg.warn and
                # then kill the old log
                logger.warning('deleting log, keeping log.gz')
                os.remove(log)
                fname = loggz
            elif os.path.exists(loggz):
                fname = loggz
            else:
                fname = log
            self.logfile = [fname, '']
            if os.path.exists(fname):
                self.checklogfile()
                if self.status != 0:
                    logger.warning('A log file exists BUT with status: %s !', self.status)

        if self.topology is None:
            logger.info(self.inputfile[0])
            raise (Exception, 'topo')
        if self.inputfile is None:
            raise (Exception, 'inp')
        if self.logfile is None:
            raise (Exception, 'log')

    def stats(self, obj):
        """ return stats on fileobj (eg bytes, or empty) """
        if obj is None:
            return ''
        elif isinstance(obj, list) and len(obj) == 2:
            if obj[1] == '':
                return str(obj[0])+': is empty'
            else:
                return str(obj[0])+':'+str(len(obj[1]))+'chars'
        elif isinstance(obj, str):
            return str(obj)+": is empty"
        else:
            logger.info(len(obj))
            raise Exception('weirdo')

    def __repr__(self):
        """ return str of self """
        out = ''
        for obj in [self.topology, self.pdbfile, self.fepfile,
                    self.inputfile, self.restartfile,
                    self.restraintfile, self.logfile, self.velocityfile,
                    self.dcdfile, self.energyfile]:
            out += self.stats(obj) + NLC
            out += self.stats(obj) + NLC
        out += 'status:' + str(self.status)
        out += 'exitcode:' + str(self.q_exitcode)
        return out

    def _parse_inputfile(self):
        ''' Parse input file and populate self with filenames '''
        files_section = False
        logger.debug(os.getcwd())
        if self.inputfile[1] is None or self.inputfile[1] == '':
            with open(self.inputfile[0], 'r') as f:
                self.inputfile[1] = f.read()
        for line in self.inputfile[1].split(NLC):
            if line.strip() == '':
                continue

            # kill comments:
            line = line.replace('#', '!').split('!')[0].strip()

            if files_section:
                if '[' in line.strip()[0]:
                    break
                else:
                    if len(line.split()) != 2:
                        continue
                    ftype, fname = line.split()
                    ftype = ftype.lower().strip()
                    if ftype == 'topology':
                        if self.topology is None:
                            self.topology = [fname, '']
                    elif ftype == 'fep':
                        if self.fepfile is None:
                            self.fepfile = [fname, '']
                    elif ftype == 'restart':
                        if self.restartfile is None:
                            self.restartfile = [fname, '']
                    elif ftype == 'restraint':
                        if self.restraintfile is None:
                            if os.path.isfile(fname):
                                self.restraintfile = [fname, '']
                            else:
                                if fname == self.restartfile[0] + "st.re":
                                    # restart + 'st.re' == restraint
                                    msg = 'RE-use restart as restraint:'
                                    msg += '-----> %s'
                                    msg += self.restartfile[0]
                                    logger.debug(msg)
                                    shutil.copy(self.restartfile[0], fname)
                                    self.restraintfile = [fname, '']
                                else:
                                    # search harddisk for similar restraint
                                    restartfile = fname[:-5]
                                    msg = 'Try: use restart as restraint:'
                                    msg += '-----> %s'
                                    msg += restartfile
                                    logger.debug(msg)
                                    if os.path.isfile(restartfile):
                                        shutil.copy(restartfile, fname)
                                    else:
                                        msg = 'cant find restraint:' + fname
                                        logger.warning(msg)
                                        # TODO: more heuristics
                                        raise (Exception, msg)
                    elif ftype == 'final' and self.velocityfile is None:
                        self.velocityfile = [fname, '']
                    elif ftype == 'trajectory' and self.dcdfile is None:
                        self.dcdfile = [fname, '']
                    elif ftype == 'energy' and self.energyfile is None:
                        self.energyfile = [fname, '']
                    else:
                        logger.warning('do not know this key here %s', ftype)
                        raise (Exception, 'do not know this key here')

            if line.lower() == '[files]':
                files_section = True

        if not files_section:
            logger.warning('Fatal: no files section in input file %s',
                           self.inputfile)
            raise (Exception, 'Fatal: no files section in input file')

    def run(self, exe):
        """ run simulation with executable exe """
        if os.path.isfile(self.logfile[0]):
            self.status = self.checklogfile()
            if self.status == 0:
                return 0

        ifname = self.inputfile[0]
        ofname = self.logfile[0]

        if len(ifname) == 1 or len(ofname) == 1:
            raise (Exception, 'WTF')

        self._deploy()

        # run q
        start = time.time()
        cmd = exe + " " + ifname
        logger.info("%s %s", hostname(), cmd)
        try:
            subprocess.check_call([exe, ifname], stdout=open(ofname, 'w'))
            self.q_exitcode = 0
        except subprocess.CalledProcessError as exitstatus:
            logger.warning('Detected a non-zero exit status!', exitstatus)
            self.q_exitcode = exitstatus.returncode

        # check logfile
        self.checklogfile()

        self.time = time.time() - start

        if self.status == 0 and self.q_exitcode == 0:
            return 0
        else:
            logger.warning('Detected status %s %s %s', self.status,
                           'and an exitcode', self.q_exitcode)
            return self.status + self.q_exitcode

    def checklogfile(self):
        """ Check Logfile """

        logger.debug("Checking log file ...")

        # TODO: do something about hot atoms.
        #       eg  a) run longer
        #           b) never read the restraint of this run
        # WARN_HOT_ATOM = 0
        self.hot_atoms = 0

        log = []

        logfile = self.logfile[0]

        if os.path.isfile(self.logfile[0]+".gz"):
            os.remove(self.logfile[0])
            self.logfile[0] = self.logfile[0]+".gz"

        try:
            if logfile[-3:] == ".gz":
                for line in gzip.open(logfile):
                    log.append(line)
                compress = False
            else:
                with open(logfile) as fil:
                    for line in fil:
                        if line in DONT_LOG_LINES_WITH:
                            continue
                        log.append(line)
                compress = True
        except IOError:
            err = "Could not open log file!"
            logger.warning(err)
            self.errMsg += err
            self.status = ERR_LOG_INEXISTENT
            return ERR_LOG_INEXISTENT

        if len(log) > 5:
            l5l = 'Last 5 lines of log file:'
            self.errMsg += NLC + l5l + NLC
            for logline in log[-5:]:
                self.errMsg += len('Last 5 Lines') * ' ' + ">" + logline
            self.errMsg += "/" + 'Last 5 Lines' + NLC + NLC

        if len(log) < 100:
            err = 'The log file is too short (less than 100 lines)!'
            logger.warning(err)
            self.errMsg += err
            self.status = ERR_LOG_TOO_SHORT
            for line in log:
                logger.info(line.strip())
            return ERR_LOG_TOO_SHORT
        else:
            # TODO: check if we have insane high energies
            for line in log[-50:]:
                allok = 0
                if NORMALTERM in line:
                    allok = 1
                    break
                elif NAN_INDICATOR in line:
                    err = "Found NaN Error: '" + NAN_INDICATOR + "'"
                    logger.warning(err)
                    self.errMsg += err
                    self.status = ERR_NAN
                    return ERR_NAN
                elif SHAKE_TERM in line:
                    err = "Found Shake Error: '" + SHAKE_TERM + "'"
                    logger.warning(err)
                    self.errMsg += err
                    self.status = ERR_SHAKE
                    return ERR_SHAKE
                elif WARNING_HOT_ATOM in line:
                    if self.hot_atoms < 1:
                        err = "Found HOT ATOM'"
                        logger.warning(err)
                        self.errMsg += err
                    self.hot_atoms += 1
            if allok != 1:
                err = 'The log file is missing ' + str(NORMALTERM) + ' string!'
                err += ' UNKNOWN ERROR '
                self.errMsg += err
                logger.warning(err)

                self.status = ERR_ABNORMAL_TERM

                return ERR_ABNORMAL_TERM

        if compress:
            # re-writing compressed logfile without rubbish lines
            gzip.open(logfile+".gz", 'wb').writelines(log)
            self.logfile[0] = logfile+".gz"
            os.remove(logfile)

        # search and kill douplicate *rest.re files
        if self.restraintfile is not None and self.restartfile is not None:
            if len(self.restraintfile) == 2 and len(self.restartfile) == 2:
                restre = self.restraintfile[0]
                restart = self.restartfile[0]
                if (len(restre) > 8 and len(restart) > 3 and
                        restre[-8:] == ".rest.re" and restart[-3:] == ".re"):
                    if os.path.isfile(restre) and os.path.isfile(restart):
                        if comparefiles(restart, restre, False):
                            os.remove(restre)

        # compress energyfile
        if self.energyfile is not None and len(self.energyfile) == 2:
            energy = self.energyfile[0]
            if os.path.isfile(energy) and energy[-3:] != ".gz":
                engz = energy + ".gz"
                with open(energy, 'rb') as fil_in:
                    with gzip.open(engz, 'wb') as fil_out:
                        shutil.copyfileobj(fil_in, fil_out)
                        os.remove(energy)

        self.status = 0
        return 0

    def _deploy(self):
        """ serialize data from memory-oject to disk (deploy to disk) """
        for data in (self.topology, self.pdbfile, self.fepfile,
                     self.inputfile, self.restraintfile, self.restartfile,
                     self.logfile, self.dcdfile,
                     self.energyfile, self.velocityfile):
            if data is None:
                continue

            if isinstance(data, str):
                continue

            fname, data = data

            if data.strip() == '':
                if WorkUnit.DEBUG:
                    logger.debug('Wont write empty file %s !', fname)
                continue

            # if os.path.isfile(fname):
            #     Msg.log('Wont overwrite existing file', fname)
            #     continue

            if isinstance(fname, str) and isinstance(data, str):
                with open(fname, 'wb') as fil:
                    fil.write(data)
                if WorkUnit.DEBUG:
                    logger.debug('Serialized: %s .', fname)
            else:
                logger.warning('This might be a problem here: type(fname) == %s %s %s !',
                               type(fname), 'and type(data) ==', type(data))
                raise (Exception, 'Expected: str() and str()')


class QdynPackage(object):

    def map_and_analyze(self, eqfil=None):
        if self.mapped is None:
            logger.debug('Mapping disabled.')
        elif self.mapped is True:
            logger.debug('is already mapped (skipping)!')
            return True
        elif self.mapped is False:
            with tools.cd(self.path):
                if eqfil is None:
                    self.mapped = analysis.main(self.map_settings, eqfil)
                else:
                    analysis.main(self.map_settings, eqfil)
        else:
            raise 'WTF'

    def parse_file(self, fname):
        """
        Read fname.
        If fname is None:
            return None
        if fname is str:
            read the file and store in memory
        if fname is [path, data]
            return [basename(path), data]
        """
        if fname is None:
            return None

        if isinstance(fname, str):
            if os.path.isfile(fname):
                with open(fname, 'r') as fil:
                    data = fil.read()
                fname = os.path.basename(fname)
                return [fname, data]
            else:
                logger.warning('Could not find %s .', fname)
                os.system('ls')
                raise 'FAILED'
        elif isinstance(fname, list) and len(fname) == 2:
            fname[0] = os.path.basename(fname[0])
            return fname
        else:
            raise (Exception, 'either str(fname) or list[fname, data])')

    def check_exe(self):
        """ check executable permissions, raises exception if not OK """
        if (isinstance(self.q_dyn5_exe, str) and
                os.path.isfile(self.q_dyn5_exe)):
            if os.access(self.q_dyn5_exe, os.X_OK):
                pass
            else:
                raise (Exception, 'executable is not executable!')
        else:
            raise (Exception, 'executable: is not file.')

    def set_exe(self, exe):
        """ set self.q_dyn5_exe to exe """
        self.q_dyn5_exe = exe

    def set_temp(self, temp):
        """ cd into temp """
        if not os.path.isdir(temp):
            raise (Exception, 'you provided a path which is not a directory!')
        else:
            os.chdir(temp)

        self.path = temp

    def __init__(self, topology, path, q_executable=None, inputfiles=None,
                 description=None, pdbfile=None, restartfile=None,
                 restraintfile=None, fepfile=None, map_settings=None):
        """
        Inititalisation of a md-simulation.
        Topology:
            string with a filename OR list [ topology-name, data ]
        Path:
            a string with the path where this unit is working in.
        Executable:
            string with a filename
        Inputfiles:
            list with inputfiles.
        Description:
            str with description
        PDBFile:
            string with a filename OR list [ pdb-name, data ]
        restartfile:
            string with a filename OR list [ restartname, data ]
        restraintfile:
            string with a filename OR list [ restraintname, data ]

        """

        self.set_temp(path)

        if q_executable is not None:
            self.set_exe(q_executable)

        self.topology = self.parse_file(topology)
        self.fepfile = self.parse_file(fepfile)

        self.pdbfile = self.parse_file(pdbfile)
        self.restartfile = self.parse_file(restartfile)
        self.restraintfile = self.parse_file(restraintfile)
        self.description = self.parse_file(description)

        self.wus = []
        self.cwu = None
        self.if_pos = 0
        self.inputfiles = []

        self.filenames = {}

        if inputfiles is not None and isinstance(inputfiles, list):
            for inp in inputfiles:
                if inp is None:
                    raise 'WTF'
                self.add_inputfile(inp)
        else:
            self.inputfiles = []

        if map_settings is None:
            self.mapped = None
        else:
            self.mapped = False
            self.map_settings = map_settings

        logger.info('Next qdyn simulation step initialized.')

    def stats(self, obj):
        if obj is None:
            return ''
        elif isinstance(obj, list) and len(obj) == 2:
            if obj[1] == '':
                return str(obj[0]) + ': is empty'
            else:
                return str(obj[0]) + ':' + str(len(obj[1])) + 'chars'
        elif isinstance(obj, str):
            return str(obj)+": is empty"
        else:
            logger.warning('weirdo %s', (len(obj)))
            raise Exception('weirdo')

    def __repr__(self):
        out = ''
        out += 'PackageFiles' + NLC
        for obj in [self.topology, self.pdbfile, self.fepfile,
                    self.restartfile, self.restraintfile, self.description]:
            if self.stats(obj) != '':
                out += ' ' + self.stats(obj) + NLC
        out += NLC
        # out += '# inputiles:' + str(self.inputfiles)
        out += 'Work-Units:' + NLC
        for each in self.wus:
            out += ' WU: ' + str(each) + NLC
        out += 'Input files:' + NLC
        for i in range(len(self.wus), len(self.inputfiles)):
            out += ' IF: ' + str(self.inputfiles[i]) + NLC
        out += '@  position:' + str(self.if_pos) + NLC
        out += '   finished:' + str(self.is_finished()) + NLC
        return out

    def get_file_with_name(self, fname):
        with tools.cd(self.path):
            if os.path.exists(fname):
                return File(fname, read=True)
        raise IOError('File not found', fname, 'in', self.path)

    def _check_eq_and_map(self):
        if '_eq' in self.inputfiles[self.if_pos][0]:
            logger.debug('is eq-file %s', (self.inputfiles[self.if_pos][0]))
            self.map_and_analyze(self.inputfiles[self.if_pos][0])
        elif '_fep' in self.inputfiles[self.if_pos][0]:
            pass
        elif '_dyn' in self.inputfiles[self.if_pos][0]:
            pass
        else:
            logger.warning('Neither temperization, nor equlilbration nor fep: %s',
                           (self.inputfiles[self.if_pos][0]))

    def compute(self):
        """
        Run Q:
        Forward trough inputfile-list to inputfile without logfile and run it.
        """
        with tools.cd(self.path):
            # Check if we're done arleady
            if self.is_finished():
                try:
                    logger.warning('Nothing to compute. %s %s', self.if_pos,
                                   self.inputfiles[self.if_pos])
                except IndexError:
                    logger.warning('Nothing to compute. %s', self.if_pos)

                # TODO: 1) automatic mapping

            else:
                # TODO: add restart-capability

                if not self.is_finished():
                    if len(self.wus) > self.if_pos:
                        old_input = self.wus[self.if_pos].inputfile[0]
                        new_input = self.inputfiles[self.if_pos]
                        if old_input == new_input:
                            if self.wus[self.if_pos].checklogfile() == 0:
                                # this WorkUnit is finished we; load next one
                                logger.warning(
                                    'this run is already done skipping %s',
                                    self.wus[self.if_pos].inputfile[0]
                                    )
                                self.if_pos += 1
                                self.compute()
                                return
                    # Generate new compute units, until one w/o logfile exists
                    while True:
                        if self.is_finished():
                            return

                        self.cwu = self.create_next_workunit()
                        if (self.cwu.status is not None and
                                self.cwu.status == 0):
                            logger.debug(
                                    'skip step %s',
                                    self.inputfiles[self.if_pos][0])
                            self._check_eq_and_map()
                            self.wus.append(self.cwu)
                            self.if_pos += 1
                            continue
                        break

                    exe = self.q_dyn5_exe
                    self.check_exe()

                    if len(self.wus) != self.cwu.unitnumber:
                        raise (Exception, 'discrepancy in input file order')

                    if self.cwu.run(exe) == 0:
                        self.wus.append(self.cwu)
                        self._check_eq_and_map()
                    else:
                        err = 'There was a problem with step: '
                        err += str(self.if_pos)
                        err += ', in inputfile'
                        err += str(self.inputfiles[self.if_pos][0])
                        err += NLC + 'The status Code was:'
                        err += str(self.cwu.status)
                        err += NLC + NLC + 'The Error Messages where: '
                        err += NLC + NLC + 'Directory' + os.getcwd()
                        err += str(self.cwu.errMsg) + NLC
                        err += 'Will Raise Exception...'
                        logger.warning(err)
                        raise (Exception, 'computation failed')

                    # increment for next step
                    self.if_pos += 1

    def _parse_inputfile(self, inputfile):
        ''' Parse input file '''
        files_section = False
        self.filenames = {}
        for line in open(inputfile, 'r'):
            if line.strip() == '':
                continue

            # kill comments:
            line = line.replace['#', '!'].split()[0].strip()

            if files_section:
                if '[' in line.strip()[0]:
                    files_section = False
                else:
                    ftype, fname = line.split()
                    ftype = ftype.lower()
                    if ftype == 'topology':
                        self.filenames[ftype] = fname
                    elif ftype == 'fepfile':
                        self.filenames[ftype] = fname
                    elif ftype == 'restart':
                        self.filenames[ftype] = fname
                    elif ftype == 'restraint':
                        self.filenames[ftype] = fname
                    elif ftype == 'final':
                        self.filenames[ftype] = fname
                    elif ftype == 'trajectory':
                        self.filenames[ftype] = fname
                    elif ftype == 'energy':
                        self.filenames[ftype] = fname
                    else:
                        raise (Exception, "Parse Input File: Unknown Filetype.")

            if line.lower() == '[files]':
                files_section = True

    def add_inputfile(self, inputfile):
        ''' Add an inputfile '''
        inputfile = self.parse_file(inputfile)
        if inputfile is None:
            raise 'WTF'
        self.inputfiles.append(inputfile)

    def create_next_workunit(self):
        ''' Load & Prepare the next input file '''

        if self.if_pos == 0:
            # initial one, may need special input
            cwu = WorkUnit(self.if_pos, self.inputfiles[self.if_pos],
                           self.topology, self.pdbfile, self.fepfile,
                           self.restraintfile, self.restartfile)

        else:
            # try to locate the restart and restraint files that might be need
            restart = None
            restraint = None
            self.parse_file(self.inputfiles[self.if_pos])
            if 'restart' in self.filenames:
                for i in range(self.if_pos):
                    old_restart = self.wus[self.if_pos-i].velocityfile[0]
                    new_restart = self.filenames['restart']
                    if old_restart == new_restart:
                        restart = self.wus[self.if_pos-i].velocityfile
                        if i != 1:
                            logger.warning('UNUSUAL RESTART')
                            logger.warning(
                                    'input: %s',
                                    self.inputfiles[self.if_pos])
                            logger.warning(
                                    'is not using restart from: %s',
                                    self.inputfiles[self.if_pos-1])
                            logger.warning(
                                    'BUT INSTEAD the one from: %s',
                                    self.inputfiles[self.if_pos-i])
            if 'restraint' in self.filenames:
                for i in range(self.if_pos):
                    old_restraint = self.wus[self.if_pos-i].velocityfile[0] + "st.re"  # NOPEP8
                    new_restraint = self.filenames['restraint']
                    if old_restraint == new_restraint:
                        restraint = self.wus[self.if_pos-i].velocityfile[:]
                        restraint[0] += "st.re"
            # TODO: multiple fep files could be taken from here as well
            cwu = WorkUnit(self.if_pos, self.inputfiles[self.if_pos],
                           self.topology, None, self.fepfile, restraint,
                           restart)
        return cwu

    def is_finished(self):
        """ Return True is simulation is finished, or False. """
        if len(self.inputfiles) <= self.if_pos:
            self.map_and_analyze()
            return True
        return False

    def progress(self):
        """ logger.info progress of md.QdynPackage """
        logger.info("%s / %s", self.if_pos, len(self.inputfiles))


class MolDynSim(object):
    """ Container for MolecularDynamicSimulation """
    def __init__(self, path, qexecutable, topology,
                 inputfiles, fepfile, restartfile,
                 restraintfile, description, pdbfile,
                 map_settings=None):
        self.pack = QdynPackage(topology, path, qexecutable, inputfiles,
                                description, pdbfile, restartfile,
                                restraintfile, fepfile, map_settings)

    def set_executable(self, exe):
        """ set executable """
        self.pack.set_exe(exe)

    def set_tempdir(self, temp):
        """ set temporary directory """
        self.pack.set_temp(temp)

    def is_finished(self):
        """ return bool """
        return self.pack.is_finished()

    def progress(self):
        """ write progress to stdout"""
        self.pack.progress()

    def __repr__(self):
        return self.pack.__repr__()

    def compute(self):
        """ run one step, return self """
        self.pack.compute()
        return self.pack
