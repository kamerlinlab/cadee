#!/usr/bin/env python

"""Ensemble Simulation with MPI4PY.

This module is responsible to coordinate Input/Output of the
ensemble simulation.

It contains instructions for IO and both,
one class for rank0 (master) and the other ranks (worker)

The master is responsible for ALL IO, to and from disk and to stdout and stderr

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""

from __future__ import print_function
from platform import node as hostname
import argparse
import os
import signal
import sys
import shutil
import time
import tarfile
import cPickle
import hashlib
import traceback

import scan

import mpi

import tools
import trajectory

if mpi.mpi:
    from mpi4py import MPI
else:
    tools.getLogger().fatal('You can not run ensemble simulation without MPI!')
    sys.exit(1)

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

logger = tools.getLogger()

logger.debug('start')

if mpi.rank == mpi.root:
    print('Rank 0: Started @', time.time(), mpi.get_info())

# This class is running a simulation


# RELEVANT FOR PYTHON - MPI Interaction:
# http://stackoverflow.com/questions/10590064/interact-with-a-mpi-binary-via-a-non-mpi-python-script

# User-Defined Constants

DEBUG = False

RAISE_EXCEPTIONS = False

# 'Dynamic' IO Adjustment, at least 8 files at the same time,
#  or mpi.size/16 (whatever is bigger)
PARALLEL_IO = 8
if int(mpi.size / 16) > PARALLEL_IO:
    PARALLEL_IO = int(mpi.size / 16)


DONT_LOG_LINES_WITH = ["--- Nonbonded pair list generation -----",
                       "---------- Timing ----------",
                       "Seconds per step (wall-clock):"]

MIN_BACKUP_INTERVAL = 600  # [s] min. sec between backups to persistent storage

NLC = '\n'

# TODO: scale parallel_io with jobsize
MD5, MTIME, SIZE = (1, 2, 3)

if DEBUG:
    RAISE_EXCEPTIONS = True
    MIN_BACKUP_INTERVAL = 60
    PARALLEL_IO = 2


def check_int_or_float(value):
    """Check if a value is castable to float or rise ArgumentTypeError"""
    try:
        return float(value)
    except:
        raise argparse.ArgumentTypeError("%s is neither integer nor float" % value)  # NOPEP8


def log_speed(elapsed, fsize, name):
    """ Write a logger.info with speed of I/O operation
    @param elapsed: elapsed time in seconds
    @param fsize:   filesize in MegaBytes (Bytes/ 1024/1024.)
    @param name: a name
    @type elapsed: int
    @type fsize: float
    @type name: str
    """
    msg = 'Backup-timing for {name}, {elapsed:5.3f}s, MB: {fsize:6.2f}'
    msg += ' Speed: {speed} MB/s'

    msg = msg.format(name=name, elapsed=elapsed,
                     fsize=fsize, speed=fsize/elapsed)
    logger.info(msg)


class Worker(object):
    """ MPI - Worker
        This is a MPI-worker. (rank>0).
    """

    def __init__(self, tempdir, a, h, force_remap):
        """
        @param tempdir: path to store temporary files
        @type tempdir: str
        @return: None
        """

        self.comm = mpi.comm
        self.rank = mpi.rank
        self.alpha = a
        self.hij = h
        self.force_remap = force_remap
        self.root = 0
        self.nanos = -1.0
        self.steps = -1
        self._executable()
        self._tempdir(tempdir, mpi.rank)
        self._md = None
        self.archive = None
        self.saved_files = {}
        self.lastbackup = time.time()
        self.alive = True

    def _executable(self):
        """ Create executable and mark it executable """
        from cadee.executables import exe
        exe = exe.which('qdyn5')
        if exe is None:
            raise 'qdyn5 not found'
        self.exe = exe

    def _next(self):
        """Receive next unit of work from Master (rank0).
        """
        if not self.alive:
            return False

        try:
            intar = self.inputarchive
        except AttributeError:
            intar = 'START'

        logger.debug('Worker send data')
        self.comm.send(intar, self.root, tag=mpi.Tags.DONE)
        logger.debug('Worker wait data')
        data = self.comm.recv(source=self.root, tag=mpi.Tags.INPUTS)
        logger.debug('Worker recd data')

        if data == 'SHUTDOWN':
            logger.debug('Worker %s got shutdown signal!', mpi.rank)
            logger.debug('Create last backup ...')
            try:
                self._store()
            except ValueError:
                logger.exception(
                    'Could not store before exiting %s', self.archive)
            logger.debug('Tell master that we exited...')
            self.comm.send('GoodBye!', self.root, tag=mpi.Tags.SHUTDOWN)
            self.alive = False
            sys.exit(0)
        else:
            logger.info('Worker reinitializing.')
            intar, outtar = data
            self.reinit(intar, outtar)
            return True

    def _tar2md(self, tarchive, map_settings):
        """
        @param tarchive: the archive that will be extracted.
        @type tarchive: str
        @return mdobj
        ::note::
        WARNING: The md-object has to be re-initialized.
        """

        logger.info('unpacking: %s', os.getcwd())
        # UNIT: make sure there is only 1 executable in this folder
        logger.debug(str(os.listdir(os.getcwd())))

        # REQUEST AND WAIT FOR IO-TICKET
        self.comm.send('', self.root, tag=mpi.Tags.IO_REQUEST)
        self.comm.recv(source=self.root, tag=mpi.Tags.IO_TICKET)
        start = time.time()
        try:
            tarfile.open(tarchive).extractall()  # TODO: UNSAVE IF TARCHIVE $@!
            fsize = os.path.getsize(tarchive) / 1024 / 1024.  # TODO:Do in scan
        # RETURN TICKET. DO NOT FORGET 'FINALLY' IS FOR E.G. CASE OF IO-ERROR
        finally:
            self.comm.send('', dest=self.root, tag=mpi.Tags.IO_FINISHED)

        elapsed = time.time() - start

        log_speed(elapsed, fsize, self.archive)

        # TODO: scan for pdbfile(s) and or description
        topology, fepfile, inputfiles, self.steps, self.nanos = scan.Scan().scan()  # NOPEP8
        pdbfile = None
        description = None

        mdobj = trajectory.MolDynSim(self.tmp, self.exe, topology, inputfiles,
                                     fepfile, None, None, description, pdbfile,
                                     map_settings)
        return mdobj

    def reinit(self, inputarchive, outputarchive):
        """
        @param tempdir: a path to store temporary files
        @param inputarchive: tarchive with inputfiles
        @param outputarchive: tarchive where results are written to
        """
        try:
            shutil.rmtree(self.tmp)
            del self._md
        except AttributeError:
            pass

        os.makedirs(self.tmp)
        os.chdir(self.tmp)
        self._executable()

        if self.alpha is not None:
            import analysis
            fname = os.path.basename(inputarchive)
            name = fname.split('_')

            if len(name) == 2:
                mutant = name[0]
                replik = int(name[1].split('.')[0])
            else:
                mutant = name[0]
                replik = 0
                logger.warning('No replik-info found. Setting to 0')

            mset = analysis.MapSettings(
                    mutant, replik, self.alpha,
                    self.hij, self.force_remap)

            self._md = self._tar2md(inputarchive, mset)
        else:
            self._md = self._tar2md(inputarchive, None)

        if outputarchive is None:
            logger.warning("NO OUTPUTARCHIVE. Appending data to inputarchive.")
            self.archive = inputarchive
        else:
            self.archive = outputarchive

        self.saved_files = {}

        # initizalize self.saved_files, so we do not add files to archive 2x
        if inputarchive == outputarchive:
            for obj in self._check_files_to_store():
                fname, mtim, md5, size = obj
                self.saved_files[MTIME][fname] = mtim
                self.saved_files[SIZE][fname] = size
                self.saved_files[MD5][fname] = md5

        self.lastbackup = time.time()

        logger.debug('Initialized on %s, in %s',
                     hostname(), os.getcwd())

    def _check_files_to_store(self):
        """
        Walks cwd and checks if files were changed and/or added.
        Return list of files to save with fn, modtime, md5hash, size
        Return: [ [path,mtim,md5,size], [path2,mtim2,md5_2,size2] ... ]
        """

        if MTIME not in self.saved_files:
            self.saved_files[MTIME] = {}
        if MD5 not in self.saved_files:
            self.saved_files[MD5] = {}
        if SIZE not in self.saved_files:
            self.saved_files[SIZE] = {}

        to_store = []

        for fname in os.listdir('.'):
            mtim = os.path.getmtime(fname)
            size = os.path.getsize(fname)

            # never backup executable
            if fname == self.exe or fname == os.path.basename(self.exe):
                continue

            # check if fn was seen before and if size or mtime changed:
            if (fname in self.saved_files[MTIME] and
                    fname in self.saved_files[MD5] and
                    fname in self.saved_files[SIZE] and
                    self.saved_files[MTIME][fname] == mtim and
                    self.saved_files[SIZE][fname] == size):
                continue

            # hashing (since expensive, only if mtime, size or hash chgd)
            md5 = hashlib.md5(open(fname, 'rb').read()).hexdigest()
            if (fname in self.saved_files[MD5] and
                    self.saved_files[MD5][fname] == md5):
                logger.debug('rehashed %s but hash didnt change!', fname)
                continue
            to_store.append([fname, mtim, md5, size])

        # TODO: pickle self.saved_files and add it to to_store

        # because logfile is used to check if a run was
        # successful, logfiles must be written last, to ensure "save restarts";
        # eg. crash during writing the data out, could cause a corrupt md run!
        #     => move logfiles to end of list
        logs = []
        other = []
        for fil in to_store:
            if fil[0][-7:] == ".log.gz":
                logs.append(fil)
            else:
                other.append(fil)

        to_store = other
        to_store.extend(logs)

        return to_store

    def _store(self):
        """ append new data to self.archive

        check if file is in self.saved_files and modification time,
              else append to archive

        WARN: THIS IS NOT THREAD SAVE:
              ITS NOT POSSIBLE TO COMPUTE & STORE SIMULTANIOUSLY"""

        to_store = self._check_files_to_store()

        if len(to_store) == 0:
            return

        # To avoid reading harddisk while having IO-Ticket:
        # http://stackoverflow.com/questions/15857792/how-to-construct-a-tarfile-object-in-memory-from-byte-buffer-in-python-3
        # REQUEST AND WAIT FOR IO-TICKET
        self.comm.send('', self.root, tag=mpi.Tags.IO_REQUEST)
        self.comm.recv(source=self.root, tag=mpi.Tags.IO_TICKET)
        try:
            start = time.time()
            fsize = 0.0
            tar = tarfile.open(self.archive, 'a')
            for obj in to_store:
                fname, mtim, md5, size = obj
                tar.add(fname)
                self.saved_files[MTIME][fname] = mtim
                self.saved_files[SIZE][fname] = size
                self.saved_files[MD5][fname] = md5
                fsize += size
            tar.close()
        # RETURN TICKET. DO NOT FORGET 'FINALLY' IS FOR E.G. CASE OF IO-ERROR
        except ValueError:
            logger.exception(
                    'Exception while opening %s for appending', self.archive)
        finally:
            self.comm.send('', self.root, tag=mpi.Tags.IO_FINISHED)
        elapsed = time.time() - start
        fsize = fsize/1024/1024.

        self.lastbackup = time.time()

        log_speed(elapsed, fsize, self.archive)

    def _tempdir(self, tempdir, rank):
        """ Create tempdir/{rank} and cd into it
        @param tempdir: temporary directory
        @param rank: rank
        @type tempdir: str
        @type rank: int
        """
        self.tmp = os.path.join(tempdir, str(rank)) + str("/")
        if not os.path.exists(self.tmp):
            os.makedirs(self.tmp)
        os.chdir(self.tmp)

    def _compute(self):
        """ Compute 1 step """
        self._md.compute()
        if (time.time() - self.lastbackup) > MIN_BACKUP_INTERVAL:
            self._store()

    def _term_handler(self, signum, frame):
        """ Signal Handler
        @param signum
        @param frame

        @raises KeyboardInterrupt
        """
        logger.warning('Signal received. %s %s', signum, frame)
        logger.info('Creating Backup')
        self._store()
        self.comm.send('GoodBye!', self.root, tag=mpi.Tags.SHUTDOWN)
        logger.warning('Backup done. Raise KeyboardInterrupt...')

        raise KeyboardInterrupt

    def run(self):
        """Compute all steps
        Compute all steps, without giving control back to caller, until
        all steps are computed or an exception is raised.

        @return: None

        ::note ::
        When a step is finished and MIN_BACKUP_INTERVAL [s] has passed,
        backup is written to disk.
        """
        if self._md is None:
            self._next()

        if not self.alive:
            return

        try:
            self._md.is_finished()
        except AttributeError as e:
            logger.debug('Attribute error %s, load next', e)
            if not self._next():
                logger.debug('Attribute Error, done!')
                return

        logger.debug('working')
        signal.signal(signal.SIGTERM, self._term_handler)
        while not self._md.is_finished():
            try:
                self._compute()
            except Exception as err:
                logger.exception('Caught Exception %s, while processing %s',
                                 err, self.archive)
                raise

        if (time.time() - self.lastbackup) > MIN_BACKUP_INTERVAL:
            self._store()

        if self._next():
            logger.debug('run done, recursive loop')
            self.run()


class Master(object):
    """MPI Rank 0:
    Distributes Inputfiles to Worker Nodes.
    Receives MPI messages with tags defined in mpi.Tags.Class
    """
    def __init__(self, tempdir, start, simpackdir, force_map=False):
        self.comm = mpi.comm
        self.tmp = tempdir + str(0) + str("/")
        self.inputlist = []
        self.numworkers = mpi.size-1
        self.io_tickets = [None]*mpi.size
        self.io_queue = []
        self.start = start
        self.simpackdir = simpackdir

        dbname = os.path.join(simpackdir, 'cadee.db')

        if os.path.exists(dbname) and force_map:
            logger.warning('Database and would be overwritten: %s', dbname)
            i = 0
            while os.path.exists(dbname + str(i)):
                i += 1
            os.rename(dbname, dbname+str(i))
            logger.info('Created Backup: %s BACKUP!', dbname+str(i))


        self.db = tools.SqlDB(
            os.path.join(simpackdir, 'cadee.db')
            )

        # TODO make sure output file does not exist!
        if not os.path.exists(self.tmp):
            os.makedirs(self.tmp)

    def _shutdown(self):
        logger.info('Preparing to end this Simulation! Syncing...')
        self.db.close()
        logger.info('Database connection closed.')
        logger.info('Removing Temporary Files...')
        import shutil
        shutil.rmtree(self.tmp)
        logger.info("DONE. Exiting")
        sys.exit()

    # TODO: Does not work for slaves
    def _term_handler(self, signum, frame):
        """
        Signal Handler for Signal Term
        """
        logger.warning('Signal received. %s %s', signum, frame)
        logger.warning('Press CTRL+C to KILL')

        if len(self.inputlist) > 0:
            logger.warning('Following Simpacks have been removed unfinished from the Queue: %s, Items: %s',
                           len(self.inputlist), str(self.inputlist))
            self.inputlist = []

        logger.warning('Waiting up to 170s for Workers to shutdown...')

        stop = time.time() + 170
        while time.time() < stop:
            self._iter()
        logger.warning('... Timeout!')
        if mpi.mpi:
            logger.warning('MPI_ABORT')
            MPI.MPI_Abort(self.comm, int(signum))
        self._shutdown()

    def _iter(self):
        sleeped = 0  # measure time w/o messaging
        # Waiting for any input
        while not self.comm.Iprobe(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG):
            if sleeped == 0.5:
                secs_since_start = round(time.time() - self.start, 1)
                logger.info('sleeping @ %s s', secs_since_start)
            # TODO: stats!
            #       observe ram usage
            #       observe free disk space
            #       observe ...
            #       make sure executed every 60? seconds
            if self.numworkers == 0:
                self._shutdown()
            # TODO: think about better use of sleeped
            sleeped += self._manage_io()

        if sleeped > 0.5:
            logger.info("slept for %s seconds", sleeped - 0.5)

        # TODO: catch mpi.send errors properly
        try:
            self._process_mpi()
        except IndexError as e:
            logger.info('IndexError happend: %s', e)

    def run(self):
        """
        Communicate with nodes.
        """
        # TODO:
        # If we scale with even more nodes, we might choose
        # to use an atomic database to manage inputs or we
        # have to split the mpi-regions.

        # signal handler, waiting for nodes to stop
        signal.signal(signal.SIGTERM, self._term_handler)

        while True:
            self._iter()

    def _process_mpi(self):
        """Process MPI Package"""
        mpistatus = MPI.Status()
        try:
            data = self.comm.recv(source=MPI.ANY_SOURCE,
                                  tag=MPI.ANY_TAG, status=mpistatus)
        except cPickle.UnpicklingError as err:
            logger.critical('ERROR! Unpickling Error: %s tag: %s, source: %s',
                            err, mpistatus.Get_tag(), mpistatus.Get_source())
            return

        tag = mpistatus.Get_tag()
        if tag == mpi.Tags.LOG:
            # print('data:', data)
            logger.handlers[0].emit(data)
        elif tag == mpi.Tags.SHUTDOWN:
            self.numworkers -= 1
            logger.info(
                'Worker %s was removed from worker-list: '
                'There are %s (out of %s) left...', 
                mpistatus.Get_source(),
                self.numworkers,
                mpi.size-1)

        elif tag == mpi.Tags.RESULTS:
            self.db.add_row(data)

        elif tag == mpi.Tags.DONE:
            logger.debug('recv mpi.Tags.DONE from %s',
                         mpistatus.Get_source())
            if data == 'START':
                pass
            else:
                # self.done.append(data)
                pass

            if len(self.inputlist) == 0:
                logger.info('send SHUTDOWN to %s', mpistatus.Get_source())
                self.comm.send('SHUTDOWN', mpistatus.Get_source(),
                               tag=mpi.Tags.INPUTS)
            else:
                data = self.inputlist.pop()
                self.comm.send(data, mpistatus.Get_source(),
                               mpi.Tags.INPUTS)

            logger.info('number jobs left on queue %s', len(self.inputlist))
            if len(self.inputlist) < 10:
                logger.debug('jobs left (list) %s', self.inputlist)
        elif tag == mpi.Tags.IO_REQUEST:
            self.io_queue.append(mpistatus.Get_source())
            logger.debug('%s into io-queue (%s)', mpistatus.Get_source(),
                         len(self.io_queue))
        elif tag == mpi.Tags.IO_FINISHED:
            self.io_tickets[mpistatus.Get_source()] = None
            ctr = self.io_tickets.count(mpi.Tags.IO_TICKET)
            logger.debug('%s release ticket. concurrency: %s',
                         mpistatus.Get_source(), ctr)
        else:
            logger.critical('got data w/ unknown tag from %s',
                            mpistatus.Get_source())
            if DEBUG:
                raise (Exception, 'unknown tag')

    def _manage_io(self):
        """
        Manage I/O - queue.

        @return:   0 if I/O managed
                   0.1 if slept

        This method is organizing together with self.io_tickets the
        input/output of the slave nodes.

        ::notes::
        If no I/O required or no tickets available, sleep 0.1s
        """
        used_tickets = self.io_tickets.count(mpi.Tags.IO_TICKET)
        if used_tickets < PARALLEL_IO and len(self.io_queue) > 0:
            worker_rank = self.io_queue.pop(0)
            self.io_tickets[worker_rank] = mpi.Tags.IO_TICKET
            if used_tickets > (PARALLEL_IO - 2):
                logger.debug("%s recv ticket. concurrencty: %s",
                             worker_rank, used_tickets+1)
            self.comm.send('', worker_rank, mpi.Tags.IO_TICKET)
        else:
            time.sleep(0.1)
            return 0.1
        return 0


def main(inputs, alpha=None, hij=None, force_map=None, simpackdir=None):
    """ Ensemble Start, Divides Work on Ranks """
    try:
        tmp = os.environ["CADEE_TMP"]
        if tmp == '':
            raise KeyError
    except KeyError:
        for tmp in ['/scratch/', '/tmp/', '/dev/shm']:
            if os.path.isdir(tmp):
                break

    tempdir = os.path.join(tmp, 'cadee')
    tempdir = os.path.join(tempdir, str(os.getpid()))

    logger.info("Working directory of rank %s: %s", mpi.rank, tempdir)

    if mpi.rank == 0:
        start = time.time()
        if simpackdir is None:
            raise Exception('Simpackdir is not defined on rank0.')
        io_rank = Master(tempdir, start, simpackdir, force_map=force_map)
        for each in inputs:
            # TODO: remove each, each
            io_rank.inputlist.append([each, each])
        try:
            io_rank.run()
        except Exception as err:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            logger.fatal('FATAL EXCEPTION: Master Died! %s', err)
            traceback.print_exception(exc_type, exc_value, exc_traceback)
            mpi.comm.Abort(1)
            raise
        logger.info("TOTALTIME: %s s", round(time.time() - start, 1))
    else:
        while True:
            try:
                Worker(tempdir, alpha, hij, force_map).run()
                break
            except KeyboardInterrupt:
                break
            except Exception as err:
                logger.exception('Rank %s raised exception %s', mpi.rank, err)
                if RAISE_EXCEPTIONS:
                    raise
                else:
                    continue
        logger.warning("COMPUTATION LOOP ENDED: %s", mpi.rank)


def priorize(inputs):
    prio0 = []
    prio1 = []
    prio2 = []
    prio3 = []
    prio9 = []
    for inp in inputs:
        if '_0' in inp:
            prio0.append(inp)
        elif '_1' in inp:
            prio1.append(inp)
        elif '_2' in inp:
            prio2.append(inp)
        elif '_3' in inp:
            prio3.append(inp)
        else:
            prio9.append(inp)
    inputs = []
    inputs.extend(prio0)
    inputs.extend(prio1)
    inputs.extend(prio2)
    inputs.extend(prio3)
    inputs.extend(prio9)
    logger.info('Prioritized')
    #for each in inputs:
    #    print(each)
    return inputs


def parse_args():
    # TODO: load defaults from somewhere
    parser = argparse.ArgumentParser('CADEE: simpack computation.')

    # Minimum Inputfiles needed
    parser.add_argument('simpackdir', action='store',
                        help='Path to folder with simpacks.')
    # Mapping
    parser.add_argument('--alpha', action='store', default=None,
                        help="Alpha to use for mapping.")

    parser.add_argument('--hij', action='store', default=None,
                        help="Hij to use for mapping.")

    parser.add_argument('--force_map', action='store_true', default=False,
                        help='forced remapping')

    args = parser.parse_args()

    simpackdir = args.simpackdir
    simpackdir = os.path.abspath(simpackdir)

    if not os.path.isdir(simpackdir):
        raise argparse.ArgumentTypeError('Not a folder:{0}'.format(simpackdir))

    if args.hij is None and args.alpha is not None:
        raise argparse.ArgumentTypeError(
                '--hij and --alpha must be either both set, or neither')

    if args.hij is not None and args.alpha is None:
        raise argparse.ArgumentTypeError(
                '--hij and --alpha must be either both set, or neither')

    if args.hij is None and args.force_map is True:
        raise argparse.ArgumentTypeError(
                '--force_map alone is invalid, you must also set --hij and --alpha')

    alpha = None
    hij = None

    if args.hij:
        hij = check_int_or_float(args.hij)

    if args.alpha:
        alpha = check_int_or_float(args.alpha)

    if not mpi.mpi:
        raise Exception('MPI not available')

    if mpi.rank == mpi.root:

        logger.info(
            'Settings: '
            'Path: %s, '
            'Alpha: %s, '
            'Hij: %s, '
            'Force mapping: %s.',
            simpackdir, alpha, hij, args.force_map)

        inputs = []
        if os.path.isdir(simpackdir):
            # we got a folder to scan & we are rank0!
            wd = os.getcwd()
            os.chdir(simpackdir)
            for fil in os.listdir(simpackdir):
                if fil[-4:] == '.tar':
                    inputs.append(os.path.abspath(fil))
                    logger.info('Add inputfile %s', fil)
            os.chdir(wd)

            inputs = priorize(inputs)

        main(inputs, alpha, hij, args.force_map, simpackdir=simpackdir)
    else:
        main(None, alpha, hij, args.force_map)

if __name__ == "__main__":
    parse_args()
