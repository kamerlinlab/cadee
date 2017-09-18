#!/usr/bin/env python

"""
Tool Collection for ensemble-Q

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function
import os
import stat
import gzip
import hashlib
import logging
import time
import sqlite3

import mpi

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

NLC = '\n'


class cd:
    """Context manager for changing the current working directory
    """
    # http://stackoverflow.com/questions/431684/how-do-i-cd-in-python/13197763#1319776

    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)


def filname(fil):
    """strip folder off path, returning filename"""
    while fil[-1] == "/":
        fil = fil[:-1]
    return os.path.split(fil)[1]


class Results(object):
    """
    Very complicated object to hold information on results.
    Feel free to simplify, aka kill subclasses.
    """
    class _EVB(object):
        def __init__(self, barr, exo):
            self.barr = barr
            self.exo = exo
            self.barr_rev = barr - exo

        def __repr__(self):
            return "dGa: {0:5.2f} dG0: {1:5.2f} dGa_rev: {2:5.2f}".format(
                self.barr, self.exo, self.barr_rev)

        def items(self):
            """ Returns List
            [ dG*, dG0, dG*_rev ]
            """
            return [self.barr, self.exo, self.barr_rev]

    class _Temp(object):
        def __init__(self, tot, free, free_solvent, free_solute):
            self.tot = tot
            self.free = free
            self.free_solvent = free_solvent
            self.free_solute = free_solute

        def __repr__(self):
            return "Temps: ttot: {0}, free: {1}, solv: {2}, solu: {3}".format(
                round(self.tot, 1), round(self.free, 1),
                round(self.free_solvent, 1), round(self.free_solute, 1))

        def items(self):
            """" Returns List
            [ Ttot, Tfree, Tfree_solvent, Tfree_solute ]
            """
            return [self.tot, self.free, self.free_solvent, self.free_solute]

    class _SumEne(object):
        def __init__(self, tot, pot, kin):
            self.tot = tot
            self.pot = pot
            self.kin = kin

        def __repr__(self):
            return "Energies: Total: {0}, Potential: {1}, Kinetic: {2}".format(
                round(self.tot, 1), round(self.pot, 1), round(self.kin, 1))

        def items(self):
            """" Returns List
            [ Total Energy, Potenial Energy, Kinetic Energy ]
            """
            return [self.tot, self.pot, self.kin]

    def __init__(self, mutant, replik, name, feptype):
        self.time = int(time.time() * 1000)
        self.mutant = mutant
        self.replik = replik
        self.name = name
        self.feptype = feptype
        self.evb = None
        self.t = None
        self.enesum = None

    def dg(self, barr, exo):
        """ Add dG* and dG0 """
        self.evb = self._EVB(barr, exo)

    def temp(self, tot, free, free_solvent, free_solute):
        """ Add Temperatures """
        self.t = self._Temp(tot, free, free_solvent, free_solute)

    def ene(self, kin, pot, tot):
        """ Add energies """
        self.enesum = self._SumEne(tot, pot, kin)

    def __repr__(self):
        ret = 'Time: {0}, Name: {1}, FepType: {2}'.format(
            self.time, self.name, self.feptype)
        for obj in self.evb, self.t, self.enesum:
            if obj is not None:
                ret += obj.__repr__()
        return ret

    def items(self):
        """ Return list of items """
        items = [self.time, self.mutant, self.replik, self.name, self.feptype]
        for obj in self.evb, self.t, self.enesum:
            if obj is not None:
                items.extend(obj.items())
        return items

    def sql(self):
        return self.items()


def report_results(results):
    """
    :param results: object containing results
    :type results: Results
    ::note::
    raises exception if mpi is enabled, but rank0 reports results
    if used outside mpi context, will not reuse database
    """
    if mpi.mpi:
        if mpi.rank != mpi.root:
            mpi.comm.send(results.items(), mpi.root, tag=mpi.Tags.RESULTS)
            return
        else:
            logger.warning('Rank0 reports results!')
            raise Exception('Rank0 reported results!')

    else:
        logger.debug('Not in MPI Session: Inefficient db access.')
        logger.info('results: %s', results.sql())
        db = SqlDB('cadee.db')
        db.add_row(results.sql())
        db.close()


class LogFileHandler(logging.StreamHandler):
    def __init__(self, logfile=None):

        self.logfile = None
        if mpi.rank is None or mpi.rank == 0:
            if logfile is not None:
                if isinstance(logfile, file):
                    self.logfile = logfile
                else:
                    try:
                        self.logfile = open(logfile, 'w')
                    except (IOError, TypeError):
                        self.logfile = None

        logging.StreamHandler.__init__(self)

        msg = 'Logger created: Logfile: %s, MPI: %s, Rank: %s'
        self.emit(msg % (self.logfile, mpi.mpi, mpi.rank))

    def emit(self, record):
        if isinstance(record, str):
            msg = record
        else:
            # print(record)
            msg = self.format(record)

        if mpi.mpi:
            if mpi.rank != 0:
                mpi.comm.send(msg, 0, tag=mpi.Tags.LOG)
                return

        print(msg, file=self.logfile)

    def flush(self):
        if self.logfile is not None:
            if mpi.mpi:
                if mpi.rank == 0:
                    self.logfile.flush()


def getLogger(name, logfile=None, level=logging.INFO):
    """ Returns a logger, if the MPI module is available, with MPI """

    # frmt = '%(asctime)s - %(name)s - %(levelname)s - %s(lineno)s - %(message)s'  # NOPEP8
    frmt = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'

    logger = logging.getLogger(name)

    if len(logger.handlers) > 0:
        return logger

    if mpi.mpi:
        formatter = logging.Formatter("{0:04d} - {1}".format(mpi.rank, frmt))  # NOPEP8
    else:
        formatter = logging.Formatter(frmt)

    lfh = LogFileHandler(logfile)
    lfh.setFormatter(formatter)
    logger.addHandler(lfh)
    logger.propagate = False
    logger.setLevel(level)

    return logger


class File(object):
    """ container for a file """
    def __init__(self, name, data=None, exe=False, exists=True,
                 ramfile=False, read=True):
        """ initalize with the filename.
        @param:
        name:       filename
        optional:
        data:       binary blob
        exe:        True, to set the executable flag
        exists:     True, if the file exists already
        ramfile:    True, to keep data in memory only
        read:       True, to read file to ram right away

        if read is True, exists must be True or Exception is raised
        """
        if read and not exists:
            raise OSError("Can't read file that does not exist.")

        if exists and not os.path.isfile(name):
            raise Exception("No such file on disk:", name)

        head, tail = os.path.split(name)

        logger.debug('head is %s tail is %s', head, tail)

        if head == '':
            self.path = os.getcwd()
        else:
            self.path = head
        if tail == '':
            raise Exception("Empty filename", name, 'in', self.path)
        self.name = tail
        self.data = data

        if self.name[-3:] == '.gz':
            self.gzipped = True
        else:
            self.gzipped = False

        self.executable = exe
        if self.data is None:
            self._inram = False
        else:
            self._inram = True
        self._ondisk = bool(os.path.isfile(self.fullpath))

        if read:
            self.update()

        if ramfile:
            self.update()
            self.path = ''
            self._ondisk = False
            self._inram = True

        if not exists:
            self._ondisk = False
            self._inram = False

        self.check_integritiy()

    def unroot(self):
        """update and remove path"""
        if self._ondisk:
            self.update()
        if not self._inram:
            raise Exception('No data in memory, cannot unroot!')
        self.path = ''
        self._ondisk = False

    def get_lines(self, comments=True, empty=False):
        """return list of non-empty lines"""
        if not os.path.isfile(self.fullpath):
            self.to_disk()

        if self.gzipped:
            loglines = gzip.open(self.fullpath).readlines()
        else:
            loglines = open(self.fullpath).readlines()

        ret = []
        for line in loglines:
            line = line.rstrip()
            if not empty:
                if line.strip() == '':
                    continue
            if not comments:
                line = line.replace('#', '!').split('!')[0]
            ret.append(line)

        return ret

    def check_integritiy(self):
        if self._inram:
            if self.data is None:
                raise Exception('No data, but in ram')
        if self._ondisk:
            if self.path == '':
                raise Exception('_ondisk=True, but no path')
            if not os.path.exists(self.fullpath):
                raise Exception('_ondisk=True, but no file in', self.fullpath)

    @property
    def fullpath(self):
        """return full path of file (path/name)"""
        return os.path.join(self.path, self.name)

    @property
    def memorized(self):
        """bool: in memory?"""
        return self._inram

    def ext(self):
        """returns what's after the last '.' in the filname"""
        return os.path.splitext(self.name)

    def update(self):
        """read from disk, update blob in ram"""
        if self._ondisk is False and not os.path.isfile(self.fullpath):
            raise Exception('Cant update; file not on_disk!', self.fullpath)
        if self.path == '' and not os.path.isfile(self.name):
            raise Exception('Cant update; file not found on disk')
        with cd(self.path):
            if not os.path.isfile(self.name):
                raise Exception('Cant update; file not found on disk')
            with open(self.name, 'rb') as fil_in:
                self.data = fil_in.read()
        self._inram = True
        self.check_integritiy()

    @property
    def sha1(self):
        if self.data is None:
            return None
        return hashlib.sha1(self.data).hexdigest()

    def to_ram(self):
        """read from disk, overwrite ram, delete from disk"""
        if self._inram:
            pass
        elif os.path.isfile(self.fullpath):
            with cd(self.path):
                self.update()
                os.remove(self.name)
        else:
            raise Exception('Fatal: File is not in ram and not ok disk!',
                            self.fullpath)
        self._inram = True
        self._ondisk = False

        self.check_integritiy()

    def to_disk(self, path=None, name=None):
        """read from ram, overwrite disk"""
        if not self._inram:
            raise Exception('Not in ram. Cant write to disk.')
        if path is None:
            path = self.path
        if name is None:
            name = self.name
        if not os.path.isdir(path):
            os.makedirs(path)
        with cd(path):
            with open(name, 'wb') as fil_out:
                fil_out.write(self.data)
            if self.executable:
                os.chmod(self.name, stat.S_IEXEC)
        self._ondisk = True

    def compress(self):
        """gzip self, writes to do disk!"""
        if self.data is None:
            return

        if self.name[-3:] == '.gz':
            return

        with cd(self.path):
            newname = self.name+".gz"
            oldname = self.name
            olddata = self.data

            try:
                self.name = newname
                gzip.open(newname, 'wb').write(self.data)
                self.update()
                self.gzipped = True
            except:
                self.name = oldname
                self.data = olddata
                raise

    def decompress(self):
        if self.data is None:
            return
        if self.name[-3:] == '.gz':
            with cd(self.path):
                newname = self.name[:-3]
                oldname = self.name
                olddata = self.data
            try:
                self.name = newname
                self.data = gzip.open(oldname, 'rb').read()
                open(newname, 'wb').write(self.data)
                self.update()
                self.gzipped = False
            except:
                self.name = oldname
                self.data = olddata
                raise

    def change_path(self, newpath):
        """read file to ram, change path, write file to disk"""
        if not os.path.isdir(newpath):
            os.makedirs(newpath)
        self.to_ram()
        self.path = newpath
        self.to_disk()

    def change_name(self, newname):
        """read file to ram, change name, write file to disk"""
        if not self._inram:
            self.to_ram()
        self.name = newname
        self.to_disk()

    def deploy(self):
        """write file to cwd/name. overwrites existing file"""
        self.write_to(os.getcwd(), self.name, overwrite=True)

    def write_to(self, path, name, overwrite=False):
        """places a copy of the file to newpathname"""
        fil = os.path.join(path, name)
        if os.path.exists(fil):
            if not overwrite:
                raise Exception('File exists!', fil)
            logger.debug('overwriting file %s', fil)
        if not self._inram:
            self.update()
        oldflg = self._ondisk
        self.to_disk(path, name)
        self._ondisk = oldflg

        self.check_integritiy()

    def info(self):
        if self.data is not None:
            size = len(self.data)
            pf = 'B '
            if size > 1024:
                size = size/1024.
                pf = 'KB'
                if size > 1024:
                    size = size/1024.
                    pf = 'MB'
                    if size > 1024:
                        size = size/1024.
                        pf = 'GB'
            sha1 = self.sha1[:8]
            size = size
        else:
            pf = ''
            size = float('nan')
            sha1 = 'N/A'
        return self.name, size, pf, sha1

    def __str__(self):
        info = self.info()
        return "{0:>15}, size: {1:4.0f}{2:>2}, sha1: {3}".format(*info)

    def __repr__(self):
        return "{0}/{1}{2}/{3}".format(*self.info())


class QmapfepSettings(object):
    """qmapfep settings container"""
    def __init__(self, hij, gps, qfep_proc, temp=300,
                 skip=250, minpts_per_bin=20, bins=20):
        """Initialize with
        @Hij        alpha
        @gps        gasphaseshift
        @qfep_proc  qfep_proc
        @temp temperature [K], default 300K"""
        self.hij = hij
        self.gps = gps
        if not isinstance(qfep_proc, File):
            self.qfep_proc = File(qfep_proc, read=True)
        else:
            self.qfep_proc = qfep_proc
        self.temp = temp
        self.skip = skip
        self.minpts_per_bin = minpts_per_bin
        self.bins = bins


class SqlDB(object):
    def __init__(self, name, interval=300):
        """Connect to database and initialize table if not exists
        :param name: path to database
        :param interval: interval (seconds) between committing changes
        :type name: str
        :type interval: int
        """
        self.conn = sqlite3.connect(name)
        self.cursor = self.conn.cursor()
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS results
        (time int, mutant text, replik int, name text, feptype text, barr_forw real, exo real, barr_back real, ttot real, tfree real, tfreesolute real, tfreesolvent real,  ene_kin real, ene_pot real, ene_tot real); ''')  # NOPEP8
        #self.cursor.execute('''CREATE VIEW IF NOT EXISTS avg as
        #SELECT avg(barr_forw), avg(exo), avg(barr_back), avg(ttot), avg(tfree), avg(ene_tot), avg(ene_pot), avg(ene_kin) FROM results; ''')  # NOPEP8
        self.commit_interval = interval
        self.commit()
        self.template = 'INSERT INTO results VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'    # NOPEP8

    def flush(self):
        self.commit()

    def commit(self):
        self.conn.commit()
        logger.info('Committed cadee.db.')
        self.last_commit = time.time()

    def add_row(self, results):
        if isinstance(results, Results):
            results = results.items()
        try:
            self.cursor.execute(self.template, results)
        except ValueError as e:
            logger = getLogger(__name__)
            logger.critical('Unable to store rows; ValueError %s', e)
            logger.critical('template: %s', self.template)
            logger.critical('results:  %s', results)
        except ProgrammingError as e:
            getLogger(__name__).exception('ProgrammingError %s', results)
        if self.last_commit + self.commit_interval < time.time() :
            self.commit()

    def close(self):
        self.commit()


logger = getLogger(__name__)
