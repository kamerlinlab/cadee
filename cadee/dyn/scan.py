#!/usr/bin/env python

"""
Scan directory and create inputlist/tree for Q.

Author: {0} ({1})

This module is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


from __future__ import print_function
import os

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

ALPHABETIC = True


class Scan():
    """ Class to scan a folder for Q input files. """
    # TODO: Better Error Handling
    SKIP_INTERMEDIATES = False
    SKIP_INTERMEDIATES = True
    INDENT = '  '

    class OutputOverWriteWarning(UserWarning):
        pass

    @staticmethod
    def read_input_file(fn):
        """
        open fn and read lines.
        return non-empty lines w/o comment, NLC.
        """
        lines = []
        for line in open(fn, 'r'):
            if line.strip() == "":
                continue
            line = line.replace('!', '#')
            line = line.split('#')[0]
            line = line.strip()
            if line == "":
                continue
            lines.append(line)
        return lines

    @staticmethod
    def get_key(line):
        """ return keyword of provided line, or '' """
        line = line.strip()
        if line == "":
            return ''
        return line.split()[0].strip().lower()

    @staticmethod
    def get_value(line):
        """ return value of provided line, or '' """
        line = line.strip()
        if line == "":
            return ''
        parts = line.split()
        if len(parts) < 2:
            return ''
        return parts[1].strip()

    @staticmethod
    def get_simtime(lines):
        """ return steps and stepsize of fn """
        steps = 0
        ss = 0.
        for line in lines:
            if Scan.get_key(line) == "stepsize":
                ss = float(Scan.get_value(line))
            elif Scan.get_key(line) == "steps":
                steps = int(Scan.get_value(line))
            if steps > 0 and ss > 0:
                break
        return steps, ss

    @staticmethod
    def get_all_io_file_names(lines):
        """locate final, restart, toplogy, restraint, fep in lines
        @input: lines
        @return: final, restart, topology, restraint and fep"""
        final = None
        restart = None
        topology = None
        restraint = None
        fep = None
        energy = None
        trajectory = None

        files = False
        for line in lines:
            if '[files]' in line.lower():
                files = True
                continue
            elif files and '[' in line:
                files = False
            if not files:
                continue
            key = Scan.get_key(line).lower().strip()
            if "final" == key:
                final = Scan.get_value(line)
            elif "restart" == key:
                restart = Scan.get_value(line)
            elif "topology" == key:
                topology = Scan.get_value(line)
            elif "restraint" == key:
                restraint = Scan.get_value(line)
            elif "fep" == key:
                fep = Scan.get_value(line)
            elif "energy" == key:
                energy = Scan.get_value(line)
            elif "trajectory" == key:
                trajectory = Scan.get_value(line)
            elif '' == key:
                continue
            else:
                print (key)
                raise (Exception, 'unknown keyword'+key)

        if topology is None:
            print(os.getcwd())
            print(lines)
            raise (Exception, 'WTF: topology keyword *MUST* be set in input!')

        return final, restart, topology, restraint, fep, energy, trajectory


    def _prepare_filelist(self):
        """
        search for files
           AND
        search and read input files ('.inp') into memory.
        """
        self.files = {}
        self.restart = {}
        self.final = {}
        self.topology = {}
        self.restraint = {}
        self.fepfile = {}
        self.energy = {}
        self.trajectory = {}
        self.stepsize = {}
        self.steps = {}

        def isint(txt):
            """Return True if @param txt, is integer"""
            try:
                int(txt)
                return True
            except TypeError:
                return False
            except ValueError:
                return False

        for fn in os.listdir(self.folder):
            if '.inp' != fn[-4:]:
                lines = []
            else:
                num = fn.split("_")[0]
                if not isint(num):
                    print('skipping', fn)
                    # TODO properly check if this file is a qdyn input file
                    #      (eg grep for [MD])
                    lines = []
                else:
                    # TODO: check if the input file is a Qdyn6 input file
                    #       (eg. is first non-comment line [MD]?
                    lines = Scan.read_input_file(fn)
                    (self.final[fn], self.restart[fn],
                            self.topology[fn], self.restraint[fn],
                            self.fepfile[fn], self.energy[fn],
                            self.trajectory[fn]) = Scan.get_all_io_file_names(lines)  # NOPEP8
                self.steps[fn], self.stepsize[fn] = self.get_simtime(lines)
            self.files[fn] = lines

    def _find_file_with_name(self, fn):
        """ Find file with name fn.

        If exists:
            return /full/path/to/file.inp
        If not try:
            return /full/path/to/another/file.inp
        if not:
            return None
        """
        if fn in self.files:
            return fn

        # could not find file with the exact name 'fn'
        # try to find file with same basename...

        fn = os.path.basename(fn)
        for key in self.files.keys():
            if fn in key:
                print ('WARNING:  located file', fn, 'with heuristic search!')
                return key

        return None

    def _find_files_with_keyword(self, keyword):
        """ Return set of files with keyword keyword"""
        files = []
        for fn in self.files.keys():
            for line in self.files[fn]:
                if Scan.get_key(line) == keyword:
                    files.append(fn)
        return set(files)

    def _get_inputfiles_restarting_from(self, re):
        files = []
        for fn in self.restart:
            if self.restart[fn] == re:
                files.append(fn)
        return set(files)

    def _walk_files_from(self, start, indent="", path=0, threads=0, stuff=[], visited=[]):    # NOPEP8
        me = "{0} {1:04d} {2} {3}".format(path, len(indent)/len(Scan.INDENT), indent, start)

        if stuff is not None:
            idx1 = len(indent)/len(Scan.INDENT)
            idx2 = path
            ok = False
            while not ok:
                ok = True
                if len(stuff) <= idx1:
                    stuff.append([])
                    ok = False
            ok = False
            while not ok:
                ok = True
                if len(stuff[idx1]) <= idx2:
                    stuff[idx1].append([])
                    ok = False
            stuff[idx1][idx2].append(start)

        if me in visited:
            return stuff
        visited.append(me)

        if not Scan.SKIP_INTERMEDIATES:
            # write a neat list of files (indented) to stdout
            print(me)

        indent += Scan.INDENT
        fn = start
        final = self.final[fn]
        if final is None:
            return stuff
        f = final[:-3]
        nxt = sorted(self._get_inputfiles_restarting_from(f + ".re"))
        for i, f in enumerate(nxt):
            ifn = f[:-4]+'.inp'
            if i > 0:
                threads += 1
            self._walk_files_from(ifn, indent, path, threads, stuff, visited)
        return stuff

    def __init__(self, workdir=None):
        if workdir is None:
            self.folder = os.getcwd()
        else:
            self.folder = workdir
            os.chdir(self.folder)

        if self.folder[-1] != '/':
            self.folder += '/'

    def _check_writeable_files(self, fn):
        """ raise Exception if any outputfile (ene, final, dcd) is written 2x"""
        for d in self.final, self.energy, self.trajectory:
            if d[fn] is not None and d[fn] in self.usedfiles:
                warnmsg = " InputFileError: Outputfile: '" + str(d[fn])
                warnmsg += "' written by '" + str(self.usedfiles[d[fn]])
                warnmsg += "' is overwritten by '" + str(fn) + "'"
                raise OutputOverWriteWarning(warnmsg)
            else:
                self.usedfiles[d[fn]] = fn

    def scan(self):
        """ Scan folder for input files. Only 1 START
        (independent input files strand) is supported.
        """


        if ALPHABETIC:
            self._prepare_filelist()
            fs = 0.
            sum_steps = 0.
            start = list(self._find_files_with_keyword('initial_temperature'))[0]
            final, restart, topology, restraint, fepfile, energy, trajectory = self.get_all_io_file_names(open(start).readlines())
            files = sorted(list(self._find_files_with_keyword('steps')))
            for fil in files:
                steps, ss = self.get_simtime(open(fil).readlines())
                if ss == 0:
                    print(fil, ss, 'WTF')
                    die()
                fs += steps/ss
                sum_steps += steps
            self.result = (topology, fepfile, files, sum_steps, fs)
            return self.result

        self._prepare_filelist()

        starts = self._find_files_with_keyword('initial_temperature')
        if len(starts) == 0:
            raise (Exception, 'NO STARTS FOUND')
        restarters = self._find_files_with_keyword('restart')
        i = 0
        for start in starts:
            if len(starts) > 1 and start in restarters:
                print ('skipping file', start)
            else:
                if i > 0:
                    print ('WARNING: MORE THAN ONE START FOUND! NOT SUPPORTED!!')
                    raise (Exception, 'MORE THAN 1 START FOUND')
                    break
                stuff = self._walk_files_from(start, "", i, 0)
                i += 1

        fd = ''
        steps = 0
        fs = 0.
        inputfiles = []
        self.usedfiles = {}
        for fn in stuff:
            fil = os.path.basename(fn[0][0].replace(self.folder, ''))
            fn = fn[0][0]
            steps += self.steps[fn]
            fs += self.steps[fn] * self.stepsize[fn]
            inputfiles.append(fil)

            self._check_writeable_files(fn)

            # Consistency tests
            if fd == '':
                # is everything in 1 single folder, no subfolder?
                fd = fn
            elif os.path.dirname(fd) != os.path.dirname(fn):
                print('WARNING: subfolders not supported!;', fd, fn)

        # TODO: Check if restart/restraint or any other files are OVERWRITTEN


        # Consistency Tests

        fv = filter(None, set(self.fepfile.values()))

        if len(fv) > 1:
            print('WARNING: more than 1 fepfile not supported')
            print(fv)
            raise (Exception, 'mulitple fepfiles')
        else:
            fepfile = fv[0]
        if len(set(self.topology.values())) > 1:
            print(set(self.topology.values()))
            print('WARNING: more than 1 topology not supported')
        else:
            topology = self.topology.values()[0]

        if fepfile is not None and not os.path.isfile(fepfile):
            print(fp)
            raise (Exception, 'fepfile does not exists')

        if not os.path.isfile(topology):
            print(ft)
            raise (Exception, 'topology does not exists')

        self.result = (topology, fepfile, inputfiles, steps, fs)
        return self.result

if __name__ == "__main__":
        import time
        import sys
        start = time.time()
        if len(sys.argv) > 1:
            s = Scan(os.path.abspath(sys.argv[1]))
        else:
            s = Scan()
        print(s.scan())
        print('timing', time.time()-start)

