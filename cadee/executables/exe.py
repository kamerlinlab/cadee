#!/usr/bin/env python
"""
This module looks up an executable.

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""

import os

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"


def which(program, enforce_subdir=False):
    """ look up if executable/program exists

        :param program: name of executable to find
        :param enforce_subdir: only search for executables in enclosed with this package
        :type program: str
        :type enforce_subdir: bool
        :return: path to executable or None
        :type return: str or None
        ::Note::
        inspired by:
        http://stackoverflow.com/questions/377017/test-if-executable-exists-in-python
    """
    def is_exe(fpath):
        return os.path.isfile(fpath) and os.access(fpath, os.X_OK)

    # find executable in folder where this script is located
    script_path = os.path.realpath(__file__)
    path = os.path.dirname(script_path)
    exe_file = os.path.join(path, program)
    if is_exe(exe_file):
        return exe_file

    # find q-executable in subfolder where this script is located
    script_path = os.path.realpath(__file__)
    path = os.path.dirname(script_path) + '/q/'
    exe_file = os.path.join(path, program)
    if is_exe(exe_file):
        return exe_file
    # special case for upper-case Qdyn etc...
    exe_file = os.path.join(path, program[0].upper()+program[1:])
    if is_exe(exe_file):
        return exe_file

    if enforce_subdir:
        return None

    fpath, fname = os.path.split(program)
    fpath = False

    if fpath:
        if is_exe(program):
            return program
    else:
        for path in os.environ["PATH"].split(os.pathsep):
            path = path.strip('"')
            exe_file = os.path.join(path, program)
            if is_exe(exe_file):
                return exe_file

    return None
