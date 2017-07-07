#!/usr/bin/env python

"""
Module for coordinating mpi

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

try:
    from mpi4py import MPI
except ImportError:
    print('mpi4py not found')

try:
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    root = 0
    size = comm.Get_size()
    mpi = True
except NameError:
    comm = 0
    rank = 0
    root = 0
    size = 0
    mpi = False
    print('MPI disabled')


class Tags(object):
    """ MPI tags """
    DONE = 1
    INPUTS = 2
    LOG = 3
    IO_TICKET = 4
    IO_REQUEST = 5
    IO_FINISHED = 6
    RESULTS = 7
    SHUTDOWN = 8


def get_info():
    ret = "MPI Info: "
    ret += " enabled: " + str(mpi) 
    ret += " rank: " + str(rank) 
    ret += " size: " + str(size)
    return ret
