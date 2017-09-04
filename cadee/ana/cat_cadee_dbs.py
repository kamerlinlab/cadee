#!/usr/bin/env python

"""
This program will concatenate two cadee.db files together

Usage: python fuse_dbs.py cadee1.db cadee2.db [cadee3.db [...]]
       Will create a concat_cadee.db file, containing both [all] database data in one file.

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""

from __future__ import print_function

import sys
import os
import sqlite3
from cadee.dyn.tools import SqlDB

# TODO: Write module to read data from sqldatabase and display in browser

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

OUTFILE = 'concat_cadee.db'


def usage():
    """Print Usage and Exit."""
    print ("Usage: ", os.path.basename(__file__), ' cadee1.db cadee2.db cadee3.db ...')
    print ("       Will create a concat_cadee.db file, containing both [all] database data in one file.")
    sys.exit(1)


def main(db_list):
    """
    param: db_list is a list of cadee.db files, i.e. ['cadee1.db', 'cadee2.db', ...]
    output: creates file OUTFILE (default: concat_cadee.db)
    error: will quit if OUTFILE exists
    """
    for dbfile in db_list:
        # connect and get values from DB
        conn = sqlite3.connect(dbfile)
        cursor = conn.cursor()
        # select only 'medium' runs
        try:
            cursor.execute("SELECT * FROM results")
        except sqlite3.DatabaseError as e:
            print("Error when accessing the database: '{}' ({})".format(
                dbfile, e))
            sys.exit(1)

        results = cursor.fetchall()
        conn.close()

        db = SqlDB(OUTFILE, 10000)

        for res in results:
            res = list(res)
            res[1] = res[1].split('_')[0]
            db.add_row(res)

        db.flush()


if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()

    if os.path.exists(OUTFILE):
        print('ERROR: Outputfile', OUTFILE, 'exists. Please remove the file and try again.')
        sys.exit(2)

    main(sys.argv)

