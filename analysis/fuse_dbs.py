#!/usr/bin/env python

"""
This program will fuse to cadee.db files together

Usage: python fuse_dbs.py cadee1.db cadee2.db
       will create a fuse.db file, containing both database data

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.
"""


import sys
import os
import sqlite3

# TODO: Write module to read data from sqldatabase and display in browser

__author__ = "Beat Amrein"
__email__ = "beat.amrein@gmail.com"

class SqlDB(object):
    def __init__(self, name, interval=10):
        """Connect to database and initialize table if not exists
        :param name: path to database
        :param interval: interval to commit changes
        :type name: str
        :type interval: int
        """
        self.conn = sqlite3.connect(name)
        self.cursor = self.conn.cursor()
        self.cursor.execute('''CREATE TABLE IF NOT EXISTS results
        (time int, mutant text, replik int, name text, feptype text, barr_forw real, exo real, barr_back real, ttot real, tfree real, tfreesolute real, tfreesolvent real,  ene_kin real, ene_pot real, ene_tot real); ''')  # NOPEP8
        self.counter = 0
        self.commit_interval = interval
        self.conn.commit()
        self.template = 'INSERT INTO results VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)'    # NOPEP8

    def flush(self):
        self.commit()

    def commit(self):
        print('flushing')
        self.conn.commit()
        self.counter = 0

    def add_row(self, results):
        try:
            self.cursor.execute(self.template, results)
        except ValueError as e:
            print('Unable to store rows; ValueError %s' % e)
            print('template: %s' % self.template)
            print('results:  %s' % results)
        self.counter += 1
        if (self.counter % self.commit_interval == 0):
            self.commit()

    def close(self):
        self.commit()


def usage():
    print ("Usage: ", os.path.basename(__file__), ' file1 file2 file3 ...')
    sys.exit(1)


if __name__ == "__main__":
    if len(sys.argv) < 2:
        usage()

    for dbfile in sys.argv[1:]:
        # connect and get values from DB
        conn = sqlite3.connect(dbfile)
        cursor = conn.cursor()
        # select only 'medium' runs
        try:
            cursor.execute("SELECT * FROM results")
        except sqlite3.DatabaseError as e:
            print("Error when accesing the database: '{}' ({})".format(
                sys.argv[1], e))
            sys.exit(1)

        results = cursor.fetchall()
        conn.close()

        db = SqlDB('fuse.db', 10000)

        for res in results:
            res = list(res)
            res[1] = res[1].split('_')[0]
            db.add_row(res)

        db.flush()
