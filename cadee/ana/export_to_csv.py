#!/usr/bin/env python

"""
The script reads extracts a database produced by CADEE
and exports the data in a CSV file for custom analysis.

usage: python extract_to_csv_medium.py cadee.pdb medium.csv
       this will create your medium.csv file (spreadsheet)

Author: {0} ({1})

This program is part of CADEE, the framework for
Computer-Aided Directed Evolution of Enzymes.

"""


from __future__ import print_function
import sys
import os
import sqlite3

__author__ = "Beat Amrein, Miha Purg"
__email__ = "beat.amrein@gmail.com, miha.purg@gmail.com"

RUN_TYPE="us"    

def main(args, what='barr_forw'):
    """
    :param args: ['cadee.db', 'output.csv']     *list of filenames*:
    :param what: 'barr_forw' (default) or 'exo' *string*:
    :return: void
    """

    if len(args) != 3:
        print("Invalid arguments", args)
        print()
        print("Usage: \n  " + os.path.basename(__file__) + " cadee.db output.csv")
        sys.exit(1)

    db=args[1]
    if not os.path.lexists(db ):
        print("File %s does not exist!" % db)
        sys.exit(1)

    outcsv=args[2]

    if os.path.exists(outcsv):
        print("File %s exists. Please remove it and retry." % outcsv)
        sys.exit(2)

    if what == 'barr_forw':
        print('Exporting Barrier')
    elif what == 'exo':
        print('Exporting deltaG')
    else:
        print("Unknown Column %s. Please use either 'barr_forw' or 'exo'", what)
        sys.exit(3)

    # connect and get values from DB
    conn = sqlite3.connect(db)
    cursor = conn.cursor()
    try:
        cursor.execute("SELECT mutant,replik,name,? FROM results WHERE feptype=?", (what, RUN_TYPE))
    except sqlite3.DatabaseError as e:
        print("Error accessing the database: '%s' (%s)" % (db, e))
        sys.exit(1)

    results = cursor.fetchall()
    conn.close()

    data = {}
    names = set()   # filenames (1100_eq, 1200_eq), set to make a list of unique values from all the mutants (in case some are missing)

    maxreplik = 0
    # get WT averages
    for res in results:
        mutant, replik, name, barr = res
        if replik > maxreplik:
            maxreplik = replik
        mutant = mutant.split('_')[0]
        if mutant not in data:
            data[mutant] = {}
        if replik not in data[mutant]:
            data[mutant][replik] = {}
        data[mutant][replik][name] = barr
        names.add(name)

    csv = []
    muts = sorted(data.keys())
    csv.append("mutant;replik;" + ";".join(muts))
    for name in sorted(list(names)):

        for replik in range(maxreplik+1):
            values = []
            for mut in muts:
                i = data[mut].get(replik, {})
                if len(i) == 0:
                    values.append("")
                    print('WARNING: empty replik', mut, replik)
                else:
                    j = i.get(name, "")
                    if j == "":
                        values.append("")
                        print('info: empty name', mut, replik, name)
                    else:
                        values.append(str(j))
            csv.append(name + ";" + str(replik) + ";" + ";".join(values))

    open(outcsv, "w").write("\n".join(csv))
    print("Success... Wrote %s..." % outcsv)

if __name__ == "__main__":
    main(sys.argv[1:])
