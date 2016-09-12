#!/usr/bin/env python

"""
The script reads extracts a database produced by CADEE
and exports the data in a CSV file for custom analysis.

usage: python extract_to_csv_us.py cadee.pdb us.csv
       this will create your us.csv file (spreadsheet)

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

if len(sys.argv) != 3:
    print("Usage: \n  " + os.path.basename(__file__) + " cadee.db output.csv")
    sys.exit(1)
elif not os.path.lexists(sys.argv[1]):
    print("File %s does not exist!" % sys.argv[1])
    sys.exit(1)


# connect and get values from DB
conn = sqlite3.connect(sys.argv[1])
cursor = conn.cursor()
try:
    cursor.execute("SELECT mutant,replik,name,barr_forw FROM results WHERE feptype=?", (RUN_TYPE,) )
except sqlite3.DatabaseError as e:
    print("Error when accesing the database: '%s' (%s)" % (sys.argv[1], e))
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

open(sys.argv[2], "w").write("\n".join(csv))
print("Success... Wrote %s..." % sys.argv[2])
