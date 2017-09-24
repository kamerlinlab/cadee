#!/usr/bin/env python2

### Original: msrescale.py
###
### matej.repic@ki.si, Oct 2013
##
### Rescale the charges obtained from QM and make them suitable for molaris.
### The script needs the amino lib entry of the residue in a separate file.
### If the absolute sum of old charges in a group exceeds the threshold
### parameter, the script prompts the user for a desired charge.

#
# Modified by yours truly to work with Q
#


import os
import sys
import time
import re
from shutil import copyfile
try:
    import argparse
except ImportError:
    import lib.argparse as argparse
from lib.common import backup_file, __version__
from qscripts_config import QScriptsConfig

parser = argparse.ArgumentParser()
parser.add_argument("lib_file", help = "Q library file (single entry only)")
parser.add_argument("-t", dest="threshold", type=float, help="Charge group charge threshold for user prompt of net charge. Default is 0.4.", default=0.4)

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

if not os.path.lexists(args.lib_file):
    print "FATAL! File %s doesn't exist." % args.lib_file
    sys.exit(1)

# qlib file parser 
def parseQlibFile(inp_file):

    section=""

    atom_dict = {}
    cg_list = []
    cg_allatoms = []  # to check for bad entries

    lines = open(inp_file, 'r').readlines()
    
    for l in lines:
        l = re.split("#|\!", l)[0].strip()
        if l == "": # ignore empty lines
            continue
        if l[0] == "{":
            resid_name = l.strip('{}')
            continue
        # if a new section begins, set the section variable 
        if l[0] == "[":    
            section = l.strip("[]")
            if section not in ( "atoms", "bonds", "impropers", "charge_groups" ):
                print "Unsupported section: %s\n" % section
                sys.exit(1)
            continue
        if section == "atoms":
            atom = l.split()
            atom_dict[ atom[1] ] = [ atom[0], atom[2], float(atom[3]) ]             # key=name, value=[ index, ff type, charge ]

        elif section == "charge_groups":
            cgrp = l.split()
            for atom_name in cgrp:
                if atom_name not in atom_dict.keys():
                    print "\nAtom %s is not defined in the [atoms] section. Fix it please." % (atom_name)
                    sys.exit(1)
            cg_list.append( cgrp )
            cg_allatoms.extend( cgrp )

    for atom_name in atom_dict.keys():
        if atom_name not in cg_allatoms:
            print "\nAtom %s is not defined in the [charge_groups] section. Fix it please." % (atom_name)
            sys.exit(1)

    return resid_name, atom_dict, "", cg_list



# rescale the charge
def rescale(_oc_, _ocs_, _ocas_, net_crg):

    # Equation provided by Ram, round to two decimals (molaris reads only two)
    # Q reads 4
    try:
        new_crg = round( _oc_ - abs(_oc_) * ( _ocs_ - net_crg ) / _ocas_ , 5 )
        return new_crg

    except ZeroDivisionError:
        return _oc_


#parse the Qlib file
rn, atom_dict, bl, cgrp_atom_list = parseQlibFile(args.lib_file)
print cgrp_atom_list

cgrp_nc_dict = {}    # net charges for charge groups

# sum charges to target for each en group
for cgrp in cgrp_atom_list :

    target = 0

    # construct a list of old charges and old absolute charges
    oc = [ atom_dict[atom][2] for atom in cgrp ]
    oca = [ abs( atom_dict[atom][2] ) for atom in cgrp ]

    # ask user for charge if the absolute sum exceeds the threshold
    if abs(sum(oc)) > args.threshold:
        print "\n%10.5f charge for group %s" % ( sum(oc), " ".join(cgrp) )
        while True:
            try:
                target = int(raw_input("Specify target charge:"))
            except ValueError:
                print "Error: Non-integer charge"
                continue
            except KeyboardInterrupt:
                print "\n\nExiting. Goodbye!"; sys.exit()
            break

    # append new charge to atom_dict
    for atom in cgrp:
        atom_dict[atom].append( rescale( atom_dict[atom][2], sum(oc), sum(oca), target ) )

    # construct lists of temporary charges and temporary absolute charges
    tc = [ round( atom_dict[atom][3], 5 ) for atom in cgrp ]
    tac = [ abs( round( atom_dict[atom][3], 5 ) ) for atom in cgrp ]

    # find largest absolute charge
    max_abs_index = int( tac.index( max(tac) ) )
    mai = max_abs_index

    # correct the excess charge
    atom_dict[cgrp[mai]][3] -= ( sum(tc) - target )

    nc = sum([ atom_dict[atom][3] for atom in cgrp ])

    cgrp_nc_dict["_".join(cgrp)] = nc


# modify the Qlib file

section=""

lib_lines = open(args.lib_file, 'r').readlines()
new_lib = []

for line in lib_lines:
    line = line.strip()
    l = re.split("#|\!", line)
    try:
        comment = "#  " + "#".join( l[1:] ).strip()
    except IndexError:
        comment = ""
    l = l[0].strip()

    if l == "": # ignore empty lines
        new_lib.append(line)
        continue
    if l[0] == "{":
        new_lib.append(line)
        new_lib.append("# Rounded and rescaled with q_rescale.py (v%s) (%s) " % (__version__, time.ctime() ))
        continue
    # if a new section begins, set the section variable 
    if l[0] == "[":    
        new_lib.append(line)
        section = l
        continue
    if section == "[atoms]":
        atom = l.split()
        index, atom_name, ff_type, old_charge = atom[0:4]
        new_charge = atom_dict[atom_name][3]

        new_line = "%5s  %-5s %-15s %10.5f   # %10s  %s" % (index, atom_name, ff_type, new_charge, old_charge, comment)
        new_lib.append(new_line)

    elif section == "[charge_groups]":
        cgrp = l.split()
        new_line = "%s      # net: %.1f" % (" ".join(cgrp), cgrp_nc_dict["_".join(cgrp)])

        new_lib.append(new_line)
    else:
        new_lib.append(line)


backup = backup_file(args.lib_file)
if backup:
    print "Backed up '%s' to '%s'" % (args.lib_file, backup)
open(args.lib_file, 'w').write( "\n".join(new_lib) )
print "The library file was successfully modified...\n"






