#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# MIT License
# 
# Copyright (c) 2016  Miha Purg <miha.purg@gmail.com>
# 
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# 
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# 
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.
#
#
#
#
# This poor excuse for a script converts FFLD output to library and parameter files for Q.
# Two arguments are required: 
# - the FFLD file created with 'ffld_server' (bundled with SCHRODINGER's Maestro)
# - the 3 letter code for the ligand
# And creates three files: .lib, .prm and .prm.chk
# 
# It also takes in an optional argument (-s): the PDB structure from which the FFLD_OUTPUT was created.
# This is used to copy the original atom names from the PDB to the lib and the parm, and to check how
# good the parameters are with regards to the structure (should be a QM optimized structure of course)
#
# Example:  
#   $SCHRODINGER/utilities/ffld_server -ipdb paraoxon.pdb -print_parameters -version 14 > paraoxon.ffld11
#   q_ffld2q.py paraoxon.ffld11 pxn -s pxn.pdb
# (creates pxn.lib and pxn.prm and pxn.prm.chk)
#
#
#

import math
import sys
import os
import time
try:
    import argparse
except ImportError:
    import lib.argparse as argparse
from lib.common import backup_file, __version__
from qscripts_config import QScriptsConfig as QScfg

parser = argparse.ArgumentParser()
parser.add_argument('ffld_output', help = 'ffld_server output')
parser.add_argument('code',  help = '3 letter code for the residue/ligand')
parser.add_argument("-s", dest="pdb", help="QM optimized PDB structure file WHICH WAS USED TO CREATE THE FFLD_OUTPUT (used to copy atom names and check parameter energies)", default=argparse.SUPPRESS)

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

args = parser.parse_args()

for k,v in vars(args).iteritems():
    if k in ['ffld_output', 'pdb'] and not os.path.lexists(v):
        print "FATAL! File %s doesn't exist." % v
        sys.exit(1)
if len(args.code) != 3:
    print "\n\nFATAL! The '3 letter code' should have 3 letters."
    sys.exit(1)

libfn = args.code + ".lib"
prmfn = args.code + ".prm"
prmchkfn = args.code + ".prm.chk"


# most of this is directly from Q/md.f90
def calc_bond(ac1, ac2, fk, r0):
    rji1 = ac2[0] - ac1[0] 
    rji2 = ac2[1] - ac1[1] 
    rji3 = ac2[2] - ac1[2] 
    r = math.sqrt( rji1**2 + rji2**2 + rji3**2 )
    dr = r-r0
    en = 0.5*fk*dr**2
    return en, r
   
def calc_angle(ac1, ac2, ac3, fk, th0):
    th0 = math.pi/180.0 * th0   # degrees to radians
    rji1 = ac1[0] - ac2[0]
    rji2 = ac1[1] - ac2[1]
    rji3 = ac1[2] - ac2[2]
    rjk1 = ac3[0] - ac2[0]
    rjk2 = ac3[1] - ac2[1]
    rjk3 = ac3[2] - ac2[2]
    
    # get the angle from the dot product equation ( A*B = |A|*|B|*cos(theta) )
    # where A and B are vectors bewteen atoms (1->2 and 2->3)
    bji2inv = 1./(rji1**2 + rji2**2 + rji3**2)
    bjk2inv = 1./(rjk1**2 + rjk2**2 + rjk3**2)
    bjiinv = math.sqrt(bji2inv)
    bjkinv = math.sqrt(bjk2inv)
    scp = (rji1*rjk1 + rji2*rjk2 + rji3*rjk3)
    scp = scp * bjiinv* bjkinv
    if scp > 1.0:
        scp =  1.0
    elif scp < -1.0:
        scp = -1.0
    theta = math.acos(scp)   # in radians
    dtheta = (theta-th0)  
    en = 0.5*fk*dtheta**2 
    return en,theta*180.0/math.pi    # kcal/mol and degrees
    
def calc_torsion(ac1, ac2, ac3, ac4, fks):
    rji1 = ac1[0] - ac2[0]
    rji2 = ac1[1] - ac2[1]
    rji3 = ac1[2] - ac2[2]
    rjk1 = ac3[0] - ac2[0]
    rjk2 = ac3[1] - ac2[1]
    rjk3 = ac3[2] - ac2[2]
    rkl1 = ac4[0] - ac3[0]
    rkl2 = ac4[1] - ac3[1]
    rkl3 = ac4[2] - ac3[2]
    
    rnj1 =  rji2*rjk3 - rji3*rjk2
    rnj2 =  rji3*rjk1 - rji1*rjk3
    rnj3 =  rji1*rjk2 - rji2*rjk1
    rnk1 = -rjk2*rkl3 + rjk3*rkl2
    rnk2 = -rjk3*rkl1 + rjk1*rkl3
    rnk3 = -rjk1*rkl2 + rjk2*rkl1

    bj2inv = 1./(rnj1**2 + rnj2**2 + rnj3**2)
    bk2inv = 1./(rnk1**2 + rnk2**2 + rnk3**2)
    bjinv = math.sqrt(bj2inv)
    bkinv = math.sqrt(bk2inv)
    
    scp = (rnj1*rnk1+rnj2*rnk2+rnj3*rnk3)*(bjinv*bkinv)
    if scp > 1.0:
        scp = 1.0
    elif scp < -1.0:
        scp = -1.0
    phi = math.acos(scp)

    if (rjk1*(rnj2*rnk3-rnj3*rnk2) + rjk2*(rnj3*rnk1-rnj1*rnk3) + rjk3*(rnj1*rnk2-rnj2*rnk1) < 0):
        phi = -phi
    
    energy = lambda x: fks[0]*(1.0+math.cos(x)) + fks[1]*(1.0-math.cos(2*x)) + fks[2]*(1.0+math.cos(3*x)) + fks[3]*(1.0+math.cos(4*x))
    # find the absolute minima of the function
    phis = range(0,185,5)
    energies = [ energy(p/180.0*math.pi) for p in phis ] 
    abs_min_i = energies.index( min(energies) )
    abs_min_phi = phis[abs_min_i]
    return energy(phi), phi*180.0/math.pi, abs_min_phi

def calc_improper(ac1, ac2, ac3, ac4, fk, phi0):
    phi0 = phi0*math.pi/180
    rji1 = ac1[0] - ac2[0]
    rji2 = ac1[1] - ac2[1]
    rji3 = ac1[2] - ac2[2]
    rjk1 = ac3[0] - ac2[0]
    rjk2 = ac3[1] - ac2[1]
    rjk3 = ac3[2] - ac2[2]
    rkl1 = ac4[0] - ac3[0]
    rkl2 = ac4[1] - ac3[1]
    rkl3 = ac4[2] - ac3[2]
    
    
# cross product to get normals to the planes
    rnj1 =  rji2*rjk3 - rji3*rjk2
    rnj2 =  rji3*rjk1 - rji1*rjk3
    rnj3 =  rji1*rjk2 - rji2*rjk1
    rnk1 = -rjk2*rkl3 + rjk3*rkl2
    rnk2 = -rjk3*rkl1 + rjk1*rkl3
    rnk3 = -rjk1*rkl2 + rjk2*rkl1
    
# doc product to get angle between normals (planes)
    bj2inv  = 1./(rnj1**2 + rnj2**2 + rnj3**2)
    bk2inv  = 1./(rnk1**2 + rnk2**2 + rnk3**2)
    bjinv = math.sqrt(bj2inv)
    bkinv = math.sqrt(bk2inv)
    
    scp = (rnj1*rnk1+rnj2*rnk2+rnj3*rnk3)*(bjinv*bkinv)
    if scp > 1.0:
        scp = 1.0
    elif scp < -1.0:
        scp = -1.0
    phi = math.acos(scp)

    if (rjk1*(rnj2*rnk3-rnj3*rnk2) + rjk2*(rnj3*rnk1-rnj1*rnk3) + rjk3*(rnj1*rnk2-rnj2*rnk1) < 0):
        phi = -phi
    
    dp = phi - phi0
    dp = dp - 2.0*math.pi*round(dp/(2.0*math.pi))
    energy = 0.5*fk*dp**2
    return energy, phi*180/math.pi


    



ATOM_MASSES = { "H": 1.0079,    
                "C": 12.011,
                "N": 14.007,
                "O": 15.999,
                "F": 18.988,
                "P": 30.974,
                "S": 32.065,
                "CL": 35.453,
                "BR": 79.904,
                "I": 126.90,
                "DU": 0 }
    
prefix = args.code.lower() + "_"
parms = { "atom_types": [], "bonds": [], "angles": [], "torsions": [], "impropers":[] }
lib = { "atoms": [], "bonds": [], "impropers": [], "charge_groups": [] }
prm_checks = []

# parse PDB file (if it was specified with -s)
pdb_atoms = []
pdb_coords = []
if hasattr(args,"pdb"):
    for line in open(args.pdb, 'r').readlines():
        if line.startswith("ATOM") or line.startswith("HETATM"):
            pdb_atoms.append( line[11:17].strip() )
            x,y,z = float(line[31:39]), float(line[39:47]), float(line[47:55])
            pdb_coords.append( (x,y,z) )
    if not pdb_atoms or not pdb_coords:
        print "\nERROR: Failed to parse the PDB file. Try running it through babel.\n"
        sys.exit(1)
else:
    print "\nWARNING: No PDB specified (-s pdb_filename), the atom names in the library file will probably not match the PDB ones.\n"

# parse the FFLD file
lookup_dict = {}   # keys are ffld atom names, values are tuples (pdb atname, (x,y,z))
section = ""
for line in open(args.ffld_output, 'r').readlines():
    l = line.strip()
    if (l == "") or ("------" in l):
        continue
    elif l == "atom   type  vdw  symbol    charge   sigma    epsilon  quality   comment":
        section = "VDW"
        continue
    elif l == "Stretch            k            r0    quality         bt        comment":
        section = "BONDS"
        continue
    elif l == "Bending                      k       theta0    quality   at  comment":
        section = "ANGLES"
        continue
    elif l == "proper Torsion                     V1      V2      V3      V4    quality  tt  comment":
        section = "TORSIONS"
        continue
    elif "improper" in l:
        section = "IMPROPERS"
        continue

    if section == "VDW":
#
#  C1      135  C1   CT      -0.0175   3.5000   0.0660 high   C: alkanes
#
        l = l.split()
        name, typ, vdw, symbol, charge, sigma, epsilon, quality, comment = l[0], l[1], l[2], l[3], float(l[4]), float(l[5]), float(l[6]), l[7], l[8:]
# get the element name from the vdw type: anything that is not a number basically (C from C1, Br from Br1)
        element = "".join( c for c in vdw if c.isalpha())
        if pdb_atoms:
            current_index = len(parms["atom_types"])
            lookup_dict[name] = ( pdb_atoms[current_index], pdb_coords[current_index] )
# get the element name from the atom name: first letter + any letter that is lowercase (C from CA, Br from Br11, H from HA2,...)
# this is not completely fool proof though
            pdb_element = "".join( c for i,c in enumerate(pdb_atoms[current_index]) if (i==0 or c.islower()) )
            if pdb_element.lower() != element.lower():
                print """ERROR: Element names '%s: %s' (ffld) and '%s: %s' (pdb) don't match.
Please check if the two-letter elements in your PDB are all upper case (should be Cl2, Br11, Ru1 instead of CL2, BR11, RU1). 
If not, the atom order might be screwed up!"""  % (name, element, pdb_atoms[current_index], pdb_element) 
                if raw_input("\nAre you sure this is ok? Are you really sure? (y/N) ") != "y":
                    print "Fine, quitting..."
                    sys.exit(1)
                else:
                    print "Ignoring..."
                    pass # ignore
            name = pdb_atoms[current_index]

        prm_name = prefix + name
        comment = " ".join(comment)
        lj_A = math.sqrt(4*epsilon*((sigma)**12))
        lj_B = math.sqrt(4*epsilon*((sigma)**6))
        lj_A_14 = lj_A/math.sqrt(2)
        lj_B_14 = lj_B/math.sqrt(2)
        element = "".join( c for c in vdw if c.isalpha() )
        try:
            mass = ATOM_MASSES[element.upper()]
        except KeyError:
            print "WARNING: Mass for element '%s' (atom '%s') not found, set it manually." % (element, name)
            mass = "<FIX>"

        s = "%-15s %10.4f %10.4f %10.4f %10.4f %10.4f %10s  ! %-6s %-5s ! %s" % \
             (prm_name, lj_A, lj_A, lj_B, lj_A_14, lj_B_14, mass, quality, symbol, comment)
        parms["atom_types"].append( s )

        ind = len(lib["atoms"]) + 1
        s = "%5d  %-5s %-15s %10.5f" % (ind, name, prm_name, charge)
        lib["atoms"].append( s )
        lib["charge_groups"].append( name )

    elif section == "BONDS":
#
#  C1      H2      340.00000    1.09000   high      140  0   CT  -HC    ==> CT  -HC
#
        l = l.split()
        at1,at2,fk,r0,quality,comment = l[0],l[1],float(l[2])*2.0,float(l[3]),l[4],l[7:]
        if pdb_atoms:
            at1,at1_coords = lookup_dict[at1]
            at2,at2_coords = lookup_dict[at2]
            en,r = calc_bond(at1_coords, at2_coords, fk, r0)
            prm_checks.append("Bond: %-5s %-5s  ###  ###    R0: %7.4f     R: %7.4f     Energy: %7.2f kcal/mol" % (at1, at2, r0, r, en))

        prm_at1 = prefix + at1
        prm_at2 = prefix + at2
        comment = " ".join(comment)
        
        s = "%-15s %-15s %10.2f %10.4f     ! %-6s ! %s" % (prm_at1, prm_at2, fk, r0, quality, comment)
        parms["bonds"].append( s )

        s = "%-5s %-5s" % (at1,at2)
        lib["bonds"].append( s )
    
    elif section == "ANGLES":
#
#  H2      C1      H3        33.00000  107.80000   high      410     0   HC  -CT  -HC    ==> HC  -CT  -HC
#
        l = l.split()
        at1,at2,at3,fk,th0,quality,comment = l[0],l[1],l[2],float(l[3])*2.0,float(l[4]),l[5],l[8:]
        if pdb_atoms:
            at1,at1_coords = lookup_dict[at1]
            at2,at2_coords = lookup_dict[at2]
            at3,at3_coords = lookup_dict[at3]
            en,th = calc_angle(at1_coords, at2_coords, at3_coords, fk, th0)
            prm_checks.append("Angle: %-5s %-5s %-5s ###  Theta0: %7.2f     Theta: %7.2f     Energy: %7.2f kcal/mol" % (at1, at2, at3, th0, th, en))

        at1 = prefix + at1
        at2 = prefix + at2
        at3 = prefix + at3
        comment = " ".join(comment)
        
        s = "%-15s %-15s %-15s %10.2f %10.4f     ! %-6s ! %s" % (at1,at2,at3,fk,th0,quality, comment)
        parms["angles"].append( s )

    elif section == "TORSIONS":
#
#  O2      P1      O3      C4        0.000   0.000   0.562   0.000    high    617   0   O2Z -P   -OS  -CT    ==> O2  -P   -OS  -CT 
#
        l = l.split()
        at1,at2,at3,at4,fks,quality,comment = l[0],l[1],l[2],l[3],l[4:8],l[8],l[11:]
        fks = [ float(fk)/2 for fk in fks ]
        if pdb_atoms:
            at1,at1_coords = lookup_dict[at1]
            at2,at2_coords = lookup_dict[at2]
            at3,at3_coords = lookup_dict[at3]
            at4,at4_coords = lookup_dict[at4]
            en,phi,minima = calc_torsion(at1_coords, at2_coords, at3_coords, at4_coords, fks)
            prm_checks.append("Torsion: %-5s %-5s %-5s %-5s   Abs.minima_phi: %7.2f     Phi: %7.2f     Energy: %7.2f kcal/mol" % (at1, at2, at3, at4, minima, phi, en))
        at1 = prefix + at1
        at2 = prefix + at2
        at3 = prefix + at3
        at4 = prefix + at4
        comment = " ".join(comment)
        
        non_zero_fks = []
        for i,fk in enumerate(fks):
            if abs(fk) > 0.000001:
                non_zero_fks.append( (i+1, fk) )
        if not non_zero_fks:
            non_zero_fks.append( (1, 0.0000) )

        path = 1.000   # no idea what this does

        for i,(rmult,fk) in enumerate(non_zero_fks):
            if rmult == 2:
                phase_shift = 180.0
            else:
                phase_shift = 0.0
            if i != len(non_zero_fks)-1: 
                rmult = -rmult
            s = "%-15s %-15s %-15s %-15s %10.4f %10.4f %10.4f %10.4f  ! %-6s ! %s" % (at1,at2,at3,at4,fk,rmult,phase_shift,path,quality,comment)
            parms["torsions"].append( s )

    elif section == "IMPROPERS":
#
#  C21     C22     C20     O19       2.200   high   aromatic atom           
#
        l = l.split()
# I have a feeling this fk should not be divided by two, but this is what everyone is using, so whatever...
        at1,at2,at3,at4,fk,quality,comment = l[0],l[1],l[2],l[3],float(l[4])/2,l[5],l[6:]
#
# It has been decided that no3 should be in second place, even though switching no2 and no3 doesn't really make a difference.
# To avoid Q doing stupid stuff later on, reorder the atoms in alphabetical order around no3 (which is on second place).
#
        a1,a3,a4 = sorted([at1,at2,at4])
        a2 = at3
        if pdb_atoms:
            a1,at1_coords = lookup_dict[a1]
            a2,at2_coords = lookup_dict[a2]
            a3,at3_coords = lookup_dict[a3]
            a4,at4_coords = lookup_dict[a4]
            en,phi = calc_improper(at1_coords, at2_coords, at3_coords, at4_coords, fk, 180.0)
            prm_checks.append("Improper: %-5s %-5s %-5s %-5s   Phi0: %7.2f     Phi: %7.2f     Energy: %7.2f kcal/mol" % (a1, a2, a3, a4, 180.0, phi, en))

        prm_at1 = prefix + a1
        prm_at2 = prefix + a2
        prm_at3 = prefix + a3
        prm_at4 = prefix + a4

        th0 = 180.0
        comment = " ".join(comment)
        
        s = "%-15s %-15s %-15s %-15s %10.4f %10.4f  ! %-6s ! %s" % (prm_at1,prm_at2,prm_at3,prm_at4,fk,th0,quality,comment)
        parms["impropers"].append( s )

        s = "%-5s %-5s %-5s %-5s" % (a1, a2, a3, a4)
        lib["impropers"].append( s )

    
# write out .prm file
print "Writing the parameter file: %s" % prmfn
backup = backup_file(prmfn)
if backup:
    print "Backed up '%s' to '%s'" % (prmfn, backup)

with open(prmfn, 'w') as prmf:
    prmf.write("# Generated with q_ffld2q.py (v%s) from '%s' (%s)\n\n" \
            % ( __version__, args.ffld_output, time.ctime() ) )
    for k in [ "atom_types", "bonds", "angles", "torsions", "impropers" ]:
# section
        prmf.write("\n[" + k + "]\n")
# parms
        for v in parms[k]:
            prmf.write(v + "\n")
        prmf.write("\n")
    
# join atoms to one line
lib["charge_groups"] = (" ".join(lib["charge_groups"]), )



# write out .lib file
print "Writing the library file: %s" % libfn
backup = backup_file(libfn)
if backup:
    print "Backed up '%s' to '%s'" % (libfn, backup)

with open(libfn, 'w') as libf:
    libf.write("{%s}\n# Generated with q_ffld2q.py (v%s) from '%s' (%s)\n\n" \
            % ( args.code.upper(), __version__, args.ffld_output, time.ctime() ) )
    for k in [ "atoms", "bonds", "impropers", "charge_groups" ]:
# section
        libf.write("[" + k + "]\n")
# parms
        for v in lib[k]:
            libf.write(v + "\n")
        libf.write("\n\n")



# write out .prmchk file
if pdb_atoms:
    print "Writing the parameter check file: %s    # Tip:  sort -n -k 11 %s" % (prmchkfn,prmchkfn)
    backup = backup_file(prmchkfn)
    if backup:
        print "Backed up '%s' to '%s'" % (prmchkfn, backup)
    open(prmchkfn, 'w').write( "\n".join(prm_checks) )
print "Done"
print "Please carefully inspect the parameters, even the 'high' quality ones, they can still suck..."
print
    

