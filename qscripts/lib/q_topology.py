#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#
# Copyright 2016    Miha Purg  <miha.purg@gmail.com>
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#
#

import re
from collections import OrderedDict as ODict

class QlibError(Exception):
    pass

class QprmError(Exception):
    pass

class Qlib(object):
    def __init__(self):
        self.residues = []   # list of _QResidue objects

    def _add_residue(self, residue):
        if residue.residue_name in [ x.residue_name for x in self.residues ]:
            raise QlibError("Duplicate library entry for residue '%s'" % residue.residue_name)
        self.residues.append(residue)


    def read_lib(self, libfile):
        with open(libfile, 'r') as lib:
            section = ""
            tmpresidue = None

            for lnumber, line in enumerate(lib.readlines()):
                lnumber += 1
                line = re.split("#|\*|\!",line, 1)[0].strip()   # remove comments 
                if line == "": continue
                if line[0] == "{":
                    if tmpresidue:
                        self._add_residue(tmpresidue)
                    tmpresidue = _QResidue(line.strip("{}"))   # get the 3 letter code and make a new object
                    continue
                if line[0] == "[":
                    section = line.strip("[]").lower()
                    continue
                if not tmpresidue or not section:
                    raise QlibError("Line #%s in LIB file '%s' is not a comment and is not inside a {XXX} section and [XXX...X] subsection:\n%s" % (lnumber, libfile, line) )
                
                if section == "atoms":
                    try:
                        atom_index,atom_name,atom_type,atom_charge = line.split()[0:4]
                        tmpresidue.atoms.append( (atom_name,atom_type,float(atom_charge)) )
                    except ValueError:
                        raise QlibError("Line #%s in LIB file '%s' couldn't be parsed (should look like this 'aindex aname atype charge ...'):\n%s" % (lnumber, libfile, line) )

                elif section == "bonds":
                    try:
                        atomnames = zip(*tmpresidue.atoms)[0]  # transpose and take first column
                        a1,a2 = line.split()
                        if a1 not in atomnames or a2 not in atomnames:
                            raise QlibError("Undefined atom(s) %s and/or %s mentioned in the bonds section of '%s', ln.%s." % (a1, a2, libfile, lnumber) )

                        tmpresidue.bonds.append( (a1,a2) )
                    except ValueError:
                        raise QlibError("Line #%s in LIB file '%s' couldn't be parsed (should look like this 'atom1 atom2'):\n%s" % (lnumber, libfile, line) )

                elif section == "impropers":
                    try:
                        atoms = line.split()[0:4]
                        atomnames = zip(*tmpresidue.atoms)[0]  # transpose and take first column
                        for atom_name in atoms:
                            if atom_name not in atomnames and atom_name not in ["+N", "-C"]: 
                                raise QlibError("Undefined atom %s mentioned in the improper section of '%s', ln.%s." % (atom_name, tmpresidue.residue_name, lnumber) )
                        tmpresidue.impropers.append( atoms )
                    except ValueError:
                        raise QlibError("Line #%s in LIB file '%s' couldn't be parsed (should look like this 'atom1 atom2 atom3 atom4'):\n%s" % (lnumber, libfile, line) )

                elif section == "connections":
                    tmpresidue.connections.append(line)

                elif section == "charge_groups":
                    cgrp = line.split()
                    atomnames = zip(*tmpresidue.atoms)[0]  # transpose and take first column
                    for atom_name in cgrp:
                        if atom_name not in atomnames:
                            raise QlibError("Undefined atom %s mentioned in the charge_groups section of '%s', ln.%s." % (atom_name, tmpresidue.residue_name, lnumber) )

                    tmpresidue.charge_groups.append( cgrp )

        self._add_residue(tmpresidue)

             

    def read_amber_lib(self, libfile):
        section = ""
        tmpresidue = None
        with open(libfile) as lib:
            for lnumber, line in enumerate(lib.readlines()):
                lnumber += 1
                line = line.strip()
                if not line: continue
                if line[:2] == "!!": continue # ignore comment
                if line[0] == "!":
                    l = line.split()[0].split(".")
                    if len(l) > 2:
                        rname = l[1]
                        if not tmpresidue:
                            tmpresidue = _QResidue(rname)
                        if tmpresidue.residue_name != rname:
                            self._add_residue(tmpresidue)
                            tmpresidue = _QResidue(rname)

                        section = l[3]
                    else: raise QlibError("Line #%s in LIB file '%s' couldn't be parsed:\n%s" % (lnumber, libfile, line) )
                    continue
        
                if not section or not tmpresidue:
                    continue
        
                if section == "atoms":
                    l = line.split()
                    try:
                        name,atype,charge=l[0].strip('"'),l[1].strip('"'),float(l[7])
                        atype = atype.replace("*", "star")
                        tmpresidue.atoms.append( (name,atype,charge) )
                    except ValueError:
                        raise QlibError("Line #%s in LIB file '%s' couldn't be parsed:\n%s" % (lnumber, libfile, line) )

                elif section == "connectivity":
                    atomnames = zip(*tmpresidue.atoms)[0]  # transpose and take first column
                    try:
                        ai1,ai2 = line.split()[:2]
                        a1 = atomnames[int(ai1)-1]
                        a2 = atomnames[int(ai2)-1]
                        if a1 not in atomnames or a2 not in atomnames:
                            raise QlibError("Undefined atom(s) %s and/or %s mentioned in the connectivity section of '%s', ln.%s." % (a1, a2, libfile, lnumber) )
                        tmpresidue.bonds.append( (a1,a2) )
                    except ValueError:
                        raise QlibError("Line #%s in LIB file '%s' couldn't be parsed:\n%s" % (lnumber, libfile, line) )

                elif section == "residueconnect":
                    atomnames = zip(*tmpresidue.atoms)[0]  # transpose and take first column
                    try:
                        ai1,ai2 = line.split()[:2]
                        a1 = atomnames[int(ai1)-1]
                        a2 = atomnames[int(ai2)-1]
                        if a1 not in atomnames or a2 not in atomnames:
                            raise QlibError("Undefined atom(s) %s and/or %s mentioned in the residueconnect section of '%s', ln.%s." % (a1, a2, libfile, lnumber) )
                        tmpresidue.connections.append( "head " + a1 )
                        tmpresidue.connections.append( "tail " + a2 )
                    except ValueError:
                        raise QlibError("Line #%s in LIB file '%s' couldn't be parsed:\n%s" % (lnumber, libfile, line) )

        self._add_residue(tmpresidue)

    def read_mol2(self, mol2file):
        pass

    def read_ffld(self, ffldfile):
        pass

    def get_string(self):
        out = ""
        for residue in self.residues:
            out += residue.get_Qlib_str()
            out += "*" + "-"*80 + "\n"
        return out



class AtomType(object):
    def __init__(self, atom_type, lj_A, lj_B, lj_r, lj_eps, mass):
        # atom_type is a string like "CA" or "OW"
        self.atom_type = atom_type
        self.atom_type_id = atom_type
        self.lj_A = float(lj_A)
        self.lj_B = float(lj_B)
        self.lj_R = float(lj_R)
        self.lj_eps = float(lj_eps)
        # Soft_core Ci (change manually in the fep) and ai 
        self.sc_Ci = 1.0
        self.sc_ai = 2.5
        self.mass = float(mass)

    def __repr__(self):
        return "AtomType (%s)" % (self.atom_type_id)

class BondType(object):
    def __init__(self, atom_types, fc, r):
        self.atom_types = atom_types   # list of atom_type strings
        self.fc = fc
        self.r = r

        # bond_id = atom_types - sorted, "CA CB" instead of "CB CA", to prevent double entries
        self.bond_id = " ".join(sorted(atom_types))

    def __repr__(self):
        return "BondType (%s)" % (self.bond_id)

class AngleType(object):
    def __init__(self, atom_types, fc, theta):
        self.atom_types = atom_types
        self.fc = fc
        self.theta = theta

        # angle_id = atom types - smaller of forward or reverse lists, "CA OS CT" not "CT OS CA"
        self.angle_id = " ".join( min(atom_types, list(reversed(atom_types))) )   # take the "smaller" of the normal and reverse lists to prevent double entries

    def __repr__(self):
        return "AngleType (%s)" % (self.angle_id)

class TorsionType(object):
    def __init__(self, atom_types, fc, nminim, phase, npaths):   # parameters == [ (fc,-1,psi,1), (fc,-2,psi,1), (fc, 3,psi,1) ]
        self.atom_types = atom_types
        self.fc = fc
        self.nminim = nminim    # -1,-2,3 or similar
        self.phase = phase
        self.npaths = npaths    # not used in oplsaa

# torsion_id = atom_types - smaller of forward or reverse lists, "CA CB CG OG1" instead of "OG1 CG CB CA", plus absolute of nminim (multiple functions for one torsion)
        self.torsion_id = " ".join( min(atom_types, list(reversed(atom_types))) ) + str(abs(int(self.nminim)))

    def __repr__(self):
        return "TorsionType (%s)" % (self.torsion_id)

class ImproperType(object):
    def __init__(self, atom_types, fc, psi):
        self.atom_types = atom_types  
        self.fc = fc
        self.psi = psi

        # improper_id = atom types - sorted, with second one fixed, "CB CG OG1 OG2" instead of "OG1 CG CB OG2" or "CB CG OG2 OG1" or ...
        center_atom = atom_types.pop(1)   # sort all but the center atom
        atom_types.sort()
        atom_types.insert(1, center_atom)
        self.improper_id = " ".join(atom_types)     # take the "smaller" of the normal and reverse lists to prevent double entries

    def __repr__(self):
        return "ImproperType (%s)" % " ".join(self.atom_types)



class Qprm(object):
    def __init__(self):
        self.options = ODict()
        self.atom_types = ODict()
        self.bonds = ODict()
        self.angles = ODict()
        self.torsions = ODict()
        self.impropers = ODict()
   
    def read_prm(self, prmfile):
        with open(prmfile, 'r') as parm:
            section=""
            for lnumber, line in enumerate(parm.readlines()):
                lnumber += 1
                line = re.split("#|\*|\!",line, 1)[0].strip()   # remove comments 
                if line == "":
                    continue
                if line[0] == "[":
                    section = line.split("]")[0].strip(" [").lower()         # it is apparently allowed to write text after the section identifier, so a simple strip isn't enough
                    continue
                if not section:
                    raise QprmError("Line %s in PARM file '%s' is not a comment and is not inside any section ([atom_types], [bonds], ...):\n%s" % (lnumber,parm_fn,line))

                if section == "atom_types":
                    parms = line.split()
                    try:
                        atom_type = parms[0]
                        lj_a, lj_b, lj_a_14, lj_b_14, mass = parms[1], parms[3], parms[4], parms[5], parms[6]
                    except IndexError:
                        print "\n\nFATAL! Could not parse line %s in [atom_types] section of parm file '%s':\n%s" % (lnumber, parm_fn, line)
                        sys.exit(1)
                    try:
                        AtomType.find_by_atomtype(atom_type)  # expecting an exception
                    except KeyError:
                        fep_types["atoms"].append( AtomType(atom_type, lj_a, lj_b, lj_a_14, lj_b_14, mass) )
                    else:
                        print "\n\nFATAL! Double VDW parameters found in parm file '%s', ln.%s: %s" % (parm_fn, lnumber, atom_type)
                        sys.exit(1)
                elif section == "bonds":
                    parms = line.split()
                    try:
                        atom_types = parms[0:2]
                        fc, r = parms[2], parms[3]
                    except IndexError:
                        print "\n\nFATAL! Could not parse line %s in [bonds] section of parm file '%s':\n%s" % (lnumber, parm_fn, line)
                        sys.exit(1)
                    try:
                        BondType.find_by_atomtypes(atom_types)  # expecting an exception
                    except KeyError:
                        fep_types["bonds"].append( BondType(atom_types, fc, r) )   # create new BondType
                    else:
                        print "\n\nFATAL! Double BOND parameters found in parm file '%s', ln.%s: %s" % (parm_fn, lnumber, " ".join(atom_types))
                        sys.exit(1)
                elif section == "angles":
                    parms = line.split()
                    try:
                        atom_types = parms[0:3]
                        fc, theta = parms[3], parms[4]
                    except IndexError:
                        print "\n\nFATAL! Could not parse line %s in [angles] section of parm file '%s':\n%s" % (lnumber, parm_fn, line)
                        sys.exit(1)
                    try:
                        AngleType.find_by_atomtypes(atom_types)  # expecting an exception
                    except KeyError:
                        fep_types["angles"].append( AngleType(atom_types, fc, theta) )  # create new AngleType
                    else:
                        print "\n\nFATAL! Double ANGLE parameters found in parm file '%s', ln.%s: %s" % (parm_fn, lnumber, " ".join(atom_types))
                        sys.exit(1)
                elif section == "torsions":
                    parms = line.split()
                    try:
                        atom_types = parms[0:4]
                        fc, rmult, psi = parms[4], parms[5], parms[6]
                    except IndexError:
                        print "\n\nFATAL! Could not parse line %s in [torsions] section of parm file '%s':\n%s" % (lnumber, parm_fn, line)
                        sys.exit(1)
        
                    if not abs(int(float(rmult))) in [1,2,3,4]:
                        print "\n\nFATAL! Found torsion with unsupported mutliplier (not +- 1,2,3 or 4) in parm file '%s', ln.%s:\n%s" % (parm_fn, lnumber, " ".join(atom_types))
                        sys.exit(1)
                    try:
                        tortypes = TorsionType.find_by_atomtypes(atom_types)  # this guy returns multiple tortypes in a dictionary or a KeyError if no torsion with this atom_types combo was defined yet
                    except KeyError:
                        fep_types["torsions"].append(TorsionType(atom_types, fc, rmult, psi))    # create new TorsionType
                    else:
                        if tortypes.has_key( abs(int(float(rmult))) ):
                            print "\n\nFATAL! Double TORSION parameters found in parm file '%s', ln.%s: %s" % (parm_fn, lnumber, " ".join(atom_types))
                            sys.exit(1)
                        else:
                            fep_types["torsions"].append(TorsionType(atom_types, fc, rmult, psi))    # create new TorsionType
                elif section == "impropers":
                    parms = line.split()
                    try:
                        atom_types = parms[0:4]
                        fc, psi = parms[4], parms[5]
                    except IndexError:
                        print "\n\nFATAL! Could not parse line %s in [torsions] section of parm file '%s':\n%s" % (lnumber, parm_fn, line)
                        sys.exit(1)
                    try:
                        ImproperType.find_by_atomtypes(atom_types)  # expecting an exception
                    except KeyError:
                        fep_types["impropers"].append( ImproperType(atom_types, fc, psi) )  # create new AngleType
                    else:
                        print "\n\nFATAL! Double IMPROPER parameters found in parm file '%s', ln.%s: %s" % (parm_fn, lnumber, " ".join(atom_types))
                        sys.exit(1)
                elif section == "options":
                    pass   # ignore this section
                else:
                    raise QprmError("Unknown section found in the parm file %s: %s" % (parm_fn, section))

    def read_amber_parm(self, parmfile):
        pass

    def read_amber_frcmod(self, frcmodfile):
        pass

    def read_ffld(self, ffldfile):
        pass

    def get_string(self):
        pass
    

class Qcoord(object):
    def __init__(self):
        pass

    def read_pdb(self, pdbfile):
        pass

    def read_mol2(self, mol2file):
        pass

    

class QTopology(object):
    def __init__(self):
        self.lib = Qlib()
        self.prm = Qprm()
        self.coord = Qcoord()





class _QResidue(object):
    def __init__(self, resname):
        self.residue_name = resname
        self.atoms = []   # [ ["CA","cx",0.567], ... ]
        self.bonds = []   # [ ("CA","CB"), ("CA", "HA1"), ... ]
        self.impropers = []   # [ ("C1","C2","C3","C4"), ... ]   # C2 is center atom
        self.connections = []   # [ "head N", "tail C" ]
        self.charge_groups = []   # [ ["C1","C4", "N1", "H1"], ["C3","H2","H3"], ... ]

    def check(self):
        """
        Checks for duplicates, missing atoms, integer charges.
        Raises an error.
        """
        return

    def get_Qlib_str(self, force=False):
        """
        Returns the Q-lib formatted string for the residue.
        """
        al,bl,il,cl = [],[],[],[]
        for i,atom in enumerate(self.atoms):
            al.append("{:>5d} {a[0]:<5s} {a[1]:<15s} {a[2]:>10.5f}".format(i+1, a=atom))
        for bond in self.bonds:
            bl.append("    {b[0]:<5s} {b[1]:<5s}".format(b=bond))
        for imp in self.impropers:
            il.append("    " + " ".join([ "{:<5s}".format(a) for a in imp ]))
        for chgr in self.charge_groups:
            cl.append("    " + " ".join([ "{:<5s}".format(a) for a in chgr ]))

        al = "\n".join(al)
        bl = "\n".join(al)
        il = "\n".join(al)
        cl = "\n".join(al)
        col = "\n".join([ "    " + conn for conn in self.connections ] )

        return """\
{{{}}}

[atoms]
{}

[bonds]
{}

[impropers]
{}

[connections]
{}

[charge_groups]
{}
""".format(self.residue_name,al,bl,il,col,cl)


