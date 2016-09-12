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
# This script can act as a standalone program with no arguments,
# or it can be imported to use the QAnalyseMaps object 
# (first argument is a list of directories which contain the qfep output fil; optional argument is the name of the output file)
# and its functions like get_dGa_mean(), get_dG0_mean(), get_summary(), get_failed(), get_analysed() and so on...
# 
# .get_analysed() returns a list of QAnalyseMap objects, which have functions to retrieve useful info like:
#
## - Reaction free energy profile normalized vs the reactants state - dG vs (E2-E1) ( get_dG_dE() )
## - Activation and reaction free energies ( get_dGa(), get_dG0() )
## - Free Energy Perturbation results - dG vs lambda ( get_dG_lambda() )
## - Average energies (Qbond,Qangle...) in each step vs frame ( get_E_lambda() )
## - dE(state1), dE(state2), LRA and REORG energies (Qbond,Qangle,...) ( get_dE1(), get_dE2(), get_lra(), get_reorg() )

import os
import subprocess
import tempfile
import re
from qscripts_config import QScriptsConfig as QScfg
from lib.common import DataContainer, np, backup_file, __version__


class QAnalyseMapError(Exception):
    pass

         



class QAnalyseMaps():
    def __init__(self, mapped_directories, qfep_out=QScfg.get("files","qfep_out"), _qanalysemaplist=None):

        # failed are those that fail to give dGa and dG0 values
        self._analysed = []     # _QAnalyseMap instances
        self._failed = []     # ( path, errorstring )    tuples
        self._exclusions = {}   # { "full_88_177": QAnalyseMaps_instance, ... }

        self._dgas = []
        self._dg0s = []

        if not isinstance(mapped_directories, list):
            mapped_directories = [mapped_directories,]
        self._mapped_directories = mapped_directories

# analyse the maps 
# if _qanalysemaplist is not passed in create new _QAnalyseMap objects
        if not _qanalysemaplist:
            qanlist = []
            for md in mapped_directories:
                    mfile=os.path.join(md, qfep_out)
                    qanlist.append( _QAnalyseMap(mfile) )
        else:
            qanlist = _qanalysemaplist

# parse and analyse the maps
        for qan in qanlist:
            try:
                qan._parse()
                dga,dg0 = qan.get_dGa(), qan.get_dG0()     
                self._analysed.append(qan)
                self._dgas.append(dga)
                self._dg0s.append(dg0)
            except Exception as e:     # catch everything 
                #raise   # enable this for debug
                self._failed.append( (qan.get_dirname(), e) )

# If _qanalysemaplist was not passed in (the object was called from outside), iterate through the analysed objects and check for exclusions
# If they exist, create new QAnalyseMaps objects from them, and add to self._exclusions dictionary
        if not _qanalysemaplist:
            exclusions = {}
            for qan in self._analysed:
                for name, excl in qan._get_exclusions().iteritems():
                    if name not in exclusions.keys():
                        exclusions[name] = []
                    exclusions[name].append( excl )
            for name, qans in exclusions.iteritems():
                 self._exclusions[name] = QAnalyseMaps(self._mapped_directories, _qanalysemaplist=qans)


    def get_exclusions(self):
        return self._exclusions

    def get_dGa_mean(self):
        return np.mean( self._dgas )

    def get_dG0_mean(self):
        return np.mean( self._dg0s )

    def get_analysed(self):
        return self._analysed    # populated in __init__

    def get_failed(self):
        return self._failed    # populated in __init__

    def get_dGa_stdev(self):
        return np.std(self._dgas, ddof=1)

    def get_dG0_stdev(self):
        return np.std(self._dg0s, ddof=1)

    def get_dGa_median(self):
        return np.median(self._dgas)
                                                             
    def get_dG0_median(self):        
        return np.median(self._dg0s)
         
    def get_extrema_lambdas(self):
# get the lambdas that are most frequent in the replicas
        min1_ls, min2_ls, max_ls = {}, {}, {}
        for an in self.get_analysed():
            a,b,c = an.get_extrema_lambdas()
            min1_ls[a] = min1_ls.get(a, 0) + 1
            max_ls[b] = max_ls.get(b, 0) + 1
            min2_ls[c] = min2_ls.get(c, 0) + 1
        min1_l = max(min1_ls, key = lambda x: min1_ls[x])
        max_l = max(max_ls, key = lambda x: max_ls[x])
        min2_l = max(min2_ls, key = lambda x: min2_ls[x])
        return (min1_l, max_l, min2_l)

    
    def get_average_GCs(self, lambda1, lambda2, first_res=1, last_res=None, cached=False, solvent=False, qmaskfile=None):

# returns averages and stds for the group contributions
        average_GCs = DataContainer( ["Residue id", "VdW(%5.4f-%5.4f)_mean" % (lambda2, lambda1), 
                                      "VdW(%5.4f-%5.4f)_stdev" % (lambda2, lambda1), 
                                      "El(%5.4f-%5.4f)_mean" % (lambda2, lambda1), 
                                      "El(%5.4f-%5.4f)_stdev" % (lambda2, lambda1) ] )

        average_GCs.comment = "Group contribution statistics at lambdas: %s, %s" % (lambda1, lambda2)
        gcs={}
        for an in self.get_analysed():
            gc = an.get_group_contributions(lambda1,lambda2,first_res=first_res,last_res=last_res, cached=cached, solvent=solvent, qmaskfile=qmaskfile)
            for rc in gc.get_rows():
                resid = rc[0]
                values = [ [val,] for val in rc[1:] ]
                if not gcs.has_key(resid):
                    gcs[resid] = values
                else:
                    for i,val in enumerate(gcs[resid]):
                        val.extend( values[i] )

# iterate through each residue 
        for resid in sorted(gcs.keys()):
            rc = gcs[resid]
            # get mean and stdev
            rc_stats = [ resid,
                         np.mean(rc[0]), np.std(rc[0]),    # vdw
                         np.mean(rc[1]), np.std(rc[1]) ]   # el 

            average_GCs.add_row(rc_stats) 
        average_GCs.comment += "\n# "
        return average_GCs


    def get_average_LRAs(self, lambda1, lambda2):
# returns averages and stds for LRAs

        an_maps = self.get_analysed()
        if len(an_maps) != 0:

            average_lras = DataContainer( [ "E_type", "(E2-E1)_10_mean", "(E2-E1)_10_std", "(E2-E1)_01_mean", "(E2-E1)_01_std", "LRA_mean", "LRA_std", "REORG_mean", "REORG_std" ] )
            average_lras.comment = " LRA statistics"

            allvals=[]   
            for an in an_maps:
                lra=an.get_lra(lambda1,lambda2)
                rows = lra.get_rows()
                for irow,row in enumerate(rows):
                    try:
                        allvals[irow].append( row )
                    except IndexError:
                        allvals.append( [row, ] )

# allvals now looks like this:
# [ 
#   [ 
#     ["EQtot", EQtot_de_st1_1, EQtot_de_st2_1, EQtot_lra_1, EQtot_reorg_1], 
#     ["EQtot", EQtot_de_st1_2, EQtot_de_st2_2, ...], ...
#   ],
#   [
#     ["EQbond", EQbond_de_st1_1, EQbond_de_st2_1, EQbond_lra_1, EQbond_reorg_1], 
#     ["EQbond", EQbond_de_st1_2, EQbond_de_st2_2, ...], ...
#   ]
# ]
# 
            for values in allvals:
                # transpose to get [ ["EQtot","EQtot"...], [ EQtot_de_st1_1, EQtot_de_st1_2,...], [EQtot_de_st2_1,EQtot_de_st2_2,...], ...]
                values = zip(*values)  
                # now they can be easily averaged and std-ed
                e_type = values[0][0]
                de_st1_mean = np.mean(values[1])
                de_st2_mean = np.mean(values[2])
                lra_mean    = np.mean(values[3])
                reo_mean    = np.mean(values[4])
                de_st1_std   = np.std(values[1])
                de_st2_std   = np.std(values[2])
                lra_std      = np.std(values[3])
                reo_std      = np.std(values[4])

                average_lras.add_row( [e_type, de_st1_mean, de_st1_std, de_st2_mean, de_st2_std, lra_mean, lra_std, reo_mean, reo_std] )

            return average_lras


    def get_summary(self):
        # This function calculates means, deviations, medians (of dGa and dG0),
        # it looks for dGa outliers (>3s) and removes them 
        # It returns a string containing all results in a neat format.

        allres = {}

        allres["n"] = len(self._dgas)
        allres["dga"] = (self.get_dGa_mean(),self.get_dGa_stdev(),self.get_dGa_median())
        allres["dg0"] = (self.get_dG0_mean(),self.get_dG0_stdev(),self.get_dG0_median())

        qams = []
        warns = []
        for qam in self._analysed:
            rp = os.path.relpath(qam.get_dirname()) + os.path.sep 
            min1,ts,min2 = qam.get_extrema_lambdas()
            qams.append( "%-40s %8.2f %8.2f      %8.4f %8.4f %8.4f" % (rp, qam.get_dGa(), qam.get_dG0(), min1, ts, min2) )
            w = qam.get_warnings()
            if w:
                warns.append( "%s:\n-- %s" % (rp, "\n-- ".join( w ) ) )
        if warns:
            allres["warnings"] = "WARNINGS: %d\n%s" %( len(warns), "\n".join(warns) )
        else:
            allres["warnings"] = "Warnings: None"

        allres["analysed"] = "\n".join( sorted( qams ) )

        if self._failed:
            errors = "\n".join( [ "%s -> %s " % (n, e) for n,e in self._failed ] )
            allres["fails"] = "FAILED: %d\n%s" % ( len(self._failed), errors )
        else:
            allres["fails"] = "Failed to analyse: None"

        return    """
---------------------------- QAnalyseMaps SUMMARY ----------------------------
Analysed with version: {version}

DIRNAMES                                     dG#      dG0          RS_l     TS_l     PS_l
{analysed}

Assuming normal distribution...

N={n}        Mean      St.dev     Median   
dG#   {dga[0]:10.2f} {dga[1]:10.2f} {dga[2]:10.2f}
dG0   {dg0[0]:10.2f} {dg0[1]:10.2f} {dg0[2]:10.2f} 

{fails}

{warnings}
------------------------------------------------------------------------------
""".format( version=__version__, **allres ) 


class _QAnalyseMap():
    def __init__(self, mappinglogfile, _logfilestring=None):

        self._mappinglogfile = mappinglogfile
        self._dirname = os.path.dirname(os.path.abspath(mappinglogfile))
        self._warnings = []
        self._exclusions = {}   # { "full_386" : _QAanalyseMap_instance, ... }

        # compiled REs
        self._PART0_RE = re.compile("(# Part 0.*?)# Part 1", re.DOTALL)
        self._PART1_RE = re.compile("(# Part 1.*?)# Part 2", re.DOTALL)
        self._PART2_RE = re.compile("(# Part 2.*?)# Part 3", re.DOTALL)
        self._PART3_RE = re.compile("(# Part 3.*?)# Part 1", re.DOTALL)  # '# Part 1' is added manually to the end of the logfile
        self._EXCL_RE  = re.compile("Calculation for system with (\w+) exclusion, residues (.*?\n)", re.DOTALL)
        # constants
        self._PART0_HEADER = "# file             state   pts   lambda    EQtot   EQbond  EQang   EQtor   EQimp    EQel   EQvdW  Eel_qq  EvdW_qq Eel_qp  EvdW_qp Eel_qw EvdW_qw Eqrstr"
        self._PART1_HEADER = "# lambda(1)      dGf sum(dGf)      dGr sum(dGr)     <dG>"
        self._PART2_HEADER = "# Lambda(1)  bin Energy gap      dGa     dGb     dGg    # pts    c1**2    c2**2"
        self._PART3_HEADER = "# bin  energy gap  <dGg> <dGg norm> pts  <c1**2> <c2**2> <r_xy>"

        self._PART0_COLUMNS = [ "Energy_file", "State", "Points", "Lambda", "EQtot", "EQbond", "EQang", "EQtor", "EQimp", "EQel", "EQvdW", "Eel_qq", "EvdW_qq", "Eel_qp", "EvdW_qp", "Eel_qw", "EvdW_qw", "Eqrstr"]
        self._PART1_COLUMNS = [ "Lambda", "dGf", "sum(dGf)", "dGr", "sum(dGr)", "dG" ]
        self._PART2_COLUMNS = [ "Lambda", "bin", "Energy_gap", "dGa", "dGb", "dGg", "Points", "c1^2", "c2^2" ] 
        self._PART3_COLUMNS = [ "Bin", "Energy_gap", "dGg", "dGg_norm", "Points",  "c1^2", "c2^2", "r_xy" ]

        # Part 0 variables
        self._part0_strings = []  # only one part0 in any output, but is a list for consistency
        self._part0_data_st1 = DataContainer( self._PART0_COLUMNS )
        self._part0_data_st2 = DataContainer( self._PART0_COLUMNS )

        # Part 1 variables 
        self._part1_strings = []   # if exclusions were enabled, more than one part 1 will exist
        self._part1_data = DataContainer( self._PART1_COLUMNS )
        
        # Part 2 variables 
        self._part2_strings = []   # if exclusions were enabled, more than one part 2 will exist
        self._part2_data = DataContainer( self._PART2_COLUMNS )

        # Part 3 variables
        self._part3_strings = []   # if exclusions were enabled, more than one part 3 will exist
        self._part3_data = DataContainer( self._PART3_COLUMNS )

        self._dga = None
        self._dg0 = None
        self._minima_bins = None
        self._maxima_bins = None
         
        # group contributions
        self._group_contrib = None
        self._group_contrib_solvent = None


        if not _logfilestring:   # reading the output of qfep
            try:
                self._mapping_output = open(mappinglogfile,'r').read()
            except IOError as e:
                raise QAnalyseMapError("Cant read the qfep output file '%s'" % mappinglogfile)
            # add "\n# Part 1" to the end of the string, so that RE works
        else:   
            self._mapping_output = _logfilestring

        self._mapping_output += "\n# Part 1"

   

    def _parse(self):
# this should be called after initializing the object

        # The criteria for an analysis to be considered ok (not failed) is to pass all the parsing and to produce dGa and dG0  (for the full system  - the first part1,part2,part3 that appear)
        # If something fails in the following commands, don't catch the exception
        # parse the shit out of the output file       
        self._parse_part0()
        self._parse_part1()
        self._parse_part2()
        self._parse_part3()
        self._get_part3_dgs()
        
        # if the logfile contains multiple entries (full + exclusions), fill the _exclusions dictionary
        if len(self._part1_strings) > 1:  # exclusions are present
            part0 = self._part0_strings[0]
            for i in range(1, len(self._part1_strings)):
                part1 = self._part1_strings[i]
                part2 = self._part2_strings[i]
                part3 = self._part3_strings[i]
                match = self._EXCL_RE.search(part1)    # something like "Calculation for system with full exclusion, residues    88  177 ", returns ("full", "  88  177\n")
                if match:
                    excl_name = "_".join(" ".join(match.groups()).split())  # convert to full_88_177
                    excl_logstring = part0 + part1 + part2 + part3
                    self._exclusions[excl_name] = _QAnalyseMap(self._mappinglogfile, _logfilestring=excl_logstring)


    def _parse_part0(self):
        # this part contains average energies from each frame
        # It only works on 2 state EVB at the moment

        self._part0_strings = self._PART0_RE.findall(self._mapping_output)
        if not self._part0_strings:
           raise QAnalyseMapError("Part 0 is missing in the mapping output.")
        part0 = self._part0_strings[0]   # take the first one (there should be only one part 0)

        lines = part0.split('\n')[1:]    # the first one is a comment
        header = lines.pop(0)    # comment with column names
        if header != self._PART0_HEADER:
            raise QAnalyseMapError("Failed to parse part0. Was expecting a comment that contains column names and looks like this:\n%s\nGot this instead:\n%s\nDid the qfep5 binary change?" % (self._PART0_HEAEDER, header) )
        
        e_lines_read = 0
        for line in lines:
            if "Could not read file header!" in line:  # fix for Q version df1658653187a6a72799a28aef6ce59906540275
                continue
            line = re.split("#|\!", line)[0].strip()
            line = re.sub("(\d)-(\d)", "\g<1> -\g<2>", line)   # fix 'sticky' values in old versions of Q
            if line == '':
                continue
            if "-->" in line:
                e_lines_read = 0
                continue
## As of version 5.10.1, group exclusions are included in q. 
# In part0 of the mapping output, additional lines are printed for energies of system with exlusions.
# To not complicate things, only the first two (original, no exlusions) are extracted here.
            if e_lines_read == 2:
                continue
## the lines have these values:
# file  state   pts   lambda    EQtot   EQbond  EQang   EQtor   EQimp    EQel   EQvdW  Eel_qq  EvdW_qq Eel_qp  EvdW_qp Eel_qw EvdW_qw Eqrstr
            cols = line.split()
            enfile = cols[0]
            state = int(cols[1])
            points = int(cols[2])
            if state > 2: 
                raise QAnalyseMapError("Detected EVB with more than 2 states. Part0 parsing will not work!")
            lamb = float(cols[3])
            if state == 2:
                lamb = 1-lamb
            energies = [ float(x) for x in cols[4:] ]
            l = [enfile, state, points, lamb, ]
            l.extend( energies )

            if state == 1:
                self._part0_data_st1.add_row( l )
            else:
                self._part0_data_st2.add_row( l )

            e_lines_read += 1
        

        if not self._part0_data_st1.get_rows() or not self._part0_data_st1.get_rows():
            raise QAnalyseMapError("Part 0 is empty in the mapping output (no rows apparently).")



    def _parse_part1(self):
        # get part 1 - FEP (dg vs lambda)
        self._part1_strings = self._PART1_RE.findall(self._mapping_output)
        if not self._part1_strings:
           raise QAnalyseMapError("Part 1 is missing in the mapping output.")
        part1 = self._part1_strings[0]  # take only the first one

        lines = part1.split('\n')[1:] # first comment is useless

        ## As of version 5.10.1, group exclusions are included in q. 
        # check for the two extra lines and remove them
        if "Calculation" in lines[1]:
            lines = lines[2:]

        header = lines.pop(0)    # comment with column names
        if header != self._PART1_HEADER:
            raise QAnalyseMapError("Failed to parse part1. Was expecting a comment that contains column names and looks like this:\n%s\nGot this instead:\n%s\nDid the qfep5 binary change?" % (self._PART1_HEADER, header) )
        
        for line in lines:
            line = re.split("#|\!", line)[0].strip()
            if line == '':
                continue
            l = [ float(x) for x in line.split() ]
            self._part1_data.add_row( l )

        if not self._part1_data.get_rows():
            raise QAnalyseMapError("Part 1 is empty in the mapping output.")




    def _parse_part2(self):

        self._part2_strings = self._PART2_RE.findall(self._mapping_output)
        if not self._part2_strings:
           raise QAnalyseMapError("Part 2 is missing in the mapping output.")

        part2 = self._part2_strings[0]  # take only the first one

        lines = part2.split('\n')[1:] # first comment is useless

        header = lines.pop(0)    # comment with column names
        if header != self._PART2_HEADER:
            raise QAnalyseMapError("Failed to parse part2. Was expecting a comment that contains column names and looks like this:\n%s\nGot this instead:\n%s\nDid the qfep5 binary change?" % (self._PART2_HEADER, header) )
        
        self._part2_data.delete_rows()

        for line in lines:
            line = re.split("#|\!", line)[0].strip()
            if line == '':
                continue
            l = [ float(x) for x in line.split() ]
            self._part2_data.add_row(l)  

        if not self._part2_data.get_rows():
            raise QAnalyseMapError("Part 2 is empty in the mapping output.")



    def _parse_part3(self):
        # get Part 3 - normalized reaction free energies vs the energy gap

        self._part3_strings = self._PART3_RE.findall(self._mapping_output)
        if not self._part3_strings:
           raise QAnalyseMapError("Part 3 is missing in the mapping output.")
        part3 = self._part3_strings[0]  # take only the first one

        lines = part3.split('\n')[1:] # first comment is useless

        coltitles = lines.pop(0)    # comment with column names
        if coltitles != self._PART3_HEADER:
            raise QAnalyseMapError("Failed to parse part3. Was expecting a comment that contains column names and looks like this:\n%s\nGot this instead:\n%s\nDid the qfep5 binary change?" % (self._PART3_HEADER, coltitles) )
        
        for line in lines:
            line = re.split("#|\!", line)[0].strip()
            if line == '':
                continue
            l = [ float(x) for x in line.split() ]
            self._part3_data.add_row(l)
            
        if not self._part3_data.get_rows():
            raise QAnalyseMapError("Part 3 is empty in the mapping output.")



    def _get_part3_dgs(self):
        # get minima and maxima without any smoothing
        # if there is more than one maxima and less or more than 2 minima, raise an exception
        # search for maxima only between 0.2*nbins and 0.8*nbins (bad sampling on the edges can raise an error)
        # also save the bins of the minima - important for LRA and GC

        cols = self._part3_data.get_columns()
        bins = cols[0]
        des  = cols[1]
        dgs  = cols[3]
        minima,maxima = [],[]
        nbins=len(bins)
        for i in range(1,nbins-1):     # from the second to the second last
            
            dg,dgnext,dgprev = dgs[i],dgs[i+1],dgs[i-1]
            if dgprev >= dg and dg < dgnext: 
                minima.append(i)
            elif dgprev <= dg and dg > dgnext and i > nbins*0.2 and i < nbins*0.8: 
                maxima.append(i)

        if len(minima) > 2 or len(maxima) > 1:
            # bad sampling, more minima and maxima than wanted
            # get the highest maxima from those found so far
            # get the absolute minima to the left and to the right of this maxima 
            # set the warning string (rough profile)
            max1 = max(maxima, key=lambda i: dgs[i])
            react = [ (dgs[i],i) for i in minima if i < max1 ] 
            prod = [ (dgs[i],i) for i in minima if i > max1 ] 
            try:
                min1 = min(react)[1]   # min() will return tuple with lowest dg
                min2 = min(prod)[1]
            except ValueError:
# multiple minima on one side, none on the other (starts/ends at the lowest point)
                raise QAnalyseMapError("Bad reaction free energy profile - reactants minima: %d, products minima: %d" % (len(react), len(prod)))

            self._warnings.append("Rough Free energy profile (%d minima and %d maxima found), look at the graphs!" % (len(minima),len(maxima)))
            maxima = [ max1, ]
            minima = [ min1, min2 ]

        if len(minima) != 2:
            raise QAnalyseMapError("Bad reaction free energy profile - %d local minima (instead of 2)" % len(minima) )
        elif len(maxima) != 1:
            raise QAnalyseMapError("Bad reaction free energy profile - %d local maxima (instead of 1)" % len(maxima) )
        self._dga = dgs[maxima[0]] - dgs[minima[0]]
        self._dg0 = dgs[minima[1]] - dgs[minima[0]]
        self._minima_bins = [ bins[mini] for mini in minima ]
        self._maxima_bins = [ bins[maxi] for maxi in maxima ]
        
        # adjust the values in part3_data so that the reactants are zero 
        for row in self._part3_data.get_rows():
            row[3] = row[3] - dgs[minima[0]]


    def _get_exclusions(self):
        return self._exclusions

    def get_warnings(self):
        if self._warnings:
            return self._warnings
        else:
            return None
        

    def get_dirname(self):
        return self._dirname


    def get_E1_lambda(self):
        e1_lambda = DataContainer( [ "Lambda", "EQtot", "EQbond", "EQang", "EQtor", "EQimp", "EQel", "EQvdW", "Eel_qq", "EvdW_qq", "Eel_qp", "EvdW_qp", "Eel_qw", "EvdW_qw", "Eqrstr" ] )
        for row in self._part0_data_st1.get_rows( columns = e1_lambda.get_column_titles() ):
            e1_lambda.add_row(row)
        return e1_lambda


    def get_E2_lambda(self):
        e2_lambda = DataContainer( [ "Lambda", "EQtot", "EQbond", "EQang", "EQtor", "EQimp", "EQel", "EQvdW", "Eel_qq", "EvdW_qq", "Eel_qp", "EvdW_qp", "Eel_qw", "EvdW_qw", "Eqrstr" ] )
        for row in self._part0_data_st2.get_rows( columns = e2_lambda.get_column_titles() ):
            e2_lambda.add_row(row)
        return e2_lambda


 #   def get_lra(self):
 #       # This only works on 2 state EVB
 #       lra = DataContainer( [ "E_type", "(E2-E1)_10", "(E2-E1)_01", "LRA", "REORG" ] )
 #        
 #       # (E2-E1)_10    (reactant state) = First row E2 - E1
 #       # (E2-E1)_01    (products state) = Last row E2 - E1
 #       # [4:] ignores the lambda value at the 0 position
 #       des_st1 = [ e2-e1 for e1,e2 in zip( self._part0_data_st1.get_rows()[0][4:],  self._part0_data_st2.get_rows()[0][4:] ) ] 
 #       des_st2 = [ e2-e1 for e1,e2 in zip( self._part0_data_st1.get_rows()[-1][4:], self._part0_data_st2.get_rows()[-1][4:]) ] 

 #       # LRA=0.5*(<E2-E1>_10+<E2-E1>_01)
 #       # REO=0.5*(<E2-E1>_10-<E2-E1>_01)
 #       es_lra = [ 0.5 * (de_st1+de_st2) for de_st1,de_st2 in zip( des_st1, des_st2 ) ] 
 #       es_reo = [ 0.5 * (de_st1-de_st2) for de_st1,de_st2 in zip( des_st1, des_st2 ) ] 
 #       
 #       e_types = self._part0_data_st1.get_column_titles()[4:]

 #       for row in zip( e_types, des_st1, des_st2, es_lra, es_reo ):
 #           lra.add_row(row)

 #       return lra



    def get_lra(self, lambda1, lambda2):
        # same as get_lra except that it calculates the values not from the first and the last energies
        # but from energies at lambda values which have the highest sampling at the bins of the reactants and products minima
        # This only works on 2 state EVB

        lra = DataContainer( [ "E_type", "(E2-E1)_10", "(E2-E1)_01", "LRA", "REORG" ] )
        
        # get the appropriate rows of energies 
        # == didn't work with some floats so I just check if the diff is smaller than 0.000001
        # [4:] ignores 'file', 'state', 'points' and 'lambda'
        # [0] in the end just takes the first (which should be the only) value in the list comprehension
        e1_st1 = [ row[4:] for row in self._part0_data_st1.get_rows() if abs(row[3] - (lambda1)) < 0.000001 ][0]
        e1_st2 = [ row[4:] for row in self._part0_data_st1.get_rows() if abs(row[3] - (lambda2)) < 0.000001 ][0]
        e2_st1 = [ row[4:] for row in self._part0_data_st2.get_rows() if abs(row[3] - (lambda1)) < 0.000001 ][0] 
        e2_st2 = [ row[4:] for row in self._part0_data_st2.get_rows() if abs(row[3] - (lambda2)) < 0.000001 ][0]

        # (E2-E1)_10    (reactant state) = First row E2 - E1
        # (E2-E1)_01    (products state) = Last row E2 - E1
        des_st1 = [ e2-e1 for e1,e2 in zip(e1_st1, e2_st1) ] 
        des_st2 = [ e2-e1 for e1,e2 in zip(e1_st2, e2_st2) ] 

        # LRA=0.5*(<E2-E1>_10+<E2-E1>_01)
        # REO=0.5*(<E2-E1>_10-<E2-E1>_01)
        es_lra = [ 0.5 * (de_st1+de_st2) for de_st1,de_st2 in zip( des_st1, des_st2 ) ] 
        es_reo = [ 0.5 * (de_st1-de_st2) for de_st1,de_st2 in zip( des_st1, des_st2 ) ] 

        lra.comment = "Lambda values: %s, %s" % (lambda1, lambda2)

        e_types = self._part0_data_st1.get_column_titles()[4:]

        for row in zip( e_types, des_st1, des_st2, es_lra, es_reo ):
            lra.add_row(row)

        return lra

    
    def get_extrema_lambdas(self):
        # takes the bin numbers which correspond to extrema from part3 and looks up the lambdas where these are most populated (part2)
        # it also looks at the files in part0 which correspond to this lambda

        # get (points,lambda) values for the appropriate bins and take only the highest
        # points==6, lambda==0, bin==1
        minima1,minima2 = self._minima_bins
        maxima = self._maxima_bins[0]
        points,minima1_lambda = max( [ (int(row[6]),row[0]) for row in self._part2_data.get_rows() if int(row[1]) == minima1 ] )
        points,minima2_lambda = max( [ (int(row[6]),row[0]) for row in self._part2_data.get_rows() if int(row[1]) == minima2 ] )
        points,maxima_lambda = max( [ (int(row[6]),row[0]) for row in self._part2_data.get_rows() if int(row[1]) == maxima ] )

        self._minima1_lambda = float(minima1_lambda)
        self._minima2_lambda = float(minima2_lambda)
        self._maxima_lambda =  float(maxima_lambda)
        return self._minima1_lambda, self._maxima_lambda, self._minima2_lambda
        

    
    def get_dG_lambda(self):
        dg_lambda = DataContainer( ["Lambda","dG","dGr","dGf"] )
        for row in self._part1_data.get_rows( columns = dg_lambda.get_column_titles() ):
            dg_lambda.add_row(row)
        return dg_lambda

    def get_Egap_lambda(self):
        de_l = DataContainer( ["Lambda", "Energy_gap"] )
        for row in self._part2_data.get_rows( columns = de_l.get_column_titles() ):
            de_l.add_row(row)
        return de_l

    def get_points_Egap(self):
        pts_egap = DataContainer( ["Energy_gap", "Points"] )
        for row in self._part3_data.get_rows( columns = pts_egap.get_column_titles() ):
            pts_egap.add_row(row)
        return pts_egap

    def get_rxy_dE(self):
        rxy_de = DataContainer( ["Energy_gap", "r_xy"] )
        for row in self._part3_data.get_rows( columns = rxy_de.get_column_titles() ):
            rxy_de.add_row(row)
        return rxy_de

    def get_dG_dE(self):
        dg_de = DataContainer( ["Energy_gap", "dGg_norm"] )
        for row in self._part3_data.get_rows( columns = dg_de.get_column_titles() ):
            dg_de.add_row(row)
        return dg_de

    def get_dGa(self):
        return self._dga

    def get_dG0(self):
        return self._dg0


    def get_group_contributions(self, lambda1, lambda2, first_res=1, last_res=None, cached=False, solvent=False, qmaskfile=None):
# lambda1 is usually reactants state, lambda2 transition state
# the energies returned will be LRA energies ( LRA=0.5*(<E2-E1>_conf1+<E2-E1>_conf2), 
# where E2 is the energy function at lambda2 and E1 at lambda1 at configurations for lambda1 and lambda2 (conf1, conf2)

# if solvent==True, ignore first_res and last_res business and calculate for water only

# qmaskfile is a file with atom indexes for which the group contribution should be calculated (space or line separated)

# check if requesting cached results from a previous calculation (lambdas have to obviously be the same)
        if cached:
            if solvent and self._group_contrib_solvent and qmaskfile == self._gc_args_solv[2] and \
                    abs(lambda1)-abs(self._gc_args_solv[0]) < 0.00001 and abs(lambda2)-abs(self._gc_args_solv[1]) < 0.00001:
                return self._group_contrib_solvent
            elif not solvent and self._group_contrib and qmaskfile == self._gc_args[4] \
                    and abs(lambda1)-abs(self._gc_args[0]) < 0.00001 and abs(lambda2)-abs(self._gc_args[1]) < 0.00001 \
                    and first_res == self._gc_args[2] and last_res == self._gc_args[3]:
                return self._group_contrib

        if solvent:
            self._gc_args_solv = (lambda1,lambda2,qmaskfile)
        else:
            self._gc_args = (lambda1,lambda2,first_res,last_res,qmaskfile)

        if (subprocess.call("type " + QScfg.get("qexec", "qcalc"), shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE) != 0):
            raise QAnalyseMapError("\nGroup contributions failed - executable '%s' not found. Check your configuration..." % QScfg.get("qexec", "qcalc"))
        
# get energy filenames from lambda values 
        lambdas=(lambda1,lambda2)
        files_en = []
        for lam in lambdas:
            files_en.append( [ enf for (enf,lamb) in self._part0_data_st1.get_rows(columns=["Energy_file", "Lambda"]) if abs(lamb-lam) < 0.000001 ][0] )
# look at all the inputs (files that end with .inp) if they contain these energy filenames
# extract topology, fep and dcd filenames
        files_dcd = ["",""]
        files_top = ["",""]   # all three should be the same of course
        files_fep = ["",""]   # all three should be the same of course
        for fn in os.listdir( self.get_dirname() ):
            if fn[-4:] == ".inp":
                fs = open(os.path.join(self.get_dirname(), fn), 'r').read() 
                for i,ex_en in enumerate(files_en):
                    if fs.find(ex_en) != -1:
                        section = ""
                        for line in fs.split("\n"):
                            l = re.split("#|\!", line)[0].strip()
                            if l == "":
                                continue
                            elif l[0] == "[":
                                section = l
                            elif section.upper() == "[FILES]":
                                l = l.split()
                                ftype = l[0].upper()
                                fname = l[1]
                                if ftype == "TRAJECTORY":
                                    files_dcd[i] = fname
                                elif ftype == "TOPOLOGY":
                                    files_top[i] = fname
                                elif ftype == "FEP":
                                    files_fep[i] = fname

# check if fep and topology filenames are the same in all three inputs (this script is trying its best to be idiot proof)
        if not all(files_fep[0] == e_fep for e_fep in files_fep) or not all(files_top[0] == e_top for e_top in files_top):
            raise QAnalyseMapError("\nGroup contributions failed - your inputs use different fep or top files (hint: start over)")
        fepfile, topfile = files_fep[0], files_top[0]

# check if the files actually exist
        for fn in [fepfile,topfile,files_dcd[0],files_dcd[1]]:
            if not os.path.exists(os.path.join(self.get_dirname(), fn)):
                raise QAnalyseMapError("\nGroup contributions failed - file '%s' not found." % fn)

       # PARSE TOP FOR NUMBER OF SOLUTE RESIDUES
        with open(os.path.join(self.get_dirname(), topfile), "r") as top:
            for line in top.readlines():
                try:
                    if "No. of residues" in line:
                        l = line.split()
                        if "No of solute residues" in line:
                            num_res_all, num_res_solute = int(l[0]), int(l[1])
                        elif solvent:   # if solvent gc calc is true and "No of solute resiudes" is not in the topology line, then break
                            raise QAnalyseMapError("Group contributions for solvent failed - there is no solvent.")
                        else:
                            num_res_solute = int(l[0])
                    elif "no. of atoms, no. of solute atoms" in line:
                        l = line.split()
                        num_atoms_all, num_atoms_solute = int(l[0]), int(l[1])
                except ValueError:
                    raise QAnalyseMapError("Group contributions failed - parsing topology file '%s' failed." % topfile)

# check if residue indexes are ok
        if solvent:
            last_res = num_res_all
            first_res = num_res_solute + 1
            if first_res > last_res:
                raise QAnalyseMapError("\nGroup contributions on solvent failed, looks like a gas phase run...")
        else:
            if last_res:
                if last_res > num_res_solute:
                    raise QAnalyseMapError("\nGroup contributions failed - only '%d' solute residues in topology." % num_residues)
            else:
                last_res = num_res_solute

            if first_res < 1:
                raise QAnalyseMapError("\nGroup contributions failed - no residues with index less than 1...")
            elif first_res > last_res:
                raise QAnalyseMapError("\nGroup contributions failed - first residue ID greater than last...")

                
        # if qmaskfile is not given:
        # PARSE FEP FOR Q ATOM NUMBERS
        if not qmaskfile:
            with open(os.path.join(self.get_dirname(), fepfile), "r") as fep:
                section=""
                q_atoms=[]
                for line in fep.readlines():
                    l = re.split("#|\!", line)[0].strip()
                    if l == "":
                        continue
                    elif l[0] == "[":
                        section = l
                    elif section == "[atoms]":
                        q_atoms.append( l.split()[1] )
                        
            q_atoms = "\n".join( [ "%s %s" % (ai,ai) for ai in q_atoms ] )
        # else check whether all are integers and within the 1-num_atoms_solute
        else:
            qm = open(qmaskfile, "r").read().split()
            if qm:
                for i,qi in enumerate(qm):
                    try:
                        j = int(qi)
                        if j > num_atoms_solute:
                            raise QAnalyseMapError("\nGroup contributions failed - in qmask file: atom index '%d' larger than no. of solute atoms..." % j)
                    except ValueError:
                        raise QAnalyseMapError("\nGroup contributions failed - in qmask file: atom index '%s' not integer..." % qi)
                    qm[i] = str(j)
            else:
                raise QAnalyseMapError("\nGroup contributions failed - failed to parse the qmask file...")
            q_atoms = [ "%s %s" % (ai,ai) for ai in qm ]
            q_atoms = "\n".join(q_atoms)


        # use qcalc5 to calculate interaction energies in given frames and lambda values
        gcs=[]
        combs = ( (files_dcd[0], lambdas[0]),   # E1_conf1
                  (files_dcd[0], lambdas[1]),   # E2_conf1
                  (files_dcd[1], lambdas[0]),   # E1_conf2
                  (files_dcd[1], lambdas[1]) )  # E2_conf2
        # example with lambdas 1.00 and 0.50 
        # fep_000_1.000.dcd, 1.00 
        # fep_000_1.000.dcd, 0.50
        # fep_025_0.500.dcd, 1.00 
        # fep_025_0.500.dcd, 0.50

        for dcdfile,lamb in combs:
            a = """{top}
{fep}
{l1:.4f} {l2:.4f}
8
{first_res}
{last_res}
{q_atoms}
end
go
{dcd}
.
""".format(top=topfile,fep=fepfile,l1=lamb,l2=1-lamb,first_res=first_res,last_res=last_res,q_atoms=q_atoms,dcd=dcdfile)

            # call qcalc5
            try:
                p = subprocess.Popen(QScfg.get("qexec", "qcalc"), stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=self.get_dirname())
            except OSError as e:
                raise QAnalyseMapError("Problem when running qcalc5: %s" % e )
            stdout,stderr = p.communicate(a)

            if stderr:
                raise QAnalyseMapError("\nGroup contributions failed - qcalc5 produced STDERR: %s" % stderr)

            # extract the results
            pos = stdout.find("        Residue     Average LJ     Average el")
            if pos == -1:
                raise QAnalyseMapError("\nGroup contributions failed - qcalc5 output could not be parsed!")
            gcl = stdout[pos:].split("\n")[1:]  # remove the first line
            gc = {}
            for line in gcl:
                l = line.split()
                if len(l) == 3:
                    resid = int(l[0])
                    vdw = float(l[1])
                    el = float(l[2])
                    gc[resid] = {"vdw": vdw, "el": el}
            gcs.append(gc)
                
        resids = sorted(gcs[0].keys())

        # do the LRA thingy    
        # LRA=0.5*(<E2-E1>_conf1+<E2-E1>_conf2)
        e2e1_st1_vdw = [ gcs[1][key]["vdw"] - gcs[0][key]["vdw"] for key in resids ]
        e2e1_st1_el  = [ gcs[1][key]["el"]  - gcs[0][key]["el"]  for key in resids ]
        e2e1_st2_vdw = [ gcs[3][key]["vdw"] - gcs[2][key]["vdw"] for key in resids ] 
        e2e1_st2_el  = [ gcs[3][key]["el"]  - gcs[2][key]["el"]  for key in resids ]

        vdw_lra = [ 0.5*(a+b) for a,b in zip(e2e1_st1_vdw, e2e1_st2_vdw) ]
        el_lra  = [ 0.5*(a+b) for a,b in zip(e2e1_st1_el,  e2e1_st2_el)  ]

        # save and return the results in the form of DataContainer objects
        if solvent:
            self._group_contrib_solvent = DataContainer( ["Solvent_id", "VdW(l=%5.4f->l=%5.4f)" % (lambda1,lambda2), "El(l=%5.4f->l=%5.4f)" % (lambda1,lambda2) ] )
            for row in zip(resids, vdw_lra, el_lra):
                self._group_contrib_solvent.add_row(row)
            return self._group_contrib_solvent
        else:
            self._group_contrib = DataContainer( ["Residue_id", "VdW(l=%5.4f->l=%5.4f)" % (lambda1,lambda2), "El(l=%5.4f->l=%5.4f)" % (lambda1,lambda2) ] )
            for row in zip(resids, vdw_lra, el_lra):
                self._group_contrib.add_row(row)
            return self._group_contrib



# If run as a standalone program (not with q_mapfep.py for instance)

if __name__ == "__main__":
    from lib import plotdata
    import sys
    try:
        from collections import OrderedDict as ODict
    except ImportError:
        import lib.OrderedDict as ODict
    try:
        import argparse
    except ImportError:
        import lib.argparse as argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("--gc", dest="gc", nargs="*", help="Enable LRA group contribution calculations. Optional arguments are first and last residue id (--gc 200 359). Use --gc_l to define l1 and l2.")
    parser.add_argument("--gc_solv", dest="gc_solv", default=False, action="store_true", help="Enable LRA group contribution calculations for the solvent. No arguments.")
    parser.add_argument("--gc_l", dest="gc_l", nargs=2, metavar=("l1", "l2"), default=(1.0, 0.5), help="Specify lambdas (state 1) at which LRA group contributions are calculated. Default is '1.0 0.5'. The script also understands 'RS', 'TS', 'PS' (most populated lambdas of energy minima and maxima). Example: '--gc_l 1.0 TS' will calculate LRA energies between states l=1.0 and l=0.52, if TS is most populated at lambda 0.52.")
    parser.add_argument("--lra_l", dest="lra_l", nargs=2, metavar=("l1", "l2"), default=(1.0, 0.0), help="Specify lambdas (state 1) at which LRA and REORG are calculated. Default is '1.0 0.0'. The script also understands 'RS', 'TS', 'PS' (most populated lambdas of energy minima and maxima). Example: '--lra_l RS PS' will calculate LRA energies between states l=0.9 and l=0.12, if RS is most populated at lambda 0.9 and PS at 0.12.")
    parser.add_argument("--gc_mask", dest="gc_mask", help="Define the atom mask for group contributions from an external file instead of directly from the FEP file.")
    parser.add_argument("--mapdirs", dest="mapdirs", nargs="+", help="Directories mapped with qfeps (default is all subdirs or current working dir)", default=[])
    parser.add_argument("--qfep_out", dest="qfep_out", help="qfep output filename (default='%s')" % QScfg.get("files", "qfep_out"), default=QScfg.get("files", "qfep_out"))
    parser.add_argument("--plots_out", dest="plots_out", help="output filename for plot data (default='%s')" % QScfg.get("files", "analysemaps_plots"), default=QScfg.get("files", "analysemaps_plots"))

    args = parser.parse_args()

    if args.gc != None:
        if len(args.gc) == 2:
            try:
                args.gc = (int(args.gc[0]), int(args.gc[1]))
            except ValueError:
                print "\nFATAL! GC arguments are not integers."
                sys.exit(1)
        elif len(args.gc) == 0:
            args.gc = (1, None)
        else:
            print "\nFATAL! GC works with zero or two (first and last residue ID) arguments."
            sys.exit(1)

    gc_l = []
    for l in args.gc_l:
        try:
            l = float(l)
            if l < 0 or l > 1: raise ValueError
        except ValueError:
            l = l.upper()
            if l not in ('RS','TS','PS'):
                print "\nFATAL! GC lambdas make no sense. Floats (0<l<1) or 'RS','TS','PS' please."
                sys.exit(1)
        gc_l.append(l)

    lra_l = []
    for l in args.lra_l:
        try:
            l = float(l)
            if l < 0 or l > 1: raise ValueError
        except ValueError:
            l = l.upper()
            if l not in ('RS','TS','PS'):
                print "\nFATAL! LRA lambdas make no sense. Floats (0<l<1) or 'RS','TS','PS' please."
                sys.exit(1)
        lra_l.append(l)


    mapdirs = args.mapdirs
    if len(args.mapdirs) == 0:
        # check for a logfile in CWD
        mapdirs = []
        ls = os.listdir( os.getcwd() )
        dirs = [ f for f in ls if os.path.isdir(f) ]
        if args.qfep_out in ls:
            print "Found '%s' in current working directory." % args.qfep_out
            mapdirs.append( os.getcwd() )
        elif len(dirs) != 0:
            for d in dirs:
                if args.qfep_out in os.listdir( d ):
                    print "Found '%s' in subdir '%s'." % (args.qfep_out,d)
                    mapdirs.append( d )
        if len(mapdirs) == 0:
            print "Couldnt find '%s' in current working directory or subdirs. Try --help." % args.qfep_out
            sys.exit(1)


    amaps = QAnalyseMaps( mapdirs, qfep_out=args.qfep_out )

    ams = amaps.get_analysed()
    print "Succesfully analysed %s" % len(ams)
    fms = amaps.get_failed()
    print "Failed to analyse %s\n" % len(fms)


# extract data, unless all failed
    if ams:
    # get the most common minima an maxima lambda values in all repeats 
    # we need these values to calculate group contributions for all repeats with the same frames
        min1_lambda, max_lambda, min2_lambda = amaps.get_extrema_lambdas()
        a = {'RS': min1_lambda, 'TS': max_lambda, 'PS': min2_lambda }
        gc_l = map(lambda x: a.get(x, x), gc_l)
        lra_l = map(lambda x: a.get(x, x), lra_l)

        # make PlotData objects and save the data 
        plots = ODict()
        plots["dgde"] = plotdata.PlotData("Free energy profile", xlabel="E1-E2  [kcal/mol]", ylabel="Free energy  [kcal/mol]")
        if args.gc:
            gc_key = "gc_el_%s%s" % (gc_l[0], gc_l[1])
            plots[gc_key] = plotdata.PlotData("LRA Group contributions (electrostatic): ddG( l=%s -> l=%s )" % (gc_l[0], gc_l[1]), xlabel="Residue index", ylabel="ddG [kcal/mol]", plot_type="bar")
            gc_vdw_key = "gc_vdw_%s%s" % (gc_l[0], gc_l[1])
            plots[gc_vdw_key] = plotdata.PlotData("LRA Group contributions (vdw): ddG( l=%s -> l=%s )" % (gc_l[0], gc_l[1]), xlabel="Residue index", ylabel="ddG [kcal/mol]", plot_type="bar")
        if args.gc_solv:
            gc_solv_key = "gc_solv_%s%s" % (gc_l[0], gc_l[1])
            plots[gc_solv_key] = plotdata.PlotData("LRA GC_el - solvent: ddG( l=%s -> l=%s )" % (gc_l[0], gc_l[1]), xlabel="Solvent molecule index", ylabel="ddG [kcal/mol]", plot_type="bar")
        if amaps.get_exclusions():
            plots["excl"] = plotdata.PlotData("Free energy profile (exclusions)", xlabel="E1-E2  [kcal/mol]", ylabel="Free energy  [kcal/mol]")

        lra_de_st1_key = "lra_de_st1_%s" % (lra_l[0])
        lra_de_st2_key = "lra_de_st2_%s" % (lra_l[1])
        lra_lra_key    = "lra_lra_%s%s" % (lra_l[0], lra_l[1])
        lra_reo_key    = "lra_reo_%s%s" % (lra_l[0], lra_l[1])
        plots[lra_de_st1_key] = plotdata.PlotData("E2-E1 (lambda=%s)" % lra_l[0], xlabel="Energy type", ylabel="Potential energy  [kcal/mol]", plot_type="bar")
        plots[lra_de_st2_key] = plotdata.PlotData("E2-E1 (lambda=%s)" % lra_l[1], xlabel="Energy type", ylabel="Potential energy  [kcal/mol]", plot_type="bar")
        plots[lra_lra_key] = plotdata.PlotData("LRA (l=%s -> l=%s)" % (lra_l[0], lra_l[1]), xlabel="Energy type", ylabel="Potential energy  [kcal/mol]", plot_type="bar")
        plots[lra_reo_key] = plotdata.PlotData("Reorganization energy (l=%s -> l=%s)" % (lra_l[0], lra_l[1]), xlabel="Energy type", ylabel="Potential energy  [kcal/mol]", plot_type="bar")
        
        plots["egapl"] = plotdata.PlotData("Sampling", xlabel="Lambda", ylabel="E1-E2  [kcal/mol]")
        plots["pts_egap"] = plotdata.PlotData("Sampling v2", xlabel="Egap", ylabel="Number of points")
        plots["dgl"] = plotdata.PlotData("dG vs Lambda", xlabel="Lambda", ylabel="Free energy  [kcal/mol]")
        plots["rxy"] = plotdata.PlotData("Reactive distance", xlabel="E1-E2  [kcal/mol]", ylabel=u"Rxy  [Å]")
        for col in ams[0].get_E1_lambda().get_column_titles()[1:]:  # get the column names from the first replica (0th is lambda)
            key = "e1l_%s" % col
            plots[key] = plotdata.PlotData("E1 vs Lambda (%s)" % col, xlabel="Lambda (state 1)", ylabel="E1 (%s)  [kcal/mol]" % col)
            key = "e2l_%s" % col
            plots[key] = plotdata.PlotData("E2 vs Lambda (%s)" % col, xlabel="Lambda (state 1)", ylabel="E2 (%s)  [kcal/mol]" % col)
    
    
        for am in sorted(ams, key=lambda x: x.get_dirname()):
    
            path = os.path.relpath(am.get_dirname())
            print "#### Working on %s ####" % path
            try:
                rxy = am.get_rxy_dE()
                columns = rxy.get_columns()
                plots["rxy"].add_subplot(path,columns[0],columns[1])  #0th is dE, 1st is rxy
            except QAnalyseMapError as e:
                print "WARNING - An exception was raised during the analysis:\n> %s" % e
            try:
                e1l = am.get_E1_lambda()
                columns = e1l.get_columns()
                for i,colname in enumerate(e1l.get_column_titles()[1:]):
                    key = "e1l_%s" % colname
                    plots[key].add_subplot(path,columns[0],columns[i+1])  #0th is lambda
            except QAnalyseMapError as e:
                print "WARNING - An exception was raised during the analysis:\n> %s" % e
    
            try:
                e2l = am.get_E2_lambda()
                columns = e2l.get_columns()
                for i,colname in enumerate(e2l.get_column_titles()[1:]):
                    key = "e2l_%s" % colname
                    plots[key].add_subplot(path,columns[0],columns[i+1])  #0th is lambda
            except QAnalyseMapError as e:
                print "WARNING - An exception was raised during the analysis:\n> %s" % e
        
            try:
                dgl = am.get_dG_lambda()
                columns = dgl.get_columns()
                plots["dgl"].add_subplot(path,columns[0],columns[1])  #0th is lambda, 1st is dG
            except QAnalyseMapError as e:
                print "WARNING - An exception was raised during the analysis:\n> %s" % e
    
            try:
                egapl = am.get_Egap_lambda()
                columns = egapl.get_columns()
                plots["egapl"].add_subplot(path,columns[0],columns[1])  #0th is lambda, 1st is Egap
            except QAnalyseMapError as e:
                print "WARNING - An exception was raised during the analysis:\n> %s" % e

            try:
                ptsegap = am.get_points_Egap()
                columns = ptsegap.get_columns()
                plots["pts_egap"].add_subplot(path,columns[0],columns[1])  #0th is egap, 1st is points
            except QAnalyseMapError as e:
                print "WARNING - An exception was raised during the analysis:\n> %s" % e
        
            try:
                dgde = am.get_dG_dE()
                columns = dgde.get_columns()
                plots["dgde"].add_subplot(path,columns[0],columns[1])  #0th is Egap, 1st is dG
                if amaps.get_exclusions():
                    plots["excl"].add_subplot(path,columns[0],columns[1])  #0th is Egap, 1st is dG
            except QAnalyseMapError as e:
                print "WARNING - An exception was raised during the analysis:\n> %s" % e
        
            try:
                # save LRA and REORG energies
                fn = os.path.join( path, "qa.LRA.dat" )
                backup = backup_file(fn)
                if backup:
                    print "Backed up '%s' to '%s'" % (fn, backup)
                open(fn,"w").write( str(am.get_lra(lra_l[0], lra_l[1])) )
                print "Wrote qa.LRA.dat"
            except QAnalyseMapError as e:
                print "WARNING - An exception was raised during the analysis:\n> %s" % e
        
            if args.gc != None:
                try:
                    # save the group contributions
                    fn = os.path.join( path, "qa.GroupCon.dat" )
                    print "Calculating group contributions at lambdas %s, %s. This might take some time..." % (gc_l[0], gc_l[1])
                    backup = backup_file(fn)
                    if backup:
                        print "Backed up '%s' to '%s'" % (fn, backup)
                    open(fn,"w").write( str( am.get_group_contributions(gc_l[0], gc_l[1], first_res=args.gc[0], last_res=args.gc[1], qmaskfile=args.gc_mask)) )
                    print "Wrote qa.GroupCon.dat"
                except QAnalyseMapError as e:
                    print "WARNING - An exception was raised during the analysis:\n> %s" % e

            if args.gc_solv != False:
                try:
                    # save the group contributions
                    fn = os.path.join( path, "qa.GroupConSolv.dat" )
                    backup = backup_file(fn)
                    if backup:
                        print "Backed up '%s' to '%s'" % (fn, backup)
                    print "Calculating group contributions for solvent at lambdas %s, %s. This might take some time..." % (gc_l[0], gc_l[1])
                    open(fn,"w").write( str( am.get_group_contributions(gc_l[0], gc_l[1], solvent=True) ) )
                    print "Wrote qa.GroupConSolv.dat"
                except QAnalyseMapError as e:
                    print "WARNING - An exception was raised during the analysis:\n> %s" % e
    

        # save the average lras
        try:
            print "Calculating the average and stdev of lras..."
            lras = amaps.get_average_LRAs(lra_l[0], lra_l[1])
            columns = lras.get_columns()
            plots[lra_de_st1_key].add_subplot("average",columns[0],columns[1],yerror=columns[2])
            plots[lra_de_st2_key].add_subplot("average",columns[0],columns[3],yerror=columns[4])
            plots[lra_lra_key].add_subplot("average",columns[0],columns[5],yerror=columns[6])
            plots[lra_reo_key].add_subplot("average",columns[0],columns[7],yerror=columns[8])
        except QAnalyseMapError as e:
            print "WARNING - An exception was raised during the analysis:\n> %s" % e
    
        if args.gc != None:
            # save the average group contributions
            try:
                print "Calculating the average and stdev of GCs..."
                gcs = amaps.get_average_GCs(gc_l[0], gc_l[1], first_res=args.gc[0], last_res=args.gc[1], cached=True, qmaskfile=args.gc_mask)
                columns = gcs.get_columns()
                plots[gc_key].add_subplot("mean(reps)",columns[0],columns[3],yerror=columns[4])  # 0th is Resid, 3rd is El_mean, 4th is El_stdev
                plots[gc_vdw_key].add_subplot("mean(reps)",columns[0],columns[1],yerror=columns[2])  # 0th is Resid, 3rd is Vdw_mean, 4th is Vdw_stdev
            except QAnalyseMapError as e:
                print "WARNING - An exception was raised during the analysis:\n> %s" % e

        if args.gc_solv != False:
            # save the average group contributions for solvent molecules
            try:
                print "Calculating the average and stdev of solvent GCs..."
                gcs = amaps.get_average_GCs(gc_l[0], gc_l[1], cached=True, solvent=True)
                columns = gcs.get_columns()
                plots[gc_solv_key].add_subplot("mean(reps)",columns[0],columns[3],yerror=columns[4])  # 0th is Resid, 3rd is El(TS-R)_mean, 4th is El(TS-R)_stdev
            except QAnalyseMapError as e:
                print "WARNING - An exception was raised during the analysis:\n> %s" % e
    
    print amaps.get_summary()

    if ams:
        exclusions = amaps.get_exclusions()
        if exclusions:
            min1,ts,min2 = amaps.get_extrema_lambdas()
            print "\nExclusions:"
            print "%-30s %10s %10s %10s %10s %10s %10s %10s" % ("", "dGa_mean", "dGa_stdev", "dG0_mean", "dG0_stdev", "RS_l", "TS_l", "PS_l")
            print "%-30s %10.2f %10.2f %10.2f %10.2f %6.2f %6.2f %6.2f" % ("None", amaps.get_dGa_mean(), amaps.get_dGa_stdev(), amaps.get_dG0_mean(), amaps.get_dG0_stdev(), min1, ts, min2)
            for name,xamaps in exclusions.iteritems():
                min1,ts,min2 = xamaps.get_extrema_lambdas()
        # print out
                print "%-30s %10.2f %10.2f %10.2f %10.2f %6.2f %6.2f %6.2f" % (name, xamaps.get_dGa_mean(), xamaps.get_dGa_stdev(), xamaps.get_dG0_mean(), xamaps.get_dG0_stdev(), min1, ts, min2)
        
        # save to plots
                for am in sorted(xamaps.get_analysed(), key=lambda x: x.get_dirname()):
                    path = os.path.relpath(am.get_dirname())
                    try:
                        dgde = am.get_dG_dE()
                        columns = dgde.get_columns()
                        plots["excl"].add_subplot(path+"|"+name,columns[0],columns[1])  #0th is Egap, 1st is dG
                    except QAnalyseMapError as e:
                        print "WARNING - An exception was raised during the analysis:\n> %s" % e
                    
        # convert plots to json and save to qam.PlotData.json
        jsonenc = plotdata.PlotDataJSONEncoder(indent=2)
        backup = backup_file(args.plots_out)
        if backup:
            print "Backed up '%s' to '%s'" % (args.plots_out, backup)
        open(args.plots_out, 'w').write(jsonenc.encode(plots))
        print "\nWrote '%s'. Use q_plot.py to visualize the plots." % (args.plots_out)


