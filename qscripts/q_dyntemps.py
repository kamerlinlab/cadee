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
# Get temperatures from Qdyn5 logfile and print out the mean, stdev and median (for all points and for last 90%)
# 

from lib.common import np
import sys


total=[]
free=[]
solute=[]
solvent=[]

with open( sys.argv[1], 'r' ) as f:
    lines = f.readlines()
    for line in lines:
     ## Temperature at step       1:         T_tot=     319.1         T_free=     319.1
     ##                              T_free_solute=     377.8 T_free_solvent=     317.7
    
        if "Temperature at" in line:
            line = line.replace("step", " ")  # fix for large step numbers
            line = line.replace(":"," ")
            a = line.split()
            total.append( float( a[4] ) )
            free.append(float(  a[6] ) )
          
        if "T_free_solute=" in line:
            a = line.split()
            solute.append( float( a[1] ) )
            solvent.append( float (a[3] ) )


# remove 10 percent
n = len(total)
st = n/10

tn10 = total[st:] 
fn10 = free[st:] 
sln10 = solute[st:] 
stn10 = solvent[st:] 

if len(solute):
    total_mean = np.mean(total)
    free_mean = np.mean(free)
    solute_mean = np.mean(solute)
    solvent_mean = np.mean(solvent)
    tn10_mean = np.mean(tn10)
    fn10_mean = np.mean(fn10)
    sln10_mean = np.mean(sln10)
    stn10_mean = np.mean(stn10)
    print "%-20s %10s %10s %10s %10s  |         %-20s %10s %10s %10s %10s" % ("All %d frames" % n, 
                                                                              "Total", 
                                                                              "Free", 
                                                                              "Solute", 
                                                                              "Solvent", 
                                                                              "First 10% removed", 
                                                                              "Total", 
                                                                              "Free", 
                                                                              "Solute", 
                                                                              "Solvent")
    print "-" * 150
    print "%-20s %10.3f %10.3f %10.3f %10.3f  |         %-20s %10.3f %10.3f %10.3f %10.3f" % ("mean", 
                                                                                              total_mean,
                                                                                              free_mean,
                                                                                              solute_mean,
                                                                                              solvent_mean, 
                                                                                              "mean", 
                                                                                              tn10_mean,
                                                                                              fn10_mean,
                                                                                              sln10_mean, 
                                                                                              stn10_mean)
  
    print "%-20s %10.3f %10.3f %10.3f %10.3f  |         %-20s %10.3f %10.3f %10.3f %10.3f" % ("median", 
                                                                                              np.median(total), 
                                                                                              np.median(free), 
                                                                                              np.median(solute), 
                                                                                              np.median(solvent), 
                                                                                              "median", 
                                                                                              np.median(tn10), 
                                                                                              np.median(fn10), 
                                                                                              np.median(sln10), 
                                                                                              np.median(stn10) )
  
    print "%-20s %10.3f %10.3f %10.3f %10.3f  |         %-20s %10.3f %10.3f %10.3f %10.3f" % ("stdev", 
                                                                                              np.std(total), 
                                                                                              np.std(free), 
                                                                                              np.std(solute), 
                                                                                              np.std(solvent), 
                                                                                              "stdev", 
                                                                                              np.std(tn10), 
                                                                                              np.std(fn10), 
                                                                                              np.std(sln10), 
                                                                                              np.std(stn10) )
  
    print "%-20s %10.3f %10.3f %10.3f %10.3f  |         %-20s %10.3f %10.3f %10.3f %10.3f" % ("max_abs_dev", 
                                                                                              max( map( lambda x: abs(x-total_mean), total) ), 
                                                                                              max( map( lambda x: abs(x-free_mean), free) ), 
                                                                                              max( map( lambda x: abs(x-solute_mean), solute) ), 
                                                                                              max( map( lambda x: abs(x-solvent_mean), solvent) ), 
                                                                                              "max_abs_dev",
                                                                                              max( map( lambda x: abs(x-tn10_mean), tn10) ), 
                                                                                              max( map( lambda x: abs(x-fn10_mean), fn10) ), 
                                                                                              max( map( lambda x: abs(x-sln10_mean), sln10) ), 
                                                                                              max( map( lambda x: abs(x-stn10_mean), stn10) ) )
elif len(total):   # gas_phase, no solvent
    print "%-20s %10s %10s  |         %-20s %10s %10s" % ("All %d frames" % n, 
                                                          "Total", 
                                                          "Free", 
                                                          "First 10% removed", 
                                                          "Total", 
                                                          "Free")
    print "-" * 150
    print "%-20s %10.3f %10.3f  |         %-20s %10.3f %10.3f" % ("mean", 
                                                                  np.mean(total), 
                                                                  np.mean(free),  
                                                                  "mean", 
                                                                  np.mean(tn10), 
                                                                  np.mean(fn10) )
  
    print "%-20s %10.3f %10.3f  |         %-20s %10.3f %10.3f" % ("median", 
                                                                  np.median(total), 
                                                                  np.median(free),  
                                                                  "median", 
                                                                  np.median(tn10), 
                                                                  np.median(fn10) )
  
    print "%-20s %10.3f %10.3f  |         %-20s %10.3f %10.3f" % ("stdev", 
                                                                  np.std(total), 
                                                                  np.std(free), 
                                                                  "stdev", 
                                                                  np.std(tn10), 
                                                                  np.std(fn10) )
  
    print "%-20s %10.3f %10.3f  |         %-20s %10.3f %10.3f" % ("max_abs_dev", 
                                                                  max( map( lambda x: abs(x-total_mean), total) ), 
                                                                  max( map( lambda x: abs(x-free_mean), free) ), 
                                                                  "max_abs_dev",
                                                                  max( map( lambda x: abs(x-tn10_mean), tn10) ), 
                                                                  max( map( lambda x: abs(x-fn10_mean), fn10) ) )
else:
    print "No data!"
