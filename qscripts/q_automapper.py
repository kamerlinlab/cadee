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
# Automapper for Q
# Run the script with no arguments first.
#
# Varies EVB parameters Hij and Gas_shift until desired (experiment or qm)
# activation and reaction free energies are obtained.
# This script should be run in the directory that contains 
# the reference reaction's replicas (gas or water), each in its own directory
# (by default all subdirectories in the current dir are used for mapping 
# or if there are no subdirectories, the current dir is mapped. 
# This can be changed with --mapdirs dir1 dir2 dir3 ... 
# 
# Initial guess values for Hij and gas_shift should be relatively close to 
# their correct values (+-50) otherwise qfep crashes. If it doesn't converge,
# change the step size (--step), number of iterations (--iter) or the threshold (--threshold).
# For information about other parameters (bins, nt, skip, temperature...) see q_mapper.py
#
#

import os
import sys
import q_mapper   # this module does the actual mapping
import q_analysemaps
from qscripts_config import QScriptsConfig as QScfg
try:
    import argparse
except ImportError:
    import lib.argparse as argparse

# some constants
ITERATIONS = 10
STEPSIZE = 10.0
THRESHOLD = 0.005

def _call_mapper(qm, hij, gs):
    try:
        qm.set_hij(hij)
        qm.set_gas_shift(gs)
        (mapped, failed) = qm.q_map()
        qanalysemaps = q_analysemaps.QAnalyseMaps(mapped)
        if not qanalysemaps.get_analysed():
            print "\nAll replicas failed to analyse... Try changing the initial guess values (probably your simulations are crap though)... Also, look at qfep output files (%s)..." % QScfg.get("files", "qfep_out")
            print "\nHere are the error messages:"
            for mapdir, error in qanalysemaps.get_failed():
                print "%s -> %s" % (os.path.relpath(mapdir), error)
            sys.exit(1)

        return qanalysemaps

    except q_mapper.QMappingError as e:
        print "%s\n\nTry changing the initial guess values and/or the step size (--step). If that doesn't work you're screwed." % e
        sys.exit(1)


def do_iteration(qmapper, hij_init, gs_init, args, means_init=None):

    if means_init:
        means0 = means_init
    else:
        # Mapping with initial values for Hij and GS
        qan_maps = _call_mapper(qmapper, hij_init, gs_init)
        means0 = ( qan_maps.get_dGa_mean(), qan_maps.get_dG0_mean() )

    print "%10s %10s %10s %10s" % ("Hij", "Gas shift", "dG#", "dG0")
    print "%10.2f %10.2f %10.2f %10.2f" % (hij_init, gs_init, means0[0], means0[1])
    
    # Step1, changing gs=gs+step_size,
    hij1 = hij_init
    gs1 = gs_init + args.step_size
    qan_maps = _call_mapper(qmapper, hij1, gs1)
    means1 = ( qan_maps.get_dGa_mean(), qan_maps.get_dG0_mean() )
    print "%10.2f %10.2f %10.2f %10.2f" % (hij1, gs1, means1[0], means1[1])
    

    # Step2, changing hij=hij+step_size,
    hij2=hij_init+args.step_size
    gs2=gs_init
    qan_maps = _call_mapper(qmapper, hij2, gs2)
    means2 = ( qan_maps.get_dGa_mean(), qan_maps.get_dG0_mean() )
    print "%10.2f %10.2f %10.2f %10.2f" % (hij2, gs2, means2[0], means2[1])
    

    # calculate optimal shifts for hij and gs
    ka1=(means1[0]-means0[0])/args.step_size
    k01=(means1[1]-means0[1])/args.step_size
 
    ka2=(means2[0]-means0[0])/args.step_size
    k02=(means2[1]-means0[1])/args.step_size
    
    delta_dga = args.ref_dga - means0[0] 
    delta_dg0 = args.ref_dg0 - means0[1]
    
    try:
        dgs=((delta_dga/ka1)-(delta_dg0*ka2/ka1/k02))/(1-k01*ka2/ka1/k02)
        dhij = (delta_dg0 - dgs*k01)/k02
    except ZeroDivisionError:
        dgs = +1
        dhij = +1
    
    # Step 3, get results
    hij3 = hij_init + dhij
    gs3 = gs_init + dgs
    qan_maps = _call_mapper(qmapper, hij3, gs3)
    means3 = ( qan_maps.get_dGa_mean(), qan_maps.get_dG0_mean() )
    print "%10.2f %10.2f %10.2f %10.2f" % (hij3, gs3, means3[0], means3[1])
    
    return (qan_maps, hij3, gs3, means3)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('ref_dga', type=float, help = 'Activation free energy (reference)')
    parser.add_argument('ref_dg0', type=float, help = 'Reaction free energy (reference)')
    parser.add_argument('init_hij', type=float, help = 'first guess Hij coupling constant')
    parser.add_argument('init_gs', type=float, help = 'first guess gas shift constant')
    parser.add_argument('--nt', dest='nthreads', type=int, help = 'Number of threads (default = %s)' % QScfg.get("mapping", "nthread"), default=QScfg.get("mapping", "nthread"))
    parser.add_argument('--bin', dest = 'bins', type=int, help = 'number of bins (default=%s)' % QScfg.get("mapping", "bin"), default=QScfg.get("mapping", "bin") )
    parser.add_argument('--skip', dest = 'skip', type=int, help = 'number of points to skip in each frame (default=%s)' % QScfg.get("mapping", "skip"), default=QScfg.get("mapping", "skip"))
    parser.add_argument('--min', dest = 'minpts_per_bin', type=int, help = 'minimum points for bin (default=%s)' % QScfg.get("mapping", "minpts_per_bin"), default=QScfg.get("mapping", "minpts_per_bin"))
    parser.add_argument('--temp', dest = 'temperature', type=float, help = 'Temperature (default=%s)' % QScfg.get("mapping", "temp"), default=QScfg.get("mapping", "temp"))
    parser.add_argument('--step', dest = 'step_size', type=float, help = 'Step size (default = %s)' % STEPSIZE, default=STEPSIZE)
    parser.add_argument('--threshold', dest = 'threshold', type=float, help = 'Convergence threshold for dG# and dG0, default = %s' % THRESHOLD, default=THRESHOLD)
    parser.add_argument('--iter', dest = 'iterations', type=int, help = 'Max number of iterations (default = %s)' % ITERATIONS, default=ITERATIONS)
    parser.add_argument('--mapdirs', nargs="+", dest = 'mapdirs', help = 'Replica directories (default=all dirs in current folder) ', default=[])
    parser.add_argument('--mapper_logfile', dest = 'mapper_logfile', help = 'q_mapper logfile name (default=%s)' % QScfg.get("files", "mapper_logfile"), default=QScfg.get("files", "mapper_logfile"))
    parser.add_argument('--nosingle', dest = 'nosingle', action='store_true', help = 'Do not run the first iteration on only 1 replica') 

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    qmap_args = { "bins"            : args.bins,
                  "skip"            : args.skip,
                  "minpts_per_bin"  : args.minpts_per_bin,
                  "temperature"     : args.temperature,
                  "mapdirectories"  : args.mapdirs,
                  "verbose"         : False,
                  "nthreads"        : args.nthreads }

    print "Attempting to fit to dG# = %s and dG0 = %s\n(number of threads = %s, stepsize = %s, threshold = %s, max iterations = %s)\n" % (args.ref_dga, args.ref_dg0, args.nthreads, args.step_size, args.threshold, args.iterations)

    

    gs_init = args.init_gs
    hij_init = args.init_hij

    # create QMapper instance with all arguments
    try:
        qmapper = q_mapper.QMapper( hij_init, gs_init, **qmap_args )
    except q_mapper.QMappingError as e:
        print "Problem initializing qmapper: %s" % str(e)
        sys.exit(1)
    mapdirs = qmapper.get_mapdirs(relative_paths=True) 
    print "These directories will be used:\n%s\n" % ( ", ".join( mapdirs ) )

    # automap with only the first replica (when we have 3 or more) to get a better init guess quickly
    if not args.nosingle and len(mapdirs) > 2:
        try:
            print "Iteration #0, using only the first folder (disable this with --nosingle)."
            single_qmap_args = dict(qmap_args)
            single_qmap_args["mapdirectories"] = (mapdirs[0],)
            single_qmapper = q_mapper.QMapper( hij_init, gs_init, **single_qmap_args )
        except q_mapper.QMappingError as e:
            print "Problem initializing qmapper: %s" % str(e)
            sys.exit(1)

        (qan_maps_final, hij_init, gs_init, means_final) = do_iteration(single_qmapper, hij_init,gs_init, args)
        print "Switching to all directories..."


    # switch to all directories and iterative mapping
    iteration = 1
    means_init = None
    while (True):
        print "Iteration #%d" % iteration
        
        (qan_maps_final, hij_final, gs_final, means_final) = do_iteration(qmapper, hij_init, gs_init, args, means_init=means_init)

        # Check for convergence
        if abs(float(means_final[0])-float(args.ref_dga))< args.threshold and abs(float(means_final[1])-float(args.ref_dg0)) < args.threshold:

            # save results summary
            logstr = "%s\n%s" % (qmapper.get_details(), qan_maps_final.get_summary() )
            open(args.mapper_logfile,'w').write( logstr )

            # print out some stuff
            print "\n\nWell done, use this on the protein:"
            print qmapper.get_input_parms()
            if qan_maps_final.get_failed() or qmapper.get_failed():
                print "WARNING: Some directories failed to map..."
            print "\nLook at %s!" % args.mapper_logfile
            break
    
        iteration += 1
        if iteration >= args.iterations:
            print "Did not converge. Try changing the step (--step), increasing number of iterations (--iter) or lowering the treshold (--threshold)"
            print
            sys.exit(1)

        means_init = means_final
        hij_init = hij_final
        gs_init = gs_final



if __name__ == "__main__":
    try:
        main()
    except KeyboardInterrupt:
        print "\nCtrl-C detected. Quitting..."
        sys.exit(1)
