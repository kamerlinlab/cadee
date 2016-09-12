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
# 2015-11-04
#
# Generate relaxation inputs for qdyn5 
#
# An example of a procedure file (only required argument) can be found in qscripts/template_examples/
#

import sys
import os
import time
import locale
import shutil
import tempfile
from lib.q_dyninp import QDynInput,QDynInputError
from qscripts_config import QScriptsConfig as QScfg
from lib.common import __version__
import q_pdbindex
import re
import random
import logging

try:
    from collections import OrderedDict as ODict
except ImportError:
    import lib.OrderedDict as ODict

class QGenrelaxError(Exception):
    pass

class SpecialFormatter(logging.Formatter):
    FORMATS = {logging.DEBUG :"DBG: %(module)s: %(lineno)d: %(message)s",
               logging.WARNING : "\nWARNING: %(message)s\n",
               logging.INFO : "# %(message)s",
               'DEFAULT' : "%(message)s"}

    def format(self, record):
        self._fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        return logging.Formatter.format(self, record)

def genrelax(relax_proc, top=None, fep=None, 
               runscript=None, pdb=None, cont=None, 
               outdir=QScfg.get("inputs", "relax_dir"), 
               logger=None):

    if not logger:
# log to console, only warnings
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.WARNING)
        handler = logging.StreamHandler(sys.stdout)
        formatter = SpecialFormatter()
        handler.setFormatter(formatter)
        logger.addHandler(handler)

    for k, v in locals().iteritems():
        if k in ['pdb', 'cont', 'relax_proc', 'fep', 'top', 'runscript', 'relax_input']:
            if v and not os.path.lexists(v):
                raise QGenrelaxError("File '%s' doesn't exist." % v)

    # constants
    PREFIX="relax_"
    DIR=os.path.join( os.getcwd(), outdir )
    if os.path.lexists(DIR):
        raise QGenrelaxError("Directory '%s' exists. Please (re)move it or set 'outdir'." % DIR)
        sys.exit(1)
    TMPDIR = tempfile.mkdtemp()
    
    header_comment = "# Generated with %s, version %s, on %s\n# CWD: %s\n# Cmdline: %s\n" % (os.path.basename(sys.argv[0]), __version__, time.ctime(), os.getcwd(), " ".join(sys.argv))
    
    # find and replace placeholders. if not PDB was given to replace them, exit
    relax_proc_str=open(relax_proc, 'r').read()
    c = q_pdbindex.findPlaceholders(relax_proc_str)
    if c and not pdb:
        raise QGenrelaxError("Found placeholders, but no PDB was given (--pdb):\n%s\n" % ", ".join(c))
    elif c:
        logger.info("These placeholders will be replaced with atom indices: " + ", ".join(c))
        try:
            relax_proc_str = q_pdbindex.convertToIndexes(relax_proc_str, pdb)
        except KeyError as e:
            raise QGenrelaxError("Failed to replace placeholders: %s" % str(e))
    
    
    # get topology and fep and others from previous input if given (--cont)
    if cont:
        if top:
            raise QGenrelaxError("'top' and 'cont' don't like each other. Difficult to continue with a different topology...")
        try:
            c = QDynInput(open(cont, 'r').read())
        except QDynInputError as e:
            raise QGenrelaxError("There is something wrong with the given input file (%s):\n%s" % (cont, str(e)))
        cont_files = c.parameters["files"]
        di = os.path.dirname( cont )
        top_fn = cont_files["topology"]
        re_fn = "cont_" + cont_files["final"]
        shutil.copy2( os.path.join(di,top_fn), TMPDIR )
        shutil.copy2( os.path.join(di,cont_files["final"]), os.path.join(TMPDIR, re_fn) )
        if fep:
            logger.warning("Using the fep file '%s', instead of the one found in the input\n" % (fep))
            fep_fn = os.path.basename(fep)
            shutil.copy2( fep, TMPDIR )
        else:
            try:
                fep_fn = cont_files["fep"]
                shutil.copy2( os.path.join(di,fep_fn), TMPDIR )
            except KeyError:
                logger.info("No FEP file found in the input")
    # or take the arguments
    else:
        if not top:
            raise QGenrelaxError("Please specify the topology file ('top') or specify 'cont' to continue from a previous relaxation.")
    
        cont_files = None
        top_fn = os.path.basename(top)
        shutil.copy2( top, TMPDIR )
        try:
            fep_fn = os.path.basename(fep)
            shutil.copy2( fep, TMPDIR )
        except AttributeError:
            logger.info("NOTE: No FEP file!")
    
    try:
        shutil.copy2( runscript, TMPDIR )
    except AttributeError:
        logger.info("No submission script was given.")
    
    
    general_inp = []
    steps_inps = [ [], ]
    script_vars = {}
    
    section = ""
    for line in relax_proc_str.split("\n"):
        # remove comments and strip whitespaces. 
        line = re.split("#|\!",line)[0].strip()
        # empty lines are useless
        if line == "":
            continue   
        # found a section
        if line[0] == "{":
            section = line.strip("{}").lower()
            continue
    
        if not section:  
            raise QGenrelaxError("Failed to parse '%s'... this line - '%s' is not inside any section:" % (relax_proc, line))
    
        if section == "script_vars":
            c = line.split()
            var = c[0]
            value = " ".join( c[1:] )
            script_vars[var] = value
        elif section == "general":
            general_inp.append(line)
        elif section == "steps":
            if "__________" in line:
                steps_inps.append( [] )
            else:
                steps_inps[-1].append( line )
    
    
    
    # check for steps with no parameters ( too many _________ lines ) and remove them
    for i in range(len(steps_inps)-1, -1, -1):
        if not steps_inps[i]:
            steps_inps.pop(i)
    
    # join lists of lines to strings and replace the placeholders
    gen_inp_s = "\n".join(general_inp) 
    for placeholder,value in script_vars.iteritems():
        gen_inp_s = gen_inp_s.replace(placeholder, value)
    
    step_inps_s = []
    for i,step_inp in enumerate(steps_inps):
        s = "\n".join(step_inp) 
        for placeholder,value in script_vars.iteritems():
            s = s.replace(placeholder, value)
        step_inps_s.append( s )
    
    # make and save the inputs
    steps = []
    overridden_prms_all = []
    step_n = 1
    inp_fns = []  # to store the filenames and use the return value
    for step_inp_s in step_inps_s:
        # create the files section
        final = "%s%03g.re" % (PREFIX, step_n)
        dcd = "%s%03g.dcd" % (PREFIX, step_n)
        files = { 'final'      : final,
                  'trajectory' : dcd,
                  'topology'   : top_fn }
        try:
            files['fep'] = fep_fn
        except NameError:
            pass
        if step_n != 1:
            prev_step = step_n - 1
            files["restart"] = "%s%03g.re" % (PREFIX, prev_step)
        elif cont_files:
            files["restart"] = re_fn
    
    
        try:
            # parse the general input 
            inp = QDynInput(gen_inp_s)
            # update the general parameters with step input, printout the overriden parms, update the files section
            overridden_prms = inp.update(step_inp_s)
            if overridden_prms:
                overridden_prms_all.append( (step_n, ", ".join(  ["%s:%s->%s" % (key,value_old,value_new) for key,(value_old,value_new) in overridden_prms.iteritems() ] )) )
    
            if "energy" in inp.parameters["intervals"]:
                files["energy"] = "%s%03g.en" % (PREFIX, step_n)
    
            inp.update(parameters = { "files": files } )
    
        except QDynInputError as e:
            raise QGenrelaxError("Problem with step no. %d:\n%s" % (step_n, str(e)))
    
        # set the random seed
        mdp = inp.parameters["md"]
        if "random_seed" in mdp and int(mdp["random_seed"]) < 1:
            rs = random.randint(1,1000000)
            inp.update(parameters = { "md": { "random_seed": rs } })
            logger.info("Generated random seed in step %d: %d" % (step_n, rs))
    
        # get the input string 
        try:
            inpstr = inp.get_string()
        except QDynInputError as e:
            raise QGenrelaxError("Error in step %d: %s\n" % (step_n, str(e)))
    
        inpfn = "%s%03g.inp" % (PREFIX, step_n)
        inp_fns.append( os.path.join(DIR,inpfn) )
        s = header_comment + inpstr
        open( os.path.join( TMPDIR, inpfn ) , 'w').write( s )
    
        steps.append(inp)
        step_n += 1
    
    try:
        shutil.copytree(TMPDIR,DIR)
    except OSError:
        raise QGenrelaxError("Cannot create directory '%s'." % DIR)
        sys.exit(1)
    # remove temporary directory
    shutil.rmtree(TMPDIR)
    logger.info("Created inputs %s%03g.inp - %s%03g.inp" % (PREFIX,1,PREFIX,len(steps)))
    
    
    
    
    # print some useful information
    
    if overridden_prms_all:
        logger.info("Overridden parameters:")
        for step_n, op in overridden_prms_all:
            logger.info("%d: %s" % (step_n, op))
    
    summary = """
Quick summary
{0:<10} {1:>5} {2:>10} {3:>10} {4:^10} {5:^10} {6:^10} {7:^30} {8:^10} {9:>10} 
""".format("Step", "T", "Stepsize", "Steps", "Seq.rest", 
"Dist.rest", "Ang.rest", "Shake", "Rand.Seed", "Data (MB)")
    locale.setlocale(locale.LC_ALL, '')
    restraints = []
    total_time = 0
    # print out how much data this run will produce, for this we need the atom count from the topology
    for line in open(os.path.join(DIR,top_fn),'r').readlines(1024):
        if "no. of atoms, no. of solute atoms" in line:
            num_atoms_all = int( line.strip().split()[0] )
            break
    
    REST_B_PER_ATOM = 48.0
    TRJ_B_PER_ATOM = 12.0
    EN_B_PER_STEP = 370.0  # very rough estimate, depends on Q version, it can double if group_contributions are calculated
    CONV_MB = 2**20
    
    qintervals = {"trj"  : ["trajectory",   100, num_atoms_all*TRJ_B_PER_ATOM ],      # q_parameter_key, q_default_value, approx_bytes_per_frame
                  "log"  : ["output",        10, 2000 ],   # 2000 is a very rough estimate of bytes_per_frame
                  "temp" : ["temperature",   10,  160 ],
                  "en"   : ["energy",        10,  EN_B_PER_STEP ],
                  "nb"   : ["non_bond",      10,   80 ] }
    total_data = { "trj": 0, "log": 0, "en": 0, "rest": 0}


    for i,step in enumerate(steps):
        nstep = i+1
        try:
    
            # get md parameters
            mdparms = step.parameters["md"]
            total_time += float(mdparms["stepsize"])*int(mdparms["steps"])
            random_seed = mdparms.get("random_seed", "")
            numsteps = int(mdparms["steps"])
    
            # get restraints
            step_rests = {"sequence_restraints" : [], "distance_restraints" : [], "angle_restraints" : []}
            for rest_type in step_rests.keys():
                for seqrest in step.parameters.get(rest_type, []):
                    if seqrest in restraints:
                        step_rests[rest_type].append( str(restraints.index(seqrest)+1) )
                    else:
                        restraints.append( seqrest )
                        step_rests[rest_type].append( str(len(restraints)) )
            seq = ",".join(step_rests["sequence_restraints"])
            dist = ",".join(step_rests["distance_restraints"])
            angle = ",".join(step_rests["angle_restraints"])
    
            # get shake parameters
            shake = []
            if mdparms.get("shake_solvent", "on") == "on":    # this is a Q default, hopefully it will not change
                shake.append("solvent")
            if mdparms.get("shake_hydrogens", "off") == "on":
                shake.append("hydrogens")
            if mdparms.get("shake_solute", "off") == "on":
                shake.append("solute")
            shake = ",".join(shake)
    
            # calculate approx amount of data
            data = {}
            for k,v in qintervals.iteritems():
                try:
                    data[k] = numsteps / int(step.parameters["intervals"][ v[0] ]) * v[2]
                except KeyError:
                    data[k] = numsteps / v[1] * v[2]  # default
                except ZeroDivisionError:
                    data[k] = 0   # no printout
                finally:
    # if energy or trajectory, check that files for output are defined, otherwise set the printout to 0
                    if v[0] in ("energy","trajectory") and not (v[0] in step.parameters["files"].keys()):
                        data[k] = 0
    
            trj_data = data["trj"]
            en_data = data["en"]
            log_data = (data["log"] + data["temp"] + data["nb"])
            rest_data = num_atoms_all * REST_B_PER_ATOM
            total_data["trj"] += trj_data
            total_data["log"] += log_data
            total_data["en"] += en_data
            total_data["rest"] += rest_data
    
            data = (trj_data + log_data + rest_data + en_data)/CONV_MB
    
    
            summary += "{0:<10} {1:>5} {2:>10} {3:>10} {4:^10} {5:^10} {6:^10} {7:^30} {8:^10} {9:>10.2f}\n".format(\
                    nstep, mdparms["temperature"], mdparms["stepsize"], locale.format('%d', numsteps, 1),\
                    seq, dist, angle, shake, random_seed, data)
        except KeyError as e:
            raise QGenrelaxError("You are missing either 'steps', 'temperature' or 'stepsize' in one of your relaxation steps. These parameters are quite important you know...")

    summary += "Restraints:\n"
    for i, rest in enumerate(restraints):
        summary += "%d: %s\n" % (i+1, rest)
    
    summary += """
Total time: {0} ps
Total wasted storage (wild approximation): \
{1:.2f} MB (trj: {2:.1f}, log: {3:.1f}, en: {4:.1f}, rest: {5:.1f}
""".format(total_time/1000.0, sum(total_data.values())/CONV_MB, total_data["trj"]/CONV_MB, 
total_data["log"]/CONV_MB, total_data["en"]/CONV_MB, total_data["rest"]/CONV_MB )

    for l in summary.split("\n"):
        logger.info(l)

    return inp_fns
    



if __name__ == "__main__":
    try:
        import argparse
    except ImportError:
        import lib.argparse as argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('relax_proc', help='input file containing the relaxation procedure')
    parser.add_argument('--top', dest='top', help='path to topology file', default=argparse.SUPPRESS)
    parser.add_argument('--fep', dest='fep', help='path to fep file (if there is one)', default=argparse.SUPPRESS)
    parser.add_argument('--rs', dest='runscript', help='submission script (for Slurm,SGE,Torque,...) ', default=argparse.SUPPRESS)
    parser.add_argument('--cont', dest='cont', help='Continue a previous relaxation, argument is the name of the last input file (e.g. "relax_012.inp")', default=argparse.SUPPRESS)
    parser.add_argument('--pdb', dest='pdb', help='PDB file created with qprep. Used to replace $RESID.ATOMNAME$ placeholders with atom indices (ex. $512.N1$ -> 5514).', default=argparse.SUPPRESS)
    parser.add_argument('--outdir', dest='outdir', help='Output directory name. Default=%s' % QScfg.get("inputs", "relax_dir"), default=QScfg.get("inputs", "relax_dir"))
    
    
    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    logger = logging.getLogger(__name__)
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler(sys.stdout)
    formatter = SpecialFormatter()
    handler.setFormatter(formatter)
    logger.addHandler(handler)
     
    kwargs = vars(args)
    kwargs["logger"] = logger

    try:
        print
        gen_inps = genrelax(**kwargs)
        #print gen_inps
    except QGenrelaxError as e:
        sys.stderr.write("\nERROR: %s\n" % str(e))


