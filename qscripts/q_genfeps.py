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
# Generate FEP inputs (this script is pure manure shrapnel, read and use at your own risk)
#
# An example of a procedure file can be found in qscripts/template_examples/
#

import sys
import os
import shutil
import time
import tempfile
from lib.q_dyninp import QDynInput,QDynInputError
from qscripts_config import QScriptsConfig as QScfg
from lib.common import __version__
import q_pdbindex
import re
import random
import copy
import logging


class QGenfepsError(Exception):
    pass


class SpecialFormatter(logging.Formatter):
    FORMATS = {logging.DEBUG :"DBG: %(module)s: %(lineno)d: %(message)s",
               logging.WARNING : "\nWARNING: %(message)s\n",
               logging.INFO : "# %(message)s",
               'DEFAULT' : "%(message)s"}

    def format(self, record):
        self._fmt = self.FORMATS.get(record.levelno, self.FORMATS['DEFAULT'])
        return logging.Formatter.format(self, record)

def genfeps(fep_proc, relax_input, restraint="inp", 
            pdb=None, fep=None, runscript=None, 
            frames=QScfg.get("inputs", "fep_frames"), 
            repeats=QScfg.get("inputs", "num_repeats"),
            fromlambda=None, 
            prefix=QScfg.get("inputs", "prefix_rep"), 
            first_frame_eq=False, logger=None):

    frames = int(frames)
    repeats = int(repeats)
    if fromlambda != None:
        fromlambda = float(fromlambda)

    if not logger:
# log to console, only warnings
        logger = logging.getLogger(__name__)
        logger.setLevel(logging.WARNING)
        handler = logging.StreamHandler(sys.stdout)
        formatter = SpecialFormatter()
        handler.setFormatter(formatter)
        logger.addHandler(handler)


    # constants
    PREFIX_EQ = "equil_"
    PREFIX_FEP = "fep_"


    # check if files exist 
    for k, v in locals().iteritems():
        if k in ['pdb', 'fep_proc', 'fep', 'top', 'runscript', 'relax_input']:
            if v and not os.path.lexists(v):
                raise QGenfepsError("File '%s' doesn't exist." % v)
    
    if restraint not in ["top", "relax", "inp"]:
        raise QGenfepsError("Argument 'restraint' has to be either 'inp', 'top' or 'relax'")
    
    
    # find and replace atom placeholders. if not PDB was given to replace them, exit
    fep_proc_str = open(fep_proc, 'r').read()
    c = q_pdbindex.findPlaceholders(fep_proc_str)
    if c and not pdb:
        raise QGenfepsError("Found placeholders, but no PDB was given (--pdb):\n%s\n" % ", ".join(c))
    elif c:
        logger.info("These placeholders will be replaced with atom indices: " + ", ".join(c))
        try:
            fep_proc_str = q_pdbindex.convertToIndexes(fep_proc_str, pdb)
        except KeyError as e:
            raise QGenfepsError("Failed to replace placeholders: %s" % str(e))
    
    
    # make a nice header comment in each input file with the details about the arguments used
    argfep = "--fep %s" % fep if fep else ""
    argfl = "--fromlambda %s" % fromlambda if fromlambda else ""
    argrs = "--rs %s" % runscript if runscript else ""
    argpdb = "--pdb %s" % pdb if pdb else ""
    header_comment = """#
# Generated with {name}, version {version}, on {date}
# CWD: {cwd}
# CMDline: {name} {fep_proc} {relax_input} --rest {restraint} --frames {frames} --repeats {repeats} {fep} {fromlambda} {rs} {pdb}
#
""".format(name=os.path.basename( __file__ ), version=__version__, date=time.ctime(), cwd=os.getcwd(),
            fep_proc=fep_proc, relax_input=relax_input, restraint=restraint, frames=frames, repeats=repeats, fep=argfep, fromlambda=argfl, rs=argrs, pdb=argpdb)
    
    
    # get topology and fep and others from the last relaxation input
    top_file_abs, fep_file_abs, re_file_abs, rest_file = None, None, None, None
    lambda_initial = None
    try:
        c = QDynInput(open(relax_input, 'r').read())
    except QDynInputError as e:
        raise QGenfepsError("There is something wrong with the given input file (%s):\n%s" % (relax_input, str(e)))
    
    
    di = os.path.dirname( relax_input )
    try:
        files = c.parameters["files"]
        lambda_initial = float(c.parameters["lambdas"].split()[0])
        top_file_abs = os.path.join(di, files["topology"])
        re_file_abs = os.path.join(di, files["final"])
        fep_file_abs = os.path.join(di, files["fep"])
        if "restraint" in files:
            rest_file = os.path.join(di, files["restraint"])
    except KeyError as e:
        raise QGenfepsError("Parsing the relaxation input file failed, keyword missing... %s" % str(e))
    
    # check if the files actually exist
    for fn,descr in [(top_file_abs,'topology'), (fep_file_abs,'fep'), (re_file_abs,'final'), (rest_file, 'restraint')]:
        if fn and not os.path.lexists(fn):
            raise QGenfepsError("When parsing the input, found this filename '%s' next to the '%s' command. Unfortunately, the file doesnt exist..." % (fn,descr))
    
    # change the FEP (if debugging your system, you might want to use an old relax and not waste 100M core hours when changing a soft core value in the fep)
    if fep: fep_file_abs = os.path.abspath(fep)
    
    # change the inital lambda (this is not recommended, the system should be properly relaxed at a particual lambda before doing FEP)
    if fromlambda != None:
        lambda_initial = float(fromlambda)
        if lambda_initial > 1.0 or lambda_initial < 0.0:
            raise QGenfepsError("Lambda value is bogus, are you on drugs?")
    
    # create lambda values, find the closest to the starting one and rearrange accordingly
    # [0.0, 0.02, 0.04, ... 0.98, 1.0]  for frames==51
    lambdas = [ float(num)/(frames-1) for num in xrange(0,frames) ]
    
    # [2,]   for lambda_initial == 0.04 (or close to 0.04) and frames==51
    l_i = [ i for i in xrange(0,frames) if (abs(lambdas[i]-lambda_initial) <= (1.0/frames)) ]
    # there should be only one 
    l_i = l_i[0]
    lambda_initial = lambdas[l_i]
    
    # [0.02, 0.0, ] for the case of lambda_initial == 0.04 and frames == 51
    forward_lambdas = list(reversed(lambdas[0:l_i]))
    # [0.06, 0.08, ..., 1.0] for the case of lambda_initial == 0.04 and frames == 51
    backward_lambdas = lambdas[l_i+1:]
    
    lambdas = [lambda_initial,] + forward_lambdas + backward_lambdas
    
    
    # print out some useful information
    logger.info("Using restart file: %s " % os.path.relpath(re_file_abs))
    logger.info("Using topology file: %s " % os.path.relpath(top_file_abs))
    logger.info("Using FEP file: %s " % os.path.relpath(fep_file_abs))
    logger.info("Starting from lambda value (state 1): %s " % lambda_initial)
    logger.info("Number of FEP frames: %s " % frames)
    
    
    # create a temporary directory to store the files that are identical in all replicas - top, fep, runscript, relax restart, restraint file (if any)
    # and copy the common files
    TMPDIR = tempfile.mkdtemp()
    top_fn = os.path.basename(top_file_abs)
    fep_fn = os.path.basename(fep_file_abs)
    relax_re_fn = "cont_" + os.path.basename(re_file_abs)
    shutil.copy2(top_file_abs, TMPDIR)
    shutil.copy2(fep_file_abs, TMPDIR)
    shutil.copy2(re_file_abs, os.path.join(TMPDIR, relax_re_fn))
    if runscript:
        shutil.copy2(runscript, TMPDIR)
    else:
        logger.info("No Q runscript given.")
    
    # handle the whole restraint coordinates crap... 
    # rest_file is the file from the relaxation input (if any)
    # rest_fn is the basename of the restraints file (either from input or relaxed.re.rest), or None if rest. to topology
    if restraint == "relax":
        logger.info("Restraining to: relaxation")
        rest_fn = "cont_" + os.path.basename(re_file_abs) + ".rest"
        shutil.copy2(re_file_abs, os.path.join(TMPDIR, rest_fn))
    elif restraint == "top":
        logger.info("Restraining to: topology")
        rest_fn = None
    else: # default, from input
        if rest_file:
            logger.info("Restraining to: %s (from input)" % os.path.relpath(rest_file))
            rest_fn = "cont_" + os.path.basename(rest_file)
            shutil.copy2(rest_file, TMPDIR)
        else:
            logger.info("Restraining to: topology (from input)")
            rest_fn = None
    
    
    
    # parse the proc file
    general_inp = []
    eq_steps_inps = [ [], ]
    fep_inp = []
    script_vars = {}
    
    section = ""
    for line in fep_proc_str.split("\n"):
        # remove comments and strip whitespaces. 
        line = re.split("#|\!", line)[0].strip()
        # empty lines are useless
        if line == "":
            continue   
        # found a section
        if line[0] == "{":
            section = line.strip("{}").lower()
            continue
    
        if not section:  
            raise QGenfepsError("Parsing the procedure file failed... This line: '%s' is not inside any section:" % line)
    
        if section == "script_vars":
            c = line.split()
            var,value = c[0], " ".join( c[1:] )
            script_vars[var] = value
        elif section == "general":
            general_inp.append(line)
        elif section == "steps_equil":
            if "__________" in line:
                eq_steps_inps.append( [] )
            else:
                eq_steps_inps[-1].append( line )
        elif section == "fep":
            fep_inp.append(line)
        else:
            raise QGenfepsError("Parsing the procedure file failed: Unsupported section: '%s'" % (section))
    
    
    # check for steps with no parameters ( too many _________ lines ) and remove them
    for i in range(len(eq_steps_inps)-1, -1, -1):
        if not eq_steps_inps[i]:
            eq_steps_inps.pop(i)
    
    # check for missing sections
    for l,n in ("general_inp", "GENERAL"), ("eq_steps_inps", "STEPS_EQUIL"), ("fep_inp", "FEP"):
        if not l:
            raise QGenfepsError("Parsing the procedure file failed: Section '%s' is missing" % (n))
    
    
    # join lists of lines to strings and replace the placeholders
    script_variables = sorted(script_vars.items(), reverse=True)
    gen_inp_s = "\n".join(general_inp) 
    fep_inp_s = "\n".join(fep_inp) 
    eq_steps_inps_s =  ["\n".join(eq_step_inp) for eq_step_inp in eq_steps_inps]
    
    for placeholder,value in script_variables:
        gen_inp_s = gen_inp_s.replace(placeholder, value)
        fep_inp_s = fep_inp_s.replace(placeholder, value)
        for eq_step_inp_s in eq_steps_inps_s:
            eq_step_inp_s = eq_step_inp_s.replace(placeholder, value)
    
    
    ####################
    # make equil. inputs
    eq_steps = []
    for step_n, eq_step_inp_s in enumerate(eq_steps_inps_s):
        # create the files section
        final = "%s%03d_%4.3f.re" % (PREFIX_EQ, step_n, lambda_initial)
        dcd = "%s%03d_%4.3f.dcd" % (PREFIX_EQ, step_n, lambda_initial)
        files = { 'final'      : final,
                  'trajectory' : dcd,
                  'topology'   : top_fn,
                  'fep'        : fep_fn }
    
        if first_frame_eq: 
            files['energy'] = "%s%03d_%4.3f.en" % (PREFIX_EQ, step_n, lambda_initial)
    
        if rest_fn:
            files["restraint"] = rest_fn
    
        if step_n != 0:
            files["restart"] = "%s%03d_%4.3f.re" % (PREFIX_EQ, step_n-1, lambda_initial)
        else:
            files["restart"] = relax_re_fn
    
        # parse the general input and update with step input and files section
        try:
            inp = QDynInput(gen_inp_s)
            inp.update(eq_step_inp_s)
            if "energy" in inp.parameters["intervals"]:
                files["energy"] = "%s%03d_%4.3f.en" % (PREFIX_EQ, step_n, lambda_initial)
            elif first_frame_eq:
                raise QGenfepsError("Argument first_frame_eq requires the energy printout defined in the intervals section of the equilibration (e.g. 'energy   10')")
    
            inp.update(parameters = { "files": files } )
            inp.update(parameters = { "lambdas": "%9.7f %9.7f" % (lambda_initial, 1-lambda_initial) })
        except QDynInputError as e:
            raise QGenfepsError("Problem with equil. step no. %d: %s" % (step_n, str(e)))
    
        # get the input string 
        try:
            inpstr = inp.get_string()
        except QDynInputError as e:
            raise QGenfepsError("Error in equil. step %d: %s" % (step_n, str(e)))
    
        # check if random seed is not defined or is fixed in the first step
        if step_n == 0:
            if repeats > 1:
                if ("random_seed" not in inp.parameters["md"] or int(inp.parameters["md"]["random_seed"]) > 0):
                    raise QGenfepsError("Fixed random seed (or restart velocities) works only with one repeat (others will be identical).\n\
    Please use 'random_seed   -1' in your first equilibration step to generate random random seeds.")
    
            elif "random_seed" not in inp.parameters["md"]:
                logger.info("No random seed in first step of equilibration, using restart velocities.")
    
                if (not rest_file and rest_fn) or (not rest_fn and rest_file) or (rest_file and (os.path.basename(rest_file) != rest_fn)):
                    logger.warning("This will not be a true continuation run! The relaxation restraint does not match yours. Use 'inp' instead of 'top' or 'relax' for the restraint.")
    
        # append the input
        eq_steps.append(inp)
    
    
    
    #################
    # make FEP inputs
    en_filenames = []
    feps = []
    
    for step_n, lam in enumerate(lambdas):
        # create the files section
        final = "%s%03d_%4.3f.re" % (PREFIX_FEP, step_n, lam)
        dcd = "%s%03d_%4.3f.dcd" % (PREFIX_FEP, step_n, lam)
        en = "%s%03d_%4.3f.en" % (PREFIX_FEP, step_n, lam)
        files = { 'final'      : final,
                  'trajectory' : dcd,
                  'topology'   : top_fn,
                  'energy'     : en,
                  'fep'        : fep_fn }
    
        # if this step is in new direction (backwards) then set the previous lambda and step to initial
        if backward_lambdas and lam == backward_lambdas[0]:
            prev_fep = feps[0]
        elif step_n == 0:
            prev_fep = eq_steps[-1]
        else:
            prev_fep = feps[-1]

        # if this flag is set, all steps that point to the first step should point to the last eq step
        if first_frame_eq:
            if step_n == 1 or (backward_lambdas and lam == backward_lambdas[0]):
                prev_fep = eq_steps[-1]

        if rest_fn:
            files["restraint"] = rest_fn
    
        files["restart"] = prev_fep.parameters["files"]["final"]
    
    
        # update the parameters and check the input
        try:
            inp = QDynInput(gen_inp_s)
            inp.update(fep_inp_s)
            if "energy" not in inp.parameters["intervals"]:
                raise QGenfepsError("FEP stage requires the energy printout defined in the intervals section (e.g. 'energy   10')")
    
            inp.update(parameters = { "files": files } )
            inp.update(parameters = { "lambdas": "%9.7f %9.7f" % (lam, 1-lam) } )
            inp.check()
        except QDynInputError as e:
            raise QGenfepsError("Error in FEP step %d: %s\n" % (step_n, str(e)))
    
        # append the input
        feps.append(inp)
    
        # add the energy filename to the list
        en_filenames.append( inp.parameters["files"]["energy"] )
    
    
    # if first_frame_eq is set add the energy file and remove the first fep frame
    if first_frame_eq:
        logger.info("Replacing the first FEP frame with the last equilibration step")
        en_filenames[0] = eq_steps[-1].parameters["files"]["energy"]
        feps.pop(0)
    

    # check random seed in fep
    if "random_seed" in feps[0].parameters["md"] and int(feps[0].parameters["md"]["random_seed"]) < 1:
        logger.warning("Generating random seeds in FEP inputs. Are you sure this is ok?")
    
    
    
    # write a file that contains the names of all energy files in proper order 
    # this file is used later by q_mapper.py
    # sort the enfiles according to lambda (1.0 -> 0.0) so that the mapping will always go from reactants to products
    enfiles_lambdas = sorted([ (enf.split("_")[-1], i) for i,enf in enumerate(en_filenames) ], reverse=True)
    en_filenames_sorted = [ en_filenames[i] for l,i in enfiles_lambdas ]
    enf = os.path.join( TMPDIR, QScfg.get("files", "en_list_fn") )
    open(enf, 'w').write("\n".join(en_filenames_sorted))
    
    
    # create directories for repeats/replicas (rep_000,rep_001,rep_002...)
    # copy everything from TMPDIR (topology, fep file, relax restart and restraint file (if any))
    # create the eq and fep inputs
    
    # first check for existing directories
    for num in xrange(0,repeats):
        rep = "%s%03d" % (prefix, num)
        if os.path.lexists(rep):
            raise QGenfepsError("Directory '%s' exists. Please (re)move it or change the prefix with --prefix." % rep)
    
    lsdir = os.listdir( TMPDIR )
    rep_dirs = []
    for num in xrange(0,repeats):
        rep = "%s%03d" % (prefix, num)
        os.mkdir(rep)
        # copy stuff from TMPDIR
        for f in lsdir:
            shutil.copy2(os.path.join(TMPDIR, f), rep)
    
        # create eq inputs
        for step_n,eq_step in enumerate(eq_steps):
    
    # check if random seed is a fixed value or not (generate random or fail)
            eqs = copy.deepcopy(eq_step)  # a copy
            if "random_seed" in eqs.parameters["md"] and int(eqs.parameters["md"]["random_seed"]) < 1:
                rs = random.randint(1,1000000)
                eqs.update(parameters = { "md": {"random_seed" : rs } })
    
            try:
                s = eqs.get_string()
            except QDynInputError as e:
                raise QGenfepsError("Error in step %d: %s" % (step_n, str(e)))
            fn = os.path.join(rep, "%s%03d_%4.3f.inp" % (PREFIX_EQ, step_n, lambda_initial))
            s = header_comment + s
            open(fn, 'w').write(s)
    
        last_eq_fn = fn
        # create FEP inputs
        for step_n,fep in enumerate(feps):
            if first_frame_eq:
                step_n += 1
    
            fs = copy.deepcopy(fep)  # a copy
            if "random_seed" in fs.parameters["md"] and int(fs.parameters["md"]["random_seed"]) < 1:
                rs = random.randint(1,1000000)
                fs.update(parameters = { "md": {"random_seed" : rs } })
    
            try:
                s = fs.get_string()
            except QDynInputError as e:
                raise QGenfepsError("Error in step %d: %s\n" % (step_n, str(e)))
            lam = lambdas[step_n]  # feps was created in lambdas iteration
            fn = os.path.join(rep, "%s%03d_%4.3f.inp" % (PREFIX_FEP, step_n, lam))
            s = header_comment + s
            open(fn, 'w').write(s)
            
        logger.info("Created inputs for repeat/replica '%s'." % rep)
        rep_dirs.append( rep )
        
        
    
    # get the amount of storage that will be wasted
    # for this we need the atom count from the topology
    for line in open(os.path.join(TMPDIR,top_fn),'r').readlines(1024):
        if "no. of atoms, no. of solute atoms" in line:
            num_atoms_all = int( line.split()[0] )
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
            # calculate approx amount of data
    for i,step in enumerate(eq_steps+feps):
        data = {}
        for k,v in qintervals.iteritems():
            try:
                data[k] = int(step.parameters["md"]["steps"]) / int(step.parameters["intervals"][ v[0] ]) * v[2]
            except KeyError:
                data[k] = int(step.parameters["md"]["steps"]) / v[1] * v[2]  # default
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
    
    logger.info("Your runs will waste approx. %.2f MB of storage. Per replica: %.2f MB (trj: %.1f, log: %.1f, en: %.1f, rest: %.1f)"  % \
                  (sum(total_data.values())/CONV_MB*repeats, 
                   sum(total_data.values())/CONV_MB, total_data["trj"]/CONV_MB, 
                   total_data["log"]/CONV_MB, 
                   total_data["en"]/CONV_MB, 
                   total_data["rest"]/CONV_MB))
    
    # remove temporary directory
    shutil.rmtree(TMPDIR)
    
    return rep_dirs







if __name__ == "__main__":

    try:
        import argparse
    except ImportError:
        import lib.argparse as argparse
    
    parser = argparse.ArgumentParser()
    parser.add_argument('fep_proc', help = 'input file containing the FEP procedure')
    parser.add_argument('relax_input', help = 'path to the input file from the last relaxation step (to extract all the relevant filenames - re,top,fep and lambda values)')
    parser.add_argument('--rest', dest='restraint', help='Sequence restraints applied to topology (top) or relaxed structure (relax). Default is whatever is in the relaxation input (inp)', default="inp")
    parser.add_argument('--rs', dest='runscript', help='shell runscript for Q', default=argparse.SUPPRESS)
    parser.add_argument('--frames', type=int, help = 'number of frames (31,51,101,...). Default=%s' % QScfg.get("inputs", "fep_frames"), default=QScfg.get("inputs", "fep_frames"))
    parser.add_argument('--repeats', type=int, help = 'number of repeats/replicas. Default=%s' % QScfg.get("inputs", "num_repeats"), default=QScfg.get("inputs", "num_repeats"))
    parser.add_argument('--fep', help = 'FEP file (default is the one in the input file)', default=argparse.SUPPRESS)
    parser.add_argument('--fromlambda', type=float, help='Starting lambda for state 1. Example: --fromlambda 0.45  will go from 0.45,0.55 in both directions, to 1.0,0.0 and 0.0,1.0. Example2: --fromlambda 0.0 will drive the reaction in reverse direction. Default is the one in the relaxation input (usually 1.0 - starting from the reactants state).', default=argparse.SUPPRESS)
    parser.add_argument('--pdb', help='PDB file created with qprep. Used to replace $RESID.ATOMNAME$ placeholders with atom indices (ex. $512.N1$ -> 5514).', default=argparse.SUPPRESS)
    parser.add_argument('--prefix', help='Prefix for repeat/replica folder names (default=%s)' % QScfg.get("inputs", "prefix_rep"), default=QScfg.get("inputs", "prefix_rep"))
    parser.add_argument('--first_frame_eq', action="store_true", help='If set, the first FEP frame will be replaced by the last equilibration step. Default is unset.', default=False)
    
    
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
        gen_dirs = genfeps(**kwargs)
        #print gen_dirs
    except QGenfepsError as e:
        sys.stderr.write("\nERROR: %s\n" % str(e))

