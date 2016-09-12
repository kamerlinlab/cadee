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
# Module (and standalone CLI script) for analysing Q dynamics output (logfiles).
# The main class is QAnalyseDyns, it accepts log filenames.
# Each logfile is parsed and _QAnalyseDyn object is returned. 
# It contains temperature, energy, distance DataContainer objects
# To access this data, use functions get_temps, get_energies, get_q_energies and get_offdiags
#

import os
import re
from qscripts_config import QScriptsConfig as QScfg
from lib.common import DataContainer, np, backup_file, __version__
try:
    from collections import OrderedDict as ODict
except ImportError:
    import lib.OrderedDict as ODict



class QAnalyseDynError(Exception):
    pass


class QAnalyseDyns(object):
    def __init__(self, logfiles, timeunit="ps", stepsize=None):
        """
        Wrapper class for QanalyseDyn for analysing a sequence of log files.
        Args:
           logfile (list):  paths/filenames of Q logfiles
           timeunit (string):  fs,ps,ns (optional, default is ps)
           stepsize (float):  in case the on in Q is 0.000 (Q printout is a work of art)

        Usage:
        
        qads = QAnalyseDyns(["fep_000.log", "fep_001.log", ...], timeunit="ps")
        
        # get average temperatures by combining all logs and skipping 10% in each one
        temps_dc = qads.get_temps(percent_skip=10)   # returns Datacontainer with all logs combined
        temps = temps_dc.get_columns()   # returns the columns
        coltitles = temps_dc.get_column_titles()       # "Time", "T_tot", "T_free", "T_free_solute", "T_free_solvent"

        for i,colt in coltitles[1:]:
            print colt, np.mean( [ x for j,x in enumerate(temps[i]) if temps[0][j] >= midpoint ] )

        """
        if not logfiles:
            raise QAnalyseDynError("No logfiles given")
        self.analysed = []
        starttime = 0
        for i,logfile in enumerate(logfiles):
            try:
                qad = _QAnalyseDyn(logfile, timeunit, stepsize=stepsize, starttime=starttime)
                self.analysed.append(qad)
            except QAnalyseDynError as e:
                raise QAnalyseDynError("%s: %s" % (logfile, str(e)))
            starttime = qad.get_endtime()
        self.n_evb_states = self.analysed[0]._evb_states
        self.en_section_keys = self.analysed[0].map_en_section.keys()
        self.qen_section_keys = self.analysed[0].map_qen_section.keys()

    
    def get_temps(self, percent_skip=0):
        """
        Get temperatures from all logfiles combined.
        Args:
           percent_skip (int):  percent of datapoints in each logfile to skip
        Returns:
           temperatures (DataContainer)
        """

        # "Time", "T_tot", "T_free", "T_free_solute", "T_free_solvent"
        temps = DataContainer( list(self.analysed[0].data_temp.get_column_titles() ) )
        for qad in self.analysed:
            rows = qad.data_temp.get_rows()
            skip = int(round(len(rows)*percent_skip/100.0))
            for x in rows[skip:]:
                temps.add_row(x)
        return temps



    def get_temp_stats(self, percent_skip=0):
        """
        Returns temperature stats in string format (used for cmdline printout)
        for all logfiles combined
        """
        tt,tf,tf_solu,tf_solv = self.get_temps(percent_skip=percent_skip).get_columns( columns=("T_tot", "T_free", "T_free_solute", "T_free_solvent") )
        tt_mean = np.mean(tt) 
        tf_mean = np.mean(tf)
        tf_solu_mean = np.mean(tf_solu)
        tf_solv_mean = np.mean(tf_solv) 

        return """Temperature stats:
{0:20s}{1:>20s}{2:>20s}{3:>20s}
{4:20s}{5:>20.2f}{6:>20.2f}{7:>20.2f}
{8:20s}{9:>20.2f}{10:>20.2f}{11:>20.2f}
{12:20s}{13:>20.2f}{14:>20.2f}{15:>20.2f}
{16:20s}{17:>20.2f}{18:>20.2f}{19:>20.2f}
""".format("","Mean", "Stdev", "Max.Abs.Dev.",\
        "T_total", tt_mean, np.std(tt), max( map( lambda x: abs(x-tt_mean), tt) ),\
        "T_free", tf_mean, np.std(tf), max( map( lambda x: abs(x-tf_mean), tf) ),\
        "T_free_solute", tf_solu_mean, np.std(tf_solu), max( map( lambda x: abs(x-tf_solu_mean), tf_solu) ),\
        "T_free_solvent", tf_solv_mean, np.std(tf_solv), max( map( lambda x: abs(x-tf_solv_mean), tf_solv) ) )



    def get_offdiags(self, percent_skip=0):
        """
        Get distances from all logfiles combined.
        Args:
           percent_skip (int):  percent of datapoints in each logfile to skip
        Returns:
           distances (dict):  example { "13_31": DataContainer, "13_18": DataContainer }
        """

        dists = DataContainer( list(self.analysed[0].data_offdiags.get_column_titles() ) )
        for qad in self.analysed:
            rows = qad.data_offdiags.get_rows()
            skip = int(round(len(rows)*percent_skip/100.0))
            for x in rows[skip:]:
                dists.add_row(x)
        return dists


    def get_energies(self, e_type, percent_skip=0):
        """
        Get energies from all logfiles combined.
        Args:
           e_type (string):  keys in QAnalyseDyn.map_en_section dictionary
           percent_skip (int):  percent of datapoints in each logfile to skip
        Returns:
           energies (DataContainer)
        """

        energies = DataContainer( list(self.analysed[0].map_en_section[ e_type ].get_column_titles() ) )

        for qad in self.analysed:
            rows = qad.map_en_section[ e_type ].get_rows()
            skip = int(round(len(rows)*percent_skip/100.0))
            for x in rows[skip:]:
                energies.add_row(x)
        return energies


    def get_q_energies(self, qe_type, evb_state, percent_skip=0):
        """
        Get Q energies from all logfiles combined.
        Args:
           qe_type (string):  keys in QAnalyseDyn.map_qen_section dictionary
           evb_state (int):  1 or 2 or 3...
           percent_skip (int):  percent of datapoints in each logfile to skip
        Returns:
           energies (DataContainer)
        """

        col_tits = list(self.analysed[0].map_qen_section[qe_type][evb_state-1].get_column_titles() )
        energies = DataContainer(col_tits)

        for qad in self.analysed:
            rows = qad.map_qen_section[qe_type][evb_state-1].get_rows()
            skip = int(round(len(rows)*percent_skip/100.0))
            for x in rows[skip:]:
                energies.add_row(x)
        return energies




class _QAnalyseDyn(object):
    def __init__(self, logfile, timeunit="ps", stepsize=None, starttime=0):
        """
        Parses a Q dynamics logfile and extracts data (temperature, energies...)
        For interfacing, use QAnalyseDyns.

        Args:
           logfile (string):  path/filename of Q logfile  
           timeunit (string):  fs,ps,ns (optional, default is ps)
           stepsize (float):  in case the one in Q is 0.000 (Q printout is a work of art)


        Usage looks like this:

        # parse
        qad = QAnalyseDyns(.....).analysed[0]

        # print out nicely formatted temperature stats
        print qad.get_temp_stats()

        # get averages for seconds half (step >= 50% of steps) of all the temperatures 
        temps = qad.data_temp.get_columns()
        coltitles = qad.data_temp.get_olumn_titles()       # "Time", "T_tot", "T_free", "T_free_solute", "T_free_solvent"
        midpoint = int(temps[0][-1])/2                 # 0 == "Time", -1 == last frame
        for i,colt in coltitles[1:]:
            print colt, np.mean( [ x for j,x in enumerate(temps[i]) if temps[0][j] >= midpoint ] )

        # get the potential energy data and just print it out
        Epot = qad.data_E_SUM.get_columns( ["Time", "Potential"] )
        print Epot

        """

        # parse the logfile:
        # first the header using RE 
        # then dynamics (_parse_dyn()) line by line using the lazy generator in 'open' (less memory consumption and faster than regular expressions)

        self._logfile = logfile
        self._starttime = starttime

        self.MAP_TIME = { "fs": 1.0, "ps": 0.001, "ns": 0.000001 }
        if timeunit not in self.MAP_TIME:
            raise QAnalyseDynError("Timeunit has to be either 'fs', 'ps' or 'ns'")
        self._timeconv = self.MAP_TIME[timeunit]
            
        self._header=""
        try:
            with open(self._logfile,'r') as lf:
                for line in lf:
                    self._header += line
                    if "Initialising dynamics" in line:
                        break
        except OSError as e:
            raise QAnalyseDynError("\nCould not read the logfile: " + str(e))

        # use RE to get some info about the simulations
        m = re.search("Build number\s*([\d\.]+)", self._header)
        if m: self._qversion = m.group(1)
        else: 
            m = re.search('QDyn version 5.06', self._header)
            if m:
                self._qversion = '5.06' 
            else:
                raise QAnalyseDynError("Not a valid Q log file or Q version is very old...")

        m = re.search("Topology file      =\s*(\S+)", self._header)
        if m: self._topfile = m.group(1)
        else: raise QAnalyseDynError("Couldn't find the topology filename!?")

        m = re.search("Number of MD steps =\s*(\d+)", self._header)
        if m: self._md_steps = int(m.group(1))
        else: raise QAnalyseDynError("Couldn't find number of steps!?")

        m = re.search("Stepsize \(fs\)    =\s*([\d\.]+)", self._header)
        if m: self._stepsize = float(m.group(1))
        else: raise QAnalyseDynError("Couldn't find the stepsize!?")

        if not stepsize:
            if abs(self._stepsize - 0.0) < 10E-8:
                raise QAnalyseDynError("Can't convert steps to time, stepsize is 0.0 in the logfile (Q sucks). Set the stepsize please.")
        else:
            if self._stepsize:
                raise QAnalyseDynError("Will not override the non-zero stepsize in the logfile...")

        m = re.search("FEP input file     =\s*(\S+)", self._header)
        if m: self._fepfile = m.group(1)
        else: self._fepfile = None

        if self._fepfile:
            m = re.search("No. of fep/evb states    =\s*(\d+)", self._header)
            if m: self._evb_states = int(m.group(1))
            else: raise QAnalyseDynError("Couldn't find the number of states!?")

        offdsection = re.search("(No. of offdiagonal \(Hij\) functions =.*?^$)", self._header, re.MULTILINE | re.DOTALL).group(1)
        offdiags = re.findall("\s+\d+\s+\d+\s+(\d+)\s+(\d+)\s+[\d\.]+\s+[\d\.]+", offdsection)


        #
        # make datacontainer variables for storing all the data

        # offdiags
        offdiags = [ "%s_%s" % (a1,a2) for a1,a2 in offdiags ]
        self._tmp_offdiags = {}
        for k in offdiags:
            self._tmp_offdiags[k] = DataContainer(["Time", "Distance"]) 

        # temperature
        self.data_temp = DataContainer(["Time", "T_tot", "T_free", "T_free_solute", "T_free_solvent"])

        # energies
        self.data_E_solute = DataContainer(["Time", "El", "VdW", "Bond", "Angle", "Torsion", "Improper"])
        self.data_E_solvent = DataContainer(["Time", "El", "VdW", "Bond", "Angle", "Torsion", "Improper"])
        self.data_E_solute_solvent = DataContainer(["Time", "El", "VdW"])
        self.data_E_LRF = DataContainer(["Time", "El"])
        self.data_E_Q_atom = DataContainer(["Time", "El", "VdW", "Bond", "Angle", "Torsion", "Improper"])
        self.data_E_restraints = DataContainer(["Time", "Total", "Fix", "Solvent_rad", "Solvent_pol", "Shell", "Solute"])
        self.data_E_SUM = DataContainer(["Time", "Total", "Potential", "Kinetic"])

        # Q energies
        q_columns1 = ("Time", "Lambda", "El", "VdW")
        q_columns2 = ("Time", "Lambda", "El", "VdW", "Bond", "Angle", "Torsion", "Improper")
        q_columns3 = ("Time", "Lambda", "Total", "Restraint")

        self.data_EQ_Q, self.data_EQ_prot, self.data_EQ_wat, self.data_EQ_surr, self.data_EQ_any, self.data_EQ_SUM = [], [], [], [], [], []
        for i in range(self._evb_states):
            self.data_EQ_Q.append( DataContainer( q_columns1 ) )
            self.data_EQ_prot.append( DataContainer( q_columns1 ) )
            self.data_EQ_wat.append( DataContainer( q_columns1 ) )
            self.data_EQ_surr.append( DataContainer( q_columns1 ) )
            self.data_EQ_any.append( DataContainer( q_columns2 ) )
            self.data_EQ_SUM.append( DataContainer( q_columns3 ) )

        # mapping of energy types (label in the output) with containers
        self.map_en_section = { "solute":           self.data_E_solute,
                                "solvent":          self.data_E_solvent,
                                "solute-solvent":   self.data_E_solute_solvent,
                                "LRF":              self.data_E_LRF,
                                "Q-atom":           self.data_E_Q_atom,
                                "SUM":              self.data_E_SUM }

        self.map_qen_section = { "Q-Q":       self.data_EQ_Q,
                                 "Q-prot":    self.data_EQ_prot,
                                 "Q-wat":     self.data_EQ_wat,
                                 "Q-surr.":    self.data_EQ_surr,
                                 "Q-any":     self.data_EQ_any,
                                 "Q-SUM":     self.data_EQ_SUM }
 
        self._parse_dyn()
        d_dcs = self._tmp_offdiags.values()
        self.data_offdiags = DataContainer( ["Time",] + self._tmp_offdiags.keys() )
        for d_row in zip(* [ d_dcs[0].get_columns([0,])[0], ]  + [ d_dc.get_columns([1,])[0] for d_dc in d_dcs ]):
            self.data_offdiags.add_row(d_row)

        self._endtime = self._md_steps * self._timeconv + self._starttime

    def get_endtime(self):
        return self._endtime


    def _parse_dyn(self):
        """
        Parses the dynamics part of the logfile (used by init)
        """
        time = self._starttime
        t_free,t_tot = None,None
        insection = False
        with open(self._logfile, 'r') as lf:
            lf.seek( len(self._header) )
            for line in lf:
                l = line.split()
                if not l: continue
                if "Initialising dynamics" in line:
                    raise QAnalyseDynError("Found more than one logfile... Don't concatenate man...")

                if t_free != None: # second line with temps
                    try:
                        tf_solute = float(l[1])
                    except: # gas phase
                        tf_solute = 0 
                    try:
                        tf_solvent = float(l[3])
                    except: # gas phase
                        tf_solvent = 0

                    self.data_temp.add_row( (time, t_tot, t_free, tf_solute, tf_solvent) )
                    t_free, t_tot = None, None

                elif "Temperature at step" in line:  # first line with temps
                    line = line.replace("step", "step ")  # fix for large step numbers
                    l = line.split()
                    step = int(l[3].strip(":"))
                    time = step * self._stepsize * self._timeconv + self._starttime
                    t_tot, t_free = float(l[5]), float(l[7])

                elif "Energy summary at step" in line or "Q-atom energies at step" in line:
                    insection=True
                    step = int(l[5])
                    time = step * self._stepsize * self._timeconv + self._starttime

                elif "FINAL  Energy summary" in line or "FINAL Q-atom energies" in line:
                    insection=True
                    time = self._md_steps * self._stepsize * self._timeconv + self._starttime

                elif "==========================================================================" in line:
                    insection=False

                if not insection:
                    continue

                # always skip the 0th step 
                if step == 0: continue

                key = l[0]
                if key in self.map_en_section:
                    self.map_en_section[ key ].add_row( [time, ] + [ float(x) for x in l[1:] ] )
                elif key in self.map_qen_section:
                    evb_index = int(l[1]) - 1 
                    self.map_qen_section[ key ][ evb_index ].add_row( [time, ] + [ float(x) for x in l[2:] ] )
                elif "dist. between" in line:
                    a1,a2,dist = l[8],l[9],float(l[11])
                    k = "%s_%s" % (a1,a2)
                    self._tmp_offdiags[k].add_row( [time, dist ] )




if __name__ == "__main__":
    import sys
    from lib import plotdata
    try:
        from collections import OrderedDict as ODict
    except ImportError:
        import lib.OrderedDict as ODict
    try:
        import argparse
    except ImportError:
        import lib.argparse as argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("logfiles", nargs="+", help="Qdyn log files")
    parser.add_argument("--plots_out", dest="plots_out", help="output filename for plot data (default='%s')" % QScfg.get("files", "analysedyn_plots"), default=QScfg.get("files", "analysedyn_plots"))
    parser.add_argument("--stepsize", dest="stepsize", help="If the stepsize in your log is 0.000, define it with this flag.", default=None)
    parser.add_argument("--timeunit", dest="timeunit", help="Which units of time should the results be in? fs, ps or ns? Default is ps.", default="ps")
    parser.add_argument("--skip", dest="skip", type=int, help="Skip percentage of data points in each log. Default=0", default=0) 

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    args = parser.parse_args()

    for logfile in args.logfiles:
        if not os.path.lexists(logfile):
            print "FATAL! File '%s' doesn't exist" % logfile
            sys.exit(1)

    try:
        qads = QAnalyseDyns(args.logfiles, timeunit=args.timeunit, stepsize=args.stepsize)
    except QAnalyseDynError as e:
        print "Error: " + str(e)
        sys.exit(1)


    print qads.get_temp_stats()

    time_label = "Time [%s]" % args.timeunit
    plots = ODict()
    plots["temp"] = plotdata.PlotData("Temperature", xlabel=time_label, ylabel="T [K]")
    plots["offdiags"] = plotdata.PlotData("Offdiagonal distances", xlabel=time_label, ylabel="Distance [A]")
    
    t_dc = qads.get_temps(percent_skip=args.skip)
    t_cs, t_cts = t_dc.get_columns(), t_dc.get_column_titles()
    for i, t_ct in enumerate(t_cts[1:]):
        plots["temp"].add_subplot(t_ct, t_cs[0], t_cs[i+1])   # 0==Time

    
    d_dc = qads.get_offdiags(percent_skip=args.skip)
    d_cs, d_cts = d_dc.get_columns(), d_dc.get_column_titles()
    for i, d_ct in enumerate(d_cts[1:]):
        plots["offdiags"].add_subplot(d_ct, d_cs[0], d_cs[i+1])   # 0==Time

    
    for k in qads.en_section_keys:
        key = "E_%s" % k
        plots[key] = plotdata.PlotData("Energy: " + k, xlabel=time_label, ylabel="Energy [kcal/mol]")
        e_dc = qads.get_energies(k)
        e_cs, e_cts = e_dc.get_columns(), e_dc.get_column_titles()
        if e_cs:
            for i, e_ct in enumerate(e_cts[1:]):
                plots[key].add_subplot(e_ct, e_cs[0], e_cs[i+1])   # 0==Time


    for k in qads.qen_section_keys:
        for evb_state in range(1,qads.n_evb_states+1):
            key = "EQ%d_%s" % (evb_state, k)
            plots[key] = plotdata.PlotData("Q Energy: %s (state %d)" % (k, evb_state), xlabel=time_label, ylabel="Energy [kcal/mol]")
            qe_dc = qads.get_q_energies(k, evb_state, percent_skip=args.skip)
            qe_cs, qe_cts = qe_dc.get_columns(), qe_dc.get_column_titles()
            if qe_cs:
                for i, qe_ct in enumerate(qe_cts[1:]):
                    plots[key].add_subplot(qe_ct, qe_cs[0], qe_cs[i+1])   # 0==Time



    jsonenc = plotdata.PlotDataJSONEncoder(indent=2)
    backup = backup_file(args.plots_out)
    if backup:
        print "Backed up '%s' to '%s'" % (args.plots_out, backup)
    open(args.plots_out, 'w').write(jsonenc.encode(plots))
    print "\nWrote '%s'. Use q_plot.py to visualize the plots." % (args.plots_out)

