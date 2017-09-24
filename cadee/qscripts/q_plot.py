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


from lib import plotdata
from qscripts_config import QScriptsConfig
try:
    import argparse
except ImportError:
    import lib.argparse as argparse
try:
    from collections import OrderedDict as ODict
except ImportError:
    import lib.OrderedDict as ODict
import os
import sys
import math
import Tkinter as Tk



class PlotApp():
    def __init__(self, parent, plotdata, plotdata_files):
        self.parent = parent
        self.plots = plotdata   # example:  { "dgde" : { 0 : QPlotData_instance (from file 1), 1 : QPlotData_instance (from file 2) }, ... }, where 0,1,... are indices of the filenames in plotdata_files
        self.plotdata_files = plotdata_files    #  [ "/home/.../pro/qa.PlotData.pickle", "/home/.../wat/qa.PlotData.pickle" ]
        self.nrows = 1
        self.ncols = 1
        self.blocked_draw = False
        self.subplot_lines = {}
        self.COLORS_ACTIVE = ("#555555","#F75D59","#1589FF", "black", "red", "blue")
        self.COLORS_INACTIVE = ("#aaaaaa","#F7bDb9","#a5a9FF", "#999999", "#FFaaaa", "#aaaaFF")
        

        self.lb1_entries = ODict()
        for plot_key, plot in self.plots.iteritems():
            self.lb1_entries[ plot.values()[0].title ] = plot_key

        self.lb1 = Tk.Listbox(self.parent, selectmode=Tk.EXTENDED, exportselection=0)

        for plot_title in self.lb1_entries.keys():
            self.lb1.insert(Tk.END, plot_title)
        
        self.lb1.pack(fill=Tk.Y, side=Tk.LEFT)
        self.lb2 = Tk.Listbox(self.parent, selectmode=Tk.EXTENDED, exportselection=0)
        self.lb2.pack(fill=Tk.Y, side=Tk.LEFT)
        
        self.figure = Figure(figsize=(5,4), dpi=100)

        
        self.canvas = FigureCanvasTkAgg(self.figure, master=self.parent)
        self.canvas.get_tk_widget().pack()
        self.canvas._tkcanvas.pack(fill=Tk.BOTH, expand=1)

        self.toolbar = NavigationToolbar2TkAgg( self.canvas, self.parent )
        self.toolbar.update()
        self.canvas._tkcanvas.pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)

        self.lb1.bind("<<ListboxSelect>>", self.on_select_lb1)
        self.lb2.bind("<<ListboxSelect>>", self.on_select_lb2)
        self.parent.bind("<Configure>", self.on_resize)


    def draw_legend(self):
        handls = []
        labls = []
        pos = "lower right"
        for i, plotdata_file in enumerate(self.plotdata_files):
            handls.append(mpatches.Patch(color=self.COLORS_ACTIVE[i]))
            labls.append("%d: %s" % (i, plotdata_file))
        self.figure.legend( handls, labls, pos, fontsize="xx-small" )


    def change_geometry(self):
    
        indices = [ int(sel) for sel in self.lb1.curselection() ]
        if not indices: return
        if len(indices) == 1:
            nc,nr = 1,1
        else:
            w,h = (float(x) for x in self.canvas.get_width_height())
        
            # get first guess layout from window size ( ncols =  width/ (height*1.5) )
            nc = math.ceil(w/(h*2.5))
            if nc > len(indices): nc = len(indices)
            nr = math.ceil(float(len(indices)) / nc)
        
            # check what the actual ratio of the plots is, and adjust accordingly
            if (w/nc) / (h/nr) > 2.5:
                nc += 1
            elif (w/nc) / (h/nr) < 0.5:  
                nc -= 1
            if not nc: nc = 1
            elif nc > len(indices): nc = len(indices)
    
            nr = math.ceil(float(len(indices)) / nc)
        
    # are they different then before?
        if nc != self.ncols or nr != self.nrows:
            self.ncols = nc
            self.nrows = nr
            return True
        return False
    


    def draw_plots(self):
    
# clear the figure 
        self.figure.clear()
# clear the subplot_lines dictionary
        self.subplot_lines = {}

# get keys for the selected plots in lb1
        plot_keys = [ self.lb1_entries[ self.lb1.get(int(index)) ] for index in self.lb1.curselection() ]

        for i,key in enumerate(plot_keys):   # example of plot_keys: [ "dgde", "egapl", "dgl", ... ]
            plt = self.figure.add_subplot(self.nrows,self.ncols,i+1)
            plots = self.plots[key]  

            for plot_number, plot in plots.iteritems():   # example of plots:   { 0: protein_plot, 1: protein2_plot, 2: water_plot }
                for subplot_label, subplot_data in plot.subplots.iteritems():   

                    if plot.plot_type == "line":
                        line, = plt.plot(subplot_data["xdata"],subplot_data["ydata"], color=self.COLORS_ACTIVE[plot_number] )
                    elif plot.plot_type == "bar":
                        width = 0.9/( len(plots) ) 
                        # color bar charts in lighter colors, they don't change anyhow
                        if isinstance(subplot_data["xdata"][0], unicode):
                            xind = range(0,len(subplot_data["xdata"]) )
                            xind = [ x - 0.45 + plot_number*(width) for x in xind ]
                            line = plt.bar(xind, subplot_data["ydata"], width = width, yerr=subplot_data["yerror"], color=self.COLORS_INACTIVE[plot_number] )
                            plt.set_xticks( xind )
                            plt.set_xticklabels( subplot_data["xdata"], rotation=70 )
                        else:
                            xind = [ x - 0.45 + plot_number*(width) for x in subplot_data["xdata"] ]
                            line = plt.bar(xind, subplot_data["ydata"], width = width, yerr=subplot_data["yerror"], color=self.COLORS_INACTIVE[plot_number] )

                    # add the line that was drawn to subplot_lines so that we can change color if lb2 selection changes
                    subplot_label = "%d/%s" % (plot_number, subplot_label)
                    if subplot_label not in self.subplot_lines.keys():
                        self.subplot_lines[subplot_label] = []
                    self.subplot_lines[subplot_label].append(line)
                
            plt.set_title(plot.title)
            plt.set_xlabel(plot.xlabel)
            plt.set_ylabel(plot.ylabel)
            
    
        if plot_keys:
            self.draw_legend()
            padding = len(self.plotdata_files) * 0.03
            self.figure.tight_layout(rect = [0,0+padding,1,1])
            self.canvas.draw()
        self.blocked_draw = False
    


    def on_select_lb1(self,event):

# remove and add all subplots to lb2 (1/rep_000,1/rep_001... 2/rep_000,2/rep_001...)
        self.lb2.delete(0, Tk.END)
# get keys for the selected plots in lb1
        plot_keys = [ self.lb1_entries[ self.lb1.get(int(index)) ] for index in self.lb1.curselection() ]

# iterate through all the selected plots (lb1)
# iterate through all the plots with the same key (protein dG_dE, water dG_dE, mutant dG_dE, ...)
# iterate through all subplot labels (rep_000, rep_001)
# if exists do not append it (subplots from different keys have the same label - protein dG_dE, protein dG_lambda, ... and should be combined)
        for i,key in enumerate(plot_keys):

            subplots_labels = []

            for plot_number, plot in self.plots[key].iteritems():
                for subplot_label in sorted(plot.subplots.keys()):
                    subplot_label = "%d/%s" % (plot_number, subplot_label)
                    if subplot_label not in subplots_labels:
                        subplots_labels.append(subplot_label)

        self.lb2.insert(0, *subplots_labels)
        self.lb2.selection_set(0, Tk.END)
    
        if not self.blocked_draw:
            self.blocked_draw = True
            self.change_geometry()
            self.draw_plots()


    def on_select_lb2(self,event):

# get selected subplots from lb2
        selected_subplots_keys = [ self.lb2.get(int(index)) for index in self.lb2.curselection() ]
        for subplot_key, subplot_line_list in self.subplot_lines.iteritems():
            for subplot_line in subplot_line_list:
                plot_number = int(subplot_key.split("/")[0])
                if subplot_key in selected_subplots_keys:
                    color = self.COLORS_ACTIVE[plot_number]
                    lw = 1.5
                else:
                    color = self.COLORS_INACTIVE[plot_number]
                    lw = 0.1
                try:
                    subplot_line.set_color(color)
                    subplot_line.set_linewidth(lw)
                except AttributeError:
                    pass   # bar chart doesn't support this 

        self.canvas.draw()

    
    def on_resize(self,event):
    
        if not self.blocked_draw:
            if self.change_geometry():
                self.blocked_draw = True   # it would otherwise redraw the whole canvas with every change in geometry
                self.parent.after(500, self.draw_plots)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("plotfiles",  nargs="+", help="qa.PlotData.pickle file(s) (up to six)")
    parser.add_argument("--export", dest="export", nargs="*", default=argparse.SUPPRESS, help="Export plots in Grace format to this directory: '%s'. Try without arguments, to see available plots." % QScriptsConfig.get("files", "plot_export_dir"))

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    args = parser.parse_args()

    if len(args.plotfiles) > 6:
        print "Currently, only six data files are supported. Deal with it..."
        sys.exit(1)
    if hasattr(args, "export") and len(args.plotfiles) > 1:
        print "Exporting works only with one plotfile at a time."
        sys.exit(1)

    allplots = ODict()    # example:  { "dgde" : { 0 : QPlotData_instance (from file 1), 1 : QPlotData_instance (from file 2) }, ... }, where 0,1,... are indices of the filenames in plotdata_files

    for pf_number, pf in enumerate(args.plotfiles):
        if not os.path.lexists(pf):
            print "File '%s' doesn't exist." % pf
            sys.exit(1)
        try:
            jsondec = plotdata.PlotDataJSONDecoder()
            plots = jsondec.decode(open(pf, 'r').read())
        except Exception as e:
            raise
            print "Could not read data file '%s'. Are you sure it is a .json file?" % pf
            sys.exit(1)
        if not isinstance(plots, ODict):
            print "Something is wrong with the data file '%s'. Aborting..." % pf
            sys.exit(1)
        for plot_id, plot in plots.iteritems():
            if not isinstance(plot, plotdata.PlotData):
                print "Something is wrong with the data file '%s'. Aborting..." % pf
                sys.exit(1)

            try:
                allplots[plot_id][pf_number] = plot
            except KeyError:
                allplots[plot_id] = ODict( [(pf_number,plot)] )

    plotdata_files = [ os.path.abspath(pf) for pf in args.plotfiles ]

    if not hasattr(args, "export"):
# run gui
        try:
            import matplotlib
            matplotlib.use('TkAgg')
            from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2TkAgg
            from matplotlib.figure import Figure
            import matplotlib.patches as mpatches
        except ImportError:
            print """Matplotlib is required for this script to work. Try something like:
            (ubuntu)    $ sudo apt-get install python-matplotlib
            (mac)       $ sudo brew install matplotlib
            (mac)       $ sudo port install py27-matplotlib
            (anything)  $ sudo pip install matplotlib
            or if you are working on a cluster, try loading a different python module..."""
            sys.exit(1)

        root = Tk.Tk()
        root.title("Q_Plot")
        app = PlotApp(root, allplots, plotdata_files)
        root.mainloop()
    else:

        exports=[]
        for ex in args.export:
            if not ex in allplots.keys() and not "all":
                print "Plot '%s' not found. Ignoring.." % (ex, )
            else: exports.append(ex)

        if not exports: # no arguments passed in
            print "\nAvailable arguments for --export:\n"
            print "\n".join( [ k for k in allplots.keys() ] )
            print "\nSpecial args: all\n"
            sys.exit(1)

        print "\nExporting the plots to ASCII Grace format (use xmgrace to plot them)\n"
        exdir = QScriptsConfig.get("files", "plot_export_dir")
        if not os.path.lexists(exdir):
            try:
                os.mkdir(exdir)
            except IOError as e:
                print "Could not create directory '%s': %s" % (exports, str(e))
                print "Quitting..."
                sys.exit(1)

        if "all" in exports: exall=True
        else: exall=False

        for plot_id, plots in allplots.iteritems():
            if (not exall) and (not plot_id in exports):  # don't save this one
                continue
            plot = plots[0]
            try:
                fn = os.path.join(exdir, "%s.agr" % plot_id)
                open(fn, 'w').write( plot.export_grace() )
                print "Wrote '%s' to %s" % (plot.title, fn)
            except IOError as e:
                print "Could not export '%s': %s" % (plot.title, str(e))
            
        

