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


try:
    from collections import OrderedDict as ODict
except ImportError:
    import lib.OrderedDict as ODict
import json
import sys


class PlotDataJSONEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, PlotData):
            return { "__type__":  "PlotData",
                     "title":     obj.title,
                     "plot_type": obj.plot_type,
                     "xlabel":    obj.xlabel,
                     "ylabel":    obj.ylabel,
                     "subplots":  obj.subplots }
        else:
            return json.JSONEncoder.default(self, obj)


class PlotDataJSONDecoder(json.JSONDecoder):
    def __init__(self):
        if sys.version_info < (2,7):
            # object_pairs_hook is supported only in version 2.7
            print "You need python 2.7 or later to run this script, sorry (it's json's fault)!"
            sys.exit(1)
        super(PlotDataJSONDecoder, self).__init__(object_pairs_hook=self.decode_plotdata)

    def decode_plotdata(self, d):
        d = ODict(d)
        if "__type__" not in d:
            return d
        t = d["__type__"]
        if t == "PlotData":
            pd = PlotData(d["title"], d["plot_type"], d["xlabel"], d["ylabel"])
            pd.subplots = d["subplots"]
            return pd
        else:
            return d


class PlotData(object):
    def  __init__(self, title, plot_type="line", 
                  xlabel=None, ylabel=None):
        self.title = title
        PLOT_TYPES = ["line", "bar"]
        if plot_type not in PLOT_TYPES:
            raise ValueError("'plot_type' %s not supported. Try one of these instead: %s" % (plot_type, ",".join(PLOT_TYPES)) )
        self.plot_type = plot_type  
        self.xlabel = xlabel
        self.ylabel = ylabel
        self.subplots = ODict()

    def add_subplot(self, label, xdata, ydata, yerror=None):
        self.subplots[label] = { "xdata": xdata, "ydata": ydata, "yerror": yerror }

    def export_grace(self):
        if self.plot_type == "line":
            typ = "xy"
        elif self.plot_type == "bar":
            typ = "bar"
        legends = self.subplots.keys()
# creates this:
# @s0 legend "rep_000"
# @s1 legend "rep_001" ...
        set_config = ""
        for i,sp in enumerate(self.subplots.keys()):
            set_config += "@s%d legend \"%s\" \n" % (i, sp)   # add legends
            if typ == "bar":
                set_config += "@s%d line type 0 \n" % (i,)   # don't show the line in bar plots

        sets = ""
        for label,sp in self.subplots.iteritems():
            if not sp["yerror"] or len(sp["yerror"]) != len(sp["xdata"]):
                yerror=[ "" for x in sp["xdata"] ]
            else:
                yerror=sp["yerror"]
                typ = typ + "dy"
            for x,y,dy in zip(sp["xdata"],sp["ydata"],yerror):
                sets += "%s %s %s\n" % (x,y,dy)
            sets += "&\n"
        
        return """#
@type {typ}
@title "{title}"
@xaxis label "{xlabel}"
@yaxis label "{ylabel}"
{set_config}
{sets}

""".format(typ=typ, title=self.title, xlabel=self.xlabel.encode("utf-8"), ylabel=self.ylabel.encode("utf-8"), set_config=set_config, sets=sets)


