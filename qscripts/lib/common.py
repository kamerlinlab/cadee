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
# Some common classes and functions

import math
import os
import shutil


__version__ = "0.1.10"


def backup_file( filename ):
    """  
    Checks if a file exists, makes a backup (#filename.1#, #filename.2#...).
    Returns the new basename as a string or empty string if the file was not found.

    Args:
        filename (string):  name of file to backup

    """
    if os.path.lexists( filename ):
        di = os.path.dirname( filename )
        fn = os.path.basename( filename )
        backup_filename = fn
        i = 1
        while os.path.lexists( os.path.join(di,backup_filename) ):
            backup_filename = "#%s.%d#" % (fn, i)
            i += 1
        shutil.copy2( filename, os.path.join(di,backup_filename) )
        return backup_filename
    return ""
    


# no need for numpy to do these basic stats
class np():
    @staticmethod
    def mean(vals):
        N = len(vals)
        if N == 0: 
            return float('nan')
        return sum(vals) * 1.0 / N

    @staticmethod
    def std(vals, ddof=1):
        N = len(vals)
        if N == 0 or N-ddof == 0: 
            return float('nan')
        mean = np.mean(vals)
        variance = map(lambda x: (x-mean)**2, vals)
        return math.sqrt( sum(variance)/(N-ddof) )

    @staticmethod
    def median(vals):
        N = len(vals)
        if N == 0:
            return float('nan')
        vals = sorted(vals)
        if N % 2 == 0: #even
            return np.mean( (vals[N/2-1], vals[N/2]) )
        else: #odd
            return vals[N/2]





class DataContainer(object):
    """
    Contains a two dimensional array of values:

    [ [ row1_column1, row1_column2, row1_column3, ...],
      [ row2_column1, row2_column2, row2_column3, ...],
      ...                                               ]

    and column titles.

    Args:
        coltitles (list): column titles

    Example of usage:
    >>> dg_de = DataContainer( ['Energy_gap', 'dG'] )
    >>> dg_de.add_row( [-300.0, 10.0 ]
    >>> rows = dg_de.get_rows( reversed(dg_de.get_column_titles()) )  # reversed rows
    >>> cols = dg_de.get_columns( columns=[0, 1] )
    """

    def __init__(self, coltitles):
        if not isinstance(coltitles, (list,tuple)): coltitles = [ coltitles, ]
        self._column_titles = list(coltitles)    
        self._rows = [] # a list containing rows of values (each row is a list with length = len(coltitles))
        self.comment = None


    def get_columns(self, columns=None):
        """
        Transposes the array and returns the columns instead of rows.

        Args:
            columns (list), optional: return only columns with these indices and/or titles

        Returns:
            list of columns (list of lists)
        """
        if not columns: columns = []
        col_inds = []
        for col in columns:
            if type(col) == int: 
                col_inds.append(col)
            else:
                col_inds.append(self._column_titles.index(str(col)))
        cols = zip(*self._rows)   # transpose
        if col_inds:
            return [ cols[i] for i in col_inds]
        else:
            return cols


    def get_rows(self, columns=None):
        """
        Returns the rows.

        Args:
            columns (list), optional: return only columns with these indices and/or titles

        Returns:
            list of rows (list of lists)
        """
        if columns:
            cols = self.get_columns(columns)
            return zip(*cols)
        else:
            return self._rows


    def get_column_titles(self):
        """
        Returns:
            list of column names (list)
        """
           
        return self._column_titles


    def add_row(self, row):
        """
        Args:
            row (list): a list of values

        Raises:
            ValueError: if number of elements in row is not equal to number of column titles
        """
        if len(row) != len(self._column_titles):
            raise ValueError("Number of elements is not equal to number of columns, in row:\n%s" % row)
        self._rows.append(list(row))


    def delete_rows(self):
        """
        Removes the rows.
        """
        self._rows = []

    def __str__(self):
        if self.comment:
            s = "#" + self.comment + "\n"
        else:
            s = ""
        for name in self._column_titles:     
            width = len(name)
            if width<10:
                width=10
            s += " {name:{width}} ".format(name=name, width=width)
        for row in self._rows:
            s += "\n"
            for i,val in enumerate(row): 
                try:
                    width=len(self._column_titles[i])
                    if width<10:
                        width=10
                except IndexError:
                    width=20
                if type(val) == float:    
                    s+=" {val:{width}.2f} ".format(val=val, width=width )
                else: 
                    s+=" {val:{width}} ".format(val=str(val), width=width)
        return s    

