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

import sys
import os
import ConfigParser


CFG_VERSION=6   # this should change when the config file format changes or something, then old config files will not work
# absolute config path
_CFG_FILE=os.path.abspath(os.path.join(os.path.dirname(__file__), "qscripts.cfg"))
_CFG_FILE_DEFAULT=os.path.abspath(os.path.join(os.path.dirname(__file__), "lib", "qscripts.cfg.default"))



class _QScriptsConfig(object):
    def __init__(self, cfgfile):

        self.cfgfile = cfgfile
        if not os.path.lexists(self.cfgfile):
            print "Configuration file '%s' not found. Please run 'qscripts_config.py'." % self.cfgfile
            sys.exit(1)

        self.config=ConfigParser.SafeConfigParser()

        try:
            self.config.read(self.cfgfile)
        except ConfigParser.ParsingError:
            print "Configuration file '%s' could not be read. Fix it or remove and run 'qscripts_config.py'." % self.cfgfile
            sys.exit(1)

        # check version of config file
        version = self.config.getint("other", "cfgversion")
        if version != CFG_VERSION:
            print "Your configuration file '%s' is outdated. Please remove it and run 'qscripts_config.py'." % self.cfgfile
            sys.exit(1)

    def get(self, section, option):
        try:
            return self.config.get(section, option)
        except (ConfigParser.NoSectionError, ConfigParser.NoOptionError) as e:
            print "Somehow your configuration file '%s' got messed up. Please fix it or remove it and rerun 'qscripts_config.py'" % self.cfgfile
            print "Details: %s" % str(e)
            sys.exit(1)

    def set(self, section, option, value):
        self.config.set(section, option, value)



def get_exec_path(name):
    #paths=[]
    #for path in os.environ["PATH"].split(os.pathsep):                       
    #    path=path.strip('"')                                                
    #    if not os.path.lexists(path):
    #        continue
    #    for ex in os.listdir(path):                                         
    #         if name in ex:                                              
    #              ex = os.path.join(path,ex)                                  
    #              if not ex in paths:
    #                  paths.append(ex)
    #
    #if not paths:
    #    print "No '%s' executable was found in your PATH. Please add the path to qfep to the config file manually." % name
    #else:
    #    print "These '%s' executables were found in your PATH. Choose the correct one or write the path.\n" % name
    #    for i,path in enumerate(paths):
    #        print "      [%d] %s" % (i, path)
    #    path=""
    #    while not path:
    #        try:
    #            inp=raw_input("? ")
    #            i=int(inp)
    #            path=paths[i]
    #        except ValueError:
    #            if inp:
    #                path=inp
    #        except IndexError:
    #            pass
    #    return path
    #return ""
    import cadee.executables.exe as exe
    path = exe.which(name)
    if path is None:
        path = exe.which(name+'5')
    if path is None:
        print('Binary not found: ', name)
        raise (Exception, 'Please install %s ', name)
    print(path)
    return path


def main():
    if os.path.lexists(_CFG_FILE):
        print "Configuration file '%s' exists. Please remove it before running this script." % _CFG_FILE
        sys.exit(1)

    QScriptsConfig = _QScriptsConfig(_CFG_FILE_DEFAULT)
    print "Creating a new configuration file...\n"

    print "Qfep executable:"
    qfep_path=get_exec_path("qfep")
    QScriptsConfig.set("qexec", "qfep", qfep_path)
    print "Qcalc executable (version which supports group contribution calculations):"
    qcalc_path=get_exec_path("qcalc")
    QScriptsConfig.set("qexec", "qcalc", qcalc_path)

    QScriptsConfig.config.write(open(_CFG_FILE, "w+"))
    print "\n\nThe following config file was created with some default values:\n%s" % _CFG_FILE


if __name__  == "__main__" or not os.path.isfile(_CFG_FILE):
    try:
        main()
    except KeyboardInterrupt:
        print "\nCtrl-C detected. Quitting..."
        sys.exit(1)
 
QScriptsConfig = _QScriptsConfig(_CFG_FILE)

