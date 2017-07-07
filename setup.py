from __future__ import print_function

from setuptools import setup

import os
import sys
from glob import glob
import cadee.executables.exe as exe

QEXES=['qdyn5', 'qprep5', 'qcalc5', 'qfep5']
qexedir=os.path.join(os.path.dirname(os.path.abspath(__file__)), 'cadee/executables/q/')

def installation_failed():
    """
    Abort the setup/installation process.
    """
    print('')
    print('CADEE installation has failed.')
    import sys
    sys.exit(1)

def q_missing(exe):
    """
    Print message about the missing Q5-executable and call installation_failed.
    """
    print()
    print('ERROR: Could not find {0}. Please install Q5 and ensure the binaries are in $PATH.'.format(exe))
    print('       or copy the binaries to {0}.'.format(qexedir))
    print()
    print('       Q5 can be obtained from {0}.'.format('http://www.icm.uu.se/cbbi/aqvist-lab/q'))
    installation_failed()

if not (sys.version_info[0] == 2 and sys.version_info[1] == 7):
    print('Need Python 2.7')
    installation_failed()

if not exe.which('qdyn5', True):
    # There are many versions of executables named qdyn5.
    # CADEE should stick to one version, so include with the cadee installation.
    
    for qexe in QEXES:
        if not exe.which(qexe, True):
            
            print('Warning: Could not find {0} in {1}'.format(qexe, qexedir))
            print('         Searching in $PATH...')
            print()
            if exe.which(qexe):
                print('         Found {0} in {1}.'.format(qexe, exe.which(qexe)))
                print()
                print('         Will now copy {0} to {1}.'.format(qexe, qexedir))
                ans = raw_input('           Proceed (y/N)?').lower()
                if ans == 'y':
                    print('         Proceed.')
                    import shutil
                    shutil.copy2(exe.which(qexe), qexedir) 
                    print('cped {0} to {1}.'.format(exe.which(qexe), qexedir)) 
                else:
                    print('         Abort.')
                    print('         Fatal: Cannot continue installation without {0}'.format(qexe))
                    q_missing(qexe)
            else:
                print('Fatal: Could not find {0} in $PATH.'.format(qexe))
    		q_missing(qexe)

if not exe.which('Scwrl4', True):
    # It's not recommended.
    print('WARNING: Could not find Scwrl4 in cadee/executables/scwrl4/Scwrl4')


if not exe.which('babel'):
    print('ERROR: Could not find babel. Please install openbabel and ensure the binaries are in $PATH.')
    print('       Using Ubuntu try: sudo apt-get install openbabel')
    installation_failed()

if not exe.which('Scwrl4'):
    print('ERROR: Scwrl4 is missing. Please install Scwrl4 and make sure the binary is in $PATH.')
    print()
    print('       Scwrl4 can be obtained from {0} (free for non-commercial use).'.format('http://dunbrack.fccc.edu/scwrl4/'))
    installation_failed()


setup(name='cadee',
      version='0.8.0',
      description='Computer Aided Directed Evolution of Enzymes',
      url='http://github.com/kamerlinlab/cadee',
      author='Beat Anton Amrein',
      author_email='beat.amrein@gmail.com',
      license='GPLv2',
      packages=['cadee', 'cadee.ana', 'cadee.dyn', 'cadee.prep', 'cadee.executables', 'cadee.qscripts' ],
      py_modules=['cadee'],
      package_data={'cadee': ['lib/*', 'qscripts/lib/*', 'qscripts/REAMDE.md', 'qscripts/LICENSE.txt', 'executables/q/q*', 'executables/scwrl4/*', 'bash_scripts/*']},
      install_requires=[
          ['mpi4py==1.3.1'],
      ],
      scripts={
            'cadee/cadee',
      },
      zip_safe=False)
