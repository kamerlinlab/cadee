from setuptools import setup

import executables.exe as exe
from glob import glob


def installation_failed():
    """
    Abort the setup/installation process.
    """
    print('')
    print('CADEE installation has failed.')
    import sys
    sys.exit(1)


def q_missing(exe):
    print('ERROR: Could not find {0}. Please install Q5 and ensure the binaries are in $PATH.'.format(exe))
    print('Q5 can be obtained from {0}.'.format('http://www.icm.uu.se/cbbi/aqvist-lab/q'))
    installation_failed()

if not exe.which('qdyn5'):
    q_missing('qdyn5')

if not exe.which('qfep5'):
    q_missing('qfep5')

if not exe.which('qcalc5'):
    q_missing('qcalc5')

if not exe.which('Scwrl4'):
    print('ERROR: Scwrl4 is missing. Please install Scwrl4 and make sure the binary is in $PATH.')
    print()
    print('       Scwrl4 can be obtained from {0}.'.format('http://dunbrack.fccc.edu/scwrl4/'))
    installation_failed()

setup(name='cadee',
      version='0.7.1',
      description='Computer Aided Directed Evolution of Enzymes',
      url='http://github.com/kamerlinlab/cadee',
      author='Beat Anton Amrein',
      author_email='beat.amrein@gmail.com',
      license='GPLv2',
      packages=['analysis', 'mutate', 'ensemble', 'qscripts', 'executables', 'executables.q', 'executables.scwrl4'],
      py_modules=['cadee'],
      scripts = [ 'bin/cadee' ],
      #data_files=[('executables/q/qdyn5', 'qdyn5')],
      install_requires=[
          ['mpi4py==1.3.1'],
      ],
      console_scripts={
        'console_scripts': [
            'cadee = cadee:main',
            'ensemble = ensemble.ensemble:main',
        ]
      },
      zip_safe=False)
