CADEE
=====

CADEE the framework of Computer-Aided Directed Enzyme Evolution.

# Installation

##Pre-Requisites:
- Python Version 2.7 
- mpi4py 1.3.1 (python module, for ensemble simulation)
  Download: https://pypi.python.org/pypi/mpi4py/1.3.1
  Load a compiler (eg intel):
    module add intel/16.0.3;
  and then install mpi4py, eg. 
    tar xf mpi4py-1.3.1.tar.gz; 
    cd mpi4py-1.3.1; 
    python setup.py install --user

##Required Executlables:
- openbabel 2.3 or newer in $PATH (for input preparation)
- Q:
  Obtain a License for Q: (Free for Non-Commercial use)
    http://www.icm.uu.se/cbbi/aqvist-lab/q/
  compile with ifort and copy the executables to:
    ./mutate/executables/scwrl4
- SCWRL4: 
  Obtain a License for SCWRL4 (Free for Non-Commercial use)
    http://dunbrack.fccc.edu/scwrl4/license/index.html
  Download and install the package to:
    ./mutate/executables/q

##Test Configuration:
./setup_check.sh
and for ensemble simulations:
./ensemble/init.sh


## Workflow 
# 1: Generate Inputs (simpacks)
= requirements: (open)babel in $PATH
1. Prerequisites:
    i)   pdbfile generated with qprep5
    ii)  fepfile for qprep5
    iii) qprep5-inputfile, you used to generate pdbfile in i)
    iv)  path to the folder with the libraries you used to generate the inputfiles.

2. Then, check if CADEE accepts your inputfiles:
    $ ./cadee.py wt.pdb wt.fep qprep5.inp ./libraries 
    Probably you were asked to adjust your qprep.inp file, implement the changes.
    (CADEE needs absolute paths to your library, for example.)

3. Prepare an alanine scan for 16 mutants:
    $ ./cadee.py wt.pdb wt.fep qprep5.inp ./libraries --alascan --nummuts 16
    This created a folder (ala_scan) with 16 subfolders (labelled 0XX_ALA, 1XX_ALA, ...) containing the topology and fepfile: mutant.top, mutant.fep, qprep5.inp

4. If you have SCWRL4 installed, you can also do arbitrary mutations, eg mutate residue 15 to glutamic acid:
    = requirements:SCWRL4 installed (see Installation Section)
    $ ./cadee.py wt.pdb wt.fep qprep5.ipn ./libraries --libmut 15:E
    This created a folder (libmut) with subfolders containing the topologies and fepfiles

   or even a saturation on position 15 
    $ ./cadee.py wt.pdb wt.fep qprep5.ipn ./libraries --libmut 15:SATURATE

    the libmut argument is very powerful, other options include:
         --libmut 137:SATURATE (20AA)
         --libmut 137:POLAR (9AA)
         --libmut 137:APOLAR (8AA)
         --libmut 137:POSITIVE (3AA)
         --libmut 137:NEGATIVE (2AA)
         --libmut 137:UNCHARGED (4AA)
         --libmut 137:SHRINK (variable)
         --libmut 137:'CGP' (3AA)
         --libmut 137:'DEKCPG' (6AA)
         -OR-
         --libmut 137:'DEKCPG' 138:'DEKCPG' (6AAx6AA=36AA)
         --libmut 137:'DEKCPG' 138:'DEKCPG' 139:SATURATE (6AAx6AAx20=720AA)

# 2: Submission of Jobs
(You might have to adjust the following instructions to fit the machine you are executing it.)
    Next you have to actually run the generated simpacks.
    see also ensemble/submit.sh for an example script

# 3: Analysis of results
    When your data is written to cadee.db, for analysis you can either generate *csv or *html files from it.
    For *.csv:
        analyse/extract_to_csv_medium.py
        analyse/extract_to_csv_us.py  
    For *.html
        analyse/analyse.py 
        
    Other tools:
        analyse/fuse_dbs.py

