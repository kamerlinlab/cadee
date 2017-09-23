CADEE 0.9
=========

CADEE the framework of Computer-Aided Directed Enzyme Evolution.

# Installation

##System Requirements
- Multicore CPU
- Storage: 6-12 GB per mutant
- *nix like OS (GNU/Linux, MacOSX)
- Python Version 2.7
- Fortran Compiler, C Compiler, OpenMPI Compiler, MPI Launcher, git, Open Babel, pip
  ```  
  sudo apt-get install gfortran gcc openmpi-bin mpich git openbabel python-pip
  ```

##External Executables:
While CADEE itself is free software, it requires non-free components.  
The MD engine that CADEE is relying on is Q5 which is free for academic users.  

For mutations other than wild-type testing and alanine-scans, CADEE requires SCWRL4,  
which is also free for Non-Commerical usage.
- Q:
  Obtain a License for Q: (Free for Non-Commercial use)
    http://www.icm.uu.se/cbbi/aqvist-lab/q/  
  Download and compile Q. Optionally install to $PATH.
- SCWRL4: 
  Obtain a License for SCWRL4: (Free for Non-Commercial use)
    http://dunbrack.fccc.edu/scwrl4/license/index.html  
  Download and install, make sure executable is in $PATH.  
  
  Note: If you are not interested in the --libmut functionality of CADEE, you may place a dummy 
  "Scwrl4" executable/script into your $PATH. This will obviously break the --libmut option.
  
##Download and Install CADEE:
We recommend using git to install CADEE, so that future releases of CADEE are easily accessible.  

Here the suggested steps:  

Download CADEE:
```
cd $HOME                
mkdir -p Downloads                                      # create a folder
cd Downloads  
git clone https://github.com/kamerlinlab/cadee cadee    # clone CADEE
```

Next, install CADEE:  
```
cd $HOME/Downloads/cadee                                # cd into cloned folder
python setup.py install --user                          # install for current user *ONLY*
```



## Workflow 
# 1: Generate Inputs (simpacks)
1. Prerequisites:
    i)   PDB file generated with qprep5
    ii)  FEP file for/of reference/wild-type
    iii) qprep5-inputfile, you used to generate PDB file in i)
    iv)  path to the folder with the libraries you used to generate the inputfiles.

2. Then, check if CADEE accepts your inputfiles:
    $ ./cadee.py wt.pdb wt.fep qprep5.inp ./libraries 
    Probably you were asked to adjust your qprep.inp file, implement the changes.
    (CADEE needs absolute paths to your library, for example.)

3. Prepare and alanine scan for 16 mutants:
    $ ./cadee.py wt.pdb wt.fep qprep5.inp ./libraries --alascan --nummuts 16
    This created a folder (ala_scan) with subfolders containing the topology and fepfile: mutant.top, mutant.fep, qprep5.inp
    It will contain 16 subfolders, labelled in scheme 0XX_ALA, 1XX_ALA.

4. If you have SCWRL4 installed, you can also do arbitraty mutations, eg mutate residue 15 to glutamic:
    = requirements:SCWRL4 installed (see Installation Section)
    $ ./cadee.py wt.pdb wt.fep qprep5.ipn ./libraries --libmut 15:E
    This created a folder (libmut) with subfolders containing the topologies and fepfiles

   or even a saturation on 15 
    $ ./cadee.py wt.pdb wt.fep qprep5.ipn ./libraries --libmut 15:SATURATE

    the libmut argument is very powerful, here other options include:
         --libmut 137:SATURATE (20AA)  
         --libmut 137:POLAR (9AA)  
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
    When your data is, you can then analyse the cadee.db to generate either *csv or *html files.
    For *.csv:
        analyse/extract_to_csv_medium.py
        analyse/extract_to_csv_us.py  
    For *.html
        analyse/analyse.py 
        
    Other tools:
        analyse/fuse_dbs.py

# Troubleshooting

##Manually installing mpi4py 1.3.1
- mpi4py 1.3.1 (CADEE will try to install mpi4py)  
  Download: https://pypi.python.org/pypi/mpi4py/1.3.1
    ```
    wget -O mpi4py-1.3.1.tar.gz https://pypi.python.org/packages/26/b4/1a9678ec113b5c654477354169131c88be3f65e665d7de7c5ef306f2f2a5/mpi4py-1.3.1.tar.gz
    ```
  If your machine supports special compilers, load one (e.g. intel):  
    ```
    module add intel/16.0.3
    ```  
  and then install mpi4py, eg. 
    ```
    tar xf mpi4py-1.3.1.tar.gz
    cd mpi4py-1.3.1
    tar xf python setup.py install --user
    ```  
