CADEE 0.9
=========

CADEE: The framework of Computer-Aided Directed Enzyme Evolution.

# Installation

##System Requirements
- Multicore CPU
- Storage: 6-12 GB per mutant
- *nix like OS (GNU/Linux, MacOSX)
- Python Version 2.7
- Fortran Compiler, C Compiler, OpenMPI Compiler, MPI Launcher, git, Open Babel, pip
  ```  
  sudo apt-get install gfortran openmpi-bin git openbabel mpich gcc python2.7 python-pip 

  ```

##External Executables:
While CADEE itself is distributed as free software (GPLv2), it requires non-free (external) executables.
These are primarily the MD engine that CADEE is relying on, __Q6__, which is free and open source software. Additionally, CADEE requires for some functionality __SCWRL4__, which is also free for academic users.  

For mutations other than wild-type testing and alanine-scans, CADEE requires SCWRL4,  
which is also free for Non-Commerical usage.
- Q:
  https://github.com/qusers/Q6
  Download and compile Q. Optionally install to $PATH.
- SCWRL4: 
  Obtain a License for SCWRL4: (Free for Non-Commercial use)
    http://dunbrack.fccc.edu/scwrl4/license/index.html  
  Download and install, make sure executable is in $PATH.  
  
  Note: If you are not interested in the --libmut functionality of CADEE, you may place a dummy 
  "Scwrl4" executable/script into your $PATH. This will obviously break the --libmut option.
- Open Babel: Open Babel is a free program that can be installed with the os. The command in the 'System Requirements'
section will install Open Babel automatically.

  
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
cadee --help                                            # display help screen
```



## Workflow 
#### 1: Generate Inputs (simpacks)
1. Prerequisites:
    i)   PDB file generated with qprep5
    ii)  FEP file for/of reference/wild-type
    iii) qprep5-inputfile, you used to generate PDB file in i)
    iv)  path to the folder with the libraries you used to generate the inputfiles.

2. Then, check if CADEE accepts your inputfiles:
    $ cadee prep wt.pdb wt.fep qprep5.inp ./libraries 
    Probably you were asked to adjust your qprep.inp file, implement the changes.
    (CADEE needs absolute paths to your library, for example.)

3. Prepare and alanine scan for 16 mutants:
    $ cadee prep wt.pdb wt.fep qprep5.inp ./libraries --alascan --nummuts 16
    This created a folder (ala_scan) with subfolders containing the topology and fepfile: mutant.top, mutant.fep, qprep5.inp
    It will contain 16 subfolders, labelled in scheme 0XX_ALA, 1XX_ALA.

4. If you have SCWRL4 installed, you can also do arbitraty mutations, eg mutate residue 15 to glutamic:
    = requirements:SCWRL4 installed (see Installation Section)
    $ cadee prep wt.pdb wt.fep qprep5.ipn ./libraries --libmut 15:E
    This created a folder (libmut) with subfolders containing the topologies and fepfiles

   or even a saturation on 15 
   ```
   cadee prep wt.pdb wt.fep qprep5.ipn ./libraries --libmut 15:SATURATE
   ```

    the libmut argument is very powerful, here other options include:

         = SINGLE POINT MUTATION =
         --libmut 137:NEGATIVE (2AA)  
         --libmut 137:UNCHARGED (4AA)  
         --libmut 137:SHRINK (variable)  
         --libmut 137:'CGP' (3AA)  
         --libmut 137:'DEKCPG' (6AA)
           
         = COMBINATORIAL MULTI POINT MUTATION =  
         --libmut 137:'DEKCPG' 138:'DEKCPG' (6AAx6AA=36AA)  
         --libmut 137:'DEKCPG' 138:'DEKCPG' 139:SATURATE (6AAx6AAx20=720AA)  

#### 2: Submission of Jobs
   You might have to adjust the following instructions to fit the machine you are executing it.
   ```
   MPI_ENVIRONMENT -n N cadee dyn /path/to/simpacks
   ```
   MPI_ENVIRONMENT and N need to be adjusted to your system, and depend on your compiler and other settings. If you use a modern laptop Computer Ubuntu 16.04 with the recommended defaults, the mpiexec.mpich -n 3 might be the right choice.

#### 3: Analysis of Results
   You can then analyse the cadee.db to generate either *csv or *html files.
   ```
   cadee ana --help
   ```

## How to Cite this Work
The development of CADEE is mainly funded by academic research grants. To help 
us fund development, we humbly ask that you cite the CADEE paper:

* CADEE: Computer-Aided Directed Evolution of Enzymes
  Amrein BA, Steffen-Munsberg F, Szeler I, Purg M, Kulkarni Y & Kamerlin SCL. 
  IUCrJ 4, 50-64 (2017)
  DOI: https://doi.org/10.1107/S2052252516018017 
