CADEE 0.9.1
===========

the framework for **C**omputer-**A**ided **D**irected **E**nzyme **E**volution.

#### Description
CADEE is a framework for computational directed evolution of enzymes. The methodology is described _en detail_ in the our publication [CADEE: Computer-Aided Directed Evolution of Enzymes](https://doi.org/10.1107/S2052252516018017). CADEE is particularly powerful for computational modeling of SN2-like reactions in an MM/EVB framework, but the approach can be used to model any enzymatic reactions.



## System Requirements
- *nix like OS (we use GNU/Linux)
- CPU: 1 or more core-weeks per mutant
- Storage: 6-12 GB per mutant
- Software: `sudo apt-get install gfortran openmpi-bin git  openbabel mpich gcc python2.7 python-pip`
 - Python Version 2.7
 - Fortran Compiler (e.g. gfortran)
 - C Compiler (e.g. gcc)
 - OpenMPI Compiler & Launcher
 - Open Babel
 - git (optional)
 - pip (optional)
 - [SLURM](https://slurm.schedmd.com/) (optional, scheduler for super computers)*

\* CADEE was written specifically to allow efficient usage of a cluster computer. It is not required to run CADEE on a cluster and for example preparation and analysis are best performed on a laptop computer. However, because enzymatic reactions have many degrees of freedom, the simulation of enzymatic reactions require considerable conformational sampling to allow for meaningful results. Therefore, we developed and tested CADEE on cluster computers which were coincidentally running SLURM. CADEE should be compatible with other schedulers but this remains untested to this day. 

# Installation

### External Executables:
Both, CADEE and the the MD engine that CADEE is relying on [**Q6**](https://doi.org/10.1016/j.softx.2017.12.001) are free and opensource software licensed under the GPL v2. However, CADEE requires for some functionality __SCWRL4__, which is proprietary altough licensed free to academic users (research purposes).

For mutations other than wild-type testing and alanine-scans, CADEE requires SCWRL4.
- Q:
  https://github.com/qusers/Q6
  Download and compile Q. Optionally install to $PATH.
- SCWRL4: 
  Obtain a License for SCWRL4: (Free for Non-Commercial use)
    http://dunbrack.fccc.edu/scwrl4/license/index.html
  Download and install, make sure executable is in $PATH.
 
  Note: If you are not interested in the --libmut functionality of CADEE, you may place a dummy 
  "Scwrl4" executable/script into your $PATH. This will obviously break the --libmut option.
- Open Babel: Open Babel is a free program that can be installed with the os. The command in the 'System Requirements' section will install Open Babel automatically.

  
### Download and Install CADEE:
We recommend using git to install CADEE, so that future releases of CADEE are easily accessible.

Here the suggested steps:



Clone into CADEE:
 
`git clone https://github.com/kamerlinlab/cadee cadee`
```

Install CADEE:

`python setup.py install --user`



## Workflow 
   See the [workflow](./workflow.md) document. 

## How to Cite this Work
The development of CADEE is mainly funded by academic research grants. To help 
us fund development, we humbly ask that you cite the CADEE paper:

* CADEE: Computer-Aided Directed Evolution of Enzymes
  Amrein BA, Steffen-Munsberg F, Szeler I, Purg M, Kulkarni Y & Kamerlin SCL. 
  IUCrJ 4, 50-64 (2017)
  DOI: [https://doi.org/10.1107/S2052252516018017](https://doi.org/10.1107/S2052252516018017) 

