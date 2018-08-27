CADEE 0.9.1: Workflow
=====================

#### Description
CADEE's workslow is imitating directed evolution. Before one can start with directed evolution iterations, a parametrized EVB reaction has to be obtained. Once a the reference reaction has been established, one iteration after the other can be applied. Usually, this can be kicked off with an alanine scan (1, 2), followed up with analysis (3). The analysis will aid the user to select hotspots and start another round in-silico mutagenesis (4).


## Workflow 
#### Prerequisites:
#####   i) PDB file generated with qprep
#####  ii) FEP file for wild-type reaction
##### iii) Qprep-inputfile used to generate PDB file (see i)
#####  iv) Path to the folder with the libraries used to generate the inputfiles.
#####   v) Adjusting inputfiles:
Next, the wildtype-reaction can be prepared:
    `cadee prep wt.pdb wt.fep qprep6.inp ./libraries` CADEE will provide information on how adjust your qprep-input file. These changes have to be implemented to continue; CADEE needs absolute paths to your library, for example.

### 1. Prepare and alanine scan for 16 mutants:
   ` cadee prep wt.pdb wt.fep qprep6.inp ./libraries --alascan --nummuts 16`
    This command created a folder (ala_scan) with subfolders containing the topology and fepfile: mutant.top, mutant.fep, qprep6.inp
    It will contain 16 subfolders, labelled in scheme 0XX_ALA, 1XX_ALA.


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

#### 4: Next iteration: If you have SCWRL4 installed, you can also do arbitraty mutations, eg mutate residue 15 to glutamic:
    = requirements:SCWRL4 installed (see Installation Section)
    `cadee prep wt.pdb wt.fep qprep5.ipn ./libraries --libmut 15:E`
    This created a folder (libmut) with subfolders containing the topologies and fepfiles

   or even a saturation on 15 
   `cadee prep wt.pdb wt.fep qprep6.inp ./libraries --libmut 15:SATURATE`

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

#### Return to [readme.md](./readme.md) 
