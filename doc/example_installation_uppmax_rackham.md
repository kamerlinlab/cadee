# Installing CADEE on the Rackham Cluster

The intel, mpi and openbabel modules are loaded:

```
module add intel/17.4 intelmpi/17.4 openbabel
```
Next Q is cloned and compiled and the executables added to $PATH
```
git clone https://github.com/qusers/qsource.git .
git checkout development/beat  # fix for 'make all'
cd qsource/src/
make all COMP=ifort
cd ../bin
export PATH="$PWD:$PATH"
```

Next, install SCWRL4 to $HOME/.local/bin/ and and check that Scwrl4 is in $PATH:
```
which Scwrl4
> ~/.local/bin/Scwrl4
```

Only if the Scwrl4 executable exists, CADEE is can be installed:
```
cd $HOME
cd Downloads
git clone https://github.com/kamerlinlab/cadee cadee
cd cadee
python setup.py install --user
```

Let's check that cadee is in $PATH:
```
which cadee
> ~/.local/bin/cadee
```

**Important**: You must load the intelmpi library to launch cadee.

```
module add intelmpi/17.4
cadee --help

or, when on a node:

srun -n 1 $(which cadee) --help
srun -n 1 $HOME/.local/bin/cadee --help
```

If you now see the following output, you have successfully installed CADEE:

```
Usage:
                   cadee [ prep(p) | dyn(d) | ana(a) | tool(t) ]

       Multi Core Tasks:
                    mpirun -n X cadee dyn
                   mpiexec -n X cadee dyn
                   X == Number of cores to use; 2+.
```