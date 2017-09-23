# Installing CADEE on the Abisko Cluster

The intel, mpi and Open Babel modules are loaded:

```
module add intel/2017b iimpi/2017b
```

On HPC2N, Open Babel is not installed. We therefor need download and compile it:
```
cd ~/pfs
wget  https://sourceforge.net/projects/openbabel/files/openbabel/2.4.1/openbabel-2.4.1.tar.gz
tar xf openbabel-2.4.1.tar.gz 
cd openbabel-2.4.1/
mkdir build
cd build
cmake ../
make
cd bin
```
Next the Open Babel binaries are added to ~/.bashrc
```
echo export PATH="$PWD:\$PATH" >> ~/.bashrc
source ~/.bashrc
```

Next Q is cloned, compiled and the executables added to $PATH:
```
git clone https://github.com/qusers/qsource.git .
cd qsource
git checkout development/beat  # fix for 'make all'
cd src
make all COMP=ifort
cd ../bin
export PATH="$PWD:$PATH"
```

This cluster does not allow accessing files in $HOME when running a job. 
We hence must adjust the '~/.local folder':

```
cd ~
mv ~/.local ~/pfs/.local
ln -s  ~/pfs/.local .local
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
module add iimpi/2017b
cadee --help
mpiexec -n 1 $(which cadee) --help

or, when on a node:

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