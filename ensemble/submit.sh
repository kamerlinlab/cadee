#!/usr/bin/env bash
#SBATCH -A snic2015-16-12
#SBATCH -n 6
#SBATCH --time=0-01:00:00

# Demo Script on how to submit ensemble run, running on 6 cores.

cd ~/cadee/ensemble
source init.sh
cd -

# will hang on abisko
mpirun -n 6 python ~/cadee/ensemble/ensemble.py ~/inputs_tarchives/tarchives --alpha 229 --hij 60

# on abisko
srun -n 6 python ~/cadee/ensemble/ensemble.py ~/inputs_tarchives/tarchives --alpha 229 --hij 60
