#!/usr/bin/env bash
#SBATCH -A snic2015-16-12       # project name
#SBATCH -n 6                    # number of cores that shall be used in parallel
#SBATCH --time=0-01:00:00       # wallclock time limit in days - hours : minutes : seconds
#SBATCH --overcommit            # needed for efficient resource usage.

# This is a demo SLURM script on how to run ensemble.py, assuming CADEE was installed to $HOME/cadee/

# TODO: check path
cd ~/cadee/ensemble
source init.sh
cd -


# TODO: check paths
mpirun -n 7 python ~/cadee/ensemble/ensemble.py ~/inputs_tarchives/tarchives --alpha 229 --hij 60

# TODO: check paths
# alternative:
srun -n 7 python ~/cadee/ensemble/ensemble.py ~/inputs_tarchives/tarchives --alpha 229 --hij 60
