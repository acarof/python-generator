#!/bin/bash --login

#PBS -N pbs-analyse-gen
#PBS -l select=serial=true:ncpus=1
#PBS -l walltime=03:00:00
#PBS -A e05-softm-blu

#module load anaconda
#source activate forgenerator

# Make sure any symbolic links are resolved to absolute path
export PBS_O_WORKDIR=$(readlink -f $PBS_O_WORKDIR)

# Change to the directory that the job was submitted from
cd $PBS_O_WORKDIR

python scripts/general.py 0