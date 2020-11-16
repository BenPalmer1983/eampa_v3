#!/bin/bash
#
#SBATCH --job-name=potfit
#SBATCH --output=jobout_fit3.txt
#SBATCH --account=readmsd02
#
#SBATCH --ntasks 20
#SBATCH --nodes 1
#SBATCH --time 2000:00
#SBATCH --mem 32GB
#
#SBATCH --get-user-env
#SBATCH --export=NONE
#
unset SLURM_EXPORT_ENV

module purge; module load bluebear
module load bear-apps/2018a
module load imkl 2019.5.281-gompi-2019b
module load Python 3.7.4-GCCcore-8.3.0
module load matplotlib 3.1.1-foss-2019b-Python-3.7.4


# Change to $PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"

mkdir /scratch/bxp912
rm -R /scratch/${USER}/*

# Set the number of threads
export OMP_NUM_THREADS=20
export PROC_COUNT=1
export PYTHONPATH=$PYTHONPATH:"/rds/homes/b/bxp912/apps/f2py_lib"

pwd > pwd.txt
ls > ls.txt
python /rds/homes/b/bxp912/apps/python/eampa.py input.in
#python eampa.py input.in
