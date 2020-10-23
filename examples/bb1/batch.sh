#!/bin/bash
#
#SBATCH --job-name=potfit
#SBATCH --output=jobout.txt
#SBATCH --account=readmsd02
#
#SBATCH --ntasks 40
#SBATCH --nodes 1
#SBATCH --time 10:00
#SBATCH --mem 120GB
#
#SBATCH --get-user-env
#SBATCH --export=NONE
#
unset SLURM_EXPORT_ENV

#!/bin/bash
#
#SBATCH --job-name=isolated
#SBATCH --output=jobout.txt
#SBATCH --account=readmsd02
#
#SBATCH --ntasks 20
#SBATCH --nodes 1
#SBATCH --time 10:00
#SBATCH --mem 32GB
#
#SBATCH --get-user-env
#SBATCH --export=NONE
#
unset SLURM_EXPORT_ENV

module purge; module load bluebear
module load bear-apps/2018a
module load iomkl/2018a
module load Python/3.6.3-iomkl-2018a
module load matplotlib/2.1.1-iomkl-2018a-Python-3.6.3

# Change to $PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"

mkdir /scratch/bxp912
rm -R /scratch/${USER}/*

# Set the number of threads to 1
export OMP_NUM_THREADS=20
export PROC_COUNT=1


# Change to $PBS_O_WORKDIR
cd "$PBS_O_WORKDIR"

mkdir /scratch/bxp912
rm -R /scratch/${USER}/*

python3 /rds/homes/b/bxp912/apps/eampa_v3/eampa.py input.in





