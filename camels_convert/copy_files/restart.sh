#!/bin/bash -l
#########################################################
#SBATCH -J CAMELS_SWIFT
#SBATCH -p cmbas
#SBATCH --mail-user=sagan.sutherland@uconn.edu
#SBATCH --mail-type=ALL
#SBATCH --constraint="skylake"
#SBATCH -o swift.log
#########################################################
#SBATCH --time=7-0
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#########################################################

module purge
module load modules
module load gcc/11.4.0
module load gsl/2.7
module load openmpi/4.0.7
module load hdf5/mpi-1.8.22
module load fftw/mpi-3.3.10

./swift --pin --cosmology --simba --threads=${SLURM_CPUS_PER_TASK} --restart params.yml
mkdir snaps
mv snapshot*.hdf5 snaps/

#########################################################
