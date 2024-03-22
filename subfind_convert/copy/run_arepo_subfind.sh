#!/bin/bash -l
#SBATCH -J CAMELS_EAGLE_CV_1_subfind
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=28
#SBATCH -o logs/job.%J.dump
#SBATCH -e logs/job.%J.err
#SBATCH -p cosma7-rp
#SBATCH -A dp276
#SBATCH --exclusive
#SBATCH --mail-type=ALL
#SBATCH --mail-user=christopher.lovell@port.ac.uk
#SBATCH --no-requeue
#SBATCH -t 05:00:00

module purge
module load intel_comp/2022.3.0 compiler mpi
module load fftw/3.3.10
module load hdf5/1.14.0
module load gsl
module load hwloc

# mpiexec -n 2 
/cosma7/data/dp004/dc-love2/codes/Arepo_subfind/Arepo /cosma7/data/dp004/dc-love2/codes/CAMELS-SWIFT/arepo_subfind_param.txt 3 90
