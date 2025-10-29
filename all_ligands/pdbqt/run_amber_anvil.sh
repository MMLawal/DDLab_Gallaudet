#!/bin/bash
#FILENAME: jobsubmit
#SBATCH -A bio230137-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --gpus-per-node=1
#SBATCH --time=48:00:00
#SBATCH -J plk1_lig1_bs1
#SBATCH -o myjob.o%j
#SBATCH -e myjob.e%j
#SBATCH -p gpu
#SBATCH --mail-user=syouremail
#SBATCH --mail-type=all

#Manage processing environment, load compilers and applications
module purge
module load modtree/gpu
module load gcc/8.4.1  openmpi/4.0.6-cu11.0.3 amber/20
module list

#Launch code
pmemd.cuda -O -i step4.0_minimization.mdin -o step4.0_minimization.log -c step3_input.rst7 -p step3_input.parm7 -r step4.0_input.rst7 -ref step3_input.rst7
pmemd.cuda -O -i step4.1_equilibration.mdin -o step4.1_equilibration.log -c step4.0_input.rst7 -p step3_input.parm7 -r step4.1_equilibration.rst7 -ref step4.0_input.rst7 -x step4.1_equilibration.nc
mpirun -np 2 pmemd.cuda.MPI -O -i prod.in -o step4_production.log -c step4.1_equilibration.rst7 -p step3_input.parm7 -r step4_production.rst7 -x step4_production.nc
