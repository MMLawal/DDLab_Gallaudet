#!/bin/bash

args=("$@")

#check the correct num of args
if [[ ${#args[@]} == 0 ]] || [[ ${#args[@]} < 3 ]]; then
  echo "create_ASMDjobfile.sh {NUM of trajectories} {Coord/RST7} {Stage Num}"
  exit
fi

num_asmd_sim=${args[0]}

if [ ! -f ${args[1]} ]; then
  echo The file ${args[1]} does not exist
else
  coord=${args[1]}
  
fi

stage=${args[2]}

cat >>_job.sh<<EOF
#!/bin/bash
#FILENAME: jobsubmit
#SBATCH -A bio230137-gpu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --gpus-per-node=2
#SBATCH --time=48:00:00
#SBATCH -J ASMD_stage$stage
#SBATCH -o myjob.o%j
#SBATCH -e myjob.e%j
#SBATCH -p gpu
#SBATCH --mail-user=monsurat.lawal@gallaudet.edu
#SBATCH --mail-type=all

#Manage processing environment, load compilers and applications
module purge
module load modtree/gpu
module load gcc/8.4.1  openmpi/4.0.6-cu11.0.3 amber/20
module list

do_parallel="mpirun -np 2 pmemd.cuda.MPI"
prmtop="prot1_plx_1_solv.prmtop"
coord="$coord"

EOF

for ((counter=1;counter<=$num_asmd_sim;counter+=1)); do
  echo SSSdo_parallel -O -i ASMD_$counter/asmd_$counter.$stage.mdin -p SSSprmtop -c SSScoord -r ASMD_$counter/ASMD_$counter.$stage.rst7 -o ASMD_$counter/ASMD_$counter.$stage.mdout -x ASMD_$counter/ASMD_$counter.$stage.nc -inf ASMD_$counter/ASMD_$counter.$stage.info >> _job.sh
done
sed 's/SSS/$/g' _job.sh > job.$stage.sh; rm -f _job.sh
    
