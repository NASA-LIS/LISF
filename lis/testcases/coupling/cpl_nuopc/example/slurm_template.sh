#!/bin/bash
#SBATCH --job-name=lis-nuopc-test
#SBATCH --ntasks=4
#SBATCH --time=00:10:00
#SBATCH --partition=<partition_names>
#SBATCH --account=<account>
#SBATCH --constraint=<list>
#SBATCH --qos=<qos>
#SBATCH --output=test.log
#SBATCH --error=test.log

# environment settings
#   TODO

# set umask
umask 022
# set limits
ulimit -t unlimited
ulimit -f unlimited
ulimit -d unlimited
ulimit -s unlimited
ulimit -c unlimited
# set core dumps
decfort_dump_flag=y

# PRE
echo "________________________________________________________________"
echo "Job ID: ${SLURM_JOB_ID}"
s_tm=`date +%s`
s_hr=`date +%H`; s_mn=`date +%M`; s_sc=`date +%S`
echo "Model Start    ${s_hr}:${s_mn}:${s_sc}"
echo "________________________________________________________________"

# RUN
mpirun -np 4 ./lis_nuopc_test.exe
exec_s=$?

# POST
echo "________________________________________________________________"
e_tm=`date +%s`
e_hr=`date +%H`; e_mn=`date +%M`; e_sc=`date +%S`
echo "Model End      ${e_hr}:${e_mn}:${e_sc}"
r_tm=$(( ${e_tm} - ${s_tm} ))
r_hr=$(( ${r_tm} / 3600 ))
r_mn=$(( ${r_tm} % 3600 / 60 ))
r_sc=$(( ${r_tm} % 60 ))
echo "Model Runtime  ${r_hr}:${r_mn}:${r_sc}"
if [ $exec_s -ne 0 ]; then
  echo "RESULT: ERROR ${exec_s}"
else
  echo "RESULT: SUCCESS"
fi
echo "________________________________________________________________"
