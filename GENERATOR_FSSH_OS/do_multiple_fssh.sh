#!/bin/bash


list_run=("RUNS-TEMP-100" "RUNS-TEMP-200")

for dir in "${list_run[@]}";
do
    cp job_cp2k_taskfarm_fssh.pbs ${dir}/
    cp cpuid.x ${dir}/
    cp task-for-fssh-taskfarm.sh ${dir}/
    cp -r topologies ${dir}/
    cd ${dir}
    pwd
    qsub job_cp2k_taskfarm_fssh.pbs &
    cd ..
done
        

