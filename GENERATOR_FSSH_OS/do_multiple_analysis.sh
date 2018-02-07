#!/bin/bash


list_run=("RUNS-TEMP-100" "RUNS-TEMP-200")

for dir in "${list_run[@]}";
do
    cp job_analyse_general.pbs ${dir}/
    cp -r scripts/ ${dir}/
    cd ${dir}
    pwd
    qsub job_analyse_general.pbs &
    cd ..
done
        

