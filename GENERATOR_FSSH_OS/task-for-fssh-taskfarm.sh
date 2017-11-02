#!/bin/bash

# DO NOT CHANGE
# AC  17 09 21

./cpuid.x $CORES_PER_EXE
core=$?
taskindex=$(( ( (node-1)*CORES_PER_NODE + core )  ))

cd run-fssh-$taskindex
# Use the .sopt for FSSH
echo -e "task $taskindex running on core $core on node $node and doing run-fssh-$taskindex"
$CP2K_PATH run.inp > run.log 
cd ..

echo -e "task $taskindex finishes on core $core on node $node"

echo "task $taskindex has completed" >> $LOGFILE
tasksfinished=`wc -l < $LOGFILE`
while [ $tasksfinished -lt $CORES_PER_JOB ]; do
  sleep 10
  tasksfinished=`wc -l < $LOGFILE`
done

exit 0


