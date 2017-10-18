#!/bin/bash

# DO NOT CHANGE
# AC  17 09 21

./cpuid.x $CORES_PER_EXE
core=$?
#repeat=6
taskindexsimple=$(( ( (node-1)*CORES_PER_NODE + core )  ))
taskindex=$(( ( (node-1)*CORES_PER_NODE + core )*repeat  ))


START_RUN=$taskindex
LAST_RUN=$(( taskindex + repeat  - 1))

for (( i=$START_RUN; i<=$LAST_RUN; i++ ))
do
   cd run-fssh-$i
   # Use the .sopt for FSSH
   echo -e "task $taskindexsimple running on core $core on node $node and doing run-fssh-$i"
   #/work/e05/e05/acarofmc/src/CP2K/nonadiabatic-cp2k/cp2k/exe/ARCHER/cp2k.sopt run.inp > run.log 
   $CP2K_PATH run.inp > run.log 
   cd ..
done

wait
echo -e "task $taskindex finishes on core $core on node $node"

#sleep 1800
sleep $sleeptimehere
#rm *restart*
#rm run-r-1.out
#rm run-mix-1.ener
#rm input-1.psf

#exit 0
