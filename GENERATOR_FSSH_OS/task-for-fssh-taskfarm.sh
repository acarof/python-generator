#!/bin/bash

# DO NOT CHANGE
# AC  17 09 21

./cpuid.x $CORES_PER_EXE
core=$?
taskindex=$(((node-1)*CORES_PER_NODE + core))

# Uncomment to debug
echo -e "task $taskindex running on core $core on node $node"

cd run-fssh-${taskindex}

# Use the .sopt for FSSH
/work/e358/e358/acarof/src/nonadiabatic-cp2k/cp2k/exe/ARCHER/cp2k.sopt run.inp > run.log &

wait
echo -e "task $taskindex finishes on core $core on node $node"

sleep 1800
#rm *restart*
#rm run-r-1.out
#rm run-mix-1.ener
#rm input-1.psf

#exit 0
