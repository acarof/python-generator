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
/work/e358/e358/acarof/src/nonadiabatic-cp2k/cp2k/exe/ARCHER/cp2k.sopt run.inp > run.log 2> run.err

exit 0
