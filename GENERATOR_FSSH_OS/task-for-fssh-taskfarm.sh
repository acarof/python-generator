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
../cp2k.sopt run.inp > run.log

exit 0
