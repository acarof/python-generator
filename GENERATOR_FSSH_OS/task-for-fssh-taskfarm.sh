#!/bin/bash

# DO NOT CHANGE
# AC  17 09 21

./cpuid.x $CORES_PER_EXE
core=$?
taskindex=$(((node-1)*CORES_PER_NODE + core))


# Uncomment to debug
echo -e "task $taskindex running on core $core on node $node"



DIR=run-fssh-${taskindex}
cd $DIR

# Use the .sopt for FSSH
#CP2K_PATH=/home/e358/e358/acarof/src/CP2K/nonadiabatic/cp2k/exe/ARCHER/cp2k.sopt
CP2K_PATH=/home/e358/e358/acarof/src/CP2K/nonadiabatic/cp2k/exe/ARCHER/
ls $CP2K_PATH

#${CP2K_PATH} run.inp > run.log

exit 0
