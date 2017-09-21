#!/bin/bash

./cpuid.x $CORES_PER_EXE
core=$?
taskindex=$(((node-1)*CORES_PER_NODE + core))

input_file=some_name_based_on_taskindex
input_parameters=some_function_of_taskindex

# Uncomment to debug
echo -e "task $taskindex running on core $core on node $node"

# Uncomment to run your executable
#my_executable $input_file $input_parameters > ${taskindex}.out 2> ${taskindex}.err &

exit 0
