#!/bin/bash

source $augeanstables/src/config/runtest_header.sh

#test_dir
test_dir=$augeanstables/src/test_files/test_bf_layer


#============================================================
#main body
#============================================================
echo ''


#test_mpi_tag
AUGEANSTABLES_PARALLEL=false
file='test_mpi_tag'
perform_test $file


#test_mpi_process
AUGEANSTABLES_PARALLEL=true
file='test_mpi_process'
change_param_input 'npx' '1'
change_param_input 'npy' '2'
change_param_input 'ntx' '10'
change_param_input 'nty' '10'
perform_test $file 2


#test_