#!/bin/bash

source $augeanstables/src/config/runtest_header.sh

#test_dir
test_dir=$augeanstables/src/test_files/test_bf_layer


#============================================================
#main body
#============================================================
echo ''

#test_bf_interior_bc_sections
file='test_mpi_process'
change_param_input 'npx' '1'
change_param_input 'npy' '2'
change_param_input 'ntx' '10'
change_param_input 'nty' '10'
perform_test $file 2