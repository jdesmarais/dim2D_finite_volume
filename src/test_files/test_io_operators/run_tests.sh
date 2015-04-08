#!/bin/bash

source $augeanstables/src/config/runtest_header.sh

#test_dir
test_dir=$augeanstables/src/test_files/test_bf_layer


#============================================================
#main body
#============================================================
AUGEANSTABLES_PARALLEL=false
change_param_input 'npx' '1'
change_param_input 'npy' '1'

echo ''

#test_io_operators
file='test_io_operators'
perform_test $file