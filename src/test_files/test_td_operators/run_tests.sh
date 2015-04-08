#!/bin/bash

source $augeanstables/src/config/runtest_header.sh

#test_dir
test_dir=$augeanstables/src/test_files/test_td_operators


#============================================================
#main body
#============================================================
AUGEANSTABLES_PARALLEL=false
change_param_input 'npx' '1'
change_param_input 'npy' '1'

echo ''

#test_td_operators
file='test_td_operators'
change_param_makefile 'pm_choice' 'simpletest_choice'
change_param_makefile 'bc_choice' 'periodic_xy_choice'
change_param_input 'bc_choice' 'periodic_xy_choice'
change_param_input 'bc_N_type_choice' 'bc_nodes_choice'
change_param_input 'bc_S_type_choice' 'bc_nodes_choice'
change_param_input 'bc_E_type_choice' 'bc_nodes_choice'
change_param_input 'bc_W_type_choice' 'bc_nodes_choice'
change_param_input 'ntx' '10'
change_param_input 'nty' '6'
change_param_input 'ne' '1'
perform_test $file