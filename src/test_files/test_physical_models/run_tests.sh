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

#test_gaussian_perturbation
file='test_gaussian_perturbation'
change_param_makefile 'pm_choice' 'dim2d_choice'
change_param_makefile 'ic_choice' 'homogeneous_liquid'
change_param_makefile 'bc_choice' 'periodic_xy_choice'
change_param_input 'bc_N_type_choice' 'bc_nodes_choice'
change_param_input 'bc_S_type_choice' 'bc_nodes_choice'
change_param_input 'bc_E_type_choice' 'bc_nodes_choice'
change_param_input 'bc_W_type_choice' 'bc_nodes_choice'
change_param_input 'x_min' "-0.5"
change_param_input 'x_max' " 0.5"
change_param_input 'y_min' "-0.5"
change_param_input 'y_max' " 0.5"
change_param_input 'ntx' "105"
change_param_input 'nty' "105"
change_param_input 'T0' "0.995d0"
change_param_input 'ic_choice' "homogeneous_liquid"
change_param_input 'ic_perturbation_ac' ".false."
change_param_input 'ne' "4"
perform_test $file