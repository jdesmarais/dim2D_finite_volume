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

#test_openbc_operators
file='test_openbc_operators'
change_param_makefile 'pm_choice' 'dim2d_choice'
change_param_makefile 'ic_choice' 'newgrdpt_test'
change_param_makefile 'bc_choice' 'hedstrom_xy_choice'
change_param_dim2d 'cv_r' '2.5d0'
change_param_input 'ntx' '6'
change_param_input 'nty' '6'
change_param_input 'ne' '4'
change_param_input 'flow_direction' 'x_direction'
change_param_input 'flow_x_side' '1.0d0'
change_param_input 'flow_y_side' '0.0d0'
change_param_input 'flow_velocity' '0.1d0'
change_param_input 'T0' '0.95d0'
change_param_input 'ic_choice' 'newgrdpt_test'
change_param_input 'obc_eigenqties_strategy' 'obc_eigenqties_lin'
change_param_input 'bc_choice' 'hedstrom_xy_choice'
change_param_input 'bc_N_type_choice' 'bc_timedev_choice'
change_param_input 'bc_S_type_choice' 'bc_timedev_choice'
change_param_input 'bc_E_type_choice' 'bc_timedev_choice'
change_param_input 'bc_W_type_choice' 'bc_timedev_choice'
perform_test $file


#test_hedstrom_xy
file='test_hedstrom_xy'
perform_test $file


#test_hedstrom_xy_anti_corner_flux
file='test_hedstrom_xy_anti_corner_flux'
perform_test $file


#test_bc_operators_openbc
file='test_bc_operators_openbc'
perform_test $file


#test_bc_operators_openbc_normal
file='test_bc_operators_openbc_normal'
perform_test $file


#test_reflection_xy
file='test_reflection_xy'
change_param_makefile 'pm_choice' 'wave2d_choice'
change_param_makefile 'ic_choice' 'peak'
change_param_input 'ntx' '6'
change_param_input 'nty' '8'
change_param_input 'ne' '3'
change_param_input 'ic_choice' 'peak'
change_param_input 'bc_choice' 'reflection_xy_choice'
change_param_input 'bc_N_type_choice' 'bc_nodes_choice'
change_param_input 'bc_S_type_choice' 'bc_nodes_choice'
change_param_input 'bc_E_type_choice' 'bc_nodes_choice'
change_param_input 'bc_W_type_choice' 'bc_nodes_choice'

perform_test $file


#test_periodic_xy
file='test_periodic_xy'
change_param_makefile 'pm_choice' 'wave2d_choice'
change_param_makefile 'ic_choice' 'peak'
change_param_input 'ntx' '8'
change_param_input 'nty' '10'
change_param_input 'ne' '3'
change_param_input 'ic_choice' 'peak'
change_param_input 'bc_choice' 'periodic_xy_choice'
change_param_input 'bc_N_type_choice' 'bc_nodes_choice'
change_param_input 'bc_S_type_choice' 'bc_nodes_choice'
change_param_input 'bc_E_type_choice' 'bc_nodes_choice'
change_param_input 'bc_W_type_choice' 'bc_nodes_choice'

perform_test $file