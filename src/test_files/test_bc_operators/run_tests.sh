#!/bin/bash

source $augeanstables/src/config/runtest_header.sh

#test_dir
test_dir=$augeanstables/src/test_files/test_bf_layer


#============================================================
#main body
#============================================================
echo ''

#test_openbc_operators
file='test_openbc_operators'
perform_test $file


#test_hedstrom_xy
file='test_hedstrom_xy'
change_param_makefile 'pm_choice' 'dim2d_choice'
change_param_makefile 'ic_choice' 'newgrdpt_test'
change_param_dim2d 'cv_r' '2.5d0'
change_param_input 'nx' '6'
change_param_input 'ny' '6'
change_param_input 'ne' '4'
change_param_input 'flow_direction' 'x_direction'
change_param_input 'flow_x_side' '1.0d0'
change_param_input 'flow_y_side' '0.0d0'
change_param_input 'flow_velocity' '0.1d0'
change_param_input 'T0' '0.95d0'
change_param_input 'ic_choice' 'newgrdpt_test'
change_param_input 'obc_eigenqties_strategy' 'obc_eigenqties_lin'
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

