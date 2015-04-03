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


#test_mpi_derived_types
AUGEANSTABLES_PARALLEL=true
file='test_mpi_derived_types'
perform_test $file 2


#test_mpi_mg_neighbours
AUGEANSTABLES_PARALLEL=true
file='test_mpi_mg_neighbours'
change_param_input 'npx' '3'
change_param_input 'npy' '2'
perform_test $file 6


#test_mpi_mg_bc_class
AUGEANSTABLES_PARALLEL=true
file='test_mpi_mg_bc'
change_param_input 'ntx' '18'
change_param_input 'nty' '16'
change_param_input 'ne'  '2'
change_param_input 'npx' '3'
change_param_input 'npy' '2'
perform_test $file 6


#test_mpi_requests
AUGEANSTABLES_PARALLEL=true
file='test_mpi_requests'
change_param_input 'ntx' '18'
change_param_input 'nty' '16'
change_param_input 'ne'  '2'
change_param_input 'npx' '3'
change_param_input 'npy' '2'
perform_test $file 6


#test_mpi_mg_ini_bc_proc
AUGEANSTABLES_PARALLEL=true
file='test_mpi_mg_ini_bc_proc'
change_param_input 'npx' '3'
change_param_input 'npy' '2'
change_param_input 'bc_choice' 'reflection_xy_choice'
perform_test $file 6


#test_mpi_mg_construct
AUGEANSTABLES_PARALLEL=true
file='test_mpi_mg_construct'
change_param_input 'ntx' '18'
change_param_input 'nty' '16'
change_param_input 'ne'  '2'
change_param_input 'npx' '3'
change_param_input 'npy' '2'
change_param_input 'bc_choice' 'periodic_xy_choice'
perform_test $file 6


#test_mpi_mg_bc_ext
AUGEANSTABLES_PARALLEL=true
file='test_mpi_mg_bc_ext'
change_param_input 'ntx' '18'
change_param_input 'nty' '16'
change_param_input 'ne'  '2'
change_param_input 'npx' '3'
change_param_input 'npy' '2'
change_param_input 'bc_choice' 'periodic_xy_choice'
perform_test $file 6



#test_mpi_mg_bc_ext
AUGEANSTABLES_PARALLEL=true
file='test_mpi_interface'
change_param_input 'ntx' '18'
change_param_input 'nty' '16'
change_param_input 'ne'  '2'
change_param_input 'npx' '3'
change_param_input 'npy' '2'
change_param_input 'bc_choice' 'reflection_xy_choice'
perform_test $file 6