#!/bin/bash

runtest=~/local/runtest/runtest.sh

#dir paths
config_dir=$augeanstables/src/config
param_dir=$augeanstables/src/parameters

#file paths
make_header=$config_dir/makefile_header.mk
param_input=$param_dir/parameters_input.f

#test_bf_interior_bc_sections
echo ''
echo 'test_bf_interior_bc_sections'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '10'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '10'
make test_bf_interior_bc_sections > /dev/null
./test_bf_interior_bc_sections
make cleanall > /dev/null
echo ''


#test_bf_layer_basic
echo ''
echo 'test_bf_layer_basic'
echo '------------------------------------------------------------'
make test_bf_layer_basic > /dev/null
./test_bf_layer_basic
make cleanall > /dev/null
echo ''


#test_bf_layer_bc_anticorner
echo ''
echo 'test_bf_layer_bc_anticorner'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'dim2d_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '4'
make test_bf_layer_bc_anticorner > /dev/null
./test_bf_layer_bc_anticorner
make cleanall > /dev/null
echo ''


#test_bf_layer_bc_checks
echo ''
echo 'test_bf_layer_bc_checks'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_N_type_choice' -v 'bc_timedev_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_S_type_choice' -v 'bc_timedev_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_E_type_choice' -v 'bc_timedev_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_W_type_choice' -v 'bc_timedev_choice'
make test_bf_layer_bc_checks > /dev/null
./test_bf_layer_bc_checks
make cleanall > /dev/null
echo ''


#test_bf_layer_bc_procedure
echo ''
echo 'test_bf_layer_bc_procedure'
echo '------------------------------------------------------------'
make test_bf_layer_bc_procedure > /dev/null
./test_bf_layer_bc_procedure
make cleanall > /dev/null
echo ''


#test_bf_layer_bc_sections
echo ''
echo 'test_bf_layer_bc_sections'
echo '------------------------------------------------------------'
make test_bf_layer_bc_sections > /dev/null
./test_bf_layer_bc_sections
make cleanall > /dev/null
echo ''


#test_bf_layer_bc_sections_overlap
echo ''
echo 'test_bf_layer_bc_sections_overlap'
echo '------------------------------------------------------------'
make test_bf_layer_bc_sections_overlap > /dev/null
./test_bf_layer_bc_sections_overlap
make cleanall > /dev/null
echo ''


#test_bf_layer_exchange
echo ''
echo 'test_bf_layer_exchange'
echo '------------------------------------------------------------'
make test_bf_layer_exchange > /dev/null
./test_bf_layer_exchange
make cleanall > /dev/null
echo ''


#test_bf_layer_extract
echo ''
echo 'test_bf_layer_extract'
echo '------------------------------------------------------------'
make test_bf_layer_extract > /dev/null
./test_bf_layer_extract
make cleanall > /dev/null
echo ''


#test_bf_layer_nf90_operators_prog
echo ''
echo 'test_bf_layer_nf90_operators_prog'
echo '------------------------------------------------------------'
make test_bf_layer_nf90_operators_prog > /dev/null
./test_bf_layer_nf90_operators_prog
make cleanall > /dev/null
echo ''


#test_bf_layer_print
echo ''
echo 'test_bf_layer_print'
echo '------------------------------------------------------------'
make test_bf_layer_print > /dev/null
./test_bf_layer_print
make cleanall > /dev/null
echo ''


#test_bf_layer_sync
echo ''
echo 'test_bf_layer_sync'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '100'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '150'
make test_bf_layer_sync > /dev/null
./test_bf_layer_sync
make cleanall > /dev/null
echo ''


#test_bf_newgrdpt_prim
echo ''
echo 'test_bf_newgrdpt_prim'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_eigenqties_strategy' -v 'obc_eigenqties_bc'
make test_bf_newgrdpt_prim > /dev/null
./test_bf_newgrdpt_prim
make cleanall > /dev/null
echo ''


#test_bf_newgrdpt_procedure
echo ''
echo 'test_bf_newgrdpt_procedure'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'newgrdpt_test'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'newgrdpt_test'
make test_bf_newgrdpt_procedure > /dev/null
./test_bf_newgrdpt_procedure
make cleanall > /dev/null
echo ''


#test_bf_newgrdpt_verification
echo ''
echo 'test_bf_newgrdpt_verification'
echo '------------------------------------------------------------'
make test_bf_newgrdpt_verification > /dev/null
./test_bf_newgrdpt_verification
make cleanall > /dev/null
echo ''


#test_bf_compute
echo ''
echo 'test_bf_compute'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '6'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '6'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '4'
make test_bf_compute > /dev/null
./test_bf_compute
make cleanall > /dev/null
echo ''


#test_bf_layer_time
echo ''
echo 'test_bf_layer_time'
echo '------------------------------------------------------------'
make test_bf_layer_time > /dev/null
./test_bf_layer_time
make cleanall > /dev/null
echo ''


#test bf_sublayer_class
echo ''
echo 'test_bf_sublayer_class'
echo '------------------------------------------------------------'
make test_bf_sublayer > /dev/null
./test_bf_sublayer
make cleanall > /dev/null
echo ''


#test bf_mainlayer_basic_class
echo ''
echo 'test_bf_mainlayer_basic_class'
echo '------------------------------------------------------------'
make test_bf_mainlayer_basic > /dev/null
./test_bf_mainlayer_basic
make cleanall > /dev/null
echo ''


#test_bf_mainlayer_print_class
echo ''
echo 'test_bf_mainlayer_print_class'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '20'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '25'
make test_bf_mainlayer_print > /dev/null
./test_bf_mainlayer_print
make cleanall > /dev/null
echo ''


#test_bf_mainlayer_sync_class
echo ''
echo 'test_bf_mainlayer_sync_class'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '20'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '25'
make test_bf_mainlayer_sync > /dev/null
./test_bf_mainlayer_sync
make cleanall > /dev/null
echo ''


#test_bf_mainlayer_class
echo ''
echo 'test_bf_mainlayer_class'
echo '------------------------------------------------------------'
make test_bf_mainlayer > /dev/null
./test_bf_mainlayer
make cleanall > /dev/null
echo ''