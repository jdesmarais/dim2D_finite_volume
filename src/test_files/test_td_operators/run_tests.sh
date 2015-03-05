#!/bin/bash

runtest=~/local/runtest/runtest.sh

#dir paths
config_dir=$augeanstables/src/config
param_dir=$augeanstables/src/parameters

#file paths
make_header=$config_dir/makefile_header.mk
param_input=$param_dir/parameters_input.f

runtest=~/local/runtest/runtest.sh

#test_td_operators.f
echo ''
echo 'test_td_operators'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'simpletest_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'periodic_xy_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_choice' -v 'periodic_xy_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_N_type_choice' -v 'bc_nodes_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_S_type_choice' -v 'bc_nodes_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_E_type_choice' -v 'bc_nodes_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_W_type_choice' -v 'bc_nodes_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '10'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '6'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '1'
make test_td_operators > /dev/null
./test_td_operators
make cleanall > /dev/null
echo ''