#!/bin/bash

runtest=~/local/runtest/runtest.sh

#dir paths
config_dir=$augeanstables/src/config
param_dir=$augeanstables/src/parameters

#file paths
make_header=$config_dir/makefile_header.mk
param_input=$param_dir/parameters_input.f

#test_dir
test_dir=$augeanstables/src/test_files/test_field


#test_field_extended_class
echo ''
echo 'test_bf_interface_time_class'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'wave2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'hedstrom_xy_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'peak'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'peak'

#generate one domain results
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '100'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '110'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '3'
$config_dir/change_parameter.sh -i $test_dir/test_bf_interface_time.f -o $test_dir/test_bf_interface_time.f -p 'generate_small_domain'  -v '.true.'
make test_field_extended > /dev/null
./test_field_extended
make cleanall > /dev/null

#compare with interior+buffer layer results
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '64'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '54'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '3'
$config_dir/change_parameter.sh -i $test_dir/test_field_extended.f -o $test_dir/test_field_extended.f -p 'generate_small_domain'  -v '.false.'
make test_field_extended > /dev/null
./test_field_extended
make cleanall > /dev/null

#remove unnecessary files
rm nodes0.out
rm timedev.out
rm nodes1st.out
echo ''
