#!/bin/bash

#dir paths
config_dir=$augeanstables/src/test_files/config
param_dir=$augeanstables/src/parameters

#file paths
make_header=$config_dir/makefile_header.mk
param_input=$param_dir/parameters_input.f

runtest=~/local/runtest/runtest.sh

#test_hedstrom_xy.f
echo ''
echo 'test_hedstrom_xy'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'wave1d_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '7'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '5'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '2'
make test_hedstrom_xy > /dev/null
./test_hedstrom_xy
make cleanall > /dev/null
echo ''


#test_hedstrom_xy_corners.f
echo ''
echo 'test_hedstrom_xy_corners'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'wave2d_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '7'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '5'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '3'
make test_hedstrom_xy_corners > /dev/null
./test_hedstrom_xy_corners
make cleanall > /dev/null
echo ''


#test_lodi_relaxation_coeff.f
echo ''
echo 'test_lodi_relaxation_coeff'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'ns2d_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'sigma_P' -v '0.25d0'
make test_lodi_relaxation_coeff > /dev/null
./test_lodi_relaxation_coeff
make cleanall > /dev/null
echo ''


#test_openbc_operators.f
echo ''
echo 'test_openbc_operators'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'simpletest_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nx' -v '7'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ny' -v '5'
make test_openbc_operators > /dev/null
./test_openbc_operators
make cleanall > /dev/null
echo ''