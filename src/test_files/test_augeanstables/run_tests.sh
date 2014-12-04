#!/bin/bash

#dir paths
config_dir=$augeanstables/src/test_files/config
param_dir=$augeanstables/src/parameters
ns2d_dir=$augeanstables/src/physical_models/ns2d

#file paths
make_header=$config_dir/makefile_header.mk
param_input=$param_dir/parameters_input.f
ns2d_param=$ns2d_dir/ns2d_parameters.f

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
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice'  -v 'ns2d_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'sigma_P'    -v '0.25d0'
$config_dir/change_parameter.sh -i $ns2d_param  -o $ns2d_param  -p 'mach_infty' -v '0.1d0'
make test_lodi_relaxation_coeff > /dev/null
./test_lodi_relaxation_coeff
make cleanall > /dev/null
echo ''


#test_openbc_operators.f
echo ''
echo 'test_openbc_operators'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'simpletest_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '7'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '5'
make test_openbc_operators > /dev/null
./test_openbc_operators
make cleanall > /dev/null
echo ''


#test_poinsot_ns2d_operators.f
echo ''
echo 'test_poinsot_ns2d_operators'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'sd_choice'  -v 'mt_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice'  -v 'ns2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice'  -v 'poinsot_xy_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx'        -v '5'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty'        -v '5'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'         -v '4'
$config_dir/change_parameter.sh -i $ns2d_param  -o $ns2d_param  -p 'Re'         -v '10.0d0'
$config_dir/change_parameter.sh -i $ns2d_param  -o $ns2d_param  -p 'Pr'         -v '1.0d0'
$config_dir/change_parameter.sh -i $ns2d_param  -o $ns2d_param  -p 'mach_infty' -v '1.0d0'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_type_N' -v 'always_outflow'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_type_S' -v 'always_inflow'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_type_E' -v 'always_outflow'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_type_W' -v 'always_outflow'
make test_poinsot_ns2d_operators > /dev/null
./test_poinsot_ns2d_operators
make cleanall > /dev/null
echo ''


#test_yoo_ns2d_edge_x.f
echo ''
echo 'test_yoo_ns2d_edge_x'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'sd_choice'      -v 'mt_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice'      -v 'ns2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice'      -v 'yoolodato_xy_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx'            -v '5'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty'            -v '5'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'             -v '4'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'sigma_P'        -v '0.25'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'flow_direction' -v 'x_direction'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice'      -v 'peak'
$config_dir/change_parameter.sh -i $ns2d_param  -o $ns2d_param  -p 'Re'         -v '10.0d0'
$config_dir/change_parameter.sh -i $ns2d_param  -o $ns2d_param  -p 'Pr'         -v '1.0d0'
$config_dir/change_parameter.sh -i $ns2d_param  -o $ns2d_param  -p 'mach_infty' -v '0.2d0'
make test_yoo_ns2d_edge_x > /dev/null
./test_yoo_ns2d_edge_x
make cleanall > /dev/null
echo ''


#test_yoo_ns2d_edge_y.f
echo ''
echo 'test_yoo_ns2d_edge_y'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'sd_choice'      -v 'mt_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice'      -v 'ns2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice'      -v 'yoolodato_xy_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx'            -v '5'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty'            -v '5'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'             -v '4'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'sigma_P'        -v '0.25'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'flow_direction' -v 'y_direction'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice'      -v 'peak'
$config_dir/change_parameter.sh -i $ns2d_param  -o $ns2d_param  -p 'Re'         -v '10.0d0'
$config_dir/change_parameter.sh -i $ns2d_param  -o $ns2d_param  -p 'Pr'         -v '1.0d0'
$config_dir/change_parameter.sh -i $ns2d_param  -o $ns2d_param  -p 'mach_infty' -v '0.2d0'
make test_yoo_ns2d_edge_y > /dev/null
./test_yoo_ns2d_edge_y
make cleanall > /dev/null
echo ''

