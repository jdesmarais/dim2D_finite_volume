#!/bin/bash

#dir paths
config_dir=$augeanstables/src/test_files/config
param_dir=$augeanstables/src/parameters
ns2d_dir=$augeanstables/src/physical_models/ns2d

#file paths
make_header=$config_dir/makefile_header.mk
param_input=$param_dir/parameters_input.f
ns2d_param=$ns2d_dir/ns2d_parameters.f


#test_field_extended_integrate.f: hedstrom b.c.
echo ''
echo 'test_field_extended_integrate: hedstrom_xy'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'sd_choice' -v 'mt_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'ns2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'peak'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'hedstrom_xy_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx'   -v '10'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty'   -v '10'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'    -v '4'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'x_min' -v '-0.1d0'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'x_max' -v '0.1d0'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'y_min' -v '-0.1d0'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'y_max' -v '0.1d0'
make test_field_extended_integrate > /dev/null
./test_field_extended_integrate
make cleanall > /dev/null
echo ''

##test_field_extended_integrate.f: hedstrom b.c. with corners
#echo ''
#echo 'test_field_extended_integrate: hedstrom_xy_corners'
#echo '------------------------------------------------------------'
#$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'sd_choice' -v 'mt_choice'
#$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'ns2d_choice'
#$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'peak'
#$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'hedstrom_xy_corners_choice'
#$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '10'
#$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '10'
#$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '4'
#$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'x_min' -v '-0.1d0'
#$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'x_max' -v '0.1d0'
#$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'y_min' -v '-0.1d0'
#$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'y_max' -v '0.1d0'
#make test_field_extended_integrate > /dev/null
#./test_field_extended_integrate
#make cleanall > /dev/null
#echo ''

#test_field_extended_integrate.f: poinsot b.c.
echo ''
echo 'test_field_extended_integrate: poinsot_xy'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'sd_choice' -v 'mt_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'ns2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'peak'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'poinsot_xy_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '10'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '10'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '4'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'x_min' -v '-0.1d0'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'x_max' -v '0.1d0'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'y_min' -v '-0.1d0'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'y_max' -v '0.1d0'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_type_N' -v 'ask_flow'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_type_S' -v 'ask_flow'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_type_E' -v 'ask_flow'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_type_W' -v 'ask_flow'
make test_field_extended_integrate > /dev/null
./test_field_extended_integrate
make cleanall > /dev/null
echo ''

#test_field_extended_integrate.f: yoo & lodato b.c.
echo ''
echo 'test_field_extended_integrate: yoolodato_xy'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'sd_choice' -v 'mt_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'ns2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'peak'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'yoolodato_xy_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '10'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '10'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '4'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'x_min' -v '-0.1d0'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'x_max' -v '0.1d0'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'y_min' -v '-0.1d0'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'y_max' -v '0.1d0'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_type_N' -v 'ask_flow'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_type_S' -v 'ask_flow'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_type_E' -v 'ask_flow'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_type_W' -v 'ask_flow'
make test_field_extended_integrate > /dev/null
./test_field_extended_integrate
make cleanall > /dev/null
echo ''