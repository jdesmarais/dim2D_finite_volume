#!/bin/bash

source $augeanstables/src/config/runtest_header.sh

#test_dir
test_dir=$augeanstables/src/test_files/test_field


#============================================================
#main body
#============================================================
echo ''

#test_field_extended_class
file='test_field_extended'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'wave2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'hedstrom_xy_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'sincos'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'sincos'

#generate one domain results
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '100'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '110'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '3'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'x_min'  -v "\-10.0d0"
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'x_max'  -v "\-0.5d0"
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'y_min'  -v "\-10.0d0"
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'y_max'  -v "11.0d0"
$config_dir/change_parameter.sh -i $test_dir/test_field_extended.f -o $test_dir/test_field_extended.f -p 'generate_small_domain'  -v '.true.'
make test_field_extended > /dev/null
./test_field_extended > /dev/null
make cleanall > /dev/null

#compare with interior+buffer layer results
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '64'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '54'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '3'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'x_min' -v "\-7.2d0"
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'x_max' -v "\-1.3d0"
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'y_min' -v "\-6.4d0"
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'y_max' -v "3.4d0"
$config_dir/change_parameter.sh -i $test_dir/test_field_extended.f -o $test_dir/test_field_extended.f -p 'generate_small_domain'  -v '.false.'
file=test_field_extended
perform_test $file

#remove unnecessary files
rm field_nodes0.out
rm field_timedev.out
rm field_nodes1st.out
rm field_nodesInt.out
rm interior_grdpts_id.nc
echo ''
