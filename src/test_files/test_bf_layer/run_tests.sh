#!/bin/bash

runtest=~/local/runtest/runtest.sh


#============================================================
#options
#============================================================
if [ -z "$1" ]
then
    title_op='-title'
    runtest_op='-d'
else
    case $1 in
	'-s')
	    title_op='-no-title'
	    runtest_op='-s'
	    ;;
    esac
fi


#============================================================
#variables
#============================================================
#dir paths
config_dir=$augeanstables/src/config
param_dir=$augeanstables/src/parameters

#file paths
make_header=$config_dir/makefile_header.mk
param_input=$param_dir/parameters_input.f

#test_dir
test_dir=$augeanstables/src/test_files/test_bf_layer



#============================================================
#functions
#============================================================
make_title(){
    case $1 in
	'-no-title') ;;
	'-title')
	    echo ''
	    echo $2
	    echo '------------------------------------------------------------'
	    ;;
    esac
}

perform_test(){
    make_title $title_op $1
    $runtest $1 $runtest_op
}



#============================================================
#main body
#============================================================
echo ''

#test_bf_interior_bc_sections
file='test_bf_interior_bc_sections'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '10'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '10'
perform_test $file


#test_bf_layer_basic
file='test_bf_layer_basic'
perform_test $file


#test_bf_layer_bc_anticorner
file='test_bf_layer_bc_anticorner'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'dim2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'bubble_transported'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '4'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'bubble_transported'
perform_test $file


#test_bf_layer_bc_checks
file='test_bf_layer_bc_checks'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_N_type_choice' -v 'bc_timedev_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_S_type_choice' -v 'bc_timedev_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_E_type_choice' -v 'bc_timedev_choice'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'bc_W_type_choice' -v 'bc_timedev_choice'
perform_test $file


#test_bf_layer_bc_procedure
file='test_bf_layer_bc_procedure'
perform_test $file


#test_bf_layer_bc_sections
file='test_bf_layer_bc_sections'
perform_test $file


#test_bf_layer_bc_sections_overlap
file=test_bf_layer_bc_sections_overlap
perform_test $file


#test_bf_layer_exchange
file=test_bf_layer_exchange
perform_test $file


#test_bf_layer_extract
file=test_bf_layer_extract
perform_test $file


#test_bf_layer_nf90_operators_prog
file=test_bf_layer_nf90_operators_prog
perform_test $file


#test_bf_layer_print
file=test_bf_layer_print
perform_test $file


#test_bf_layer_sync
file=test_bf_layer_sync
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '100'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '150'
perform_test $file


#test_bf_newgrdpt_prim
file=test_bf_newgrdpt_prim
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'obc_eigenqties_strategy' -v 'obc_eigenqties_bc'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'newgrdpt_test'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'newgrdpt_test'
perform_test $file


#test_bf_newgrdpt_procedure
file=test_bf_newgrdpt_procedure
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'newgrdpt_test'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'newgrdpt_test'
perform_test $file


#test_bf_newgrdpt_verification
file=test_bf_newgrdpt_verification
perform_test $file


#test_bf_compute_basic
file=test_bf_compute_basic
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '6'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '6'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '4'
perform_test $file


#test_bf_compute_time
file=test_bf_compute_time
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '6'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '6'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '4'
perform_test $file


#test_bf_layer_time
file=test_bf_layer_time
perform_test $file


#test_bf_sublayer_class
file=test_bf_sublayer
perform_test $file


#test_bf_mainlayer_basic_class
file=test_bf_mainlayer_basic
perform_test $file


#test_bf_mainlayer_print_class
file=test_bf_mainlayer_print
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '20'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '25'
perform_test $file


#test_bf_mainlayer_sync_class
file=test_bf_mainlayer_sync
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '20'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '25'
perform_test $file


#test_bf_mainlayer_time_class
file=test_bf_mainlayer_time
perform_test $file


#test_mainlayer_interface_basic_class
file=test_mainlayer_interface_basic
perform_test $file


#test_mainlayer_interface_sync_class
file=test_mainlayer_interface_sync
perform_test $file


#test_mainlayer_interface_dyn_class
file=test_mainlayer_interface_dyn
perform_test $file


#test_bf_interface_basic_class
file=test_bf_interface_basic
perform_test $file


#test_bf_interface_print_class
file=test_bf_interface_print
perform_test $file


#test_bf_interface_sync_class
file=test_bf_interface_sync
perform_test $file


#test_bf_interface_dyn_class
file=test_bf_interface_dyn
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '30'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '35'
perform_test $file


#test_bf_mainlayer_bc_sections_module
file=test_bf_mainlayer_bc_sections
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '30'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '35'
perform_test $file


#test_bf_interface_time_class
file=test_bf_interface_time

$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'wave2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'hedstrom_xy_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'peak'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'peak'

#generate one domain results
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '100'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '110'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '3'
$config_dir/change_parameter.sh -i $test_dir/test_bf_interface_time.f -o $test_dir/test_bf_interface_time.f -p 'generate_small_domain'  -v '.true.'
make test_bf_interface_time > /dev/null
./test_bf_interface_time >/dev/null
make cleanall > /dev/null

#compare with interior+buffer layer results
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '64'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '54'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '3'
$config_dir/change_parameter.sh -i $test_dir/test_bf_interface_time.f -o $test_dir/test_bf_interface_time.f -p 'generate_small_domain'  -v '.false.'
perform_test $file

#remove unnecessary files
rm nodes0.out
rm timedev.out
rm nodes1st.out


#test_bf_newgrdpt_dispatch_module
file=test_bf_newgrdpt_dispatch
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'dim2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'hedstrom_xy_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'newgrdpt_test'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'newgrdpt_test'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '64'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '54'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '4'
perform_test $file


#test_bf_compute_newgrdpt_class
file=test_bf_compute_newgrdpt
perform_test $file


#test_bf_layer_newgrdpt_class
file=test_bf_layer_newgrdpt
perform_test $file


#test_mainlayer_interface_newgrdpt
file=test_mainlayer_interface_newgrdpt
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'dim2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'hedstrom_xy_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'newgrdpt_test'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'newgrdpt_test'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '64'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '54'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '4'
perform_test $file


#test_bf_bc_interior_pt_crenel
file=test_bf_bc_interior_pt_crenel
perform_test $file


#test_bf_bc_pt_crenel
file=test_bf_bc_pt_crenel
perform_test $file


#test_bf_layer_grdpts_id_update
file=test_bf_layer_grdpts_id_update
perform_test $file


#test_mainlayer_interface_grdpts_id_update
file=test_mainlayer_interface_grdpts_id_update
perform_test $file


#test_bf_increase_coords_module
file=test_bf_increase_coords
perform_test $file


#test_bf_mainlayer_class
file=test_bf_mainlayer
perform_test $file


#test_bf_interface_coords_class
file=test_bf_interface_coords
perform_test $file


#test_icr_path_class
file=test_icr_path
perform_test $file


#test_icr_path_chain_class
file=test_icr_path_chain
perform_test $file


#test_icr_interface_class
file=test_icr_interface
perform_test $file


#test_bf_sorting_module
file=test_bf_sorting
perform_test $file


#test_bf_layer_icr_class
file=test_bf_layer_icr
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'wave2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'hedstrom_xy_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'peak'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'peak'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '64'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '54'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '3'
perform_test $file


#test_bf_interface_icr_class
file=test_bf_interface_icr
perform_test $file


#test_bf_decrease_module
file=test_bf_decrease
perform_test $file


#test_dcr_chain_class
file=test_dcr_chain
perform_test $file


#test_dcr_list_class
file=test_dcr_list
perform_test $file


#test_bf_layer_class
file=test_bf_layer
perform_test $file


#test_dcr_interface_class
file=test_dcr_interface
perform_test $file


#test_bf_interface_class
file=test_bf_interface
perform_test $file
