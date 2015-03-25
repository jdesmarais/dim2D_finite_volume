#!/bin/bash

runtest=~/local/runtest/runtest.sh

#dir paths
config_dir=$augeanstables/src/config
param_dir=$augeanstables/src/parameters

#file paths
make_header=$config_dir/makefile_header.mk
param_input=$param_dir/parameters_input.f

#test_dir
test_dir=$augeanstables/src/test_files/test_bf_layer

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
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'bubble_transported'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '4'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'bubble_transported'
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
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'newgrdpt_test'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'newgrdpt_test'
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


#test_bf_compute_basic
echo ''
echo 'test_bf_compute_basic'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '6'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '6'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '4'
make test_bf_compute_basic > /dev/null
./test_bf_compute_basic
make cleanall > /dev/null
echo ''


#test_bf_compute_time
echo ''
echo 'test_bf_compute_time'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '6'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '6'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne' -v '4'
make test_bf_compute_time > /dev/null
./test_bf_compute_time
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


#test_bf_mainlayer_time_class
echo ''
echo 'test_bf_mainlayer_time_class'
echo '------------------------------------------------------------'
make test_bf_mainlayer_time > /dev/null
./test_bf_mainlayer_time
make cleanall > /dev/null
echo ''


#test_mainlayer_interface_basic_class
echo ''
echo 'test_mainlayer_interface_basic_class'
echo '------------------------------------------------------------'
make test_mainlayer_interface_basic > /dev/null
./test_mainlayer_interface_basic
make cleanall > /dev/null
echo ''


#test_mainlayer_interface_sync_class
echo ''
echo 'test_mainlayer_interface_sync_class'
echo '------------------------------------------------------------'
make test_mainlayer_interface_sync > /dev/null
./test_mainlayer_interface_sync
make cleanall > /dev/null
echo ''


#test_mainlayer_interface_dyn_class
echo ''
echo 'test_mainlayer_interface_dyn_class'
echo '------------------------------------------------------------'
make test_mainlayer_interface_dyn > /dev/null
./test_mainlayer_interface_dyn
make cleanall > /dev/null
echo ''


#test_bf_interface_basic_class
echo ''
echo 'test_bf_interface_basic_class'
echo '------------------------------------------------------------'
make test_bf_interface_basic > /dev/null
./test_bf_interface_basic
make cleanall > /dev/null
echo ''


#test_bf_interface_print_class
echo ''
echo 'test_bf_interface_print_class'
echo '------------------------------------------------------------'
make test_bf_interface_print > /dev/null
./test_bf_interface_print
make cleanall > /dev/null
echo ''


#test_bf_interface_sync_class
echo ''
echo 'test_bf_interface_sync_class'
echo '------------------------------------------------------------'
make test_bf_interface_sync > /dev/null
./test_bf_interface_sync
make cleanall > /dev/null
echo ''


#test_bf_interface_dyn_class
echo ''
echo 'test_bf_interface_dyn_class'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '30'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '35'
make test_bf_interface_dyn > /dev/null
./test_bf_interface_dyn
make cleanall > /dev/null
echo ''


#test_bf_mainlayer_bc_sections_module
echo ''
echo 'test_bf_mainlayer_bc_sections_module'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '30'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '35'
make test_bf_mainlayer_bc_sections > /dev/null
./test_bf_mainlayer_bc_sections
make cleanall > /dev/null
echo ''


#test_bf_interface_time_class
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
make test_bf_interface_time > /dev/null
./test_bf_interface_time
make cleanall > /dev/null

#compare with interior+buffer layer results
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '64'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '54'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '3'
$config_dir/change_parameter.sh -i $test_dir/test_bf_interface_time.f -o $test_dir/test_bf_interface_time.f -p 'generate_small_domain'  -v '.false.'
make test_bf_interface_time > /dev/null
./test_bf_interface_time
make cleanall > /dev/null

#remove unnecessary files
rm nodes0.out
rm timedev.out
rm nodes1st.out
echo ''


#test_bf_newgrdpt_dispatch_module
echo ''
echo 'test_bf_newgrdpt_dispatch_module'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'dim2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'hedstrom_xy_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'newgrdpt_test'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'newgrdpt_test'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '64'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '54'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '4'
make test_bf_newgrdpt_dispatch > /dev/null
./test_bf_newgrdpt_dispatch
make cleanall > /dev/null
echo ''


#test_bf_compute_newgrdpt_class
echo ''
echo 'test_bf_compute_newgrdpt_class'
echo '------------------------------------------------------------'
make test_bf_compute_newgrdpt > /dev/null
./test_bf_compute_newgrdpt
make cleanall > /dev/null
echo ''


#test_bf_layer_newgrdpt_class
echo ''
echo 'test_bf_layer_newgrdpt_class'
echo '------------------------------------------------------------'
make test_bf_layer_newgrdpt > /dev/null
./test_bf_layer_newgrdpt
make cleanall > /dev/null
echo ''


#test_mainlayer_interface_newgrdpt
echo ''
echo 'test_mainlayer_interface_newgrdpt_class'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'dim2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'hedstrom_xy_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'newgrdpt_test'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'newgrdpt_test'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '64'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '54'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '4'
make test_mainlayer_interface_newgrdpt > /dev/null
./test_mainlayer_interface_newgrdpt
make cleanall > /dev/null
echo ''


#test_bf_bc_interior_pt_crenel
echo ''
echo 'test_bf_bc_interior_pt_crenel'
echo '------------------------------------------------------------'
make test_bf_bc_interior_pt_crenel > /dev/null
./test_bf_bc_interior_pt_crenel
make cleanall > /dev/null
echo ''


#test_bf_bc_pt_crenel
echo ''
echo 'test_bf_bc_pt_crenel'
echo '------------------------------------------------------------'
make test_bf_bc_pt_crenel > /dev/null
./test_bf_bc_pt_crenel
make cleanall > /dev/null
echo ''


#test_bf_layer_grdpts_id_update
echo ''
echo 'test_bf_layer_grdpts_id_update'
echo '------------------------------------------------------------'
make test_bf_layer_grdpts_id_update > /dev/null
./test_bf_layer_grdpts_id_update
make cleanall > /dev/null
echo ''


#test_mainlayer_interface_grdpts_id_update
echo ''
echo 'test_mainlayer_interface_grdpts_id_update'
echo '------------------------------------------------------------'
make test_mainlayer_interface_grdpts_id_update > /dev/null
./test_mainlayer_interface_grdpts_id_update
make cleanall > /dev/null
echo ''


#test_bf_increase_coords_module
echo ''
echo 'test_bf_increase_coords_module'
echo '------------------------------------------------------------'
make test_bf_increase_coords > /dev/null
./test_bf_increase_coords
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


#test_bf_interface_coords_class
echo ''
echo 'test_bf_interface_coords_class'
echo '------------------------------------------------------------'
make test_bf_interface_coords > /dev/null
./test_bf_interface_coords
make cleanall > /dev/null
echo ''


#test_icr_path_class
echo ''
echo 'test_icr_path_class'
echo '------------------------------------------------------------'
make test_icr_path > /dev/null
./test_icr_path
make cleanall > /dev/null
echo ''


#test_icr_path_chain_class
echo ''
echo 'test_icr_path_chain_class'
echo '------------------------------------------------------------'
make test_icr_path_chain > /dev/null
./test_icr_path_chain
make cleanall > /dev/null
echo ''


#test_icr_interface_class
echo ''
echo 'test_icr_interface_class'
echo '------------------------------------------------------------'
make test_icr_interface > /dev/null
./test_icr_interface
make cleanall > /dev/null
echo ''


#test_bf_sorting_module
echo ''
echo 'test_bf_sorting_module'
echo '------------------------------------------------------------'
make test_bf_sorting_module > /dev/null
./test_bf_sorting_module
make cleanall > /dev/null
echo ''


#test_bf_layer_icr_class
echo ''
echo 'test_bf_layer_icr_class'
echo '------------------------------------------------------------'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'pm_choice' -v 'wave2d_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'bc_choice' -v 'hedstrom_xy_choice'
$config_dir/change_parameter.sh -i $make_header -o $make_header -p 'ic_choice' -v 'peak'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ic_choice' -v 'peak'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ntx' -v '64'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'nty' -v '54'
$config_dir/change_parameter.sh -i $param_input -o $param_input -p 'ne'  -v '3'
make test_bf_layer_icr > /dev/null
./test_bf_layer_icr
make cleanall > /dev/null
echo ''


#test_bf_interface_icr_class
echo ''
echo 'test_bf_interface_icr_class'
echo '------------------------------------------------------------'
make test_bf_interface_icr > /dev/null
./test_bf_interface_icr
make cleanall > /dev/null
echo ''


#test_bf_remove_module
echo ''
echo 'test_bf_remove_module'
echo '------------------------------------------------------------'

make test_bf_remove_module > /dev/null
./test_bf_remove_module
make cleanall > /dev/null
echo ''