#!/bin/bash

source $augeanstables/src/config/runtest_header.sh

#test_dir
test_dir=$augeanstables/src/test_files/test_bf_layer


#============================================================
#main body
#============================================================
AUGEANSTABLES_PARALLEL=false
change_param_input 'npx' '1'
change_param_input 'npy' '1'

change_param_input 'bc_N_type_choice' 'bc_timedev_choice'
change_param_input 'bc_S_type_choice' 'bc_timedev_choice'
change_param_input 'bc_E_type_choice' 'bc_timedev_choice'
change_param_input 'bc_W_type_choice' 'bc_timedev_choice'

change_param_input 'bc_choice' 'hedstrom_xy_choice'
change_param_makefile 'bc_choice' 'hedstrom_xy_choice'

echo ''

#test_bf_interior_bc_sections
file='test_bf_interior_bc_sections'
change_param_input 'ntx' '10'
change_param_input 'nty' '10'
perform_test $file


#test_bf_layer_basic
file='test_bf_layer_basic'
perform_test $file


#test_bf_layer_bc_fluxes
file='test_bf_layer_bc_fluxes'
change_param_makefile 'pm_choice' 'dim2d_choice'
change_param_makefile 'ic_choice' 'bubble_transported'
change_param_input 'ne' '4'
change_param_input 'ic_choice' 'bubble_transported'
perform_test $file


#test_bf_layer_bc_checks
file='test_bf_layer_bc_checks'
change_param_input 'bc_N_type_choice' 'bc_timedev_choice'
change_param_input 'bc_S_type_choice' 'bc_timedev_choice'
change_param_input 'bc_E_type_choice' 'bc_timedev_choice'
change_param_input 'bc_W_type_choice' 'bc_timedev_choice'
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
change_param_input 'ntx' '100'
change_param_input 'nty' '150'
perform_test $file


#test_bf_newgrdpt_prim
file=test_bf_newgrdpt_prim
change_param_input 'flow_direction' 'x_direction'
change_param_input 'flow_x_side' '1.0d0'
change_param_input 'flow_y_side' '1.0d0'
change_param_input 'flow_velocity' '0.1d0'
change_param_input 'T0' '0.95d0'
change_param_input 'ic_choice' 'newgrdpt_test'
change_param_input 'obc_eigenqties_strategy' 'obc_eigenqties_bc'
change_param_makefile 'ic_choice' 'newgrdpt_test'

perform_test $file


#test_bf_newgrdpt_procedure
file=test_bf_newgrdpt_procedure
change_param_makefile 'ic_choice' 'newgrdpt_test'
change_param_input 'ic_choice' 'newgrdpt_test'
perform_test $file


#test_bf_newgrdpt_verification
file=test_bf_newgrdpt_verification
perform_test $file


#test_bf_compute_basic
file=test_bf_compute_basic
change_param_input 'ntx' '6'
change_param_input 'nty' '6'
change_param_input 'ne' '4'
perform_test $file


#test_bf_compute_time
file=test_bf_compute_time
change_param_input 'ntx' '6'
change_param_input 'nty' '6'
change_param_input 'ne' '4'
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
change_param_input 'ntx' '20'
change_param_input 'nty' '25'
perform_test $file


#test_bf_mainlayer_sync_class
file=test_bf_mainlayer_sync
change_param_input 'ntx' '20'
change_param_input 'nty' '25'
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


#test_mainlayer_interface_time_class
file=test_mainlayer_interface_time
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
change_param_input 'ntx' '30'
change_param_input 'nty' '35'
perform_test $file


#test_bf_mainlayer_bc_sections_module
file=test_bf_mainlayer_bc_sections
change_param_input 'ntx' '30'
change_param_input 'nty' '35'
perform_test $file


#test_bf_interface_time_class
file=test_bf_interface_time

change_param_makefile 'pm_choice' 'wave2d_choice'
change_param_makefile 'bc_choice' 'hedstrom_xy_choice'
change_param_makefile 'ic_choice' 'peak'
change_param_input 'ic_choice' 'peak'

#generate one domain results
change_param_input 'ntx' '100'
change_param_input 'nty' '110'
change_param_input 'ne'  '3'
change_param $test_dir/test_bf_interface_time.f 'generate_small_domain' ".true."
make test_bf_interface_time > /dev/null
./test_bf_interface_time >/dev/null
make cleanall > /dev/null

#compare with interior+buffer layer results
change_param_input 'ntx' '64'
change_param_input 'nty' '54'
change_param_input 'ne'  '3'
change_param $test_dir/test_bf_interface_time.f 'generate_small_domain' ".false."
perform_test $file

#remove unnecessary files
rm nodes0.out
rm timedev.out
rm nodes1st.out


#test_bf_newgrdpt_dispatch_module
file=test_bf_newgrdpt_dispatch
change_param_makefile 'pm_choice' 'dim2d_choice'
change_param_makefile 'bc_choice' 'hedstrom_xy_choice'
change_param_makefile 'ic_choice' 'newgrdpt_test'
change_param_input 'ic_choice' 'newgrdpt_test'
change_param_input 'ntx' '64'
change_param_input 'nty' '54'
change_param_input 'ne'  '4'
perform_test $file


#test_bf_compute_newgrdpt_class
file=test_bf_compute_newgrdpt
perform_test $file


#test_bf_layer_newgrdpt_class
file=test_bf_layer_newgrdpt
perform_test $file


#test_mainlayer_interface_newgrdpt
file=test_mainlayer_interface_newgrdpt
change_param_makefile 'pm_choice' 'dim2d_choice'
change_param_makefile 'bc_choice' 'hedstrom_xy_choice'
change_param_makefile 'ic_choice' 'newgrdpt_test'
change_param_input 'ic_choice' 'newgrdpt_test'
change_param_input 'ntx' '64'
change_param_input 'nty' '54'
change_param_input 'ne'  '4'
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
change_param_makefile 'pm_choice' 'wave2d_choice'
change_param_makefile 'bc_choice' 'hedstrom_xy_choice'
change_param_makefile 'ic_choice' 'peak'
change_param_input 'ic_choice' 'peak'
change_param_input 'ntx' '64'
change_param_input 'nty' '54'
change_param_input 'ne'  '3'
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
file=test_bf_interface_dcr
perform_test $file
