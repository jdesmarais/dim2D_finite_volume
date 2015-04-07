#!/bin/bash

source $augeanstables/src/config/runtest_header.sh

#test_dir
test_dir=$augeanstables/src/test_files/test_field


#============================================================
#main body
#============================================================
AUGEANSTABLES_PARALLEL=false
change_param_input 'npx' '1'
change_param_input 'npy' '1'

echo ''

###test_field_extended_class
###------------------------------------------------------------
##file='test_field_extended'
##change_param_makefile 'pm_choice' 'wave2d_choice'
##change_param_makefile 'bc_choice' 'hedstrom_xy_choice'
##change_param_makefile 'ic_choice' 'sincos'
##change_param_input 'ic_choice' 'sincos'
##
##
###generate one domain results
##change_param_input 'ntx' '100'
##change_param_input 'nty' '110'
##change_param_input 'ne'  '3'
##change_param_input 'x_min' "\-10.0d0"
##change_param_input 'x_max' "\-0.5d0"
##change_param_input 'y_min' "\-10.0d0"
##change_param_input 'y_max' "11.0d0"
##change_param $test_dir/test_field_extended.f 'generate_small_domain' '.true.'
##
##make test_field_extended > /dev/null
##./test_field_extended > /dev/null
##make cleanall > /dev/null
##
##
###compare with interior+buffer layer results
##change_param_input 'ntx'  '64'
##change_param_input 'nty'  '54'
##change_param_input 'ne'   '3'
##change_param_input 'x_min'  "\-7.2d0"
##change_param_input 'x_max'  "\-1.3d0"
##change_param_input 'y_min'  "\-6.4d0"
##change_param_input 'y_max'  "3.4d0"
##change_param_input 'debug_adapt_computational_domain'  ".false."
##change_param $test_dir/test_field_extended.f 'generate_small_domain' '.false.'
##
##
##file=test_field_extended
##perform_test $file
##
##
###remove unnecessary files
##rm field_nodes0.out
##rm field_timedev.out
##rm field_nodes1st.out
##rm field_nodesInt.out
##rm interior_grdpts_id.nc
##echo ''


#test field_par
#------------------------------------------------------------
AUGEANSTABLES_PARALLEL=true
file='test_field_par'
#change_param_makefile 'pm_choice' 'wave2d_choice'
#
#change_param_input 'bc_N_type_choice' 'bc_nodes_choice'
#change_param_input 'bc_S_type_choice' 'bc_nodes_choice'
#change_param_input 'bc_E_type_choice' 'bc_nodes_choice'
#change_param_input 'bc_W_type_choice' 'bc_nodes_choice'
#
#change_param_input 'bc_choice' 'reflection_xy_choice'
#change_param_makefile 'bc_choice' 'reflection_xy_choice'
#
#change_param_makefile 'ic_choice' 'sincos'
#change_param_input 'ic_choice' 'sincos'
#
#
##generate one domain results
#change_param_input 'ntx' '100'
#change_param_input 'nty' '110'
#change_param_input 'npx' '1'
#change_param_input 'npy' '1'
#change_param_input 'ne'  '3'
#change_param_input 'x_min' "\-10.0d0"
#change_param_input 'x_max' "\-0.5d0"
#change_param_input 'y_min' "\-10.0d0"
#change_param_input 'y_max' "11.0d0"
#change_param $test_dir/test_field_par.f 'generate_small_domain' '.true.'
#
#make test_field_par #> /dev/null
#./test_field_par #> /dev/null
#make cleanall #> /dev/null


#compare with interior+buffer layer results
AUGEANSTABLES_PARALLEL=true
change_param_input 'ntx'  '104'
change_param_input 'nty'  '114'
change_param_input 'npx' '2'
change_param_input 'npy' '2'
change_param_input 'ne'   '3'
change_param $test_dir/test_field_par.f 'generate_small_domain' '.false.'


file=test_field_par
perform_test $file 4


##remove unnecessary files
#rm field_nodes0.out
#rm field_timedev.out
#rm field_nodes1st.out
#rm field_nodesInt.out
#echo ''