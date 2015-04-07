#!/bin/bash

source $augeanstables/src/config/runtest_header.sh

#test_dir
test_dir=$augeanstables/src/test_files/test_geometry_updates


#============================================================
#functions
#============================================================
generate_geometry_update_case(){

    #replace the case_id in the fortran file
    change_param_makefile 'pm_choice' 'wave2d_choice'
    change_param_input 'ne' '3'
    change_param $test_dir/test_geometry_update.f 'main_case_choice' $1
    change_param $test_dir/test_geometry_update.f 'sub_case_choice' $2

    #create the netcdf files
    file='test_geometry_update'
    perform_test $file
    echo 'netcdf files generated'
    
    #create the frames using visit
    plot_domain "-g -d $test_dir --x_min=-35 --x_max=35 --y_min=-35 --y_max=35"
    echo 'movie generated'

    #clean the netcdf files created
    mkdir visit_grdpts_id/nc_files
    mv *.nc ./visit_grdpts_id/nc_files

    #rename the directory for the output
    case $1 in
	'spot_transported')
	    case $2 in
		'1') mv visit_grdpts_id spot_transported_E ;;
		'2') mv visit_grdpts_id spot_transported_W ;;
		'3') mv visit_grdpts_id spot_transported_N ;;
		'4') mv visit_grdpts_id spot_transported_S ;;
		'5') mv visit_grdpts_id spot_transported_NE ;;
		'6') mv visit_grdpts_id spot_transported_NW ;;
		'7') mv visit_grdpts_id spot_transported_SE ;;
		'8') mv visit_grdpts_id spot_transported_SW ;;
	    esac
	    ;;

	'bubble_expansion')
	    mv visit_grdpts_id bubble_expansion
	    ;;
    esac
}


#============================================================
#main body
#============================================================
AUGEANSTABLES_PARALLEL=false
change_param_input 'npx' '1'
change_param_input 'npy' '1'

echo ''

# spot transported cases
generate_geometry_update_case spot_transported 1
generate_geometry_update_case spot_transported 2
generate_geometry_update_case spot_transported 3
generate_geometry_update_case spot_transported 4
generate_geometry_update_case spot_transported 5
generate_geometry_update_case spot_transported 6
generate_geometry_update_case spot_transported 7
generate_geometry_update_case spot_transported 8


#bubble expansion
generate_geometry_update_case bubble_expansion 1

echo ''
