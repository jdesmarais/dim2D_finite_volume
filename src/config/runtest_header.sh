
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
param_dim2d=$augeanstables/src/physical_models/dim2d/dim2d_parameters.f


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

change_param(){
    $config_dir/change_parameter.sh -i $1 -o $1 -p $2 -v $3
}


change_param_makefile(){
    change_param $make_header $1 $2
}


change_param_input(){
    change_param $param_input $1 $2
}

comment_param_dim2d(){
    tmp=temp.f
    sed -e s/"$1"/"!$1"/g < $param_dim2d > $tmp
    mv $tmp $param_dim2d
}

uncomment_param_dim2d(){
    tmp=temp.f
    sed -e s/"!$1"/"$1"/g < $param_dim2d > $tmp
    mv $tmp $param_dim2d
}

use_test_cv_r_dim2d(){

    comment_param_dim2d 'real(rkind), parameter :: cv_r = dim2d_M'
    uncomment_param_dim2d 'real(rkind), parameter :: cv_r = 2.5d0'

}

use_normal_cv_r_dim2d(){

    uncomment_param_dim2d 'real(rkind), parameter :: cv_r = dim2d_M'
    comment_param_dim2d 'real(rkind), parameter :: cv_r = 2.5d0'

}

use_test_param_dim2d(){

    echo '******************************************************'
    echo '*WARNING: using test parameters in dim2d_parameters.f*'
    echo '******************************************************'

    use_test_cv_r_dim2d

    comment_param_dim2d 'real(rkind), parameter :: viscous_r = dim2d_nu'
    comment_param_dim2d 'real(rkind), parameter :: Re = rho_c'
    comment_param_dim2d 'real(rkind), parameter :: We = (length_c'
    comment_param_dim2d 'real(rkind), parameter :: Pr = dim2d_mu'
    comment_param_dim2d 'real(rkind), parameter :: gravity = gravity_amp'

    uncomment_param_dim2d 'real(rkind), parameter :: viscous_r = -1.5d0'
    uncomment_param_dim2d 'real(rkind), parameter :: Re = 5.0d0'
    uncomment_param_dim2d 'real(rkind), parameter :: We = 10.0d0'
    uncomment_param_dim2d 'real(rkind), parameter :: Pr = 20.0d0'
    uncomment_param_dim2d 'real(rkind), parameter :: gravity = 9.81d0'

}


use_normal_param_dim2d(){

    echo '********************************************************'
    echo '*WARNING: using normal parameters in dim2d_parameters.f*'
    echo '********************************************************'

    use_normal_cv_r_dim2d

    uncomment_param_dim2d 'real(rkind), parameter :: viscous_r = dim2d_nu'
    uncomment_param_dim2d 'real(rkind), parameter :: Re = rho_c'
    uncomment_param_dim2d 'real(rkind), parameter :: We = (length_c'
    uncomment_param_dim2d 'real(rkind), parameter :: Pr = dim2d_mu'
    uncomment_param_dim2d 'real(rkind), parameter :: gravity = gravity_amp'

    comment_param_dim2d 'real(rkind), parameter :: viscous_r = -1.5d0'
    comment_param_dim2d 'real(rkind), parameter :: Re = 5.0d0'
    comment_param_dim2d 'real(rkind), parameter :: We = 10.0d0'
    comment_param_dim2d 'real(rkind), parameter :: Pr = 20.0d0'
    comment_param_dim2d 'real(rkind), parameter :: gravity = 9.81d0'

}


perform_test(){
    make_title $title_op $1
    if [ -z "$2" ]
    then
	$runtest $1 $runtest_op
    else
	$runtest $1 $runtest_op $2
    fi
}


plot_domain(){
    python -i $augeanstables/scripts_py/scripts_visit/plot_domain.py $1
}
