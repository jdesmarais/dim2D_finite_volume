
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

change_param_dim2d(){
    change_param $param_dim2d $1 $2
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
