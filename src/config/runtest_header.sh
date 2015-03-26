
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
