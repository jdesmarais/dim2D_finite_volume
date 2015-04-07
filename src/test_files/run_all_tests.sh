#!/bin/bash


runtest_folder(){
    echo '----------------------------------------------------------------------'
    echo $1
    echo '----------------------------------------------------------------------'
    cd $1
    chmod 777 *.sh
    ./run_tests.sh -s
    echo ''
}


# tests for the extension of the interior domain
folder=$augeanstables/src/test_files/test_bf_layer
runtest_folder $folder


# tests for the boundary conditions
folder=$augeanstables/src/test_files/test_bc_operators
runtest_folder $folder


# tests for the time discretization operators
folder=$augeanstables/src/test_files/test_td_operators
runtest_folder $folder


# tests for the i/o operators
folder=$augeanstables/src/test_files/test_io_operators
runtest_folder $folder


# test for the mpi_operators
folder=$augeanstables/src/test_files/test_mpi_operators
runtest_folder $folder


# tests for the field
folder=$augeanstables/src/test_files/test_field
runtest_folder $folder
