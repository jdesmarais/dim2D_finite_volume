#!/bin/bash

#clean directory from previous compilation
#and data files
make cleanall
rm *.dat

#compile the test case for the interface
make test_bf_interface_dcr_update_prog2

#run the test case
./test_bf_interface_dcr_update_prog2

#run the visualization using python
folder_path=$(pwd)
pushd ../../python_files
python make_movie_bf_layer.py -f $folder_path
