#!/bin/bash

#clean directory from previous compilation
#and data files
make cleanall
rm *.dat

#compile the test case for the interface
make test_interface_abstract_prog

#run the test case
./test_interface_abstract_prog

#run the visualization using python
folder_path=$(pwd)
pushd ../../python_files
python plot_test_interface_abstract.py -f $folder_path
