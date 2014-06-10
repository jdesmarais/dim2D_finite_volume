#!/bin/bash

#clean directory from previous compilation
#and data files
make cleanall
rm *.dat

#compile the test case for the interface
make test_bf_interface_icr_update_prog

#run the test case
./test_bf_interface_icr_update_prog

#run the visualization using python
folder_path=$(pwd)
pushd ../../python_files
python plot_bf_interface_icr_update.py -f $folder_path