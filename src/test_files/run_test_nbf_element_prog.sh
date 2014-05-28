#!/bin/bash

#clean directory from previous compilation
#and data files
make cleanall 1>/dev/null 2>&1
rm *.dat 1>/dev/null 2>&1


#clean directory from previous compilation
#and data files
make cleanall 1> /dev/null 2>&1
rm *.dat 1> /dev/null 2>&1

#compile the test case for buffer layer
make test_nbf_element_prog

#run the test case
./test_nbf_element_prog

#run the visualization using python
folder_path=$(pwd)
pushd ../../python_files
python plot_test_nbf_element.py -f $folder_path
