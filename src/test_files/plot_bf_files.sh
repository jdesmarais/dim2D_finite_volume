#!/bin/bash


#compile the test case for buffer layer
make test_bf_layer_abstract_prog

#run the test case
./test_bf_layer_abstract_prog

#run the visualization using python
folder_path=$(pwd)
pushd ../../python_files
python plot_bf_layers.py -f $folder_path
