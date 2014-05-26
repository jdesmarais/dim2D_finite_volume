#!/bin/bash

#cleaning
make cleanall
rm *.dat

#compilation
make test_bf_layer_update_grdpts_id_prog

#run
./test_bf_layer_update_grdpts_id_prog

#print
folder_path=$(pwd)
pushd ../../python_files
python plot_test_bf_layer_update_grdpts_id.py -f $folder_path