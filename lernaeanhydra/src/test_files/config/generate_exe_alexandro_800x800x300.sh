#!/bin/bash


# serial version: 800x800x300
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p npx -v 1
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p npy -v 1
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p x_min -v -4
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p x_max -v 4
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p y_min -v -4
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p y_max -v 4
./config.py -i inputs_comparison_alexandro.txt -c
cp ../sim_dim2d /projects/fluiddyn/comparison_alexandro/800x800x300/run_serial
echo 'serial ready: 800x800x300'


# 2x2 version: 800x800x300
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p npx -v 2
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p npy -v 2
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p x_min -v -4
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p x_max -v 4
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p y_min -v -4
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p y_max -v 4
./config.py -i inputs_comparison_alexandro.txt -c
cp ../sim_dim2d /projects/fluiddyn/comparison_alexandro/800x800x300/run_2x2
echo '2x2 ready: 800x800x300'


# 4x4 version: 800x800x300
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p npx -v 4
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p npy -v 4
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p x_min -v -4
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p x_max -v 4
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p y_min -v -4
./change_parameter.sh -i inputs_comparison_alexandro.txt -o inputs_comparison_alexandro.txt -p y_max -v 4
./config.py -i inputs_comparison_alexandro.txt -c
cp ../sim_dim2d /projects/fluiddyn/comparison_alexandro/800x800x300/run_4x4
echo '4x4 ready: 800x800x300'