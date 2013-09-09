#!/bin/bash


# serial version
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 1
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 1
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d /scratch-local/jdesmara/run_serial
echo 'serial ready'


# 1x2 version
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 1
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 2
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_1x2 /scratch-local/jdesmara/run_1x2
echo '1x2 ready'


# 1x3 version
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 1
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 3
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_1x3 /scratch-local/jdesmara/run_1x3
echo '1x3 ready'


# 2x2 version
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 2
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 2
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_2x2 /scratch-local/jdesmara/run_2x2
echo '2x2 ready'


# 2x3 version
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 2
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 3
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_2x3 /scratch-local/jdesmara/run_2x3
echo '2x3 ready'


# 3x3 version
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 3
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 3
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_3x3 /scratch-local/jdesmara/run_3x3
echo '3x3 ready'


# 3x4 version
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 3
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 4
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_3x4 /scratch-local/jdesmara/run_3x4
echo '3x4 ready'


# 4x4 version
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 4
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 4
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_4x4 /scratch-local/jdesmara/run_4x4
echo '4x4 ready'
