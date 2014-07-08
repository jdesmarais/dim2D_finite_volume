#!/bin/bash


# serial version: 1200x1200
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 1
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 1
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_min -v -0.6
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_max -v 0.6
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_min -v -0.6
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_max -v 0.6
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d /scratch-local/jdesmara/scaling/run_serial
echo 'serial ready: 1200x1200'


# 2x1 version: 2400x1200
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 2
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 1
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_min -v -1.2
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_max -v 1.2
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_min -v -0.6
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_max -v 0.6
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_2x1 /scratch-local/jdesmara/scaling/run_2x1
echo '2x1 ready: 2400x1200'


# 3x1 version: 3600x1200
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 3
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 1
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_min -v -1.8
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_max -v 1.8
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_min -v -0.6
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_max -v 0.6
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_3x1 /scratch-local/jdesmara/scaling/run_3x1
echo '3x1 ready: 3600x1200'


# 2x2 version: 2400x2400
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 2
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 2
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_min -v -1.2
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_max -v 1.2
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_min -v -1.2
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_max -v 1.2
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_2x2 /scratch-local/jdesmara/scaling/run_2x2
echo '2x2 ready: 2400x2400'


# 3x2 version: 3600x2400
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 3
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 2
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_min -v -1.8
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_max -v 1.8
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_min -v -1.2
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_max -v 1.2
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_3x2 /scratch-local/jdesmara/scaling/run_3x2
echo '3x2 ready: 3600x2400'


# 3x3 version: 3600x3600
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 3
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 3
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_min -v -1.8
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_max -v 1.8
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_min -v -1.8
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_max -v 1.8
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_3x3 /scratch-local/jdesmara/scaling/run_3x3
echo '3x3 ready: 3600x3600'


# 4x3 version: 4800x3600
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 4
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 3
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_min -v -2.4
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_max -v 2.4
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_min -v -1.8
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_max -v 1.8
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_4x3 /scratch-local/jdesmara/scaling/run_4x3
echo '4x3 ready: 4800x3600'


# 4x4 version: 4800x4800
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npx -v 4
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p npy -v 4
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_min -v -2.4
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p x_max -v 2.4
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_min -v -2.4
./change_parameter.sh -i inputs_drop_retractation.txt -o inputs_drop_retractation.txt -p y_max -v 2.4
./config.py -i inputs_drop_retractation.txt -c
cp ../sim_dim2d_4x4 /scratch-local/jdesmara/scaling/run_4x4
echo '4x4 ready: 4800x4800'
