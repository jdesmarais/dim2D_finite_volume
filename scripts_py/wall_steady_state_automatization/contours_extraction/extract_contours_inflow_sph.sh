#!/bin/bash

source $augeanstables/scripts_py/wall_steady_state_automatization/contours_extraction/library_extract_contours.sh


#============================================================
# extraction of the contours for the nucleation simulations
#============================================================
T=0.95
v=0.2
hca=0.0

# contact angle study at constant velocity
for v in 0.3 # 0.2 0.3 0.4
do
    for ca in 67.5 #22.5 45.0 67.5 90.0 112.5 135.0
    do
	generate_inflow_sph_contours $T $ca $v $hca [350,450,2] [-0.15,0.40] [0,0.250]
    done
done
