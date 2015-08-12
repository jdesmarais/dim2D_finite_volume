#!/bin/bash

source $augeanstables/scripts_py/wall_steady_state_automatization/contours_extraction/library_extract_contours.sh


#============================================================
# extraction of the contours for the nucleation simulations
#============================================================
T=0.95
v=0.2
hca=0.0

# contact angle study at constant velocity
for v in 0.5
do
    for ca in 90.0
    do
	generate_linear_inflow_sph_contours $T $ca $v $hca [0,502,2] [-0.15,0.80] [0,0.25]
    done
done