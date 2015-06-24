#!/bin/bash

source $augeanstables/scripts_py/wall_steady_state_automatization/contours_extraction/library_extract_contours.sh


#============================================================
# extraction of the contours for the nucleation simulations
#============================================================
T=0.95
v=0.1
hca=0.0

# contact angle study at constant velocity
for ca in 22.5 #22.5 45.0 67.5 90.0 112.5 135.0
do
    generate_inflow_sph_contours $T $ca $v $hca [0,505,5] [-0.15,0.15] [0,0.150]
done
