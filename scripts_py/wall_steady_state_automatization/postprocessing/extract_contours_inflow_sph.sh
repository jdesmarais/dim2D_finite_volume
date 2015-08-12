#!/bin/bash

source $augeanstables/scripts_py/wall_steady_state_automatization/contours_extraction/library_extract_contours.sh


#============================================================
# extraction of the contours for the nucleation simulations
#============================================================
T=0.95
v=0.2
hca=0.0

# contact angle study at constant velocity
for v in 0.1 0.2 0.3 0.4 0.5
do
    for ca in 22.5 45.0 67.5 #90.0 112.5 135.0
    do
	generate_inflow_sph_contours $T $ca $v $hca [0,500,2] [-0.15,0.80] [0,0.25]
    done
done


#for v in 0.1 0.2 0.3 0.4 #0.1 0.2 0.3 0.4 0.5
#do
#    for ca in 90.0 112.5 135.0 #135.0 #90.0 112.5 135.0
#    do
#	generate_inflow_sph_contours $T $ca $v $hca [0,500,2] [-0.15,0.80] [0,0.25]
#    done
#done
