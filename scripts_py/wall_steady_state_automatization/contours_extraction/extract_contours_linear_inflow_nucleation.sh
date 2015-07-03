#!/bin/bash

source $augeanstables/scripts_py/wall_steady_state_automatization/contours_extraction/library_extract_contours.sh


#============================================================
# extraction of the contours for the nucleation simulations
#============================================================

## flow inflow study at constant contact angle and heat flux
T=0.95
ca=22.5
sh=-0.04
hca=0.0

for sh in -0.04 #-0.05 -0.06 -0.08 -0.1
do
    for ca in 22.5
    do
	for v in 0.3 0.4 0.5 #0.1 0.2 0.3 0.4 0.5
	do
	    options="$T $ca $sh $v $hca [0,500,1] [-0.1,0.20] [0,0.125]"
	    generate_linear_inflow_nucleation_contours $options
	done
    done
done
