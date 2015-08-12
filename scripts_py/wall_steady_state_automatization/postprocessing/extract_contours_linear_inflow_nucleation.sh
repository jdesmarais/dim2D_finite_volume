#!/bin/bash

source $augeanstables/scripts_py/wall_steady_state_automatization/contours_extraction/library_extract_contours.sh


#============================================================
# extraction of the contours for the nucleation simulations
#============================================================

## flow inflow study at constant contact angle and heat flux
T=0.95
#ca=22.5
#fh=0.04
hca=0.0


steps=[0,550,1]
xlimits=[-0.1,0.425]
ylimits=[0,0.125]

for fh in 0.04
do
    for ca in 22.5
    do
	for v in 0.0 #0.5 #0.1 0.2 0.3 0.4 0.5
	do
	    options="$T $ca $fh $v $hca $steps $xlimits $ylimits"
	    generate_linear_inflow_nucleation_contours $options
	    #generate_inflow_nucleation_contours $options
	done
    done
done

#for fh in 0.05 0.06 0.08 0.1
#do
#    for ca in 22.5
#    do
#	for v in 0.2 0.4
#	do
#	    options="$T $ca $fh $v $hca $steps $xlimits $ylimits"
#	    generate_linear_inflow_nucleation_contours $options
#	    #generate_inflow_nucleation_contours $options
#	done
#    done
#done
