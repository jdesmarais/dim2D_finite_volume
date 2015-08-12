#!/bin/bash

source $augeanstables/scripts_py/wall_steady_state_automatization/contours_extraction/library_extract_contours.sh


#============================================================
# extraction of the contours for the nucleation simulations
#============================================================
T=0.95

## heat flux study at constant contact angle
#ca=90.0
#for fh in 0.04 0.06 0.08 0.1
#do
#    generate_nucleation_contours $T $ca $fh [0,215,1] [-0.15,0.15] [0,0.125]
#done

# contact angle study at constant heat flux
fh=0.02
data=[0,215,1]
data=[0,170,1]
xlimits=[-0.15,0.15]
ylimits=[0,0.125]

for ca in 135.0 #22.5 45.0 67.5 90.0 112.5 135.0 #67.5 90.0 112.5 ##
do
    generate_nucleation_contours $T $ca $fh $data $xlimits $ylimits
done
