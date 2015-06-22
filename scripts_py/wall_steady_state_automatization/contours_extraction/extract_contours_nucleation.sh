#!/bin/bash

source $augeanstables/scripts_py/wall_steady_state_automatization/contours_extraction/library_extract_contours.sh


#============================================================
# extraction of the contours for the nucleation simulations
#============================================================
T=0.95

## heat flux study at constant contact angle
#ca=90.0
#for sh in -0.04 -0.06 -0.08 -0.1
#do
#    generate_non_st_contours $T $ca $sh [0,200,1] [-0.15,0.15] [0,0.125]
#done

# contact angle study at constant heat flux
sh=-0.02
for ca in 135.0 #22.5 45.0 67.5 90.0 112.5 135.0
do
    generate_non_st_contours $T $ca $sh [0,180,1] [-0.15,0.15] [0,0.125]
done
