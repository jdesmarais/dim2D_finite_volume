#!/bin/bash

source $augeanstables/scripts_py/wall_steady_state_automatization/contours_extraction/library_extract_contours.sh


#============================================================
# extraction of the contours for the steady state simulations
#============================================================
T=0.999

# 22.5 45.0 67.5 90.0
for ca in 22.5 45.0 67.5 90.0
do
    generate_st_contours $T $ca [0,1010,10] [-0.4,0.4] [0,0.3]
done

# 112.5 135.0
for ca in 112.5 135.0
do
    generate_st_contours $T $ca [0,2424,24] [-0.4,0.4] [0,0.3]
done

