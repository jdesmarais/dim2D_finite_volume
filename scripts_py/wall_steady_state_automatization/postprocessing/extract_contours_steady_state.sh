#!/bin/bash

source $augeanstables/scripts_py/wall_steady_state_automatization/contours_extraction/library_extract_contours.sh


#============================================================
# extraction of the contours for the steady state simulations
#============================================================
T=0.95

# 22.5 45.0 67.5 90.0
for ca in 22.5 45.0 67.5 90.0 112.5 135.0
do
    generate_st_contours $T $ca [0,125,1] [-0.15,0.15] [0,0.15] 3
done

## 112.5 135.0
#for ca in 112.5 135.0
#do
#    generate_st_contours $T $ca [0,5300,50] [-0.4,0.4] [0,0.3]
#done

## 22.5 45.0 67.5 90.0
#ca=135.0
#generate_st_contours $T $ca [0,5300,50] [-0.4,0.4] [0,0.3]