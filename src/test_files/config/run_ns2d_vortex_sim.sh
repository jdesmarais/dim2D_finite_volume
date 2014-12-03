#!/bin/bash

#./default_inputs/wave1d_hedstrom_x_reflection_y.txt
#./default_inputs/wave2d_hedstrom_x_reflection_y.txt
#./default_inputs/wave2d_hedstrom_xy.txt
#./default_inputs/wave2d_hedstrom_xy_corners.txt
#./default_inputs/ns2d_vortex_hedstrom_detailled.txt
#./default_inputs/ns2d_vortex_hedstrom_corners_detailled.txt
#./default_inputs/ns2d_peak_poinsot.txt
#./default_inputs/ns2d_vortex_poinsot.txt
#./default_inputs/ns2d_vortex_poinsot_detailled.txt
#./default_inputs/ns2d_peak_yoolodato.txt
#./default_inputs/ns2d_vortex_yoolodato.txt
#./default_inputs/ns2d_vortex_yoolodato_detailled.txt

#settings
INPUT=./default_inputs/ns2d_vortex_yoolodato_detailled.txt
DEST_FOLDER=~/projects/ns2d_vortex_yoolodato

#compile code
./config.py -i $INPUT -c

#move the executable to the corresponding project folder
#mv ../sim_dim2d $DEST_FOLDER
#
##move to project folder
#cd $DEST_FOLDER
#rm *.nc
#./sim_dim2d
#
##plot results
#exit

