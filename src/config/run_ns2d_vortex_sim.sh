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

main_dir=$augeanstables/src/main

#settings
#INPUT=./default_inputs/dim2d_bubble_transported_hedstrom_xy.txt #_detailled.txt
#DEST_FOLDER=~/projects/dim2d_bubble_trs_hedstrom_xy_bf
INPUT=./default_inputs/wave2d/wave2d_hedstrom_xy_intermittentforcing.txt #_detailled.txt
DEST_FOLDER=~/projects/test_wave2d_bf
EXE=sim_dim2d_bf

#compile code
#./config.py -i $INPUT -c
./config.py -i $INPUT -c -b #buffer layers

#move executable to the project folder
mv $main_dir/$EXE $DEST_FOLDER

#run simulation
qsub scripts_pbs/run_sim_dim2d_bf.job

#watch simulation
watch qstat

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

