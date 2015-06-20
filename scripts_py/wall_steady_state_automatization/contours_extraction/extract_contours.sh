#!/bin/bash

# directory for the simulation at (T,ca)
get_st_dir(){
    projectDir=/home/jdesmarais/projects
    dir="$projectDir/dim2d_"$1"_ca"$2"_vap"
    echo "$dir"
}

# generate the contours for steady state simulations
generate_st_contours(){
    pythonDir=/home/jdesmarais/Code/augeanstables/scripts_py

    cd $pythonDir
    cd wall_steady_state_automatization/contours_extraction

    inputDir=$( get_st_dir $1 $2)
    ca=$2

    ./extract_bubble_contours.py -i $inputDir -c $ca -t $3 -p -g > cur_st_contours.out 2>&1
    #./extract_bubble_contours.py -i $inputDir -c $ca -t $3 -s > cur_st_contours.out 2>&1

    echo "generate_st_contours( "$T" , "$ca" ) : done"
}


# directory for the simulation at (T,ca,sh)
get_non_st_dir(){
    projectDir=/home/jdesmarais/projects
    dir="$projectDir/dim2d_"$1"_ca"$2"_vap_sh"$3
    echo "$dir"
}


# generate the contours for the non-steady state simulations
generate_contours(){
    pythonDir=/home/jdesmarais/Code/augeanstables/scripts_py

    cd $pythonDir
    cd wall_steady_state_automatization/contours_extraction

    inputDir=$( get_non_st_dir $1 $2 $3 )
    ca=$2

    ./extract_bubble_contours.py -i $inputDir -c $ca -t $4 -p -g > cur_contours.out 2>&1

    echo "generate_contours( "$1" , "$2" , "$3" ) : done"
}


#============================================================
# extraction of the contours for the steady state simulations
#============================================================
T=0.999

#22.5 45.0 67.5 90.0 112.5 130.0
for ca in 22.5 45.0 67.5 90.0
do
    generate_st_contours $T $ca [0,1010,10]
done


for ca in 112.5 135.0
do
    generate_st_contours $T $ca [0,2424,24]
done


#============================================================
# extraction of the contours for the nucleation simulations
#============================================================
T=0.95

# heat flux study at constant contact angle
ca=90.0
for sh in -0.04 -0.06 -0.08 -0.1
do
    generate_non_st_contours $T $ca $sh [0,200,2]
done

# contact angle study at constant heat flux
sh=-0.02
for ca in 22.5 45.0 67.5 90.0 112.5 135.0
do
    generate_non_st_contours $T $ca $sh [0,200,2]
done

