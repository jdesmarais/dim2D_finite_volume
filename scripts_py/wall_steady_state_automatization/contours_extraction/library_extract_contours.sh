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

    ./extract_bubble_contours.py -i $inputDir -c $ca -t $3 -p -g -w > cur_st_contours.out 2>&1
    #./extract_bubble_contours.py -i $inputDir -c $ca -t $3 -s #> cur_st_contours.out 2>&1

    mv cur_st_contours.out $inputDir/contours

    echo "generate_st_contours( "$T" , "$ca" ) : done"
}


# directory for the simulation at (T,ca,sh)
get_non_st_dir(){
    projectDir=/home/jdesmarais/projects
    dir="$projectDir/dim2d_"$1"_ca"$2"_vap_sh"$3
    echo "$dir"
}


# generate the contours for the non-steady state simulations
generate_non_st_contours(){
    pythonDir=/home/jdesmarais/Code/augeanstables/scripts_py

    cd $pythonDir
    cd wall_steady_state_automatization/contours_extraction

    inputDir=$( get_non_st_dir $1 $2 $3 )
    ca=$2

    #./extract_bubble_contours.py -i $inputDir -c $ca -t $4 -p -g -w #>cur_contours.out 2>&1
    ./extract_bubble_contours.py -i $inputDir -c $ca -t $4 -p -x $5 -y $6 -s #> cur_contours.out 2>&1

    #mv cur_contours.out $inputDir/contours

    echo "=================================================="
    echo "generate_contours( "$1" , "$2" , "$3" ) : done"
    echo "=================================================="
    echo ''
}