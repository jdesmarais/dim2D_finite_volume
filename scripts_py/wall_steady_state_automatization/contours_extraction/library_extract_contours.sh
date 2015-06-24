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

    echo "=================================================="
    echo "generate_st_contours( "$T" , "$ca" ) : done"
    echo "=================================================="
    echo ''
}


# directory for the nucleation simulation on
# uniform surface at (T,ca,sh)
get_nucleation_dir(){
    projectDir=/home/jdesmarais/projects
    dir="$projectDir/dim2d_"$1"_ca"$2"_vap_sh"$3
    echo "$dir"
}


# generate the contours for the uniform surface nucleation simulations
generate_nucleation_contours(){
    pythonDir=/home/jdesmarais/Code/augeanstables/scripts_py

    cd $pythonDir
    cd wall_steady_state_automatization/contours_extraction

    inputDir=$( get_nuclation_dir $1 $2 $3 )
    ca=$2

    #./extract_bubble_contours.py -i $inputDir -c $ca -t $4 -p -g -w #>cur_contours.out 2>&1
    ./extract_bubble_contours.py -i $inputDir -c $ca -t $4 -p -x $5 -y $6 -s #> cur_contours.out 2>&1

    #mv cur_contours.out $inputDir/contours

    echo "=================================================="
    echo "generate_contours( "$1" , "$2" , "$3" ) : done"
    echo "=================================================="
    echo ''
}


# directory for the nucleation simulation on
# uniform surface at (T,ca,sh)
get_inflow_sph_dir(){
    projectDir=/home/jdesmarais/projects
    dir="$projectDir/dim2d_"$1"_ca"$2"_vap_v"$3"_hca"$4"_sph"
    echo "$dir"
}


# generate the contours for the uniform surface nucleation simulations
generate_inflow_sph_contours(){
    pythonDir=/home/jdesmarais/Code/augeanstables/scripts_py

    cd $pythonDir
    cd wall_steady_state_automatization/contours_extraction

    T=$1
    ca=$2
    v=$3
    hca=$4

    inputDir=$( get_inflow_sph_dir $T $ca $v $hca)


    ./extract_bubble_contours.py -i $inputDir -c $ca -t $5 -p -g -w -x $6 -y $7 #>cur_contours.out 2>&1
    #./extract_bubble_contours.py -i $inputDir -c $ca -t $5 -p -x $6 -y $7 -s #> cur_contours.out 2>&1

    #mv cur_contours.out $inputDir/contours

    echo "=================================================="
    echo "generate_contours( "$T" , "$ca" , "$v", "$hca" ) : done"
    echo "=================================================="
    echo ''
}