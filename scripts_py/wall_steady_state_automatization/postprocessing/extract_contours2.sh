#!/bin/sh

get_inflow_sph_dir(){
    projectDir=/home/jdesmarais/projects
    dir="$projectDir/dim2d_"$1"_ca"$2"_vap_v"$3"_hca"$4"_sph"
    echo "$dir"
}

generate_inflow_sph_contours(){
    pythonDir=/home/jdesmarais/Code/augeanstables/scripts_py

    cd $pythonDir
    cd wall_steady_state_automatization/contours_extraction

    T=$1
    ca=$2
    v=$3
    hca=$4

    inputDir=$( get_inflow_sph_dir $T $ca $v $hca)

    
    #options="-i $inputDir -c $ca -t $5 -x $6 -y $7 -p -l max_gradient -g -w"
    options="-i $inputDir -c $ca -t $5 -x $6 -y $7 -p -l max_gradient -s"

    ./extract_bubble_contours.py  $options #>cur_contours.out 2>&1


    echo "=================================================="
    echo "generate_contours( "$T" , "$ca" , "$v", "$hca" ) : done"
    echo "=================================================="
    echo ''
}


T=0.95
v=0.1
ca=135.0
hca=0.0

generate_inflow_sph_contours $T $ca $v $hca [0,500,2] [-0.15,0.375] [0,0.25]
