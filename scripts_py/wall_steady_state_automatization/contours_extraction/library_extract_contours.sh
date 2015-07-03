#!/bin/bash

pythonDir=/home/jdesmarais/Code/augeanstables/scripts_py

# directory for the simulation at (T,ca)
get_st_dir(){
    projectDir=/home/jdesmarais/projects
    dir="$projectDir/dim2d_"$1"_ca"$2"_vap"
    echo "$dir"
}

# generate the contours for steady state simulations
generate_st_contours(){

    cd $pythonDir
    cd wall_steady_state_automatization/contours_extraction

    inputDir=$( get_st_dir $1 $2)
    ca=$2

    options="-i $inputDir -c $ca -t $3 -p -g -r -w -l max_gradient"
    #options="-i $inputDir -c $ca -t $3 -r -l max_gradient -s -x $4 -y $5"

    ./extract_bubble_contours.py $options #> cur_st_contours.out 2>&1

    #mv cur_st_contours.out $inputDir/contours

    echo "=================================================="
    echo "generate_st_contours( "$T" , "$ca" ) : done"
    echo "=================================================="
    echo ''
}


# directory for the nucleation simulation on
# uniform surface at (T,ca,sh)
get_nucleation_dir(){
    projectDir=/home/jdesmarais/projects/
    projectDir+=20150624_dim2d_0.95_ca22.5-135.0_vap_sh-0.02-0.1/
    dir="$projectDir/dim2d_"$1"_ca"$2"_vap_sh"$3
    echo "$dir"
}


# generate the contours for the uniform surface nucleation simulations
generate_nucleation_contours(){

    cd $pythonDir
    cd wall_steady_state_automatization/contours_extraction

    inputDir=$( get_nucleation_dir $1 $2 $3 )
    ca=$2

    #options="-i $inputDir -c $ca -t $4 -x $5 -y $6 -p -l wall_max_gradient -g -w -r"
    options="-i $inputDir -c $ca -t $4 -x $5 -y $6 -p -l wall_max_gradient -s -r"

    ./extract_bubble_contours.py $options #-i $inputDir -c $ca -t $4 -p -g -w #>cur_contours.out 2>&1
    #./extract_bubble_contours.py -i $inputDir -c $ca -t $4 -p -x $5 -y $6 -s #> cur_contours.out 2>&1

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


# directory for the nucleation simulation on
# uniform surface at (T,ca,sh)
get_linear_inflow_sph_dir(){
    projectDir=/home/jdesmarais/projects
    dir="$projectDir/dim2d_"$1"_ca"$2"_vap_vl"$3"_hca"$4"_sph"
    echo "$dir"
}


# generate the contours for the uniform surface nucleation simulations
generate_linear_inflow_sph_contours(){

    cd $pythonDir
    cd wall_steady_state_automatization/contours_extraction

    T=$1
    ca=$2
    v=$3
    hca=$4

    inputDir=$( get_linear_inflow_sph_dir $T $ca $v $hca)

    
    options="-i $inputDir -c $ca -t $5 -x $6 -y $7 -p -l max_gradient -g -w"
    #options="-i $inputDir -c $ca -t $5 -x $6 -y $7 -p -l max_gradient -s"

    ./extract_bubble_contours.py  $options #>cur_contours.out 2>&1


    echo "=================================================="
    echo "generate_contours( "$T" , "$ca" , "$v", "$hca" ) : done"
    echo "=================================================="
    echo ''
}


# directory for the nucleation simulation on
# uniform surface at (T,ca,sh)
get_inflow_nucleation_dir(){
    projectDir=/home/jdesmarais/projects
    dir="$projectDir/dim2d_"$1"_ca"$2"_vap_sh"$3"_v"$4"_hca"$5
    echo "$dir"
}


# generate the contours for the uniform surface nucleation simulations
generate_inflow_nucleation_contours(){

    cd $pythonDir
    cd wall_steady_state_automatization/contours_extraction

    T=$1
    ca=$2
    sh=$3
    v=$4
    hca=$5

    inputDir=$( get_inflow_nucleation_dir $T $ca $sh $v $hca)

    
    options="-i $inputDir -c $ca -t $6 -x $7 -y $8 -p -l max_gradient -g -w"
    #options="-i $inputDir -c $ca -t $6 -x $7 -y $8 -p -l max_gradient -s"

    ./extract_bubble_contours.py  $options #>cur_contours.out 2>&1


    echo "======================================================="
    echo "generate_contours( "$T" , "$ca" , "$sh", "$v", "$hca" ) : done"
    echo "======================================================="
    echo ''
}


# directory for the nucleation simulation on
# uniform surface at (T,ca,sh)
get_linear_inflow_nucleation_dir(){
    projectDir=/home/jdesmarais/projects
    dir="$projectDir/dim2d_"$1"_ca"$2"_vap_sh"$3"_vl"$4"_hca"$5
    echo "$dir"
}


# generate the contours for the uniform surface nucleation simulations
generate_linear_inflow_nucleation_contours(){

    cd $pythonDir
    cd wall_steady_state_automatization/contours_extraction

    T=$1
    ca=$2
    sh=$3
    v=$4
    hca=$5

    inputDir=$( get_linear_inflow_nucleation_dir $T $ca $sh $v $hca)

    
    #options="-i $inputDir -c $ca -t $6 -x $7 -y $8 -p -l max_gradient -g -w"
    options="-i $inputDir -c $ca -t $6 -x $7 -y $8 -p -l max_gradient -s -e"

    ./extract_bubble_contours.py  $options #>cur_contours.out 2>&1


    echo "======================================================="
    echo "generate_contours( "$T" , "$ca" , "$sh", "$v", "$hca" ) : done"
    echo "======================================================="
    echo ''
}
