#!/bin/bash

source $augeanstables/scripts_py/wall_steady_state_automatization/postprocessing/library_extract_contours.sh

true="true"
false="false"


steadyState_Homogeneous=$false    # 1) steady state on homogeneous surface

nucleation_Homogeneous=$false     # 2) nucleation on homogenenous surface

detachment_ParabolicInflow=$false # 3) detachment induced by parabolic inflow

detachment_LinearInflow=$true     # 4) detachment induced by linear inflow

nucleation_ParabolicInflow=$false # 5) nucleation parabolic velocity profile

nucleation_LinearInflow=$false    # 6) nucleation linear velocity profile


#============================================================
# 1) steady state on homogeneous surface
#============================================================
if [ "$steadyState_Homogeneous" = "$true" ]
then

    T=0.95

    for ca in 22.5 45.0 67.5 90.0 112.5 135.0
    do
	generate_st_contours $T $ca [0,125,1] [-0.15,0.15] [0,0.15] 3
    done

fi


#============================================================
# 2) nucleation on homogeneous surface
#============================================================
if [ "$nucleation_Homogeneous" = "$true" ]
then

    T=0.95

    # constant contact angle, varying heat flux
    ca=90.0
    
    for fh in 0.04 0.06 0.08 0.1
    do
	options="$T $ca $fh [0,215,1] [-0.15,0.15] [0,0.125]"
	generate_nucleation_contours $options
    done

    # constant flux, varying contact angle
    fh=0.02
    data=[0,170,1]
    xlimits=[-0.15,0.15]
    ylimits=[0,0.125]
    
    for ca in 22.5 45.0 67.5 90.0 112.5 135.0
    do
    	options="$T $ca $fh $data $xlimits $ylimits"
    	generate_nucleation_contours $options
    done
    
fi


#============================================================
# 3) detachment with incoming parabolic velocity profile
#============================================================
if [ "$detachment_ParabolicInflow" = "$true" ]
then

    T=0.95
    hca=0.0

    for v in 0.1 0.2 0.3 0.4 0.5
    do
	for ca in 22.5 45.0 67.5 90.0 112.5 135.0
	do
	    options="$T $ca $v $hca [0,500,2] [-0.15,0.80] [0,0.25]"
	    generate_inflow_sph_contours $options
	done
    done

fi


#============================================================
# 4) detachment with incoming linear velocity profile
#============================================================
if [ "$detachment_LinearInflow" = "$true" ]
then

    T=0.95
    hca=0.0
    nbContours=7
    xFigsize=12
    yFigsize=2

    for v in 0.4 #0.1 0.2 0.3 0.4 0.5
    do
	for ca in 135.0 #22.5 45.0 67.5 90.0 112.5 135.0
	do
	    options="$T $ca $v $hca [0,502,2] [-0.15,0.80] [0,0.2] $nbContours $xFigsize $yFigsize"
	    generate_linear_inflow_sph_contours $options
	done
    done

fi


#============================================================
# 5) nucleation with incoming parabolic velocity profile
#============================================================
if [ "$nucleation_ParabolicInflow" = "$true" ]
then

    T=0.95
    ca=22.5
    hca=0.0
    xlimits=[-0.1,0.425]
    ylimits=[0,0.125]

    # constant heat flux, varying maximum velocity
    #fh=0.04    
    #
    #for v in 0.1 0.2 0.3 0.4 0.5
    #do
    #	options="$T $ca $fh $v $hca [0,500,1] $xlimits $ylimits"
    #	generate_inflow_nucleation_contours $options
    #done

    # constant maximum velocity, varying heat flux
    v=0.4
    
    for fh in 0.06 #0.05 0.06 0.08 0.1
    do
    	options="$T $ca $fh $v $hca [0,500,1] $xlimits $ylimits"
    	generate_inflow_nucleation_contours $options
    done
    
fi


#============================================================
# 6) nucleation with incoming linear velocity profile
#============================================================
if [ "$nucleation_LinearInflow" = "$true" ]
then

    T=0.95
    ca=22.5
    hca=0.0
    xlimits=[-0.1,0.425]
    ylimits=[0,0.135]
    tlimits=[0,325,1]

    # constant heat flux, varying maximum velocity
    #fh=0.04
    #
    #for v in 0.1 0.2 0.3 0.4 0.5
    #do
    #	options="$T $ca $fh $v $hca $tlimits $xlimits $ylimits"
    #	generate_linear_inflow_nucleation_contours $options
    #done

    # constant maximum velocity, varying heat flux
    v=0.4

    for fh in 0.04 0.05 0.06 0.08
    do
	options="$T $ca $fh $v $hca $tlimits $xlimits $ylimits"
	generate_linear_inflow_nucleation_contours $options
    done
    
fi
