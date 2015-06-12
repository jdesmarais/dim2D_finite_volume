#!/usr/bin/python

'''
@description
create one input file to run the DIM2D simulation for 
the collapse of a bubble at a wall:
inputs_wall_st.txt
'''

import os
import sys
import inspect
import getopt
import subprocess
import shlex
import shutil

debug=True


#add the python files from sm_lg_automatization
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],"../sm_lg_domain_automatization")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from automatization_csts import (nb_pts_in_interface_default,
                                 ratio_bubble_interface_default,
                                 CFL_constant_default)

from automatization_wall_st_csts import (ratio_eq_length_domain_default,
                                         ratio_drop_length_domain_default,
                                         total_nb_files_default)

from library_sm_lg_inputs import (get_we,
                                  get_interface_length,
                                  get_bubble_diameter,
                                  get_interface_space_step,
                                  get_bubble_diameter,
                                  get_cv_r,
                                  get_max_speed_of_sound,
                                  get_dt_max,
                                  get_detail_print)

from library_wall_st_inputs import get_wall_domain_extent

from create_sm_lg_inputs import (get_parameter,
                                 create_inputFile)


# determine the inputs modified in the input file
def get_inputsToBeModified(
    steady_state_ac,
    temperature,
    micro_contact_angle,
    phase_at_center,
    gravity_ac,
    gravity_amp,
    nb_pts_in_interface,
    ratio_bubble_interface,
    CFL_constant,
    total_nb_files):

    '''
    @description
    determine the inputs to be modified in the template.txt
    template input file to run the DIM2D simulation for a drop/bubble
    collapse on a wall
    '''

    # extract length_c, dim2d_a, dim2d_b, dim2d_M, dim2d_cv, dim2d_R
    # and dim2d_K from the dim2d_parameters.f fortran file
    dim2dParamPath = os.path.join(os.getenv('augeanstables'),
                                  'src',
                                  'physical_models',
                                  'dim2d',
                                  'dim2d_parameters.f')

    if(os.path.isfile(dim2dParamPath)):
        length_c  = float(get_parameter('length_c', dim2dParamPath))
        dim2d_a   = float(get_parameter( 'dim2d_a', dim2dParamPath))
        dim2d_b   = float(get_parameter( 'dim2d_b', dim2dParamPath))
        dim2d_M   = float(get_parameter( 'dim2d_M', dim2dParamPath))
        dim2d_cv  = float(get_parameter('dim2d_cv', dim2dParamPath))
        dim2d_R   = float(get_parameter( 'dim2d_R', dim2dParamPath))
        dim2d_K   = float(get_parameter( 'dim2d_K', dim2dParamPath))

    else:
        sys.exit('*** '+dim2dParamPath+' does not exist***')


    # compute the Weber number
    we = get_we(length_c, dim2d_a, dim2d_b, dim2d_M, dim2d_K)
    if(debug): print 'we: ', we

    # compute the interface length from the temperature
    interface_lgh = get_interface_length(we,temperature)
    if(debug): print 'interface_length: ', interface_lgh

    # compute the bubble diameter from the interface length
    bubble_diameter = get_bubble_diameter(interface_lgh,
                                          2.0)
    if(debug): print 'bubble_diameter: ', bubble_diameter
    

    # compute the x_max of the domain
    x_max     = 4.0*bubble_diameter
    if(debug): print 'x_max: ', x_max

    # compute the y_max of the domain
    y_max     = x_max
    if(debug): print 'y_max: ', y_max

    # compute the maximum space step from the interface length
    dx_max = get_interface_space_step(interface_lgh,
                                      nb_pts_in_interface)
    if(debug): print 'dx_max: ', dx_max

    # compute the extent of the domain as a matrix
    # x_min : domain_extent[0][0]
    # x_max : domain_extent[1][0]
    # y_min : domain_extent[0][1]
    # y_max : domain_extent[1][1]
    domain_extent = get_wall_domain_extent(x_max,y_max,dx_max)
    if(debug): print 'domain_extent: ', domain_extent

    # compute the reduced heat capacity
    cv_r = get_cv_r(dim2d_M,dim2d_cv,dim2d_R)

    # compute the maximum speed of sound in the flow
    speed_of_sound = get_max_speed_of_sound(temperature,cv_r)

    # compute the maximum time step ensuring numerical stability
    flow_velocity = 0.0
    speed_max     = speed_of_sound
    if(debug): print 'speed_of_sound: ', speed_of_sound

    dt_max        = get_dt_max(dx_max,speed_max,CFL_constant,precision_c=6)
    if(debug): print 'dt_max: ', dt_max

    # determine the maximum simulation time
    simulation_time = 6.0

    # determine the detail print
    detail_print = get_detail_print(total_nb_files,
                                    simulation_time,
                                    dt_max)
    if(debug): print 'detail_print: ', detail_print


    # gather the inputs to be modified in a dictionnary
    inputsToBeModified = {
        'detail_print'                : detail_print,
        'dt'                          : dt_max,
        't_max'                       : simulation_time,
        'steady_state_ac'             : steady_state_ac,
        'dx'                          : dx_max,
        'x_min'                       : domain_extent[0][0],
        'x_max'                       : domain_extent[1][0],
        'dy'                          : dx_max,
        'y_min'                       : domain_extent[0][1],
        'y_max'                       : domain_extent[1][1],
        'flow_velocity'               : flow_velocity,
        'temperature'                 : temperature,
        'micro_contact_angle'         : micro_contact_angle,
        'phase_at_center'             : phase_at_center,
        'ratio_bubble_interface'      : ratio_bubble_interface,
        'gravity_ac'                  : gravity_ac,
        'gravity_amp'                 : gravity_amp}

    return inputsToBeModified


# create the input file for the wall DIM2D simulation
def create_wall_bubblecollapse_inputs(
    steady_state_ac,
    temperature,
    micro_contact_angle,
    phase_at_center,
    gravity_ac,
    gravity_amp,
    model_input,
    inputs_wall_modified      = 'inputs_wall.txt',
    nb_pts_in_interface       = nb_pts_in_interface_default,
    ratio_bubble_interface    = 2.0,
    CFL_constant              = CFL_constant_default,
    total_nb_files            = total_nb_files_default):
    '''
    @description:
    create the small and large domain inputs for the simulation    
    '''

    # determine the inputs to be modified in
    # template.txt to run the simulation on
    # small and large scale domains
    inputs_wall = get_inputsToBeModified(
        steady_state_ac,
        temperature,
        micro_contact_angle,
        phase_at_center,
        gravity_ac,
        gravity_amp,
        nb_pts_in_interface,
        ratio_bubble_interface,
        CFL_constant,
        total_nb_files)
    

    # create the input file for the small
    # domain simulation
    create_inputFile(inputs_wall,
                     model_input,
                     inputs_wall_modified)

    return
 
