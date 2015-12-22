#!/usr/bin/python

"""
Create the input file used to configure a simulation with two bubbles
transported by the flow. The distance between the two bubbles transported,
the temperature and the mean flow velocity will be varied
"""

import os
import sys
import inspect

#add the python files from ../sm_lg_automatization
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],
    "../sm_lg_domain_automatization")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from automatization_csts import\
    dct_distance_default,\
    md_threshold_ac_default,\
    md_threshold_default,\
    nb_pts_inside_interface_required,\
    ratio_bubble_interface_default,\
    CFL_constant_default,\
    ratio_interface_influence_default,\
    total_nb_files_default

from library_sm_lg_inputs import\
    get_we,\
    get_interface_length,\
    get_interface_space_step,\
    get_cv_r,\
    get_max_speed_of_sound,\
    get_dt_max,\
    get_detail_print

from create_sm_lg_inputs import\
    get_parameter,\
    create_inputFile

# libraries from this folder
from library_bb_tr_inputs import\
    get_domain_sizes,\
    get_domain_extent,\
    get_simulation_time
    

# compute the inputs needed to run the simulation
# with two bubbles transported by the flow
def compute_inputsToBeModified(temperature,
                               flow_velocity,
                               nb_pts_in_interface,
                               ratio_interface_separation,
                               ratio_bubble_interface,
                               ratio_interface_influence,
                               CFL_constant,
                               total_nb_files,
                               dct_distance,
                               md_threshold_ac,
                               md_threshold,
                               flow_direction='x'):

    """Determine the inputs to be modified in the input.txt
    file for the simulation with two bubbles transported by the
    mean flow

    Args:
        temperature (double) : mean flow temperature
        flow_velocity (double) : mean flow velocity
        nb_pts_in_interface (int) : number of grid-points needed
            to resolve the interface profile
        ratio_interface_separation (double) : length (expressed as
            a fraction of the width of the interface) separating the
           two bubbles
        ratio_bubble_interface (double) : ratio between bubble
           diameter and interface width
        ratio_interface_influence (double) : length threshold
           (expressed as a fraction of the width of the interface)
           above which the interface is not supposed to interact
           with the border
        CFL_constant (double) : CFL threshold (0-1)
        total_nb_files (int) : total nb of files written
        dct_distance (int) : distance between the detectors and
           the boundary expressed as a number of gridpoints
        md_threshold_ac (int) : activate the increase of the computational
           domain with a mass density threshold
        md_threshold (double) : mass density threshold to activate the
           increase of the computational domain (between 0 and 1)

    Returns:
        inputsToBeModified (dict) :
           - dict[key] : name of the input to be modified
           - dict[value] : value of the input to be modified
    """

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


    # compute the interface length from the temperature
    interface_length = get_interface_length(we,temperature)

    
    # compute the bubble diameter
    bubble_diameter = interface_length*ratio_bubble_interface

    
    # compute the extent of the domain
    [Lx,Ly] = get_domain_sizes(bubble_diameter,
                               interface_length,
                               ratio_interface_influence,
                               ratio_interface_separation)


    # compute the maximum space step from the interface length
    dx_max = get_interface_space_step(interface_length,
                                      nb_pts_in_interface)


    # get the extent of the domain
    # x_min : domain_extent[0][0]
    # x_max : domain_extent[1][0]
    # y_min : domain_extent[0][1]
    # y_max : domain_extent[1][1]
    domain_extent = get_domain_extent(Lx,Ly,dx_max,dx_max)

    
    # compute the reduced heat capacity
    cv_r = get_cv_r(dim2d_M,dim2d_cv,dim2d_R)


    # compute the maximum speed of sound in the flow
    speed_of_sound = get_max_speed_of_sound(temperature,cv_r)


    # compute the maximum time step ensuring numerical stability
    speed_max = speed_of_sound + abs(flow_velocity)
    dt_max    = get_dt_max(dx_max,speed_max,CFL_constant)


    # determine the simulation time needed to let the two bubbles
    # leave the computational domain
    simulation_time = get_simulation_time(Lx,
                                          bubble_diameter,
                                          ratio_interface_separation*interface_length,
                                          ratio_interface_influence*interface_length,
                                          flow_velocity)


    # determine the detail print
    detail_print = get_detail_print(total_nb_files,
                                    simulation_time,
                                    dt_max)


    # gather the inputs to be modified in a dictionnary
    inputsToBeModified = {
        'detail_print'                : detail_print,
        'dt'                          : dt_max,
        't_max'                       : simulation_time,
        'dx'                          : dx_max,
        'dy'                          : dx_max,
        'flow_velocity'               : flow_velocity,
        'temperature'                 : temperature,
        'li_separation'               : ratio_interface_separation,
        'openbc_detector_distance'    : dct_distance,
        'openbc_md_threshold_ac'      : md_threshold_ac,
        'openbc_md_threshold'         : md_threshold}

    if(flow_direction=='x'):

        inputsToBeModified.update({'x_min'          : domain_extent[0][0],
                                   'x_max'          : domain_extent[1][0],
                                   'y_min'          : domain_extent[0][1],
                                   'y_max'          : domain_extent[1][1],
                                   'flow_direction' : 'E'})
    else:

        if(flow_direction!='y'):
            sys.exit('create_bb_tr_input: '+
                     'compute_inputsToBeModified: '+
                     'flow_direction not recognized')
            sys.exit(2)

        inputsToBeModified.update({'x_min'          : domain_extent[0][1],
                                   'x_max'          : domain_extent[1][1],
                                   'y_min'          : domain_extent[0][0],
                                   'y_max'          : domain_extent[1][0],
                                   'flow_direction' : 'N'})

    return inputsToBeModified
    
    
# create the input file for the simulation with two bubbles
# transported by the flow
def create_bb_tr_inputs(temperature,
                        flow_velocity,
                        ratio_interface_separation,
                        model_input,
                        flow_direction            = 'x',
                        output_file               = 'inputs.txt',
                        dct_distance              = dct_distance_default,
                        md_threshold_ac           = md_threshold_ac_default,
                        md_threshold              = md_threshold_default,
                        nb_pts_in_interface       = nb_pts_inside_interface_required,
                        ratio_bubble_interface    = ratio_bubble_interface_default,
                        CFL_constant              = CFL_constant_default,
                        ratio_interface_influence = ratio_interface_influence_default,
                        total_nb_files            = total_nb_files_default):
    """Create the text input file to run the simulation with
    two bubbles transported by the mean flow

    Args:
        temperature (double) : mean flow temperature
        flow_velocity (double) : mean flow velocity
        ratio_interface_separation (double) : length (expressed as
            a fraction of the width of the interface) separating the
           two bubbles
        model_input (string) : path for the input file whose input
            parameters are modified for the simulation of two bubbles
            transported by the flow
        flow_direction (string) : flow direction (either 'x' or 'y')
        output_file (string) : path for the output file
        dct_distance (int) : distance between the detectors and
           the boundary expressed as a number of gridpoints
        md_threshold_ac (int) : activate the increase of the computational
           domain with a mass density threshold
        md_threshold (double) : mass density threshold to activate the
           increase of the computational domain (between 0 and 1)
        nb_pts_in_interface (int) : number of grid-points needed
            to resolve the interface profile
        ratio_bubble_interface (double) : ratio between bubble
           diameter and interface width
        CFL_constant (double) : CFL threshold (0-1)
        ratio_interface_influence (double) : length threshold
           (expressed as a fraction of the width of the interface)
           above which the interface is not supposed to interact
           with the border
        total_nb_files (int) : total nb of files written
        
    Returns:
        None
    """

    # determine the inputs to be modified in
    # template.txt to run the simulation
    inputs = compute_inputsToBeModified(temperature,
                                        flow_velocity,
                                        nb_pts_in_interface,
                                        ratio_interface_separation,
                                        ratio_bubble_interface,
                                        ratio_interface_influence,
                                        CFL_constant,
                                        total_nb_files,
                                        dct_distance,
                                        md_threshold_ac,
                                        md_threshold,
                                        flow_direction=flow_direction)
    

    # create the input file for the simulation
    create_inputFile(inputs,
                     model_input,
                     output_file)

    return


if __name__=='__main__':
    
    temperature = 0.99
    flow_velocity = 0.1
    ratio_interface_separation = 3.0
    model_input = os.path.join(os.getenv('augeanstables'),
                               'src',
                               'config',
                               'default_inputs',
                               'dim2d',
                               'dim2d_bubbles_transported_hedstrom_xy.txt')

    create_bb_tr_inputs(temperature,
                        flow_velocity,
                        ratio_interface_separation,
                        model_input,
                        flow_direction='y')
