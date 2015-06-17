#!/usr/bin/python

'''
@description
script used to generate the main results of the
bubble heat exchanges closed to a wall

'''

import os
import sys
import inspect

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],"../sm_lg_domain_automatization")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from library_wall_nonst_results import generate_wall_nonst_results

from library_sm_lg_inputs import (get_we,
                                  get_interface_length)

from create_sm_lg_inputs import get_parameter


def extract_interface_length(dim2dParamPath,temperature):

    if(os.path.isfile(dim2dParamPath)):
        length_c  = float(get_parameter('length_c', dim2dParamPath))
        dim2d_a   = float(get_parameter( 'dim2d_a', dim2dParamPath))
        dim2d_b   = float(get_parameter( 'dim2d_b', dim2dParamPath))
        dim2d_M   = float(get_parameter( 'dim2d_M', dim2dParamPath))
        dim2d_K   = float(get_parameter( 'dim2d_K', dim2dParamPath))

    else:
        sys.exit('*** '+dim2dParamPath+' does not exist***')

    # compute the Weber number
    we = get_we(length_c, dim2d_a, dim2d_b, dim2d_M, dim2d_K)

    # compute the interface length from the temperature
    interface_lgh = get_interface_length(we,temperature)

    return interface_lgh


if __name__=="__main__":

    #file with the DIM2d parameters
    dim2dParamPath = os.path.join(os.getenv('augeanstables'),
                                  'src',
                                  'physical_models',
                                  'dim2d',
                                  'dim2d_parameters.f')

    # main directory where the simulations
    # are saved
    mainDir='/home/jdesmarais/projects'


    # input used as template
    #------------------------------------------------------------
    maindir_input = os.path.join(os.getenv('augeanstables'),
                                 'src','config',
                                 'default_inputs','dim2d')

    model_input = os.path.join(
        maindir_input,
        'dim2d_bubble_nucleation_at_wall.txt')

    conductionHeatStudy = False
    sourceHeatStudy     = False
    contactAngleStudy   = False
    flowVelocityStudy   = True


    #1) conduction heat study
    if(conductionHeatStudy):

        simulationDuration  = 100
        steady_state_ac     = 0
        temperature         = 0.95
        contact_angle       = 90.0
        phase_at_center     = 'vapor'
        ratio               = 2.0
        gravity_ac          = 0
        gravity             = 0.000
        heat_source_choice  = 'gaussian_heat_source'
        max_heat_flux_array = [0.0001,0.001,0.01]
        heat_source_center  = 0.0
        heat_source_variance = 2.0*extract_interface_length(dim2dParamPath,temperature)
        total_nb_files      = 500

    
        for max_heat_flux in max_heat_flux_array:
            
            PBSnameRun = 'dim2d_'+str(temperature)+'_fh'+str(max_heat_flux)

            [destDir, nameRun] =\
                \
                generate_wall_nonst_results(
                mainDir,
                model_input,
                PBSnameRun,
                simulationDuration,
                steady_state_ac           = steady_state_ac,
                temperature               = temperature,
                micro_contact_angle       = contact_angle,
                phase_at_center           = phase_at_center,
                ratio_bubble_interface    = ratio,
                gravity_ac                = gravity_ac,
                gravity_amp               = gravity,
                wall_heat_source_choice   = heat_source_choice,
                wall_maximum_heat_flux    = max_heat_flux,
                wall_heat_source_center   = heat_source_center,
                wall_heat_source_variance = heat_source_variance,
                total_nb_files            = total_nb_files)


    #2) source heat study: constant contact angle
    if(sourceHeatStudy):

        simulationDuration  = 100
        steady_state_ac     = 0
        temperature         = 0.95
        contact_angle       = 90.0
        phase_at_center     = 'vapor'
        ratio               = 2.0
        gravity_ac          = 0
        gravity             = 0.000
        heat_source_choice  = 'gaussian_heat_source'
        max_heat_flux_array = [-0.02,-0.04,-0.06,-0.08,-0.1]
        heat_source_center  = 0.0
        heat_source_variance = 2.0*extract_interface_length(dim2dParamPath,temperature)
        total_nb_files      = 500

    
        for max_heat_flux in max_heat_flux_array:
            
            PBSnameRun = 'dim2d_'+str(temperature)+'_sh'+str(max_heat_flux)

            [destDir, nameRun] =\
                \
                generate_wall_nonst_results(
                mainDir,
                model_input,
                PBSnameRun,
                simulationDuration,
                steady_state_ac                 = steady_state_ac,
                temperature                     = temperature,
                micro_contact_angle             = contact_angle,
                phase_at_center                 = phase_at_center,
                ratio_bubble_interface          = ratio,
                gravity_ac                      = gravity_ac,
                gravity_amp                     = gravity,
                wall_extra_heat_source_choice   = heat_source_choice,
                wall_maximum_extra_heat_flux    = max_heat_flux,
                wall_extra_heat_source_center   = heat_source_center,
                wall_extra_heat_source_variance = heat_source_variance,
                total_nb_files                  = total_nb_files)


    #3) contact angles study for fixed heat flux
    if(contactAngleStudy):

        simulationDuration   = 100
        steady_state_ac      = 0
        temperature          = 0.95
        contact_angle_array  = [112.5,130.0] #[22.5,45.0,67.5,112.5,130.0]
        phase_at_center      = 'vapor'
        ratio                = 2.0
        gravity_ac           = 0
        gravity              = 0.000
        heat_source_choice   = 'gaussian_heat_source'
        max_heat_flux        = -0.02
        heat_source_center   = 0.0
        heat_source_variance = 2.0*extract_interface_length(dim2dParamPath,temperature)
        total_nb_files       = 500

    
        for contact_angle in contact_angle_array:
            
            PBSnameRun = 'dim2d_'+str(temperature)+'_sh'+str(max_heat_flux)+'_ca'+str(contact_angle)

            [destDir, nameRun] =\
                \
                generate_wall_nonst_results(
                mainDir,
                model_input,
                PBSnameRun,
                simulationDuration,
                steady_state_ac                 = steady_state_ac,
                temperature                     = temperature,
                micro_contact_angle             = contact_angle,
                phase_at_center                 = phase_at_center,
                ratio_bubble_interface          = ratio,
                gravity_ac                      = gravity_ac,
                gravity_amp                     = gravity,
                wall_extra_heat_source_choice   = heat_source_choice,
                wall_maximum_extra_heat_flux    = max_heat_flux,
                wall_extra_heat_source_center   = heat_source_center,
                wall_extra_heat_source_variance = heat_source_variance,
                total_nb_files                  = total_nb_files)


    #4) flow velocity study for the bubble
    #   nucleation
    if(flowVelocityStudy):

        simulationDuration   = 100
        steady_state_ac      = 0
        temperature          = 0.95
        contact_angle_array  = [45.0] #[22.5,45.0,67.5,112.5,130.0]
        phase_at_center      = 'vapor'
        flow_velocity_array  = [0.1] #[0.05, 0.1, 0.15, 0.2]
        ratio                = 2.0
        gravity_ac           = 0
        gravity              = 0.000
        heat_source_choice   = 'gaussian_heat_source'
        max_heat_flux        = -0.02
        heat_source_center   = 0.0
        heat_source_variance = 2.0*extract_interface_length(dim2dParamPath,temperature)
        total_nb_files       = 500

    
        for contact_angle in contact_angle_array:

            for flow_velocity in flow_velocity_array:
            
                PBSnameRun =\
                    'dim2d_'+str(temperature)+\
                    '_sh'+str(max_heat_flux)+\
                    '_ca'+str(contact_angle)+\
                    '_v'+str(flow_velocity)

                [destDir, nameRun] =\
                    \
                    generate_wall_nonst_results(
                    mainDir,
                    model_input,
                    PBSnameRun,
                    simulationDuration,
                    steady_state_ac                 = steady_state_ac,
                    temperature                     = temperature,
                    micro_contact_angle             = contact_angle,
                    flow_velocity                   = flow_velocity,
                    phase_at_center                 = phase_at_center,
                    ratio_bubble_interface          = ratio,
                    gravity_ac                      = gravity_ac,
                    gravity_amp                     = gravity,
                    wall_extra_heat_source_choice   = heat_source_choice,
                    wall_maximum_extra_heat_flux    = max_heat_flux,
                    wall_extra_heat_source_center   = heat_source_center,
                    wall_extra_heat_source_variance = heat_source_variance,
                    total_nb_files                  = total_nb_files)
