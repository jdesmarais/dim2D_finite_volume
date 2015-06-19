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
    '''
    @description: extract the interface length from the 
                  parameters in dim2d_parameters and the
                  temperature for the initial conditions
    '''

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


def get_heater_length(Li):
    '''
    @description: determine the length of the heater from
                  the size of the interface length
                  (as the radius of the bubble is set to 2.0Li)
    '''
    
    return 4.0*Li


def get_heater_variation_angle_length(Li):
    '''
    @description: determine the characteristic length over which
                  the contact angle varies at the edge of the
                  heater
    '''

    return 0.5*Li


def generate_wall_nonst_results_uniform_surface(
    mainDir,
    model_input,
    PBSnameRun,
    simulationDuration,
    steady_state_ac,
    temperature,
    flow_velocity,
    phase_at_center,
    ratio_bubble_interface,
    gravity_ac,
    gravity_amp,
    wall_surface_type,
    wall_micro_contact_angle,
    wall_heater_center,
    wall_heater_length,
    wall_heat_source_choice,
    wall_maximum_heat_flux,
    total_nb_files,
    spherical_cap):
    '''
    @description: generate the results for a surface
                  with uniform contact angle
    '''

    [destDir, nameRun] =\
        \
        generate_wall_nonst_results(
        mainDir,
        model_input,
        PBSnameRun,
        simulationDuration,
        steady_state_ac               = steady_state_ac,
        temperature                   = temperature,
        flow_velocity                 = flow_velocity,
        phase_at_center               = phase_at_center,
        ratio_bubble_interface        = ratio_bubble_interface,
        gravity_ac                    = gravity_ac,
        gravity_amp                   = gravity_amp,
        wall_surface_type             = wall_surface_type,
        wall_micro_contact_angle      = contact_angle,
        wall_heater_center            = wall_heater_center,
        wall_heater_length            = wall_heater_length,
        wall_extra_heat_source_choice = wall_heat_source_choice,
        wall_maximum_extra_heat_flux  = max_heat_flux,
        total_nb_files                = total_nb_files,
        spherical_cap                 = spherical_cap)


def generate_wall_nonst_results_surface_with_heaters(
    mainDir,                           
    model_input,                       
    PBSnameRun,                        
    simulationDuration,                
    steady_state_ac,
    temperature,
    flow_velocity,
    phase_at_center,
    ratio_bubble_interface,
    gravity_ac,
    gravity_amp,
    wall_surface_type,
    wall_heater_center,
    wall_heater_length,
    wall_heater_variation_angle_length,
    wall_heater_micro_contact_angle,
    wall_micro_contact_angle,
    wall_heat_source_choice,
    wall_maximum_heat_flux,
    total_nb_files,
    spherical_cap):
    '''
    @description: generate the results for a surface with heaters
    '''

    [destDir, nameRun] =\
        \
        generate_wall_nonst_results(
        mainDir,
        model_input,
        PBSnameRun,
        simulationDuration,
        steady_state_ac                    = steady_state_ac,
        temperature                        = temperature,
        flow_velocity                      = flow_velocity,
        phase_at_center                    = phase_at_center,
        ratio_bubble_interface             = ratio_bubble_interface,
        gravity_ac                         = gravity_ac,
        gravity_amp                        = gravity_amp,
        wall_surface_type                  = wall_surface_type,
        wall_heater_center                 = wall_heater_center,
        wall_heater_length                 = wall_heater_length,
        wall_heater_variation_angle_length = wall_heater_variation_angle_length,
        wall_heater_micro_contact_angle    = wall_heater_micro_contact_angle,
        wall_micro_contact_angle           = wall_micro_contact_angle,
        wall_extra_heat_source_choice      = wall_heat_source_choice,
        wall_maximum_extra_heat_flux       = wall_maximum_heat_flux,
        total_nb_files                     = total_nb_files,
        spherical_cap                      = spherical_cap)


if __name__=="__main__":

    #file with the DIM2d parameters
    #------------------------------------------------------------
    dim2dParamPath = os.path.join(os.getenv('augeanstables'),
                                  'src',
                                  'physical_models',
                                  'dim2d',
                                  'dim2d_parameters.f')


    # main directory where the simulations are saved
    #------------------------------------------------------------
    mainDir='/home/jdesmarais/projects'


    # input used as template
    #------------------------------------------------------------
    maindir_input = os.path.join(os.getenv('augeanstables'),
                                 'src','config',
                                 'default_inputs','dim2d')

    model_input = os.path.join(maindir_input,
                               'dim2d_bubble_nucleation_at_wall.txt')


    # choice of the results generated
    #------------------------------------------------------------
    uniformNucleation_conductionHeatStudy = False
    uniformNucleation_sourceHeatStudy     = False
    uniformNucleation_contactAngleStudy   = False
    uniformNucleation_flowVelocityStudy   = False

    heaterNucleation_sourceHeatStudy   = False
    heaterNucleation_contactAngleStudy = False
    heaterNucleation_flowVelocityStudy = False

    uniformSphericalC_flowVelocityStudy = False
    heaterSphericalC_flowVelocityStudy  = True

    
    # default parameters for the generation of results
    #------------------------------------------------------------
    # wall_micro_contact_angle_uniform_nucleation
    #        when studying the influence of heat flux on
    #        nucleation, the contact angle on the surface
    #        is kept constant
    #
    # wall_maximum_heat_flux_nucleation
    #        when studying the influence of the contact
    #        angle on the nucleation, the heat flux is
    #        kept constant
    #
    # wall_micro_contact_angle_with_heaters
    #        when studying the heat flux, the contact
    #        angle on the heater, the flow velocity...,
    #        the contact angle, away from the heater, is
    #        kept constant
    #------------------------------------------------------------
    wall_micro_contact_angle_uniform_nucleation =  90.0 
    wall_maximum_heat_flux_nucleation           = -0.02
    wall_micro_contact_angle_with_heaters       =  0.0


    #============================================================
    #1) uniform surface: source heat study for one contact angle
    #============================================================
    if(uniformNucleation_conductionHeatStudy):


        simulationDuration  = 100
        steady_state_ac     = 0
        temperature         = 0.95
        flow_velocity       = 0.0

        Li = extract_interface_length(dim2dParamPath,temperature)

        phase_at_center          = 'vapor'
        ratio                    = 2.0
        gravity_ac               = 0
        gravity                  = 0.000

        wall_surface_type        = 'uniform_surface'

        wall_micro_contact_angle = wall_micro_contact_angle_uniform_nucleation

        wall_heater_center       = 0.0
        wall_heater_length       = get_heater_length(Li) #4.0*Li

        wall_heat_source_choice  = 'gaussian_heat_source'
        wall_max_heat_flux_array = [0.0001,0.001,0.01]

        total_nb_files           = 500

    
        for max_heat_flux in wall_max_heat_flux_array:
            
            PBSnameRun = 'dim2d_'+str(temperature)+'_fh'+str(max_heat_flux)

            [destDir, nameRun] =\
                \
                generate_wall_nonst_results(
                mainDir,
                model_input,
                PBSnameRun,
                simulationDuration,
                steady_state_ac          = steady_state_ac,
                temperature              = temperature,
                phase_at_center          = phase_at_center,
                ratio_bubble_interface   = ratio,
                gravity_ac               = gravity_ac,
                gravity_amp              = gravity,
                wall_surface_type        = wall_surface_type,
                wall_micro_contact_angle = contact_angle,
                wall_heater_center       = wall_heater_center,
                wall_heater_length       = wall_heater_length,
                wall_heat_source_choice  = wall_heat_source_choice,
                wall_maximum_heat_flux   = max_heat_flux,
                total_nb_files           = total_nb_files)
            

    #============================================================
    #2) source heat study for one contact angle
    #============================================================
    if(uniformNucleation_sourceHeatStudy or
       heaterNucleation_sourceHeatStudy):

        simulationDuration  = 100
        steady_state_ac     = 0
        temperature         = 0.95
        flow_velocity       = 0.0

        Li = extract_interface_length(dim2dParamPath,temperature)

        phase_at_center          = 'vapor'
        ratio_bubble_interface   = 2.0
        gravity_ac               = 0
        gravity_amp              = 0.000

        wall_heater_center       = 0.0
        wall_heater_length       = get_heater_length(Li)
        wall_heat_source_choice  = 'gaussian_heat_source'
        wall_max_heat_flux_array = [-0.02] #[-0.02,-0.04,-0.06,-0.08,-0.1]
        total_nb_files           = 500

        spherical_cap = False

        #-------------------------------------------------------------
        #2.1) uniform surface: source heat study for one contact angle
        #-------------------------------------------------------------
        if(uniformNucleation_sourceHeatStudy):
            wall_surface_type        = 'uniform_surface'
            wall_micro_contact_angle = wall_micro_contact_angle_uniform_nucleation

            for max_heat_flux in wall_max_heat_flux_array:
            
                PBSnameRun = 'dim2d_'+str(temperature)+'_sh'+str(max_heat_flux)

                generate_wall_nonst_results_uniform_surface(
                    mainDir,
                    model_input,
                    PBSnameRun,
                    simulationDuration,
                    steady_state_ac,
                    temperature,
                    flow_velocity,
                    phase_at_center,
                    ratio_bubble_interface,
                    gravity_ac,
                    gravity_amp,
                    wall_surface_type,
                    wall_micro_contact_angle,
                    wall_heater_center,
                    wall_heater_length,
                    wall_heat_source_choice,
                    max_heat_flux,
                    total_nb_files,
                    spherical_cap)


        #-------------------------------------------------------------
        #2.2) surface with heaters: source heat study for one contact angle
        #-------------------------------------------------------------
        if(heaterNucleation_sourceHeatStudy):
            wall_surface_type                  = 'surface_with_heaters'
            wall_heater_variation_angle_length = get_heater_variation_angle_length(Li)
            wall_heater_micro_contact_angle    = wall_micro_contact_angle_uniform_nucleation
            wall_micro_contact_angle           = wall_micro_contact_angle_with_heaters

            for max_heat_flux in wall_max_heat_flux_array:
            
                PBSnameRun = 'dim2d_'+str(temperature)+'_sh'+str(max_heat_flux)

                generate_wall_nonst_results_surface_with_heaters(
                    mainDir,                           
                    model_input,                       
                    PBSnameRun,                        
                    simulationDuration,                
                    steady_state_ac,
                    temperature,
                    flow_velocity,
                    phase_at_center,
                    ratio_bubble_interface,
                    gravity_ac,
                    gravity_amp,
                    wall_surface_type,
                    wall_heater_center,
                    wall_heater_length,
                    wall_heater_variation_angle_length,
                    wall_heater_micro_contact_angle,
                    wall_micro_contact_angle,
                    wall_heat_source_choice,
                    max_heat_flux,
                    total_nb_files,
                    spherical_cap)

                
    #============================================================
    #3) contact angle study for fixed heat flux
    #============================================================
    if(uniformNucleation_contactAngleStudy or
       heaterNucleation_contactAngleStudy):

        simulationDuration  = 100
        steady_state_ac     = 0
        temperature         = 0.95
        flow_velocity       = 0.0

        Li = extract_interface_length(dim2dParamPath,temperature)

        phase_at_center          = 'vapor'
        ratio_bubble_interface   = 2.0
        gravity_ac               = 0
        gravity_amp              = 0.000

        wall_micro_contact_angle = [112.5,130.0] #[22.5,45.0,67.5,112.5,130.0]
        wall_heater_center       = 0.0
        wall_heater_length       = get_heater_length(Li)
        wall_heat_source_choice  = 'gaussian_heat_source'
        wall_maximum_heat_flux   = wall_maximum_heat_flux_nucleation
        total_nb_files           = 500

        spherical_cap = False

    
        #-------------------------------------------------------------
        #3.1) uniform surface: contact angle study for one heat flux
        #-------------------------------------------------------------
        if(uniformNucleation_contactAngleStudy):

            wall_surface_type = 'uniform_surface'

            for contact_angle in contact_angle_array:
            
                PBSnameRun = 'dim2d_'+str(temperature)+'_sh'+str(max_heat_flux)+'_ca'+str(contact_angle)

                generate_wall_nonst_results_uniform_surface(
                    mainDir,
                    model_input,
                    PBSnameRun,
                    simulationDuration,
                    steady_state_ac,
                    temperature,
                    flow_velocity,
                    phase_at_center,
                    ratio_bubble_interface,
                    gravity_ac,
                    gravity_amp,
                    wall_surface_type,
                    contact_angle,
                    wall_heater_center,
                    wall_heater_length,
                    wall_heat_source_choice,
                    wall_maximum_heat_flux,
                    total_nb_files,
                    spherical_cap)

        #-------------------------------------------------------------
        #3.2) surface with heaters: contact angle study for one heat flux
        #-------------------------------------------------------------
        if(heaterNucleation_contactAngleStudy):

            wall_surface_type                  = 'surface_with_heaters'
            wall_heater_variation_angle_length = get_heater_variation_angle_length(Li)
            wall_micro_contact_angle           = wall_micro_contact_angle_with_heaters

            for contact_angle in contact_angle_array:
            
                PBSnameRun = 'dim2d_'+str(temperature)+\
                             '_ca'+str(wall_micro_contact_angle)+\
                             '_hca'+str(contact_angle)

                generate_wall_nonst_results_surface_with_heaters(
                    mainDir,                           
                    model_input,                       
                    PBSnameRun,                        
                    simulationDuration,                
                    steady_state_ac,
                    temperature,
                    flow_velocity,
                    phase_at_center,
                    ratio_bubble_interface,
                    gravity_ac,
                    gravity_amp,
                    wall_surface_type,
                    wall_heater_center,
                    wall_heater_length,
                    wall_heater_variation_angle_length,
                    contact_angle,
                    wall_micro_contact_angle,
                    wall_heat_source_choice,
                    wall_maximum_heat_flux,
                    total_nb_files,
                    spherical_cap)


    #============================================================
    #4) flow velocity study for the bubble nucleation
    #============================================================
    if(uniformNucleation_flowVelocityStudy or
       heaterNucleation_flowVelocityStudy):

        simulationDuration   = 100
        steady_state_ac      = 0
        temperature          = 0.95

        Li = extract_interface_length(dim2dParamPath,temperature)

        contact_angle_array  = [45.0] #[22.5,45.0,67.5,112.5,130.0]
        phase_at_center      = 'vapor'
        flow_velocity_array  = [0.1] #[0.05, 0.1, 0.15, 0.2]

        ratio_bubble_interface = 2.0
        gravity_ac             = 0
        gravity_amp            = 0.000

        wall_micro_contact_angle = [112.5,130.0] #[22.5,45.0,67.5,112.5,130.0]
        wall_heater_center       = 0.0
        wall_heater_length       = get_heater_length(Li)
        wall_heat_source_choice  = 'gaussian_heat_source'
        max_heat_flux_array      = [-0.04,-0.06,-0.08,-0.1]

        total_nb_files           = 500
        spherical_cap            = False


        #-------------------------------------------------------------
        #4.1) uniform surface: flow velocity study
        #-------------------------------------------------------------
        if(uniformNucleation_flowVelocityStudy):

            wall_surface_type = 'uniform_surface'

            for max_heat_flux in max_heat_flux_array:

                for contact_angle in contact_angle_array:

                    for flow_velocity in flow_velocity_array:
            
                        PBSnameRun =\
                            'dim2d_'+str(temperature)+\
                            '_sh'+str(max_heat_flux)+\
                            '_ca'+str(contact_angle)+\
                            '_v'+str(flow_velocity)                        

                        generate_wall_nonst_results_uniform_surface(
                            mainDir,
                            model_input,
                            PBSnameRun,
                            simulationDuration,
                            steady_state_ac,
                            temperature,
                            flow_velocity,
                            phase_at_center,
                            ratio_bubble_interface,
                            gravity_ac,
                            gravity_amp,
                            wall_surface_type,
                            wall_micro_contact_angle,
                            wall_heater_center,
                            wall_heater_length,
                            wall_heat_source_choice,
                            wall_maximum_heat_flux,
                            total_nb_files,
                            spherical_cap)


        #-------------------------------------------------------------
        #4.2) surface with heaters: flow velocity study
        #-------------------------------------------------------------
        if(heaterNucleation_flowVelocityStudy):

            wall_surface_type                  = 'surface_with_heaters'
            wall_heater_variation_angle_length = get_heater_variation_angle_length(Li)
            wall_micro_contact_angle           = wall_micro_contact_angle_with_heaters


            for max_heat_flux in max_heat_flux_array:

                for contact_angle in contact_angle_array:

                    for flow_velocity in flow_velocity_array:
            
                        PBSnameRun =\
                            'dim2d_'+str(temperature)+\
                            '_sh'+str(max_heat_flux)+\
                            '_ca'+str(wall_micro_contact_angle)+\
                            '_v'+str(flow_velocity)+\
                            '_hca'+str(contact_angle)
                        
                        generate_wall_nonst_results_surface_with_heaters(
                            mainDir,
                            model_input,
                            PBSnameRun,
                            simulationDuration,
                            steady_state_ac,
                            temperature,
                            flow_velocity,
                            phase_at_center,
                            ratio_bubble_interface,
                            gravity_ac,
                            gravity_amp,
                            wall_surface_type,
                            wall_heater_center,
                            wall_heater_length,
                            wall_heater_variation_angle_length,
                            contact_angle,
                            wall_micro_contact_angle,
                            wall_heat_source_choice,
                            wall_maximum_heat_flux,
                            total_nb_files,
                            spherical_cap)
                            
                            
    #============================================================
    #5) flow velocity study with a spherical cap approximation
    #   for the initial bubble: no heat flux
    #============================================================
    if(uniformSphericalC_flowVelocityStudy or
       heaterSphericalC_flowVelocityStudy):

        simulationDuration     = 100
        steady_state_ac        = 0
        temperature            = 0.95

        Li = extract_interface_length(dim2dParamPath,temperature)

        contact_angle_array    = [45.0] #[22.5,45.0,67.5,112.5,130.0]
        phase_at_center        = 'vapor'
        flow_velocity_array    = [0.1] #[0.05, 0.1, 0.15, 0.2]
        ratio_bubble_interface = 2.0
        gravity_ac             = 0
        gravity_amp            = 0.000

        wall_heat_source_choice = 'no_heat_source'
        max_heat_flux_array     = [0.0]

        wall_heater_center     = 0.0
        wall_heater_length     = get_heater_length(Li)

        total_nb_files         = 500
        spherical_cap          = True


        #-------------------------------------------------------------
        #5.1) uniform surface: flow velocity study (spherical_cap)
        #-------------------------------------------------------------
        if(uniformSphericalC_flowVelocityStudy):

            for max_heat_flux in max_heat_flux_array:        

                for contact_angle in contact_angle_array:

                    for flow_velocity in flow_velocity_array:

                        PBSnameRun =\
                            'dim2d_'+str(temperature)+\
                            '_ca'+str(contact_angle)

                        if(max_heat_flux!=0.0):
                            PBSnameRun += '_sh'+str(max_heat_flux)
                            
                        if(flow_velocity!=0.0):
                            PBSnameRun += '_v'+str(flow_velocity)
                                
                        if(spherical_cap):
                            PBSnameRun += '_sph'
                                    
                        generate_wall_nonst_results_uniform_surface(
                            mainDir,
                            model_input,
                            PBSnameRun,
                            simulationDuration,
                            steady_state_ac,
                            temperature,
                            flow_velocity,
                            phase_at_center,
                            ratio_bubble_interface,
                            gravity_ac,
                            gravity_amp,
                            wall_surface_type,
                            wall_micro_contact_angle,
                            wall_heater_center,
                            wall_heater_length,
                            wall_heat_source_choice,
                            max_heat_flux,
                            total_nb_files,
                            spherical_cap)

        #-------------------------------------------------------------
        #5.2) surface with heaters: flow velocity study (spherical cap)
        #-------------------------------------------------------------
        if(heaterSphericalC_flowVelocityStudy):

            wall_surface_type                  = 'surface_with_heaters'
            wall_heater_variation_angle_length = get_heater_variation_angle_length(Li)
            wall_micro_contact_angle           = wall_micro_contact_angle_with_heaters

            for max_heat_flux in max_heat_flux_array:

                for contact_angle in contact_angle_array:

                    for flow_velocity in flow_velocity_array:
            
                        PBSnameRun =\
                            'dim2d_'+str(temperature)+\
                            '_ca'+str(contact_angle)

                        if(max_heat_flux!=0.0):
                            PBSnameRun += '_sh'+str(max_heat_flux)
                            
                        if(flow_velocity!=0.0):
                            PBSnameRun += '_v'+str(flow_velocity)
                                
                        PBSnameRun+='_hca'+str(contact_angle)

                        if(spherical_cap):
                            PBSnameRun += '_sph'
                        
                        generate_wall_nonst_results_surface_with_heaters(
                            mainDir,
                            model_input,
                            PBSnameRun,
                            simulationDuration,
                            steady_state_ac,
                            temperature,
                            flow_velocity,
                            phase_at_center,
                            ratio_bubble_interface,
                            gravity_ac,
                            gravity_amp,
                            wall_surface_type,
                            wall_heater_center,
                            wall_heater_length,
                            wall_heater_variation_angle_length,
                            contact_angle,
                            wall_micro_contact_angle,
                            wall_heat_source_choice,
                            max_heat_flux,
                            total_nb_files,
                            spherical_cap)


