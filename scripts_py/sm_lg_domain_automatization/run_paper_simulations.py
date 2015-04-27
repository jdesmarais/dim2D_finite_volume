#!/usr/bin/python

'''
@description
script used to generate the main results of the
study of the impact of the DIM open b.c in 2D

the results are obtained by comparing simulations
on small domains with an equivalent simulation on
a large domain where the perturbations do not have
time to reach the edges and to come back

1) temperature studyT
\f$ T \in [0.95,0.99,0.995,0.999]\f$ at \f$v=0.1\f$

2) flow mean velocity study
\f$ v \in [0.05, 0.1, 0.25, 0.5]\f$ at \f$T=0.99\f$

3) threshold study
threshold \in [0.0001, 0.001, 0.01, 0.1, 0.2, 0.3]
for \f$ T \in [0.95,0.99,0.995,0.999] \f$
and \f$ v \in [0.05, 0.1, 0.25, 0.5]  \f$

'''

import os

from library_sm_lg_results import generate_sm_lg_results


if __name__=="__main__":


    # main directory where the simulations
    # are saved
    mainDir='/home/jdesmarais/projects'


    # input used as template
    #------------------------------------------------------------
    maindir_input = os.path.join(os.getenv('augeanstables'),
                                 'src','config',
                                 'default_inputs','dim2d')
    ##for small domain simulations
    model_input=os.path.join(maindir_input,
        'dim2d_bubble_transported_hedstrom_xy.txt')

    #for large domain simulations
    #model_input=os.path.join(maindir_input,
    #    'dim2d_bubble_transported_periodic.txt')
    

    large_domain_run          = False
    small_domain_run          = True
    nb_tiles_option           = [8,8]

    temperatureStudy          = True
    velocityStudy             = False
    thresholdTemperatureStudy = False
    thresholdVelocityStudy    = False
    icPerturbationStudy       = False
    bcPerturbationStudy_T0    = False
    bcPerturbationStudy_vx0   = False
    bcPerturbationStudy_vy0   = False


    #1) temperature study
    if(temperatureStudy):

        temperature_array = [0.95,0.99,0.995,0.999]
        flow_velocity     = 0.1
        md_threshold_ac   = 0
        
    
        for temperature in temperature_array:
            
            [destDir,
             nameRun_sm_domain,
             nameRun_lg_domain] =\
             \
             generate_sm_lg_results(mainDir,
                                    temperature,
                                    flow_velocity,
                                    model_input,
                                    bf_layer_option=True,
                                    small_domain_run=small_domain_run,
                                    large_domain_run=large_domain_run,
                                    nb_tiles_option=nb_tiles_option)
             
    
    #2) meanflow velocity study
    if(velocityStudy):

        temperature         = 0.99
        flow_velocity_array = [0.05,0.25,0.5]
        md_threshold_ac     = 0
        
        for flow_velocity in flow_velocity_array:
            
            [destDir,
             nameRun_sm_domain,
             nameRun_lg_domain] =\
             \
             generate_sm_lg_results(mainDir,
                                    temperature,
                                    flow_velocity,
                                    model_input,
                                    bf_layer_option=True,
                                    small_domain_run=small_domain_run,
                                    large_domain_run=large_domain_run,
                                    nb_tiles_option=nb_tiles_option)
             

    #3) threshold study
    if(thresholdTemperatureStudy):

        temperature_array  = [0.95,0.99,0.995,0.999]
        flow_velocity      = 0.1
        md_threshold_array = [0.01,0.05, 0.1, 0.2, 0.3]
        md_threshold_ac    = 1
        
        for md_threshold in md_threshold_array:
            for temperature in temperature_array:
                
                [destDir,
                 nameRun_sm_domain,
                 nameRun_lg_domain] =\
                 \
                 generate_sm_lg_results(mainDir,
                                        temperature,
                                        flow_velocity,
                                        model_input,
                                        bf_layer_option=True,
                                        small_domain_run=small_domain_run,
                                        large_domain_run=False,
                                        nb_tiles_option=nb_tiles_option,
                                        md_threshold_ac=md_threshold_ac,
                                        md_threshold=md_threshold)
    
    if(thresholdVelocityStudy):
        temperature         = 0.99
        flow_velocity_array = [0.05,0.25,0.5]
        md_threshold_array  = [0.01, 0.05, 0.1, 0.2, 0.3]
        md_threshold_ac     = 1
    
        for md_threshold in md_threshold_array:
            for flow_velocity in flow_velocity_array:
                
                [destDir,
                 nameRun_sm_domain,
                 nameRun_lg_domain] =\
                 \
                 generate_sm_lg_results(mainDir,
                                        temperature,
                                        flow_velocity,
                                        model_input,
                                        bf_layer_option=True,
                                        small_domain_run=small_domain_run,
                                        large_domain_run=False,
                                        nb_tiles_option=nb_tiles_option,
                                        md_threshold_ac=md_threshold_ac,
                                        md_threshold=md_threshold)


    #4) perturbation study
    if(icPerturbationStudy):

        temperature_array     = [0.95,0.99,0.995,0.999]
        flow_velocity         = 0.1
        ic_perturbation_array = [0.5] #[0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1]
        ic_perturbation_ac    = 1
        
        for ic_perturbation in ic_perturbation_array:
            for temperature in temperature_array:

                [destDir,
                 nameRun_sm_domain,
                 nameRun_lg_domain] =\
                 \
                 generate_sm_lg_results(mainDir,
                                        temperature,
                                        flow_velocity,
                                        model_input,
                                        bf_layer_option=True,
                                        small_domain_run=small_domain_run,
                                        large_domain_run=False,
                                        nb_tiles_option=nb_tiles_option,
                                        ic_perturbation_ac=ic_perturbation_ac,
                                        ic_perturbation_amp=ic_perturbation)
    
    if(bcPerturbationStudy_T0):

        temperature_array     = [0.95,0.99,0.995,0.999]
        flow_velocity         = 0.1
        bc_perturbation_array = [0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05]
        bc_perturbation_ac    = 1
        
        for bc_perturbation in bc_perturbation_array:
            for temperature in temperature_array:

                bc_perturbation_amp = temperature*bc_perturbation

                [destDir,
                 nameRun_sm_domain,
                 nameRun_lg_domain] =\
                 \
                 generate_sm_lg_results(mainDir,
                                        temperature,
                                        flow_velocity,
                                        model_input,
                                        bf_layer_option=True,
                                        small_domain_run=small_domain_run,
                                        large_domain_run=False,
                                        nb_tiles_option=nb_tiles_option,
                                        bc_perturbation_T0_ac=bc_perturbation_ac,
                                        bc_perturbation_T0_amp=bc_perturbation_amp)

    if(bcPerturbationStudy_vx0):

        temperature           = [0.95,0.99,0.995,0.999]
        flow_velocity         = 0.1
        bc_perturbation_array = [0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05]
        bc_perturbation_ac    = 1
        
        for bc_perturbation in bc_perturbation_array:
            for temperature in temperature_array:

                bc_perturbation_amp = bc_perturbation

                [destDir,
                 nameRun_sm_domain,
                 nameRun_lg_domain] =\
                 \
                 generate_sm_lg_results(mainDir,
                                        temperature,
                                        flow_velocity,
                                        model_input,
                                        bf_layer_option=True,
                                        small_domain_run=small_domain_run,
                                        large_domain_run=False,
                                        nb_tiles_option=nb_tiles_option,
                                        bc_perturbation_vx0_ac=bc_perturbation_ac,
                                        bc_perturbation_vx0_amp=bc_perturbation_amp)

    if(bcPerturbationStudy_vy0):

        temperature           = [0.95,0.99,0.995,0.999]
        flow_velocity         = 0.1
        bc_perturbation_array = [0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05]
        bc_perturbation_ac    = 1
        
        for bc_perturbation in bc_perturbation_array:
            for temperature in temperature_array:

                bc_perturbation_amp = flow_velocity*bc_perturbation

                [destDir,
                 nameRun_sm_domain,
                 nameRun_lg_domain] =\
                 \
                 generate_sm_lg_results(mainDir,
                                        temperature,
                                        flow_velocity,
                                        model_input,
                                        bf_layer_option=True,
                                        small_domain_run=small_domain_run,
                                        large_domain_run=False,
                                        nb_tiles_option=nb_tiles_option,
                                        bc_perturbation_vy0_ac=bc_perturbation_ac,
                                        bc_perturbation_vy0_amp=bc_perturbation_amp)

