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
    model_input=os.path.join(
        os.getenv('augeanstables'),'src','test_files','config','default_inputs','dim2d',
        'dim2d_bubble_transported_hedstrom_xy.txt')


    temperatureStudy          = True
    velocityStudy             = True
    thresholdTemperatureStudy = True
    thresholdVelocityStudy    = True


    #1) temperature study
    if(temperatureStudy):

        temperature_array = [0.95,0.99,0.995,0.999]
        flow_velocity     = 0.1
        md_threshold_ac   = 0
        large_domain_run  = False
        
    
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
                                    large_domain_run=large_domain_run)
             
    
    #2) flow mean velocity study
    if(velocityStudy):

        temperature         = 0.99
        flow_velocity_array = [0.05,0.25,0.5]
        md_threshold_ac     = 0
        large_domain_run    = False
        
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
                                    large_domain_run=large_domain_run)
             

    #3) threshold study
    if(thresholdTemperatureStudy):

        temperature_array  = [0.95,0.99,0.995,0.999]
        flow_velocity      = 0.1
        md_threshold_array = [0.01, 0.05, 0.1, 0.2, 0.3]
        md_threshold_ac    = 1
        large_domain_run   = False
        
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
                                        large_domain_run=large_domain_run,
                                        md_threshold_ac=md_threshold_ac,
                                        md_threshold=md_threshold)
    
    if(thresholdVelocityStudy):
        temperature         = 0.99
        flow_velocity_array = [0.05,0.25,0.5]
        md_threshold_array  = [0.01, 0.05, 0.1, 0.2, 0.3]
        md_threshold_ac     = 1
        large_domain_run    = False
    
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
                                        large_domain_run=large_domain_run,
                                        md_threshold_ac=md_threshold_ac,
                                        md_threshold=md_threshold)
    
    

    
    

    
        
