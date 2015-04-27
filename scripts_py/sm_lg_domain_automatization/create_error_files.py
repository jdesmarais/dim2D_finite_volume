#!/usr/bin/python

'''
@description
script used to generate the error files of the
study of the impact of the DIM open b.c in 2D

the results are obtained by comparing simulations
on small domains with an equivalent simulation on
a large domain where the perturbations do not have
time to reach the edges and to come back

1) temperature study
\f$ T \in [0.95,0.99,0.995,0.999]\f$ at \f$v=0.1\f$

2) flow mean velocity study
\f$ v \in [0.05, 0.1, 0.25, 0.5]\f$ at \f$T=0.99\f$

3) threshold study
threshold \in [0.0001, 0.001, 0.01, 0.1, 0.2, 0.3]
for \f$ T \in [0.95,0.99,0.995,0.999] \f$
and \f$ v \in [0.05, 0.1, 0.25, 0.5]  \f$

'''

import os

from library_sm_lg_results import get_simulation_dir
from library_sm_lg_error   import generate_error_files


def get_paths_lg_simulations(main_lg_dirs,
                             temperature_array,
                             flow_velocity_array,
                             temperature_dict,
                             flow_velocity_dict):
    '''
    @description
    create the paths to the results of the large domain
    '''

    lg_dirs = []
    lg_dirs.append([])
    for temperature in temperature_array:
        for flow_velocity in flow_velocity_array:

            dir_name = get_simulation_dir(temperature,flow_velocity,0)
            lg_dir   = os.path.join(main_lg_dirs,dir_name,'lg_domain','tr_files')
            
            lg_dirs[temperature_dict[str(temperature)]].append(lg_dir)

        lg_dirs.append([])

    return lg_dirs


def generate_simulation_error_files(
    main_sm_dirs,
    main_lg_dirs,
    temperature_dict,
    flow_velocity_dict,
    lg_dirs,
    temperature_array,
    flow_velocity_array,
    md_threshold_array=[0],
    ic_perturbation_array=[0],
    bc_perturbation_T0_array=[0],
    bc_perturbation_vx0_array=[0],
    bc_perturbation_vy0_array=[0]
    ):
    '''
    @description
    generate the error files for the results
    '''

    for bc_perturbation_T0 in bc_perturbation_T0_array:
        for bc_perturbation_vx0 in bc_perturbation_vx0_array:
            for bc_perturbation_vy0 in bc_perturbation_vy0_array:
                for ic_perturbation in ic_perturbation_array:
                    for md_threshold in md_threshold_array:
                        for flow_velocity in flow_velocity_array:
                            for temperature in temperature_array:
            
                                #print the parameters for which the error files are generated
                                print 'temperature        : '+str(temperature)
                                print 'velocity           : '+str(flow_velocity)
                                print 'threshold          : '+str(md_threshold)
                                print 'ic_perturbation    : '+str(ic_perturbation)
                                print 'bc_perturbation_vy0: '+str(bc_perturbation_vy0)
                                print 'bc_perturbation_vx0: '+str(bc_perturbation_vx0)
                                print 'bc_perturbation_T0 : '+str(bc_perturbation_T0)
                                print '------------------------------------------------------------'

                                #small domain simulation results
                                dir_name = get_simulation_dir(temperature,
                                                              flow_velocity,
                                                              md_threshold,
                                                              ic_perturbation_amp=ic_perturbation,
                                                              bc_perturbation_T0_amp=bc_perturbation_T0*temperature,
                                                              bc_perturbation_vx0_amp=bc_perturbation_vx0,
                                                              bc_perturbation_vy0_amp=bc_perturbation_vy0*flow_velocity)

                                dataDir_sm_domain = os.path.join(main_sm_dirs,dir_name,'sm_domain')

                                #large domain simulation results
                                T_i = temperature_dict[str(temperature)]
                                v_i = flow_velocity_dict[str(flow_velocity)]
                                dataDir_lg_domain = lg_dirs[T_i][v_i]

                                #check whether the sm_domain path exists
                                sm_domain_exists = os.path.isdir(dataDir_sm_domain)
                                if(not sm_domain_exists):
                                    print 'sm_domain_path: ', dataDir_sm_domain
                                    print '**** path does not exist ****'
                
                                #check whether the lg_domain path exists
                                lg_domain_exists = os.path.isdir(dataDir_lg_domain)
                                if(not lg_domain_exists):
                                    print 'lg_domain_path: ', dataDir_lg_domain
                                    print '**** path does not exist ****'

                                #generate the error files by comparing the
                                #two simulations
                                if(sm_domain_exists and lg_domain_exists):

                                    generate_error_files(
                                        dataDir_sm_domain,
                                        dataDir_lg_domain)


if __name__=="__main__":


    # main directory where the simulations are saved
    #------------------------------------------------------------
    mainDir = os.path.join(os.getenv('HOME'),
                           'projects')

    main_sm_dirs = os.path.join(mainDir,
                                '20150424_dim2d_bb_trans_cv_r3.5_search4_over2_lin_crenel')

    main_lg_dirs = os.path.join(mainDir,
                                '20150422_dim2d_bb_trans_cv_r3.5_lg')


    #paths for the large domain simulations
    #------------------------------------------------------------
    flow_velocity_dict = {'0.05':0, '0.1':1, '0.25':2,  '0.5':3}
    temperature_dict   = {'0.95':0,'0.99':1,'0.995':2,'0.999':3}


    #create paths for large domani simulation results
    #------------------------------------------------------------
    lg_dirs = get_paths_lg_simulations(main_lg_dirs,
                                       [0.95,0.99,0.995,0.999],
                                       [0.05,0.1,0.25,0.5],
                                       temperature_dict,
                                       flow_velocity_dict)
    


    # options for the creation of the error files
    #------------------------------------------------------------
    # - temperatureStudy          : create error files for the
    #                               temperature study
    # - velocityStudy             : create error files for the
    #                               velocity study
    # - thresholdTemperatureStudy : create error files for the
    #                               threshold study on temperature
    # - thresholdVelocityStudy    : create error files for the
    #                               threshold study on velocity
    #------------------------------------------------------------
    temperatureStudy          = False
    velocityStudy             = False
    thresholdTemperatureStudy = False
    thresholdVelocityStudy    = False
    icPerturbationStudy       = True
    bcPerturbationStudy_T0    = False
    bcPerturbationStudy_vx0   = False
    bcPerturbationStudy_vy0   = False


    #1) temperature study
    if(temperatureStudy):

        temperature_array   = [0.99,0.995] #[0.95,0.99,0.995,0.999]
        flow_velocity_array = [0.1]
        
        generate_simulation_error_files(
            main_sm_dirs,
            main_lg_dirs,
            temperature_dict,
            flow_velocity_dict,
            lg_dirs,
            temperature_array,
            flow_velocity_array)

    
    #2) flow mean velocity study
    if(velocityStudy):

        temperature_array   = [0.99]
        flow_velocity_array = [0.05,0.25,0.5]
        
        generate_simulation_error_files(
            main_sm_dirs,
            main_lg_dirs,
            temperature_dict,
            flow_velocity_dict,
            lg_dirs,
            temperature_array,
            flow_velocity_array)


    #3) threshold studies
    if(thresholdTemperatureStudy):

        temperature_array   = [0.95,0.99,0.995,0.999]
        flow_velocity_array = [0.1]
        md_threshold_array  = [0.05, 0.1, 0.2, 0.3]
        
        
        generate_simulation_error_files(
            main_sm_dirs,
            main_lg_dirs,
            temperature_dict,
            flow_velocity_dict,
            lg_dirs,
            temperature_array,
            flow_velocity_array,
            md_threshold_array=md_threshold_array)


    if(thresholdVelocityStudy):
        
        temperature_array   = [0.99]
        flow_velocity_array = [0.05] #[0.05,0.1,0.25,0.5]
        md_threshold_array  = [0.05, 0.1, 0.2, 0.3]
        
        
        generate_simulation_error_files(
            main_sm_dirs,
            main_lg_dirs,
            temperature_dict,
            flow_velocity_dict,
            lg_dirs,
            temperature_array,
            flow_velocity_array,
            md_threshold_array=md_threshold_array)

    
    #4) perturbation studies
    if(icPerturbationStudy):

        temperature_array     = [0.95,0.99,0.995] #[0.95,0.99,0.995,0.999]
        flow_velocity_array   = [0.1]
        ic_perturbation_array = [0.00001] #[0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05,0.1]

        generate_simulation_error_files(
            main_sm_dirs,
            main_lg_dirs,
            temperature_dict,
            flow_velocity_dict,
            lg_dirs,
            temperature_array,
            flow_velocity_array,
            ic_perturbation_array=ic_perturbation_array)

        
    if(bcPerturbationStudy_T0):

        temperature_array     = [0.95,0.99,0.995,0.999]
        flow_velocity_array   = [0.1]
        bc_perturbation_array = [0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05]

        generate_simulation_error_files(
            main_sm_dirs,
            main_lg_dirs,
            temperature_dict,
            flow_velocity_dict,
            lg_dirs,
            temperature_array,
            flow_velocity_array,
            bc_perturbation_T0_array=bc_perturbation_array)

     
    if(bcPerturbationStudy_vx0):

        temperature_array     = [0.95,0.99,0.995,0.999]
        flow_velocity_array   = [0.1]
        bc_perturbation_array = [0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05]

        generate_simulation_error_files(
            main_sm_dirs,
            main_lg_dirs,
            temperature_dict,
            flow_velocity_dict,
            lg_dirs,
            temperature_array,
            flow_velocity_array,
            bc_perturbation_vx0_array=bc_perturbation_array)

    
    if(bcPerturbationStudy_vy0):

        temperature_array     = [0.95,0.99,0.995,0.999]
        flow_velocity_array   = [0.1]
        bc_perturbation_array = [0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05]

        generate_simulation_error_files(
            main_sm_dirs,
            main_lg_dirs,
            temperature_dict,
            flow_velocity_dict,
            lg_dirs,
            temperature_array,
            flow_velocity_array,
            bc_perturbation_vy0_array=bc_perturbation_array)        
