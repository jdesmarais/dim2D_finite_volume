#!/usr/bin/python

'''
@description
useful functions to generate the executables and
the pbs scripts to run the wall steady state simulations
'''


import os
import sys
import inspect
import subprocess
import math
import string
import shutil
import shlex
import time


#add the python files from ../sm_lg_automatization
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],
    "../sm_lg_domain_automatization")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from automatization_csts import bc_size

from create_sm_lg_inputs import (get_parameter,
                                 create_sm_lg_inputs)

from library_sm_lg_results import (generate_exe,
                                   estimate_simulation_duration,
                                   estimate_wall_time,
                                   create_pbs_script,
                                   get_name_run,
                                   run_simulation)

from create_wall_st_inputs import create_wall_st_inputs


# get the name of the folder for the simulation
def get_simulation_dir(temperature,
                       micro_contact_angle,
                       phase_at_center,
                       collapse_ratio=2.0,
                       gravity_amp=0):
    '''
    @description
    get the name of the folder where the simulation results
    are saved
    '''

    simDir =\
        'dim2d_'+\
        str(temperature)+'_'+\
        'ca'+str(micro_contact_angle)+'_'+\
        phase_at_center[0:3]

    if(gravity_amp!=0):
        simDir+='_g'+str(gravity_amp)

    if(collapse_ratio!=2.0):
        simDir+='_ra'+str(collapse_ratio)

    return simDir


# create the simulation
def create_simulation(destDir,
                      inputPath):
    '''
    @description
    create the executable corresponding to the inputPath, save it
    to the destDir, create the PBS script to run the simulation
    and save it in the destDir
    '''

    #0) check that destDir exists
    if(not os.path.isdir(destDir)):
        print 'library_wall_st_results'
        print 'create_simulation'
        sys.exit('***directory '+destDir+' does not exist***')

    
    #1) create the executable corresponding to the inputPath
    exePath = generate_exe(inputPath)


    #2) move the executable to destDir
    if(not os.path.isfile(exePath)):
        print 'library_wall_st_results'
        print 'create_simulation'
        sys.exit('***exe file'+exePath+' does not exist***')

    newExePath = destDir
    newExePath+='/'+os.path.basename(exePath)
    shutil.move(exePath,newExePath)


    #3) create the PBS script
    temperature         = float(get_parameter('temperature',inputPath))
    micro_contact_angle = float(get_parameter('micro_contact_angle',inputPath))
    nameRun = get_name_run(temperature,micro_contact_angle)
    simulation_duration = 4.0*60.0*60.0 #estimate_simulation_duration(inputPath)
    walltime = estimate_wall_time(simulation_duration,
                                  safety_ratio=6.0)

    pbsScriptPath = destDir+'/run_sim.job'

    create_pbs_script(
        pbsScriptPath,
        newExePath,
        nameRun=nameRun,
        walltime=walltime)

    return [pbsScriptPath,nameRun]


# create and run the simulation for test
def generate_wall_st_results(mainDir,
                             steady_state_ac,
                             temperature,
                             micro_contact_angle,
                             phase_at_center,
                             model_input,
                             gravity_ac=0,
                             gravity_amp=0):
    '''
    @description
    create the directory to save the simulation results,
    generate the input file needed to run the wall steady
    state simulations,
    generate the executable,
    create the PBS scripts file, and
    run the simulation
    '''
    
    #1) test whether 'mainDir' can be used as a reference
    #   directory where the directory to save the simulation
    #   is created
    if(not os.path.isdir(mainDir)):
        print 'library_wall_st_results'
        print 'generate_wall_st_results'
        sys.exit('*** '+mainDir+' is not a directory***')


    #2) create the directory to save the simulation    
    destDir = get_simulation_dir(temperature,
                                 micro_contact_angle,
                                 phase_at_center,
                                 gravity_amp=gravity_amp)
    destDir = os.path.join(mainDir,destDir)

    # if there is already an existing directory, the function
    # throws an error
    if(os.path.isdir(destDir)):
        print 'library_wall_st_results'
        print 'generate_wall_st_results'
        sys.exit('*** '+destDir+' already exists***')
    os.mkdir(destDir)
    
    
    #3) create the inputs for the simulation
    inputPath = 'inputs_wall.txt'

    # remove old input files
    if(os.path.isfile(inputPath)):
        os.remove(inputPath)

    #create input
    create_wall_st_inputs(steady_state_ac,
                          temperature,
                          micro_contact_angle,
                          phase_at_center,
                          gravity_ac,
                          gravity_amp,
                          model_input,
                          inputs_wall_modified = inputPath)

    #4) create dir, generate executable,
    #   create PBS script file
    [pbsScriptPath,nameRun] = create_simulation(destDir,
                                                inputPath)
        

    #5) run the simulation
    run_simulation(pbsScriptPath)

    return [destDir,nameRun]


if __name__ == "__main__":

    #test: generate_exe
    inputPath = os.getenv('augeanstables')+'/src/test_files/config/'
    inputPath+='default_inputs/dim2d/dim2d_bubble_transported_hedstrom_xy.txt'
    ##exePath   = generate_exe(inputPath,bf_layer_option=True)
    ##print 'exePath: ', exePath
    #
    ##test: estimate_simulation_duration
    #print 'estimate_simulation_duration: ', estimate_simulation_duration(inputPath)
    #
    ##test: estimate_wall_time
    #print 'estimate_wall_time: ', estimate_wall_time(90*3600.0+2*60+3)
    #
    #
    ##test: get_name_run
    #print 'get_name_run: ', get_name_run(0.95,0.1)
    #
    #
    ##test: create_pbs_script
    #outputPath = os.getenv('augeanstables')
    #outputPath+='/scripts_py/sm_lg_domain_automatization'
    #outputPath+='/script_pbs.job'
    #
    #exePath = '/home/jdesmarais/projects'
    #exePath+='/dim2d_bubble_trs_hedstrom_xy_bf'
    #exePath+='/test_verify_x_symmetry'
    #
    #create_pbs_script(outputPath,
    #                  exePath,
    #                  nameRun='dim2d_test',
    #                  walltime='01:00:00',
    #                  nodes=1,
    #                  ppn=1,
    #                  email='j.desmarais@tue.nl')
    #
    ##test: create_simulation
    #destDir = os.getenv('augeanstables')
    #destDir+= '/scripts_py/sm_lg_domain_automatization/sim_dir'
    #[pbsScriptPath,nameRun] = create_simulation(destDir,
    #                                            inputPath,
    #                                            bf_layer_option=True)
    #print 'pbsScriptPath: ', pbsScriptPath
    #
    #
    ##test: run_simulation
    #run_simulation(pbsScriptPath)
    #
    #
    ##test: get_job_status
    #get_job_status('dim2d_0.995_0.1')
    #
    #
    ##test: check_simulation_status
    #check_simulation_status([nameRun],sleepTime=2)


    #test: create_and_run_test
    #mainDir='/home/jdesmarais/Code/augeanstables/scripts_py/sm_lg_domain_automatization'

    mainDir=os.path.join(os.getenv('HOME'),'projects')

    steady_state_ac     = 1
    temperature         = 0.999
    micro_contact_angle = 50.0
    phase_at_center     = 'liquid'

    model_input=os.path.join(os.getenv('augeanstables'),
                             'src','config',
                             'default_inputs',
                             'dim2d',
                             'dim2d_bubble_next_to_wall.txt')

    [destDir,nameRun] = generate_wall_st_results(
        mainDir,
        steady_state_ac,
        temperature,
        micro_contact_angle,
        phase_at_center,
        model_input,
        gravity_ac=0,
        gravity_amp=0)
