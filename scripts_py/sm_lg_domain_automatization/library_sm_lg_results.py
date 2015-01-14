#!/usr/bin/python

'''
@description
useful functions to generate the executables and the pbs scripts to
run the small and large domain simulations
'''


import os
import sys
import subprocess
import math
import string
import shutil
import shlex
import time

from automatization_csts import bc_size

from create_sm_lg_inputs import (get_parameter,
                                 create_sm_lg_inputs)


# generate exe
def generate_exe(inputPath,bf_layer_option=False):
    '''
    @description:
    generate an optimized executable for the simulation
    with the inputs stored in 'inputPath'
    if the bf_layer_option is enabled, the simulation is
    runned with the potential use of buffer layers
    '''

    #1) check if the inputFile exists
    if(os.path.isfile(inputPath)):


        #2) set the environment variable $AUGEANSTABLES_PROFILE
        #   to false to ask for an optimized executable
        os.putenv('AUGEANSTABLES_PROFILE','false')


        #3) create the executable depending on the buffer layer
        #   option
        cmd = os.getenv('augeanstables')+"/src/test_files/config/config.py"
        cmd+= " -i "+str(inputPath)
        cmd+= " -c "
        if(bf_layer_option):
            cmd+=" -b "
        subprocess.call(cmd, shell=True)


        #4) get the executable name
        exePath = os.getenv('augeanstables')+"/src/test_files/"
        if(bf_layer_option):
            exePath += "sim_dim2d_bf"
        else:
            exePath += "sim_dim2d"


        #5) verify that the executable has been created
        if(os.path.isfile(exePath)):

            return exePath

        else:
            print 'library_sm_lg_results'
            print 'generate_exe'
            print 'error when compiling executable'
            sys.exit('*** '+exePath+' does not exist***')        

    else:
        print 'library_sm_lg_results'
        print 'generate_exe'
        sys.exit('*** '+inputPath+' does not exist***')


# estimate the total simulation time
def estimate_simulation_duration(inputPath):
    '''
    @description:
    estimate the simulation duration based on the
    grid size in x- and y-directions and the total
    number of time steps
    '''

    #1) verify that the file 'inputPath' exists
    if(os.path.isfile(inputPath)):


        #2) extract the following parameters
        #   [dt,t_max,dx,x_min,x_max,dy,y_min,y_max]
        #   from 'inputPath'
        dt    = float(get_parameter('dt',inputPath))
        t_max = float(get_parameter('t_max',inputPath))
        dx    = float(get_parameter('dx',inputPath))
        x_min = float(get_parameter('x_min',inputPath))
        x_max = float(get_parameter('x_max',inputPath))
        dy    = float(get_parameter('dy',inputPath))
        y_min = float(get_parameter('y_min',inputPath))
        y_max = float(get_parameter('y_max',inputPath))


        #3) determine nx,ny,nt, the number of space steps
        #   in the x- and y-directions and the total number
        #   of time steps in the simulation
        nx = math.ceil((x_max-x_min)/dx) + 2*bc_size
        ny = math.ceil((y_max-y_min)/dy) + 2*bc_size
        nt = math.ceil(t_max/dt)


        #4) we know that for 104*104 grid size and 4000 time steps
        #   the simulation lasts around 5.0min
        #   for each time step, the RK scheme takes 3 steps
        #   so the average time spent computing each grid point at
        #   each intermediate time step is (in seconds):
        dt = 5.0*60.0/(3.0*4000.0*104*104)
        
        
        #5) for [nx*ny] grid size and 3*nt intermediate time steps
        #   the simulation will approximately lasts:
        time = dt*nx*ny*3*nt

        return time

    else:
        print 'library_sm_lg_results'
        print 'estimate_simulation_duration'
        sys.exit('*** '+inputPath+' does not exist***')
        
    
# convert the simulation duration in wall time
def estimate_wall_time(simulation_time,safety_ratio=1.0):
    '''
    @description:
    convert the simulation duration in wall time for the
    pbs script
    '''

    #use a safety margin
    wall_time = simulation_time*safety_ratio

    #estimate the number of hours
    hours = math.floor(wall_time/3600.0)

    #estimate the number of minutes
    minutes = math.ceil((wall_time-hours*3600.0)/60.0)

    #return the wall time
    return "%02d" % hours+':'+"%02d" % minutes + ':00'


# create name run for pbs script
def get_name_run(temperature,flow_velocity):
    '''
    @description:
    determine the name for the simulation in the pbs
    script
    '''

    return 'dim2d_'+str(temperature).strip()+'_'+str(flow_velocity).strip()


# create the pbs script to run the simulation
def create_pbs_script(
    outputPath,
    exePath,
    nameRun='dim2d',
    walltime='01:00:00',
    nodes=1,
    ppn=1,
    email='j.desmarais@tue.nl'):
    '''
    @description:
    create the pbs script needed to run the simulation
    '''

    #1) if the output file already exists, the
    #   previous file is removed
    if(os.path.isfile(outputPath)):
        os.remove(outputPath)

    
    #2) check whether the executable exists
    if(not os.path.isfile(exePath)):
        print 'library_sm_lg_results'
        print 'create_pbs_script'
        sys.exit('***exe file '+exePath+' does not exist***')

    #3) open the output file and write the PBS script
    f = open(outputPath,'w')
    
    f.write('#!/bin/bash\n')
    f.write('#\n')
    f.write('#PBS -l walltime='+walltime+'\n')                 #wall time
    f.write('#PBS -l nodes='+str(nodes)+':ppn='+str(ppn)+'\n') #nb of nodes + nb of processors per node
    f.write('#PBS -j oe \n')                                   #join error messages to output log
    f.write('#PBS -N '+nameRun+'\n')                           #name for the run as it appears with qstat
    f.write('#PBS -M <'+email+'>\n')                           #email
    f.write('#PBS -m e\n')                                     #send email at the end
    f.write('\n')
    f.write('cd '+os.path.dirname(exePath)+'\n')
    f.write('./'+os.path.splitext(os.path.basename(exePath))[0]+' 1>sim_dim2d.out 2>sim_dim2d.err \n')

    f.close()


# create the simulation
def create_simulation(destDir,
                      inputPath,
                      bf_layer_option=False,
                      nameRun_suffix=''):
    '''
    @description
    create the executable corresponding to the inputPath, save it
    to the destDir, create the PBS script to run the simulation
    and save it in the destDir
    '''

    #0) check that destDir exists
    if(not os.path.isdir(destDir)):
        print 'library_sm_lg_results'
        print 'create_simulation'
        sys.exit('***directory '+destDir+' does not exist***')

    
    #1) create the executable corresponding to the inputPath
    exePath = generate_exe(inputPath,
                           bf_layer_option=bf_layer_option)


    #2) move the executable to destDir
    if(not os.path.isfile(exePath)):
        print 'library_sm_lg_results'
        print 'create_simulation'
        sys.exit('***exe file'+exePath+' does not exist***')

    newExePath = destDir
    newExePath+='/'+os.path.basename(exePath)
    shutil.move(exePath,newExePath)


    #3) create the PBS script
    temperature   = float(get_parameter('temperature',inputPath))
    flow_velocity = float(get_parameter('flow_velocity',inputPath))
    nameRun = get_name_run(temperature,flow_velocity)
    nameRun+= nameRun_suffix
    simulation_duration = estimate_simulation_duration(inputPath)
    walltime = estimate_wall_time(simulation_duration,safety_ratio=1.5)

    pbsScriptPath = destDir+'/run_sim.job'

    create_pbs_script(
        pbsScriptPath,
        newExePath,
        nameRun=nameRun,
        walltime=walltime)

    return [pbsScriptPath,nameRun]


# submit the PBS script file to the job queue
def run_simulation(pbsScriptPath):
    '''
    @description
    submit the PBS script file to the job queue
    '''
    
    if(not os.path.isfile(pbsScriptPath)):
        print 'library_sm_lg_results'
        print 'run_simulation'
        sys.exit('***the PBS script '+pbsScriptPath+' does not exist***')

    cmd = 'qsub '+pbsScriptPath
    subprocess.call(cmd, shell=True)


#get job status
def get_job_status(nameRun):
    '''
    @description
    get the status of a job
    '''

    # ask for the status of all processes
    cmd = "qstat"
    args = shlex.split(cmd)
    p_qstat = subprocess.Popen(args,
                               stdout=subprocess.PIPE)

    # ask for the line countaining our job whose name
    # is 'nameRun'
    cmd = "grep "+nameRun
    args = shlex.split(cmd)
    p_grep = subprocess.Popen(args,
                              stdin=p_qstat.stdout,
                              stdout=subprocess.PIPE)

    p_qstat.stdout.close()
    output_grep = p_grep.communicate()[0]

    
    # if there is no line matching, it means that the job
    # is not running anymore, so the status is 'completed'='C'
    if(output_grep==''):
        status='C'

    # otherwise the status is the fifth column in the line
    # split removes the blank spaces and 4 is for the 5th
    # element in the array
    else:
        status = string.split(output_grep)[4]

    return status


# check the simulation status: exit the function when all jobs are completed
def check_simulation_status(nameRuns,sleepTime=300,detailled=False):
    '''
    @description
    exit the function when all jobs are finished
    '''
    
    notAlljobsCompleted = True

    while notAlljobsCompleted:

        time.sleep(sleepTime)

        allJobsCompleted = True

        for nameRun in nameRuns:
            
            statusRun     = get_job_status(nameRun)
            allJobsCompleted = allJobsCompleted and (statusRun=='C')
            
        notAlljobsCompleted=not allJobsCompleted


# create and run the simulation for test
def generate_sm_lg_results(mainDir,
                           temperature,
                           flow_velocity,
                           model_input,
                           bf_layer_option=False):
    '''
    @description
    create the directory to save the simulation results (small and large
    domains), generate the input files needed to run the small and large
    domain simulations, generate the executables for each, create the PBS
    scripts files for each and run the simulations
    '''
    
    #1) test whether 'mainDir' can be used as a directory
    #   to save the simulations
    if(not os.path.isdir(mainDir)):
        print 'library_sm_lg_results'
        print 'create_and_run_test'
        sys.exit('*** '+mainDir+' is not a directory***')


    #2) create the directory to save the simulations
    destDir = mainDir
    destDir+= '/dim2d_'+str(temperature)+'_'+str(flow_velocity)
    
    # if there is already an existing directory, the function
    # throws an error
    if(os.path.isdir(destDir)):
        print 'library_sm_lg_results'
        print 'create_and_run_test'
        sys.exit('*** '+destDir+' already exists***')
    os.mkdir(destDir)
    
    
    #3) create the inputs for the small and the large
    #   domain simulations
    inputPath_sm_domain = 'inputs_sm_domain.txt'
    inputPath_lg_domain = 'inputs_lg_domain.txt'

    # remove old input files
    if(os.path.isfile(inputPath_sm_domain)):
        os.remove(inputPath_sm_domain)
    if(os.path.isfile(inputPath_lg_domain)):
        os.remove(inputPath_lg_domain)

    #create inputs
    create_sm_lg_inputs(temperature,
                        flow_velocity,
                        model_input,
                        sm_domain=inputPath_sm_domain,
                        lg_domain=inputPath_lg_domain)


    #4) create dir, generate executable, create PBS
    #   script file for small domain
    destDir_sm_domain = destDir+'/sm_domain'
    os.mkdir(destDir_sm_domain)
    
    [pbsScriptPath_sm_domain,nameRun_sm_domain] = create_simulation(destDir_sm_domain,
                                                                    inputPath_sm_domain,
                                                                    bf_layer_option=bf_layer_option,
                                                                    nameRun_suffix='_sm')
     

    #5) create dir, generate executable, create PBS
    #   script file for large domain
    destDir_lg_domain = destDir+'/lg_domain'
    os.mkdir(destDir_lg_domain)

    [pbsScriptPath_lg_domain,nameRun_lg_domain] = create_simulation(destDir_lg_domain,
                                                                    inputPath_lg_domain,
                                                                    bf_layer_option=bf_layer_option,
                                                                    nameRun_suffix='_lg')
    

    #6) run the small domain simulation
    run_simulation(pbsScriptPath_sm_domain)

    
    #7) run the large domain simulation
    run_simulation(pbsScriptPath_lg_domain)


    return [destDir,nameRun_sm_domain,nameRun_lg_domain]


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
    mainDir='/home/jdesmarais/projects'

    temperature=0.99
    flow_velocity=0.1

    model_input=os.getenv('augeanstables')+'/src/test_files/config/'
    model_input+='default_inputs/dim2d/dim2d_bubble_transported_hedstrom_xy.txt'

    [destDir,nameRun_sm_domain,nameRun_lg_domain] = generate_sm_lg_results(mainDir,
                                                                           temperature,
                                                                           flow_velocity,
                                                                           model_input,
                                                                           bf_layer_option=True)
