#!/usr/bin/python

"""
Run the simulations studying how the computational domain is extended
when two bubbles are transported by the flow: whether the domain extension
handles both or whether each bubble leaves the domain separately
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

#add the python files from ../wall_steady_state_automatization
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],
    "../wall_steady_state_automatization")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

# load modules
from library_sm_lg_results import run_simulation

from library_wall_nonst_results import create_simulation

from create_bb_tr_input import create_bb_tr_inputs


if __name__=="__main__":


    # main directory where the simulations are saved
    #------------------------------------------------------------
    mainDir='/home/jdesmarais/projects'
    
    
    # input used as template
    #------------------------------------------------------------
    maindir_input = os.path.join(os.getenv('augeanstables'),
                                 'src','config',
                                 'default_inputs','dim2d')

    model_input = os.path.join(maindir_input,
                               'dim2d_bubbles_transported_hedstrom_xy.txt')


    # set of parameters studied
    #------------------------------------------------------------
    temperature = 0.999
    flow_velocity = 0.1
    flow_direction ='y'
    ratios_interface_separation = [1.5,2.0,4.0]
    dct_distance    = 8
    md_threshold_ac = 1
    md_threshold    = 0.2
    adapt_domain    = True


    # run the simulations corresponding to the parameters studied
    #------------------------------------------------------------
    for ratio_interface_separation in ratios_interface_separation:

        # name for the input file created
        inputPath = 'inputs.txt'


        # create the input file for the simulation
        create_bb_tr_inputs(temperature,
                            flow_velocity,
                            ratio_interface_separation,
                            model_input,
                            flow_direction = flow_direction,
                            output_file    = inputPath)
        

        # choose the directory for the simulation outputs
        destDir =\
            'dim2d_'+\
            str(temperature)+\
            '_'+str(flow_velocity)+\
            '_sep'+str(ratio_interface_separation)

        destDir = os.path.join(mainDir,destDir)

        
        # if there is already an existing directory, the function
        # throws an error
        if(os.path.isdir(destDir)):
            raise OSError('library_wall_nonst_results',
                          'generate_wall_nonst_results',
                          '*** '+destDir+' already exists***')
        os.mkdir(destDir)
            
        
        # choose the name for the run in PBS script
        PBSnameRun = 'dim2d_sep'+str(ratio_interface_separation)+'_'+str(temperature)


        # create the PBS script to run the simulation
        [pbsScriptPath,nameRun] = create_simulation(destDir,
                                                    inputPath,
                                                    PBSnameRun,
                                                    adapt_domain=adapt_domain)


        # run the simulation
        run_simulation(pbsScriptPath)
        
