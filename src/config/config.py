#!/usr/bin/python

import os
import sys
import getopt
import subprocess
import shlex
import string

import argparse

from library_config import\
    get_flow_config,\
    get_bc_config
    


#root of the code
augeanstablesPath= os.getenv('augeanstables')

# main configuration directories
configPath       = os.path.join(augeanstablesPath,'src','config')
paramPath        = os.path.join(augeanstablesPath,'src','parameters')

#sh script files used to modify the code by the configuration
changeParameterPath= os.path.join(configPath,'change_parameter.sh')
getParameterPath   = os.path.join(configPath,'get_parameter.sh')

#files modified by the code configuration
paramCstPath       = os.path.join(paramPath,'parameters_constant.f')
paramInputPath     = os.path.join(paramPath,'parameters_input.f')
makefileHeaderPath = os.path.join(configPath,'makefile_header.mk')

#folder where the executables are compiled
exeDir = os.path.join(augeanstablesPath,'src','main')


# analyze the options passed to the program
def parse_argv(argv):
    '''
    @description:
    parse the program arguments: get the input file
    '''

    #definition of the potential options
    parser = argparse.ArgumentParser(description='Configuration of the source code')

    parser.add_argument("-c", "--compile", action="store_true",
                        help="ask to compile the code")

    parser.add_argument("-i", "--input", default='None', type=str, nargs=1,
                        help="input file describing the test case")

    parser.add_argument("-b", "--buffer", action="store_true",
                        help="enable the use of domain extension")

    parser.add_argument("-p", "--parallel", default=[1,1], type=int, nargs=2,
                        help="number of tiles for parallel run")

    args = parser.parse_args()


    # process the options
    ## input file
    inputFile= args.input[0]
    if(not os.path.isfile(inputFile)):
        print 'input file supplied does not exist'
        print 'input:', inputFile
        sys.exit(2)
    
    ## compilation
    compileCode=args.compile

    ## buffer layer
    compileCodeBuffer=args.buffer
    
    ## parallel
    nbTiles = args.parallel


    # print the configuration
    print ''
    print 'input file      : ', inputFile.replace(configPath,".")
    print 'code_compilation: ', compileCode
    print 'domain_extension: ', compileCodeBuffer
    print 'number of tiles : ', args.parallel
    print ''

    return [inputFile,compileCode,compileCodeBuffer,nbTiles]


# read the SHA reference number of the last commit to 
# the git repository
def read_commit():
    '''
    @description:
    determine the program version using the last
    commit SHA reference number
    '''
    cmd="git log -1 --format=%H"
    args = shlex.split(cmd)
    output = subprocess.Popen(args,stdout=subprocess.PIPE).communicate()[0]
    
    return output.strip()


# set the SHA reference number of the last commit to the
# git repository as a parameter in parameters_constant.f
def set_commit(file_path):
    '''
    @description:
    set the parameter 'commit' in the file to the last commit
    SHA reference number
    '''

    commit_ID = read_commit()
    commit_ID = commit_ID

    cmd=changeParameterPath
    cmd+=" -i "+str(file_path)
    cmd+=" -o "+str(file_path)
    cmd+=" -p "+'commit'
    cmd+=" -v "+commit_ID
    cmd+=" -q"
    subprocess.call(cmd, shell=True)

    print 'update parameters_constant.f for commit '+commit_ID
    

# convert an integer into a logical string
def int_to_logical_str(input_int):
    '''
    @description:
    convert an integer into a logical string
    '''    

    if(input_int==1):
        input_loc = '.true.'
    else:
        input_loc = '.false.'
    
    return input_loc


# read the parameters saved in the input text file
def read_inputs(filename, inputs_needed):
    '''
    @description:
    read the input file
    '''    
    # define a dictionnary to store the inputs needed
    # then, if one needs the input read 'x_min', one simply
    # use : inputs_read['x_min']
    inputs_read={}

    for input_param in inputs_needed:
        cmd=getParameterPath
        cmd+=" -i "+str(filename)
        cmd+=" -p "+input_param
        args = shlex.split(cmd)
        output = subprocess.Popen(args,stdout=subprocess.PIPE).communicate()[0]

        
        # the input asked should be converted to a float
        try :
            inputs_read[input_param] = float(output)
            
        # if it does not work it is because this input in a chain
        # of characters and therefore it should remain as such
        except ValueError:
            output=output.replace('\r','')
            output=output.replace('\n','')
            inputs_read[input_param] = output

            # if the input was not found in the input file
            if(output==''):
                print '['+str(input_param)+'] not found in input file'
                sys.exit(2)
            
    return inputs_read


# compute the number of grid points in the computational domain
# considering the borders, the space step and the extent of the
# boundary layer
def compute_n(x_min,x_max,dx,bc_size):
    '''
    @description:
    compute the number of gridpoints for the tile
    such that the space step coresponds to dx
    '''
    n = int(round((x_max-x_min)/dx))+2*bc_size+1
    return n


# compute the number of grid points in the computational domain
# considering the borders, the space step and the extent of the
# boundary layer for a parallel computation
def compute_n_par(npx,x_min,x_max,dx,bc_size):
    '''
    @description:
    compute the number of gridpoints for the tile
    such that the space step coresponds to dx
    '''

    n_total         = int(round((x_max-x_min)/dx))+1           #total number of grid points w/o b.c.
    n_per_proc      = int(round(float(n_total)*1./float(npx))) #number of grdpts per processor w/o the b.c.
    n_per_proc_w_bc = n_per_proc + 2*bc_size                   #number of grdpts per processor w/ the b.c.

    
    if(not n_per_proc*npx==n_total):
        print '****total number of grdpts different b/w sequential and parallel****'

    return n_per_proc_w_bc


# compute the total extent of the computational domain
# depending on the number of processors and the number
# of space steps
def compute_ntx_and_nty(npx,npy,x_min,x_max,dx,y_min,y_max,dy,bc_size):
    '''
    @description:
    compute the total extent of the computational domain
    depending on the number of processors and the number
    of space steps
    '''
    if(npx==1):
        nx = compute_n(x_min,x_max,dx,bc_size)
    else:
        nx = compute_n_par(npx,x_min,x_max,dx,bc_size)
    
    if(npy==1):
        ny = compute_n(y_min,y_max,dy,bc_size)
    else:
        ny = compute_n_par(npy,y_min,y_max,dy,bc_size)

    ntx=int(round(npx*nx))
    nty=int(round(npy*ny))

    return [ntx,nty]


# turn the parameters read from the input file into
# inputs saved in the parameters_input.f file of the
# program
def compute_code_inputs(inputFileName,nbTiles):
    '''
    @description
    compute all the inputs needed by the code
    '''

    # codes for ic_choice, bc_choice, bc_type_choice,
    # gravity_choice: these codes are defined in
    # parameters_constant.f
    pm_code      = ['simpletest_choice',
                    'wave1d_choice',
                    'wave2d_choice',
                    'ns2d_choice',
                    'dim2d_choice']

    wave2d_ic_code = ['non-existing',
                      'peak',
                      'non-existing',
                      'non-existing',
                      'non-existing',
                      'negative_spot']

    ns2d_ic_code = ['steady_state',
                    'peak',
                    'vortex',
                    'sym_x',
                    'sym_y']

    dim2d_ic_code= ['steady_state',
                    'drop_retraction',
                    'bubble_ascending',
                    'homogeneous_liquid',
                    'drop_collision',
                    'phase_separation',
                    'bubble_transported',
                    'bubble_next_to_wall',
                    'bubble_collapse',
                    'bubble_nucleation',
                    'bubble_spherical_cap',
                    'newgrdpt_test',
                    'bubbles_transported']

    bc_code      = ['periodic_xy_choice',
                    'reflection_xy_choice',
                    'hedstrom_xy_choice',
                    'poinsot_xy_choice',
                    'yoolodato_xy_choice',
                    'wall_xy_choice',
                    'wall_S_reflection_choice',
                    'wall_S_open_choice',
                    'half_wall_S_open_choice']

    surface_type_code = ['uniform_surface',
                         'surface_with_heaters']

    wave_forcing_code = ['no_wave_forcing',
                         'oscillatory_forcing',
                         'intermittent_oscillatory_forcing',
                         'moving_oscillatory_forcing']

    # read the input file
    inputs_needed=['detail_print',
                   'dt','t_max','steady_state_ac',
                   'dx','x_min','x_max',
                   'dy','y_min','y_max',
                   'pm_choice',
                   'bc_choice',
                   'openbc_detector_distance',
                   'openbc_md_threshold_ac',
                   'openbc_md_threshold',
                   'openbc_perturbation_T0_ac',
                   'openbc_perturbation_vx0_ac',
                   'openbc_perturbation_vy0_ac',
                   'openbc_perturbation_T0_amp',
                   'openbc_perturbation_vx0_amp',
                   'openbc_perturbation_vy0_amp',
                   'wall_surface_type',
                   'wall_micro_contact_angle',
                   'wall_heater_center',
                   'wall_heater_length',
                   'wall_heater_variation_angle_length',
                   'wall_heater_micro_contact_angle',
                   'wall_heat_source_choice',
                   'wall_maximum_heat_flux',
                   'wall_extra_heat_source_choice',
                   'wall_maximum_extra_heat_flux',
                   'ic_choice',
                   'flow_direction',
                   'flow_velocity',
                   'flow_profile',
                   'temperature',
                   'phase_at_center',
                   'ratio_bubble_interface',
                   'ic_perturbation_ac',
                   'ic_perturbation_amp',
                   'li_perturbation_ac',
                   'li_perturbation_amp',
                   'li_separation',
                   'dim2d_lowTemperature',
                   'gravity_ac',
                   'gravity_amp',
                   'wave_forcing']

    inputs=read_inputs(inputFileName, inputs_needed)
    
    inputs_computed = {}

    #------------------------------------------------------------
    # division of the computational domain into tiles
    #------------------------------------------------------------
    # get the number of processors
    inputs['npx']=nbTiles[0]
    inputs['npy']=nbTiles[1]

    # compute the ntx and nty determining the
    # extent of the computational domain
    bc_size=2
    [inputs_computed['ntx'],
     inputs_computed['nty']]=compute_ntx_and_nty(
        inputs['npx'], inputs['npy'],
        inputs['x_min'],inputs['x_max'],inputs['dx'],
        inputs['y_min'],inputs['y_max'],inputs['dy'],
        bc_size)


    #------------------------------------------------------------
    # duration of the simulation
    #------------------------------------------------------------
    # determine whether it is a steady state simulation
    inputs['steady_state_ac'] = int_to_logical_str(
        int(inputs['steady_state_ac']))


    #------------------------------------------------------------
    # type of physical model
    #------------------------------------------------------------
    # compute the pm_choice
    inputs['pm_choice'] = pm_code[int(inputs['pm_choice'])]
    pm_choice = inputs['pm_choice']

    # compute the ne
    if(pm_choice=='simpletest_choice'):
        inputs_computed['ne'] = 1

    if(pm_choice=='wave1d_choice'):
        inputs_computed['ne'] = 2

    if(pm_choice=='wave2d_choice'):
        inputs_computed['ne'] = 3

    if(pm_choice=='ns2d_choice' or pm_choice=='dim2d_choice'):
        inputs_computed['ne'] = 4


    #------------------------------------------------------------
    # type of initial conditions
    #------------------------------------------------------------
    # compute the ic_choice
    if(pm_choice=='wave2d_choice'):
        inputs['ic_choice'] = wave2d_ic_code[int(inputs['ic_choice'])]

    elif(pm_choice=='ns2d_choice'):
        inputs['ic_choice'] = ns2d_ic_code[int(inputs['ic_choice'])]

    elif(pm_choice=='dim2d_choice'):
        inputs['ic_choice'] = dim2d_ic_code[int(inputs['ic_choice'])]

    else:
        inputs['ic_choice'] = ns2d_ic_code[0]

    # compute the ic_perturbation_ac
    inputs['ic_perturbation_ac'] = int_to_logical_str(
        int(inputs['ic_perturbation_ac']))

    # compute the ic_perturbation_ac
    inputs['li_perturbation_ac'] = int_to_logical_str(
        int(inputs['li_perturbation_ac']))

    # determine the flow parameters
    get_flow_config(inputs['flow_direction'],inputs_computed)

    
    #------------------------------------------------------------
    # type of boundary conditions
    #------------------------------------------------------------
    # compute the bc_choice    
    inputs['bc_choice'] = bc_code[int(inputs['bc_choice'])]
    bc_choice = inputs['bc_choice']
    
    get_bc_config(bc_choice,inputs_computed)


    #------------------------------------------------------------
    # perturbation of open boundary conditions
    #------------------------------------------------------------
    # compute the openbc_md_threshold_ac
    inputs['openbc_md_threshold_ac'] = int_to_logical_str(
        int(inputs['openbc_md_threshold_ac']))

    # determine the openbc_detector_distance
    inputs['openbc_detector_distance'] =\
        int(inputs['openbc_detector_distance'])

    # compute the openbc_perturbation_T0_ac
    inputs['openbc_perturbation_T0_ac'] = int_to_logical_str(
        int(inputs['openbc_perturbation_T0_ac']))

    # compute the openbc_perturbation_vx0_ac
    inputs['openbc_perturbation_vx0_ac'] = int_to_logical_str(
        int(inputs['openbc_perturbation_vx0_ac']))

    # compute the openbc_perturbation_vy0_ac
    inputs['openbc_perturbation_vy0_ac'] = int_to_logical_str(
        int(inputs['openbc_perturbation_vy0_ac']))


    #------------------------------------------------------------
    # body forces
    #------------------------------------------------------------
    # compute the wave_forcing_choice
    inputs['wave_forcing']   = wave_forcing_code[int(inputs['wave_forcing'])]

    # compute the dim2d_lowTemperature
    inputs['dim2d_lowTemperature'] = int_to_logical_str(
        int(inputs['dim2d_lowTemperature']))

    # compute gravity_ac
    inputs['gravity_ac'] = int_to_logical_str(
        int(inputs['gravity_ac']))
    #------------------------------------------------------------


    return [inputs,inputs_computed]


# update the 'parameters_input.f' file with the inputs
# of the simulation
def update_parameters_inputs(file_path,
                             inputs,
                             inputs_computed):
    '''
    @description
    update the constants defined in the 'parameters_input'
    file
    '''
    
    # change the constant that do not require a special
    # output treatment (integer,character...)
    constants_changed1={
        'npx'                              : inputs['npx'],
        'npy'                              : inputs['npy'],
        'ntx'                              : inputs_computed['ntx'],
        'nty'                              : inputs_computed['nty'],
        'ne'                               : inputs_computed['ne'],
        'flow_profile'                     : inputs['flow_profile'],
        'pm_choice'                        : inputs['pm_choice'],
        'ic_choice'                        : inputs['ic_choice'],
        'phase_at_center'                  : inputs['phase_at_center'],
        'ic_perturbation_ac'               : inputs['ic_perturbation_ac'],
        'li_perturbation_ac'               : inputs['li_perturbation_ac'],
        'bc_choice'                        : inputs['bc_choice'],
        'bc_N_choice'                      : inputs_computed['bc_N_choice'],
        'bc_S_choice'                      : inputs_computed['bc_S_choice'],
        'bc_E_choice'                      : inputs_computed['bc_E_choice'],
        'bc_W_choice'                      : inputs_computed['bc_W_choice'],
        'bc_NW_choice'                     : inputs_computed['bc_NW_choice'],
        'bc_NE_choice'                     : inputs_computed['bc_NE_choice'],
        'bc_SW_choice'                     : inputs_computed['bc_SW_choice'],
        'bc_SE_choice'                     : inputs_computed['bc_SE_choice'],
        'bc_order1'                        : inputs_computed['bc_order1'],
        'bc_order2'                        : inputs_computed['bc_order2'],
        'bc_order3'                        : inputs_computed['bc_order3'],
        'bc_order4'                        : inputs_computed['bc_order4'],
        'bc_order5'                        : inputs_computed['bc_order5'],
        'bc_order6'                        : inputs_computed['bc_order6'],
        'bc_order7'                        : inputs_computed['bc_order7'],
        'bc_order8'                        : inputs_computed['bc_order8'],
        'bc_N_type_choice'                 : inputs_computed['bc_N_type_choice'],
        'bc_S_type_choice'                 : inputs_computed['bc_S_type_choice'],
        'bc_E_type_choice'                 : inputs_computed['bc_E_type_choice'],
        'bc_W_type_choice'                 : inputs_computed['bc_W_type_choice'],
        'bc_NW_type_choice'                : inputs_computed['bc_NW_type_choice'],
        'bc_NE_type_choice'                : inputs_computed['bc_NE_type_choice'],
        'bc_SW_type_choice'                : inputs_computed['bc_SW_type_choice'],
        'bc_SE_type_choice'                : inputs_computed['bc_SE_type_choice'],
        'wall_surface_type'                : inputs['wall_surface_type'],
        'wall_heat_source_choice'          : inputs['wall_heat_source_choice'],
        'wall_extra_heat_source_choice'    : inputs['wall_extra_heat_source_choice'],
        'obc_dct_distance'                 : inputs['openbc_detector_distance'],
        'obc_perturbation_T0_ac'           : inputs['openbc_perturbation_T0_ac'],
        'obc_perturbation_vx0_ac'          : inputs['openbc_perturbation_vx0_ac'],
        'obc_perturbation_vy0_ac'          : inputs['openbc_perturbation_vy0_ac'],
        'adapt_N_choice'                   : inputs_computed['adapt_N_choice'],
        'adapt_S_choice'                   : inputs_computed['adapt_S_choice'],
        'adapt_E_choice'                   : inputs_computed['adapt_E_choice'],
        'adapt_W_choice'                   : inputs_computed['adapt_W_choice'],
        'bf_openbc_md_threshold_ac'        : inputs['openbc_md_threshold_ac'],
        'gravity_ac'                       : inputs['gravity_ac'],
        'wave_forcing'                     : inputs['wave_forcing'],
        'flow_direction'                   : inputs_computed['flow_direction'],
        'dim2d_lowTemperature'             : inputs['dim2d_lowTemperature'],
        'debug_adapt_computational_domain' : inputs['adapt_computational_domain'],
        'steady_state_simulation'          : inputs['steady_state_ac']}

    for key, value  in constants_changed1.items():

        cmd=changeParameterPath
        cmd+=" -i "+str(file_path)
        cmd+=" -o "+str(file_path)
        cmd+=" -p "+key
        cmd+=" -v "+str(value)
        subprocess.call(cmd, shell=True)


    # change the constant that do require a special
    # output treatment (output format)
    constants_changed2={
        'x_min'                              : inputs['x_min'],
        'x_max'                              : inputs['x_max'],
        'y_min'                              : inputs['y_min'],
        'y_max'                              : inputs['y_max'],
        't_max'                              : inputs['t_max'],
        'dt'                                 : inputs['dt'],
        'detail_print'                       : inputs['detail_print'],
        'flow_x_side'                        : inputs_computed['flow_x_side'],
        'flow_y_side'                        : inputs_computed['flow_y_side'],
        'flow_velocity'                      : inputs['flow_velocity'],
        'T0'                                 : inputs['temperature'],
        'ic_perturbation_amp'                : inputs['ic_perturbation_amp'],
        'li_perturbation_amp'                : inputs['li_perturbation_amp'],
        'li_separation'                      : inputs['li_separation'],
        'gravity_amp'                        : inputs['gravity_amp'],
        'ratio_bubble_interface'             : inputs['ratio_bubble_interface'],
        'wall_micro_contact_angle'           : inputs['wall_micro_contact_angle'],
        'wall_heater_center'                 : inputs['wall_heater_center'],
        'wall_heater_length'                 : inputs['wall_heater_length'],
        'wall_heater_variation_angle_length' : inputs['wall_heater_variation_angle_length'],
        'wall_heater_micro_contact_angle'    : inputs['wall_heater_micro_contact_angle'],
        'wall_maximum_heat_flux'             : inputs['wall_maximum_heat_flux'],
        'wall_maximum_extra_heat_flux'       : inputs['wall_maximum_extra_heat_flux'],
        'obc_perturbation_T0_amp'            : inputs['openbc_perturbation_T0_amp'],
        'obc_perturbation_vx0_amp'           : inputs['openbc_perturbation_vx0_amp'],
        'obc_perturbation_vy0_amp'           : inputs['openbc_perturbation_vy0_amp'],
        'bf_openbc_md_threshold'             : inputs['openbc_md_threshold']}

    for key, value in constants_changed2.items():

        cmd=changeParameterPath
        cmd+=" -i "+str(file_path)
        cmd+=" -o "+str(file_path)
        cmd+=" -p "+key
        cmd+=" -v "+"%10.10fd0"%value
        subprocess.call(cmd, shell=True)

    print 'update parameters_input.f'


# update the makefile with the path to the folders
# needed for the simulation
def update_makefile(file_path,inputs):
    '''
    @description
    update the folder for the compilation
    of the boundary conditions
    '''

    # define the constants changed in the file
    constants_changed={
        'pm_choice':inputs['pm_choice'],
        'ic_choice':inputs['ic_choice'],
        'bc_choice':inputs['bc_choice']}


    # change the constant in the file
    for key, value in constants_changed.items():

        cmd=changeParameterPath
        cmd+=" -i "+str(file_path)
        cmd+=" -o "+str(file_path)
        cmd+=" -p "+key
        cmd+=" -v "+value
        subprocess.call(cmd, shell=True)

    print 'update makefile_header.mk'


def compile_code(inputs,compileCodeBuffer):
    '''
    @description
    compile the code using the previous configuration
    '''

    #name of the executable file obtained
    fname='sim_dim2d'
    fname+='_'+str(inputs['npx'])+'x'+str(inputs['npy'])

    #commands for compilation
    cmd_serial  ='cd '+exeDir+' && make cleanall && make sim_dim2d && make clean'

    cmd_serial_bf = 'cd '+exeDir+' && make cleanall && make sim_dim2d_bf && make clean'

    cmd_parallel='cd '+exeDir+' && make cleanall && make sim_dim2d_par && make clean'
    cmd_parallel+=' && mv sim_dim2d_par '+fname

    
    #serial compilation
    if(inputs['npx']*inputs['npy']==1):

        #with domain adaptation: buffer layers
        if(compileCodeBuffer):
            cmd      = cmd_serial_bf
            name_exe = 'sim_dim2d_bf'

        #without domain adaptation
        else:
            cmd      = cmd_serial
            name_exe = 'sim_dim2d'

    #parallel compilation
    else:
        cmd      = cmd_parallel
        name_exe = fname
        
    subprocess.call(cmd, shell=True)


    if(os.path.isfile(os.path.join(exeDir,name_exe))):

        print ''
        print 'executable ready: '+name_exe+' in '+exeDir
        print ''

    else:

        print ''
        print 'compilation error'
        print ''


if __name__ == "__main__":

    # define the paths for the files modified by the
    # configuration
    param_path           = paramInputPath
    makefile_path        = makefileHeaderPath
    param_cst_path       = paramCstPath


    # parse the program arguments
    [inputFileName,
     compileCode,
     compileCodeBuffer,
     nbTiles]=parse_argv(sys.argv[1:])


    # compute the code inputs
    [inputs,inputs_computed]=compute_code_inputs(inputFileName,nbTiles)

    
    #check whether the computational domain is
    #adapted or not
    if(compileCodeBuffer):
        inputs['adapt_computational_domain']='.true.'
    else:
        inputs['adapt_computational_domain']='.false.'


    # replace the inputs in the 'parameters_input' file
    update_parameters_inputs(param_path,
                             inputs,
                             inputs_computed)


    # replace the inputs in the 'makefile'
    update_makefile(makefile_path,inputs)


    # replace the commit SHA number in the
    # 'parameters_constant'
    set_commit(param_cst_path)


    # print the major results
    print ''
    print '(ntx,nty,ne)',\
        inputs_computed['ntx'],\
        inputs_computed['nty'],\
        inputs_computed['ne']


    # print the end of the configuration
    print ''
    print 'end of configuration'
    print ''

    
    # compile the code
    if(compileCode):
        compile_code(inputs,compileCodeBuffer)
        
