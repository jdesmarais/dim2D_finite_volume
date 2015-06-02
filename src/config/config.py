#!/usr/bin/python

import os
import sys
import getopt
import subprocess
import shlex
import string

import argparse


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

        # convert the parameter read into float
        # except for the parameter 'flow_direction'
        # which must remain of character type
        if(input_param!='flow_direction'):
            try :
                inputs_read[input_param]=float(output)
            except ValueError:
                print str(input_param)+' not found in input file'
                sys.exit(2)
                
            
        else:
            output=output.replace('\r','')
            output=output.replace('\n','')
            inputs_read[input_param]=output

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
                    'bubble_next_to_wall']

    bc_code      = ['periodic_xy_choice',
                    'reflection_xy_choice',
                    'wall_xy_choice',
                    'wall_x_reflection_y_choice',
                    'hedstrom_xy_choice',
                    'hedstrom_xy_corners_choice',
                    'hedstrom_x_reflection_y_choice',
                    'poinsot_xy_choice',
                    'yoolodato_xy_choice']

    bc_type_code = ['bc_nodes_choice',
                    'bc_fluxes_choice',
                    'bc_timedev_choice',
                    'bc_flux_and_node_choice']

    gravity_code = ['no_gravity_choice',
                    'earth_gravity_choice']

    wave_forcing_code = ['no_wave_forcing',
                         'oscillatory_forcing',
                         'intermittent_oscillatory_forcing',
                         'moving_oscillatory_forcing']

    #in order to set the correct direction of the flow
    #(N,S,E,W,NE,NW,SE,SW) from the inputs.txt, three
    #parameters are initialized in parameters_input.f
    #[1]: the direction of the flow is either horizontal
    #(x_direction), vertical (y_direction) or diagonal
    #(xy_direction)
    #[2]: whether the flow is right or left (1.0/-1.0)
    #[3]: whether the flow is upward or downwards (1.0/-1.0)
    flow_direction_code = {
        'N':  [ 'y_direction', 1.0, 1.0],
        'S':  [ 'y_direction', 1.0,-1.0],
        'E':  [ 'x_direction', 1.0, 1.0],
        'W':  [ 'x_direction',-1.0, 1.0],
        'NE': ['xy_direction', 1.0, 1.0],
        'NW': ['xy_direction',-1.0, 1.0],
        'SE': ['xy_direction', 1.0,-1.0],
        'SW': ['xy_direction',-1.0,-1.0]}

    # read the input file
    inputs_needed=['x_min','x_max','dx',
                   'y_min','y_max','dy',
                   'dt','t_max','detail_print',
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
                   'flow_direction',
                   'flow_velocity',
                   'temperature',
                   'ic_choice',
                   'ic_perturbation_ac',
                   'ic_perturbation_amp',
                   'li_perturbation_ac',
                   'li_perturbation_amp',
                   'gravity_choice',
                   'wave_forcing',
                   'dim2d_lowTemperature']
    inputs=read_inputs(inputFileName, inputs_needed)
    

    # update the type of the inputs
    inputs['npx']=nbTiles[0]
    inputs['npy']=nbTiles[1]


    # compute the ntx and nty determining the
    # extent of the computational domain
    bc_size=2
    [ntx,nty]=compute_ntx_and_nty(
        inputs['npx'], inputs['npy'],
        inputs['x_min'],inputs['x_max'],inputs['dx'],
        inputs['y_min'],inputs['y_max'],inputs['dy'],
        bc_size)

    # compute the pm_choice
    pm_choice = pm_code[int(inputs['pm_choice'])]

    # compute the ne
    if(pm_choice=='simpletest_choice'):
        ne = 1
    if(pm_choice=='wave1d_choice'):
        ne = 2
    if(pm_choice=='wave2d_choice'):
        ne = 3
    if(pm_choice=='ns2d_choice'):
        ne = 4
    if(pm_choice=='dim2d_choice'):
        ne = 4

    # compute the ic_choice
    if(pm_choice=='wave2d_choice'):
        ic_choice = wave2d_ic_code[int(inputs['ic_choice'])]

    elif(pm_choice=='ns2d_choice'):
        ic_choice = ns2d_ic_code[int(inputs['ic_choice'])]

    elif(pm_choice=='dim2d_choice'):
        ic_choice = dim2d_ic_code[int(inputs['ic_choice'])]

    else:
        ic_choice = ns2d_ic_code[0]

    # compute the ic_perturbation_ac
    ic_perturbation_ac = int_to_logical_str(
        int(inputs['ic_perturbation_ac']))

    # compute the ic_perturbation_ac
    li_perturbation_ac = int_to_logical_str(
        int(inputs['li_perturbation_ac']))

    # determine the flow parameters
    flow_direction = flow_direction_code[inputs['flow_direction']][0]
    flow_x_side    = flow_direction_code[inputs['flow_direction']][1]
    flow_y_side    = flow_direction_code[inputs['flow_direction']][2]

    
    # compute the bc_choice    
    bc_choice = bc_code[int(inputs['bc_choice'])]
    
    # compute the bc_type_choice
    if(bc_choice=='periodic_xy_choice' or
       bc_choice=='reflection_xy_choice'):

        bc_N_type_choice = bc_type_code[0]
        bc_S_type_choice = bc_type_code[0]
        bc_E_type_choice = bc_type_code[0]
        bc_W_type_choice = bc_type_code[0]
        

    if(bc_choice=='wall_xy_choice'):

        bc_N_type_choice = bc_type_code[3]
        bc_S_type_choice = bc_type_code[3]
        bc_E_type_choice = bc_type_code[3]
        bc_W_type_choice = bc_type_code[3]

       
    if(bc_choice=='wall_x_reflection_y_choice'):

        bc_N_type_choice = bc_type_code[3]
        bc_S_type_choice = bc_type_code[3]
        bc_E_type_choice = bc_type_code[0]
        bc_W_type_choice = bc_type_code[3]


    if(bc_choice=='hedstrom_xy_choice' or
       bc_choice=='hedstrom_xy_corners_choice' or
       bc_choice=='poinsot_xy_choice' or
       bc_choice=='yoolodato_xy_choice'):

        bc_N_type_choice = bc_type_code[2]
        bc_S_type_choice = bc_type_code[2]
        bc_E_type_choice = bc_type_code[2]
        bc_W_type_choice = bc_type_code[2]


    if(bc_choice=='hedstrom_x_reflection_y_choice'):

        bc_N_type_choice = bc_type_code[0]
        bc_S_type_choice = bc_type_code[0]
        bc_E_type_choice = bc_type_code[2]
        bc_W_type_choice = bc_type_code[2]


    # compute the openbc_md_threshold_ac
    openbc_md_threshold_ac = int_to_logical_str(
        int(inputs['openbc_md_threshold_ac']))

    # determine the openbc_detector_distance
    inputs['openbc_detector_distance'] =\
        int(inputs['openbc_detector_distance'])

    # compute the openbc_perturbation_T0_ac
    openbc_perturbation_T0_ac = int_to_logical_str(
        int(inputs['openbc_perturbation_T0_ac']))

    # compute the openbc_perturbation_vx0_ac
    openbc_perturbation_vx0_ac = int_to_logical_str(
        int(inputs['openbc_perturbation_vx0_ac']))

    # compute the openbc_perturbation_vy0_ac
    openbc_perturbation_vy0_ac = int_to_logical_str(
        int(inputs['openbc_perturbation_vy0_ac']))

    # compute the gravity_choice
    gravity_choice = gravity_code[int(inputs['gravity_choice'])]
    wave_forcing   = wave_forcing_code[int(inputs['wave_forcing'])]

    # compute the dim2d_lowTemperature
    dim2d_lowTemperature = int_to_logical_str(
        int(inputs['dim2d_lowTemperature']))

    return [inputs,
            ntx,nty,ne,
            pm_choice,
            ic_choice,
            ic_perturbation_ac,
            li_perturbation_ac,
            bc_choice,
            openbc_md_threshold_ac,
            openbc_perturbation_T0_ac,
            openbc_perturbation_vx0_ac,
            openbc_perturbation_vy0_ac,
            bc_N_type_choice,
            bc_S_type_choice,
            bc_E_type_choice,
            bc_W_type_choice,
            gravity_choice,
            wave_forcing,
            flow_direction,
            flow_x_side,
            flow_y_side,
            dim2d_lowTemperature]


# update the 'parameters_input.f' file with the inputs
# of the simulation
def update_parameters_inputs(file_path,inputs,ntx,nty,ne,
                             pm_choice,
                             ic_choice,
                             ic_perturbation_ac,
                             li_perturbation_ac,
                             bc_choice,
                             openbc_md_threshold_ac,
                             openbc_perturbation_T0_ac,
                             openbc_perturbation_vx0_ac,
                             openbc_perturbation_vy0_ac,
                             bc_N_type_choice,
                             bc_S_type_choice,
                             bc_E_type_choice,
                             bc_W_type_choice,
                             gravity_choice,
                             wave_forcing,
                             flow_direction,
                             flow_x_side,
                             flow_y_side,
                             dim2d_lowTemperature):
    '''
    @description
    update the constants defined in the 'parameters_input'
    file
    '''
    
    # change the constant that do not require a special
    # output treatment (integer,character...)
    constants_changed1={
        'npx':inputs['npx'],
        'npy':inputs['npy'],
        'ntx':ntx,
        'nty':nty,
        'ne':ne,
        'pm_choice':pm_choice,
        'ic_choice':ic_choice,
        'ic_perturbation_ac':ic_perturbation_ac,
        'li_perturbation_ac':li_perturbation_ac,
        'bc_choice':bc_choice,
        'bf_openbc_md_threshold_ac':openbc_md_threshold_ac,
        'obc_dct_distance':inputs['openbc_detector_distance'],
        'obc_perturbation_T0_ac':openbc_perturbation_T0_ac,
        'obc_perturbation_vx0_ac':openbc_perturbation_vx0_ac,
        'obc_perturbation_vy0_ac':openbc_perturbation_vy0_ac,
        'bc_N_type_choice':bc_N_type_choice,
        'bc_S_type_choice':bc_S_type_choice,
        'bc_E_type_choice':bc_E_type_choice,
        'bc_W_type_choice':bc_W_type_choice,
        'gravity_choice':gravity_choice,
        'wave_forcing':wave_forcing,
        'flow_direction':flow_direction,
        'dim2d_lowTemperature':dim2d_lowTemperature}

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
        'x_min'                   : inputs['x_min'],
        'x_max'                   : inputs['x_max'],
        'y_min'                   : inputs['y_min'],
        'y_max'                   : inputs['y_max'],
        't_max'                   : inputs['t_max'],
        'dt'                      : inputs['dt'],
        'detail_print'            : inputs['detail_print'],
        'flow_x_side'             : flow_x_side,
        'flow_y_side'             : flow_y_side,
        'flow_velocity'           : inputs['flow_velocity'],
        'T0'                      : inputs['temperature'],
        'ic_perturbation_amp'     : inputs['ic_perturbation_amp'],
        'li_perturbation_amp'     : inputs['li_perturbation_amp'],
        'bf_openbc_md_threshold'  : inputs['openbc_md_threshold'],
        'obc_perturbation_T0_amp' : inputs['openbc_perturbation_T0_amp'],
        'obc_perturbation_vx0_amp': inputs['openbc_perturbation_vx0_amp'],
        'obc_perturbation_vy0_amp': inputs['openbc_perturbation_vy0_amp']}

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
def update_makefile(file_path,bc_choice):
    '''
    @description
    update the folder for the compilation
    of the boundary conditions
    '''

    # define the constants changed in the file
    constants_changed={
        'pm_choice':pm_choice,
        'ic_choice':ic_choice,
        'bc_choice':bc_choice}


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
    [inputs,ntx,nty,ne,
     pm_choice,
     ic_choice,
     ic_perturbation_ac,
     li_perturbation_ac,
     bc_choice,
     openbc_md_threshold_ac,
     openbc_perturbation_T0_ac,
     openbc_perturbation_vx0_ac,
     openbc_perturbation_vy0_ac,
     bc_N_type_choice,
     bc_S_type_choice,
     bc_E_type_choice,
     bc_W_type_choice,
     gravity_choice,
     wave_forcing,
     flow_direction,
     flow_x_side,
     flow_y_side,
     dim2d_lowTemperature]=compute_code_inputs(inputFileName,nbTiles)


    # replace the inputs in the 'parameters_input' file
    update_parameters_inputs(param_path,
                             inputs,
                             ntx,nty,ne,
                             pm_choice,
                             ic_choice,
                             ic_perturbation_ac,
                             li_perturbation_ac,
                             bc_choice,
                             openbc_md_threshold_ac,
                             openbc_perturbation_T0_ac,
                             openbc_perturbation_vx0_ac,
                             openbc_perturbation_vy0_ac,
                             bc_N_type_choice,
                             bc_S_type_choice,
                             bc_E_type_choice,
                             bc_W_type_choice,
                             gravity_choice,
                             wave_forcing,
                             flow_direction,
                             flow_x_side,
                             flow_y_side,
                             dim2d_lowTemperature)


    # replace the inputs in the 'makefile'
    update_makefile(makefile_path,bc_choice)


    # replace the commit SHA number in the
    # 'parameters_constant'
    set_commit(param_cst_path)


    # print the major results
    print ''
    print '(ntx,nty,ne)', ntx,nty,ne


    # print the end of the configuration
    print ''
    print 'end of configuration'
    print ''

    
    # compile the code
    if(compileCode):
        compile_code(inputs,compileCodeBuffer)
        
