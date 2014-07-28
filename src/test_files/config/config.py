#!/usr/bin/python

import sys
import getopt
import subprocess
import shlex


def display_help():
    '''
    @description:
    display program help
    '''
    print ''
    print 'configure the augeanstables code with the input'
    print './config -i <input_file>'
    print 'configure the augeanstables code and compile'
    print './config -i <input_file> -c'
    print ''


def parse_argv(argv):
    '''
    @description:
    parse the program arguments: get the input file
    '''
    
    #< store the options and the arguments in opts, args
    try:
        opts, args = getopt.getopt(argv,"hi:c", ["help","input="])
    except getopt.GetoptError:
        display_help()
        sys.exit(2)

    if(len(opts)==0):
        display_help()
        sys.exit(2)

    compileCode=False

    for opt, arg in opts:

        if opt == '-h':
            display_help()
            sys.exit(2)

        elif opt in ("-i", "--input"):
            inputFile = arg
            print 'input file: ', arg

        elif opt in ("-c"):
            compileCode=True

    return [inputFile,compileCode]


def read_inputs(filename, inputs_needed):
    '''
    @description:
    read the input file
    '''    
    #< define a dictionnary to store the inputs needed
    #> then, if one needs the input read 'x_min', one simply
    #> use : inputs_read['x_min']
    inputs_read={}

    for input_param in inputs_needed:
        cmd="./get_parameter.sh -i "+str(filename)+" -p "+input_param
        args = shlex.split(cmd)
        output = subprocess.Popen(args,stdout=subprocess.PIPE).communicate()[0]
        inputs_read[input_param]=float(output)

    return inputs_read


def compute_n(x_min,x_max,dx,bc_size):
    '''
    @description:
    compute the number of gridpoints for the tile
    such that the space step coresponds to dx
    '''
    n = int((x_max-x_min)/dx+2*bc_size)
    return n


def compute_n_par(npx,x_min,x_max,dx,bc_size):
    '''
    @description:
    compute the number of gridpoints for the tile
    such that the space step coresponds to dx
    '''
    n = int(1./npx*( (x_max-x_min)/dx ) +2*bc_size)
    return n


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

    ntx=int(npx*nx)
    nty=int(npy*ny)

    return [ntx,nty]


def compute_code_inputs(inputFileName):
    '''
    @description
    compute all the inputs needed by the code
    '''

    #< codes for ic_choice, bc_choice, bc_type_choice,
    #> gravity_choice: these codes are defined in
    #> parameters_constant.f
    pm_code      = ['simpletest_choice',
                    'wave2d_choice',
                    'dim2d_choice']
    
    ic_code      = ['steady_state',
                    'drop_retraction',
                    'bubble_ascending',
                    'homogeneous_liquid',
                    'drop_collision',
                    'phase_separation']

    bc_code      = ['periodic_xy_choice',
                    'reflection_xy_choice',
                    'wall_xy_choice',
                    'wall_x_reflection_y_choice']

    bc_type_code = ['bc_nodes_choice',
                    'bc_fluxes_choice']

    gravity_code = ['no_gravity_choice',
                    'earth_gravity_choice']


    #< read the input file
    inputs_needed=['x_min','x_max','dx',
                   'y_min','y_max','dy',
                   'dt','t_max','detail_print',
                   'npx', 'npy',
                   'pm_choice',
                   'bc_choice',
                   'ic_choice',
                   'gravity_choice']
    inputs=read_inputs(inputFileName, inputs_needed)
    

    #< update the type of the inputs
    inputs['npx']=int(inputs['npx'])
    inputs['npy']=int(inputs['npy'])


    #< compute the ntx and nty determining the
    #> extent of the computational domain
    bc_size=2
    [ntx,nty]=compute_ntx_and_nty(
        inputs['npx'], inputs['npy'],
        inputs['x_min'],inputs['x_max'],inputs['dx'],
        inputs['y_min'],inputs['y_max'],inputs['dy'],
        bc_size)

    #< compute the pm_choice
    pm_choice = pm_code[int(inputs['pm_choice'])]

    #< compute the ne
    if(pm_choice=='simpletest_choice'):
        ne = 1
    if(pm_choice=='wave2d_choice'):
        ne = 3
    if(pm_choice=='dim2d_choice'):
        ne = 4

    #< compute the ic_choice
    ic_choice = ic_code[int(inputs['ic_choice'])]
    
    #< compute the bc_choice    
    bc_choice = bc_code[int(inputs['bc_choice'])]
    
    #< compute the bc_type_choice
    if(bc_choice=='periodic_xy_choice' or
       bc_choice=='reflection_xy_choice'):

        bcx_type_choice = bc_type_code[0]
        bcy_type_choice = bc_type_code[0]

    if(bc_choice=='wall_xy_choice' or
       bc_choice=='wall_x_reflection_y_choice'):

        bcx_type_choice = bc_type_code[1]
        bcy_type_choice = bc_type_code[1]
    
    #< compute the gravity_choice
    gravity_choice = gravity_code[int(inputs['gravity_choice'])]


    return [inputs,ntx,nty,ne,
            pm_choice,
            ic_choice,
            bc_choice,
            bcx_type_choice,bcy_type_choice,
            gravity_choice]


def update_parameters_inputs(file_path,inputs,ntx,nty,ne,
                             pm_choice,
                             ic_choice,
                             bc_choice,
                             bcx_type_choice,bcy_type_choice,
                             gravity_choice):
    '''
    @description
    update the constants defined in the 'parameters_input'
    file
    '''
    
    #< change the constant that do not require a special
    #> output treatment (double,real,integer...)
    constants_changed1={
        'npx':inputs['npx'],
        'npy':inputs['npy'],
        'ntx':ntx,
        'nty':nty,
        'ne':ne,
        'pm_choice':pm_choice,
        'ic_choice':ic_choice,
        'bc_choice':bc_choice,
        'bcx_type_choice':bcx_type_choice,
        'bcy_type_choice':bcy_type_choice,
        'gravity_choice':gravity_choice}

    for key, value  in constants_changed1.items():

        cmd="./change_parameter.sh"
        cmd+=" -i "+str(file_path)
        cmd+=" -o "+str(file_path)
        cmd+=" -p "+key
        cmd+=" -v "+str(value)
        subprocess.call(cmd, shell=True)


    #< change the constant that do require a special
    #> output treatment (output format)
    constants_changed2={
        'x_min':inputs['x_min'],
        'x_max':inputs['x_max'],
        'y_min':inputs['y_min'],
        'y_max':inputs['y_max'],
        't_max':inputs['t_max'],
        'dt':inputs['dt'],
        'detail_print':inputs['detail_print']}

    for key, value in constants_changed2.items():

        cmd="./change_parameter.sh"
        cmd+=" -i "+str(file_path)
        cmd+=" -o "+str(file_path)
        cmd+=" -p "+key
        cmd+=" -v "+"%10.10fd0"%value
        subprocess.call(cmd, shell=True)    

    print 'update ', file_path        


def update_makefile(file_path,bc_choice):
    '''
    @description
    update the folder for the compilation
    of the boundary conditions
    '''

    #< define the constants changed in the file
    constants_changed={'bc_choice':bc_choice}


    #< change the constant in the file
    for key, value in constants_changed.items():

        cmd="./change_parameter.sh"
        cmd+=" -i "+str(file_path)
        cmd+=" -o "+str(file_path)
        cmd+=" -p "+key
        cmd+=" -v "+value
        subprocess.call(cmd, shell=True)

    print 'update ', file_path


def compile_code(inputs):
    '''
    @description
    compile the code using the previous configuration
    '''

    fname='sim_dim2d'
    fname+='_'+str(inputs['npx'])+'x'+str(inputs['npy'])

    cmd_serial  ='cd .. && make cleanall && make sim_dim2d && make clean'
    cmd_parallel='cd .. && make cleanall && make sim_dim2d_par && make clean'
    cmd_parallel+=' && mv sim_dim2d_par '+fname

    if(inputs['npx']*inputs['npy']==1):
        cmd      = cmd_serial
        name_exe = 'sim_dim2d'
    else:
        cmd      = cmd_parallel
        name_exe = fname
        
    subprocess.call(cmd, shell=True)

    print ''
    print 'executable ready: '+name_exe
    print ''

if __name__ == "__main__":

    #< define the paths for the files modified by the
    #> configuration
    sim_paths={}
    sim_paths['serial']  = '../sim_dim2d.f'
    sim_paths['parallel']= '../sim_dim2d_par.f'

    param_path           = '../../parameters/parameters_input.f'
    makefile_path        = './makefile_header.mk'


    #< parse the program arguments
    [inputFileName,compileCode]=parse_argv(sys.argv[1:])


    #< compute the code inputs
    [inputs,ntx,nty,ne,
     pm_choice,
     ic_choice,
     bc_choice,
     bcx_type_choice,bcy_type_choice,
     gravity_choice]=compute_code_inputs(inputFileName)


    #< replace the inputs in the 'parameters_input' file
    update_parameters_inputs(param_path,inputs,ntx,nty,ne,
                             pm_choice,
                             ic_choice,
                             bc_choice,
                             bcx_type_choice,bcx_type_choice,
                             gravity_choice)


    #< replace the inputs in the 'makefile'
    update_makefile(makefile_path,bc_choice)


    #< print the major results
    print '(ntx,nty)', ntx,nty
    print '(ne)', ne


    #< print the end of the configuration
    print ''
    print 'end of configuration'
    print ''
    
    #< compile the code
    if(compileCode):
        compile_code(inputs)
        
