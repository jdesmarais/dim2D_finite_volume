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
    print 'configure the lernaeanhydra_opt4 code with the input'
    print './config -i <input_file>'
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
        ny = compute_n_par(npy,x_min,y_max,dy,bc_size)

    ntx=int(npx*nx)
    nty=int(npy*ny)

    return [ntx,nty]


def compute_code_inputs(inputFileName):
    '''
    @description
    compute all the inputs needed by the code
    '''
    #< read the input file
    inputs_needed=[\
        'x_min','x_max','dx',\
        'y_min','y_max','dy',\
        'dt','t_max','detail_print',\
        'npx', 'npy',\
        'bc_choice', 'ic_choice']
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

    
    #< compute the bc_choice
    boundary_code = ['periodic_xy_choice', 'reflection_xy_choice']
    bc_choice = boundary_code[int(inputs['bc_choice'])]

    return [inputs,ntx,nty,bc_choice]


def update_parameters_inputs(file_path,inputs,ntx,nty,bc_choice):
    '''
    @description
    update the constants defined in the 'parameters_input'
    file
    '''
    
    constants_changed={
        'npx':inputs['npx'],
        'npy':inputs['npy'],
        'ntx':ntx,
        'nty':nty,
        'bc_choice':bc_choice}

    for key, value  in constants_changed.items():

        cmd="./change_parameter.sh"
        cmd+=" -i "+str(file_path)
        cmd+=" -o "+str(file_path)
        cmd+=" -p "+key
        cmd+=" -v "+str(value)
        subprocess.call(cmd, shell=True)

    print 'update ', file_path        


def update_sim_dim2d(sim_paths,inputs):
    '''
    @description
    '''

    #< choose the tile modified
    if(inputs['npx']*inputs['npy']>1):
        file_path=sim_paths['parallel']
    else:
        file_path=sim_paths['serial']


    #< define the constants changed in the file
    constants_changed={
        'x_min':inputs['x_min'],
        'x_max':inputs['x_max'],
        'y_min':inputs['y_min'],
        'y_max':inputs['y_max'],
        't_max':inputs['t_max'],
        'dt':inputs['dt'],
        'detail_print':inputs['detail_print']}    


    #< change the constants in the file
    for key, value in constants_changed.items():

        cmd="./change_parameter.sh"
        cmd+=" -i "+str(file_path)
        cmd+=" -o "+str(file_path)
        cmd+=" -p "+key
        cmd+=" -v "+"%fd0"%value
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

    cmd_serial  ='cd .. && make sim_dim2d'
    cmd_parallel='cd .. && make sim_dim2d_par && mv sim_dim2d_par '+fname

    if(inputs['npx']*inputs['npy']==1):
        cmd=cmd_serial
    else:
        cmd=cmd_parallel
        
    subprocess.call(cmd, shell=True)

    print ''
    print 'executable ready'
    print ''

if __name__ == "__main__":

    #< define the paths for the files modified by the
    #> configuration
    sim_paths={}
    sim_paths['serial']  = '../sim_dim2d.f'
    sim_paths['parallel']= '../sim_dim2d_par.f'

    param_path           = '../../parameters/parameters_input.f'
    makefile_path        = '../makefile'


    #< parse the program arguments
    [inputFileName,compileCode]=parse_argv(sys.argv[1:])


    #< compute the code inputs
    [inputs,ntx,nty,bc_choice]=compute_code_inputs(inputFileName)


    #< replace the inputs in the 'parameters_input' file
    update_parameters_inputs(param_path,inputs,ntx,nty,bc_choice)


    #< replace the inputs in the 'sim_dim2d' or 'sim_dim2d_par'
    update_sim_dim2d(sim_paths,inputs)


    #< replace the inputs in the 'makefile'
    update_makefile(makefile_path,bc_choice)  


    #< print the major results
    print '(ntx,nty)', ntx,nty


    #< print the end of the configuration
    print ''
    print 'end of configuration'
    print ''
    
    #< compile the code
    if(compileCode):
        compile_code(inputs)
        
    