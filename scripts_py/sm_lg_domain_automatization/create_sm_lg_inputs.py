#!/usr/bin/python

'''
@description
create two input files to run the DIM2D simulation on two domains:
inputs_sm_domain.txt:  inputs for the small domain simulation
inputs_lg_domain.txt:  inputs for the large domain simulation
'''

import os
import sys
import getopt
import subprocess
import shlex
import shutil

from automatization_csts import (md_threshold_ac_default,
                                 md_threshold_default,
                                 ic_perturbation_ac_default,
                                 ic_perturbation_amp_default,
                                 bc_perturbation_T0_ac_default,
                                 bc_perturbation_T0_amp_default,
                                 bc_perturbation_vx0_ac_default,
                                 bc_perturbation_vx0_amp_default,
                                 bc_perturbation_vy0_ac_default,
                                 bc_perturbation_vy0_amp_default,
                                 nb_pts_in_interface_default,
                                 ratio_bubble_interface_default,
                                 CFL_constant_default,
                                 ratio_interface_influence_default,
                                 total_nb_files_default)

from library_sm_lg_inputs import (get_we,
                                  get_interface_length,
                                  get_bubble_diameter,
                                  get_domain_length,
                                  get_interface_space_step,
                                  get_small_domain_extent,
                                  get_cv_r,
                                  get_max_speed_of_sound,
                                  get_dt_max,
                                  get_simulation_time,
                                  get_detail_print,
                                  get_large_domain_extent)
                                       

changeParameterPath = os.path.join(os.getenv('augeanstables'),
                                   'src',
                                   'config',
                                   'change_parameter.sh')
getParameterPath    = os.path.join(os.getenv('augeanstables'),
                                   'src',
                                   'config',
                                   'get_parameter.sh')
                                   
# display the help for the program
def display_help():
    '''
    @description:
    display program help
    '''
    print ''
    print 'configure two inputs files to run DIM2D on a small and'
    print 'a large domain'
    print ''
    print 'options:'
    print '--------'
    print '-h (--help)             : display this help'
    print '-T (--temperature=)     : temperature'
    print '-v (--flow_velocity=)   : flow_velocity'
    print '-i (--model_input=)     : input file used as template'
    print '-s (--sm_domain=)       : filename for the small domain input'
    print '-l (--lg_domain=)       : filename for the large domain input'
    print '--md_threshold_ac=      : activate or not the mass density threshold'
    print '--md_threshold=         : value set for the mass density threshold'
    print ''
    print 'example:'
    print '---------'
    print 'configure input files at T=0.95 and flow_velocity=0.1 using'
    print 'template.txt as a template for the input file. The default'
    print 'outputs are inputs_sm_domain.txt and inputs_lg_domain.txt'
    print './config -i template.txt - T 0.95 -v 0.1'
    print ''


#parse the input options
def parse_argv(argv):
    '''
    @description:
    parse the program options
    '''

    #default values
    smDomainInputDefault = 'inputs_sm_domain.txt'
    lgDomainInputDefault = 'inputs_lg_domain.txt'


    # store the options and the arguments
    # in opts, args
    try:
        opts, args = getopt.getopt(argv,
                                   "hT:v:i:s:l:",
                                   ["help",
                                    "temperature=",
                                    "flow_velocity=",
                                    "model_input=",
                                    "sm_domain=",
                                    "lg_domain=",
                                    "md_threshold_ac=",
                                    "md_threshold="])
    except getopt.GetoptError:
        display_help()
        sys.exit(2)

    if(len(opts)==0):
        display_help()
        sys.exit(2)


    # initialization of the boolean controlling if
    # all the needed inputs are provided
    temperatureProvided  = False
    flowVelocityProvided = False
    modelInputProvided   = False

    md_threshold_ac = 0
    md_threshold    = 0.0

    smDomainInput = smDomainInputDefault
    lgDomainInput = lgDomainInputDefault


    # analyse the program options
    for opt, arg in opts:

        if opt in ("-h", "--help"):
            display_help()
            sys.exit(2)

        elif opt in ("-T", "--temperature"):
            temperature = float(arg)

            if(temperature>0 and temperature<1):
                temperatureProvided=True
            else:
                print 'the temperature should be in ]0,1['
            
        elif opt in ("-v","--flow_velocity"):
            flow_velocity = float(arg)

            if(flow_velocity>0 and flow_velocity<1):
                flowVelocityProvided=True
            else:
                print 'the flow velocity should be in ]0,1['

        elif opt in ("-i","--model_input"):
            modelInputPath = arg

            if(os.path.isfile(modelInputPath)):
                modelInputProvided=True
            else:
                print modelInputPath+' does not exist'

        elif opt in ("-s", "--sm_domain"):
            smDomainInput = arg

        elif opt in ("-l", "--lg_domain"):
            lgDomainInput = arg

        elif opt in ("--md_threshold"):
            md_threshold = float(arg)

        elif opt in ("--md_threshold_ac"):
            md_threshold_ac = int(arg)        

    inputsProvided = temperatureProvided
    inputsProvided = inputsProvided and flowVelocityProvided
    inputsProvided = inputsProvided and modelInputProvided

    if(not inputsProvided):
        display_help()
        sys.exit('***some inputs were not provided***')

    else:
        inputs = {'temperature'    : temperature,
                  'flow_velocity'  : flow_velocity,
                  'model_input'    : modelInputPath,
                  'sm_domain'      : smDomainInput,
                  'lg_domain'      : lgDomainInput,
                  'md_threshold_ac': md_threshold_ac,
                  'md_threshold'   : md_threshold}
                    
        return inputs


# extract 'value' for parameter 'param' from the text
# file 'filePath' if the definition of 'value' is as
# follow 'param = value'
def get_parameter(param,filePath):
    '''
    @description:
    extract 'value' for parameter 'param' from the text
    file 'filePath' if the definition of 'value' is as
    follow 'param = value'
    '''

    # check whether the file exists
    if(os.path.isfile(filePath)):

        # check whether the 'get_parameter.sh'
        # script exists
        if(os.path.isfile(getParameterPath)):
               
           # create process extracting the parameter
           # using the 'get_parameter.sh' script
           cmd=getParameterPath+' -p '+param+' -i '+filePath
           args = shlex.split(cmd)
           output = subprocess.Popen(args,stdout=subprocess.PIPE).communicate()[0]

        else:

            print '***get_parameter.sh script does not exist***'
            output = ''
        
    else:

        print '***'+filePath+' does not exist***'
        output = ''

    return output.strip()
    

# compute the inputs that are modified in the template.txt
# template input file to run the DIM2D simulation on a small
# and on a large domains
def get_inputsToBeModified(temperature,
                           flow_velocity,
                           md_threshold_ac,
                           md_threshold,
                           ic_perturbation_ac,
                           ic_perturbation_amp,
                           bc_perturbation_T0_ac,
                           bc_perturbation_T0_amp,
                           bc_perturbation_vx0_ac,
                           bc_perturbation_vx0_amp,
                           bc_perturbation_vy0_ac,
                           bc_perturbation_vy0_amp,
                           nb_pts_in_interface,
                           ratio_bubble_interface,
                           CFL_constant,
                           ratio_interface_influence,
                           total_nb_files):
    '''
    @description:
    compute the inputs that are modified in the template.txt
    template input file to run the DIM2D simulation on a small
    and on a large domains
    '''

    # extract length_c, dim2d_a, dim2d_b, dim2d_M, dim2d_cv, dim2d_R
    # and dim2d_K from the dim2d_parameters.f fortran file
    dim2dParamPath = os.path.join(os.getenv('augeanstables'),
                                  'src',
                                  'physical_models',
                                  'dim2d',
                                  'dim2d_parameters.f')

    if(os.path.isfile(dim2dParamPath)):
        length_c  = float(get_parameter('length_c', dim2dParamPath))
        dim2d_a   = float(get_parameter( 'dim2d_a', dim2dParamPath))
        dim2d_b   = float(get_parameter( 'dim2d_b', dim2dParamPath))
        dim2d_M   = float(get_parameter( 'dim2d_M', dim2dParamPath))
        dim2d_cv  = float(get_parameter('dim2d_cv', dim2dParamPath))
        dim2d_R   = float(get_parameter( 'dim2d_R', dim2dParamPath))
        dim2d_K   = float(get_parameter( 'dim2d_K', dim2dParamPath))

    else:
        sys.exit('*** '+dim2dParamPath+' does not exist***')


    # compute the Weber number
    we = get_we(length_c, dim2d_a, dim2d_b, dim2d_M, dim2d_K)


    # compute the interface length from the temperature
    interface_lgh = get_interface_length(we,temperature)


    # compute the bubble diameter from the interface length
    bubble_diameter = get_bubble_diameter(interface_lgh,
                                          ratio_bubble_interface)

    # compute the domain length from the bubble diameter
    domain_length = get_domain_length(bubble_diameter,
                                      interface_lgh,
                                      ratio_interface_influence)


    # compute the maximum space step from the interface length
    dx_max = get_interface_space_step(interface_lgh,
                                      nb_pts_in_interface)


    # compute the extent of the small domain as a matrix:
    # x_min : small_domain_extent[0][0]
    # x_max : small_domain_extent[1][0]
    # y_min : small_domain_extent[0][1]
    # y_max : small_domain_extent[1][1]
    small_domain_extent = get_small_domain_extent(domain_length,
                                                  dx_max)


    # compute the reduced heat capacity
    cv_r = get_cv_r(dim2d_M,dim2d_cv,dim2d_R)


    # compute the maximum speed of sound in the flow
    speed_of_sound = get_max_speed_of_sound(temperature,cv_r)


    # compute the maximum time step ensuring numerical stability
    speed_max = speed_of_sound + abs(flow_velocity)
    dt_max    = get_dt_max(dx_max,speed_max,CFL_constant)


    # determine the simulation time needed to let the bubble
    # leave the computational domain
    simulation_time = get_simulation_time(domain_length,
                                          bubble_diameter,
                                          interface_lgh,
                                          ratio_interface_influence,
                                          flow_velocity)


    # determine the detail print
    detail_print = get_detail_print(total_nb_files,
                                    simulation_time,
                                    dt_max)


    # determine the large domain extent
    large_domain_extent = get_large_domain_extent(small_domain_extent,
                                                  dx_max,dx_max,
                                                  flow_velocity,
                                                  speed_of_sound,
                                                  simulation_time)


    # gather the inputs to be modified in dictionnaries
    inputsToBeModified_sm_domain = {
        'detail_print'                : detail_print,
        'dt'                          : dt_max,
        't_max'                       : simulation_time,
        'dx'                          : dx_max,
        'x_min'                       : small_domain_extent[0][0],
        'x_max'                       : small_domain_extent[1][0],
        'dy'                          : dx_max,
        'y_min'                       : small_domain_extent[0][1],
        'y_max'                       : small_domain_extent[1][1],
        'flow_velocity'               : flow_velocity,
        'temperature'                 : temperature,
        'openbc_md_threshold_ac'      : md_threshold_ac,
        'openbc_md_threshold'         : md_threshold,
        'ic_perturbation_ac'          : ic_perturbation_ac,
        'ic_perturbation_amp'         : ic_perturbation_amp,
        'openbc_perturbation_T0_ac'   : bc_perturbation_T0_ac,
        'openbc_perturbation_T0_amp'  : bc_perturbation_T0_amp,
        'openbc_perturbation_vx0_ac'  : bc_perturbation_vx0_ac,
        'openbc_perturbation_vx0_amp' : bc_perturbation_vx0_amp,
        'openbc_perturbation_vy0_ac'  : bc_perturbation_vy0_ac,
        'openbc_perturbation_vy0_amp' : bc_perturbation_vy0_amp}

    inputsToBeModified_lg_domain = {
        'detail_print'                : detail_print,
        'dt'                          : dt_max,
        't_max'                       : simulation_time,
        'dx'                          : dx_max,
        'x_min'                       : large_domain_extent[0][0],
        'x_max'                       : large_domain_extent[1][0],
        'dy'                          : dx_max,
        'y_min'                       : large_domain_extent[0][1],
        'y_max'                       : large_domain_extent[1][1],
        'flow_velocity'               : flow_velocity,
        'temperature'                 : temperature,
        'openbc_md_threshold_ac'      : 0,
        'openbc_md_threshold'         : 0.0,
        'ic_perturbation_ac'          : 0,
        'ic_perturbation_amp'         : 0.0,
        'openbc_perturbation_T0_ac'   : 0,
        'openbc_perturbation_T0_amp'  : 0.0,
        'openbc_perturbation_vx0_ac'  : 0,
        'openbc_perturbation_vx0_amp' : 0.0,
        'openbc_perturbation_vy0_ac'  : 0,
        'openbc_perturbation_vy0_amp' : 0.0}

        
    return [inputsToBeModified_sm_domain,
            inputsToBeModified_lg_domain]


# modify the parameters in the template.txt file to create 
# a new input file
def create_inputFile(paramModified,
                     templateFilePath,
                     newInputFilePath):
    '''
    @description:
    modify the parameters in the template.txt file to
    create a new input file
    '''    
    

    #1) verify that the templateFilePath exists
    if(os.path.isfile(templateFilePath)):


        #2) make a copy of the template file for the output
        #   file newInputFilePath
        if(os.path.isfile(newInputFilePath)):
            os.remove(newInputFilePath)
        
        shutil.copyfile(templateFilePath,newInputFilePath)

        
        #3) read the parameters to be modified in paramModified
        #   and replace them in newInputFilePath
        for key, value  in paramModified.items():

            cmd=os.path.join(changeParameterPath)
            cmd+=" -i "+str(newInputFilePath)
            cmd+=" -o "+str(newInputFilePath)
            cmd+=" -p "+key
            cmd+=" -v "+str(value)
            subprocess.call(cmd, shell=True)

    else:
        sys.exit('***the file '+templateFilePath+' does not exist***')
        

    return


#create the small and large domain inputs for the simulation
def create_sm_lg_inputs(temperature,
                        flow_velocity,
                        model_input,
                        sm_domain                 = 'inputs_sm_domain.txt',
                        lg_domain                 = 'inputs_lg_domain.txt',
                        md_threshold_ac           = md_threshold_ac_default,
                        md_threshold              = md_threshold_default,
                        ic_perturbation_ac        = ic_perturbation_ac_default,
                        ic_perturbation_amp       = ic_perturbation_amp_default,
                        bc_perturbation_T0_ac     = bc_perturbation_T0_ac_default,
                        bc_perturbation_T0_amp    = bc_perturbation_T0_amp_default,
                        bc_perturbation_vx0_ac    = bc_perturbation_vx0_ac_default,
                        bc_perturbation_vx0_amp   = bc_perturbation_vx0_amp_default,
                        bc_perturbation_vy0_ac    = bc_perturbation_vy0_ac_default,
                        bc_perturbation_vy0_amp   = bc_perturbation_vy0_amp_default,
                        nb_pts_in_interface       = nb_pts_in_interface_default,
                        ratio_bubble_interface    = ratio_bubble_interface_default,
                        CFL_constant              = CFL_constant_default,
                        ratio_interface_influence = ratio_interface_influence_default,
                        total_nb_files            = total_nb_files_default):
    '''
    @description:
    create the small and large domain inputs for the simulation    
    '''

    # determine the inputs to be modified in
    # template.txt to run the simulation on
    # small and large scale domains
    [inputs_sm_domain, inputs_lg_domain] = get_inputsToBeModified(
        temperature,
        flow_velocity,
        md_threshold_ac,
        md_threshold,
        ic_perturbation_ac,
        ic_perturbation_amp,
        bc_perturbation_T0_ac,
        bc_perturbation_T0_amp,
        bc_perturbation_vx0_ac,
        bc_perturbation_vx0_amp,
        bc_perturbation_vy0_ac,
        bc_perturbation_vy0_amp,
        nb_pts_in_interface,
        ratio_bubble_interface,
        CFL_constant,
        ratio_interface_influence,
        total_nb_files)
    

    # create the input file for the small
    # domain simulation
    create_inputFile(inputs_sm_domain,
                     model_input,
                     sm_domain)


    # create the input file for the large
    # domain simulation
    create_inputFile(inputs_lg_domain,
                     model_input,
                     lg_domain)

    return


if __name__=="__main__":


    # extract the command line inputs
    #inputs = parse_argv(sys.argv[1:])

    inputs = {}
    
    inputs['temperature']     = 0.99
    inputs['flow_velocity']   = 0.05
    inputs['model_input']     = os.path.join(os.getenv('augeanstables'),
                                             'src','config','default_inputs','dim2d',
                                             'dim2d_bubble_transported_hedstrom_xy.txt')

    inputs['sm_domain']           = 'inputs_sm_domain.txt'
    inputs['lg_domain']           = 'inputs_lg_domain.txt'
    inputs['md_threshold_ac']     = 1
    inputs['md_threshold']        = 0.2
    inputs['ic_perturbation_ac']  = 0
    inputs['ic_perturbation_amp'] = 0.0

    inputs['bc_perturbation_T0']      = 0
    inputs['bc_perturbation_T0_amp']  = 0.0
    inputs['bc_perturbation_vx0']     = 0
    inputs['bc_perturbation_vx0_amp'] = 0.0
    inputs['bc_perturbation_vy0']     = 0
    inputs['bc_perturbation_vy0_amp'] = 0.0


    # create the inputs
    create_sm_lg_inputs(inputs['temperature'],
                        inputs['flow_velocity'],
                        inputs['model_input'],
                        sm_domain=inputs['sm_domain'],
                        lg_domain=inputs['lg_domain'],
                        md_threshold_ac=inputs['md_threshold_ac'],
                        md_threshold=inputs['md_threshold'],
                        ic_perturbation_ac=inputs['ic_perturbation_ac'],
                        ic_perturbation_amp=inputs['ic_perturbation_amp'],
                        bc_perturbation_T0_ac=inputs['bc_perturbation_T0'],
                        bc_perturbation_T0_amp=inputs['bc_perturbation_T0_amp'],
                        bc_perturbation_vx0_ac=inputs['bc_perturbation_vx0'],
                        bc_perturbation_vx0_amp=inputs['bc_perturbation_vx0_amp'],
                        bc_perturbation_vy0_ac=inputs['bc_perturbation_vy0'],
                        bc_perturbation_vy0_amp=inputs['bc_perturbation_vy0_amp'])

    
    
