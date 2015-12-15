#!/usr/bin/python

"""
Create the input file used to configure a simulation with two bubbles
transported by the flow. The distance between the two bubbles transported,
the temperature and the mean flow velocity will be varied
"""

#add the python files from ../sm_lg_automatization
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],
    "../sm_lg_domain_automatization")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)


from automatization_csts import\
    nb_pts_inside_interface_required

from create_sm_lg_inputs import\
    get_parameter


from library_bb_tr_inputs import\
    

# 
def compute_inputsToBeModified(temperature,
                               flow_velocity,
                               ratio_interface_separation,
                               ratio_bubble_interface,
                               ratio_interface_influence,
                               CFL_constant,
                               total_nb_files):

    """Determine the inputs to be modified in the input.txt
    file for the simulation with two bubbles transported by the
    mean flow

    Args:
        temperature (double) : mean flow temperature
        flow_velocity (double) : mean flow velocity
        ratio_interface_separation (double) : length (expressed as
            a fraction of the width of the interface) separating the
           two bubbles
        ratio_bubble_interface (double) : ratio between bubble
           diameter and interface width
        ratio_interface_influence (double) : length threshold
           (expressed as a fraction of the width of the interface)
           above which the interface is not supposed to interact
           with the border
        CFL_constant (double) : CFL threshold (0-1)
        total_nb_files (int) : total nb of files written

    Returns:
        inputsToBeModified (dict) :
           - dict[key] : name of the input to be modified
           - dict[value] : value of the input to be modified
    """

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
    Li = get_interface_length(we,temperature)

    
    # compute the bubble diameter
    bubble_diameter = Li*ratio_bubble_interface

    
    # compute the extent of the domain
    
    
    
