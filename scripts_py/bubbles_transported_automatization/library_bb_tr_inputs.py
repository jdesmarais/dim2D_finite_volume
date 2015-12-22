#!/usr/bin/python

'''
@description
useful functions to compute the input parameters for simulations with
two bubbles transported by the flow
'''

import os
import sys  #for system fcts
import inspect
import math #for mathematical fcts


#add the python files from ../sm_lg_automatization
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],
    "../sm_lg_domain_automatization")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)


from automatization_csts import\
    nb_pts_inside_interface_required


# get the total size of the domain
def get_domain_sizes(bubble_diameter,
                     interface_length,
                     ratio_interface_influence,
                     ratio_interface_separation):
    """Determine the length of the domain in the x-
    and y- directions

    Args:
        bubble_diameter (double) : diameter of the bubble
        interface_length (double) : width of the interface
        ratio_interface_influence (double) : length threshold
           (expressed as a fraction of the width of the interface)
           above which the interface is not supposed to interact
           with the border
        ratio_interface_separation (double) : length (expressed as
           a fraction of the width of the interface) separating the
           two bubbles

    Returns:
        [Lx,Ly]: size of the domain in the x- and y-directions
           - Lx (double): length of the domain in the x-direction
           - Ly (double): length of the domain in the y-direction

    """
    
    Lx = 2.0*bubble_diameter+\
        (ratio_interface_separation+2.0*ratio_interface_influence)*interface_length

    Ly = bubble_diameter+2.0*ratio_interface_influence*interface_length

    return [Lx,Ly]


# compute the extent of the simulation domain
def get_domain_extent(Lx,Ly,dx,dy):
    """Determine the extent of the simulation domain
    [[x_min,x_max],[y_min,y_max]]

    Args:
        Lx (double): length of the domain in the x-direction
        Ly (double): length of the domain in the y-direction
        dx (double): grid spacing in the x-direction
        dy (double): grid spacing in the y-direction

    Returns:
        domain_extent (2x2 double list)
           - domain_extent[0][0] : x_min
           - domain_extent[1][0] : x_max
           - domain_extent[0][1] : y_min
           - domain_extent[1][1] : y_max

    """

    nx_half = math.ceil(Lx/(2.0*dx))
    ny_half = math.ceil(Ly/(2.0*dy))
    
    domain_extent = [[0 for x in range(2)] for x in range(2)]

    domain_extent[0][0] = - nx_half*dx
    domain_extent[1][0] = + nx_half*dx
    domain_extent[0][1] = - ny_half*dy
    domain_extent[1][1] = + ny_half*dy
    
    return domain_extent


# compute the t_max of the simulation
def get_simulation_time(Lx,
                        bubble_diameter,
                        bubble_separation_lgh,
                        interface_influence_lgh,
                        flow_velocity):
    """Determine the total duration of the simulation
    to let the two bubbles leave the computational domain

    Args:
        Lx (double): length of the domain in the direction of the flow
        bubble_diameter (double) : diameter of the bubbles in the domain
        bubble_separation_lgh (double) : separation length between the
           bubbles
        interface_influence_lgh (double) : threshold length above which 
           the interface is not interacting with a object at this length 
        flow_velocity (double) : mean flow velocity

    Returns:
        simulation time t_max
    """

    total_length = 0.5*(bubble_separation_lgh+Lx) +\
        bubble_diameter + 2.0*interface_influence_lgh

    simulation_time = float(math.ceil(
            total_length/flow_velocity*
            10**4))/(10**4)

    return simulation_time


