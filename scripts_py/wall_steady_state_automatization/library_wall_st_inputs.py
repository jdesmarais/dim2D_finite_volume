#!/usr/bin/python

'''
@description
create an input file to run the DIM2D simulation with a wall
inputs_wall.txt:  inputs for the wall simulation
'''

import os
import sys
import getopt
import subprocess
import shlex
import shutil

from math import pi, sqrt, sin, cos


def get_equilibrium_length(micro_contact_angle,initial_radius):
    '''
    @description
    compute the extent of the droplet when the equilibrium
    is reached. The initial shape is half a disc. At equilibrium, 
    a disc-cut approximation is used where the micro-contact angle
    is reached at the borders
    '''
    
    contact_angle_vap = pi*(1.0-micro_contact_angle/180.0)

    length = initial_radius*sqrt(pi/2.0)*\
        sin(contact_angle_vap)/\
        sqrt(contact_angle_vap-sin(2.0*contact_angle_vap)/2.0)
        

    return length


def get_x_max(initial_radius,equilibrium_length,ratio_eq_length_domain):
    '''
    @description
    compute the extent of the droplet when the equilibrium
    is reached. The initial shape is half a disc. At equilibrium, 
    a disc-cut approximation is used where the micro-contact angle
    is reached at the borders
    '''

    x_max = max(4.0*initial_radius,equilibrium_length+2.0*initial_radius)

    return x_max


def get_y_max(radius,ratio_drop_length_domain):
    '''
    @description
    compute the domain height as the radius of the droplet
    with a safety margin
    '''

    y_max = radius*ratio_drop_length_domain

    return y_max


def get_wall_domain_extent(x_max,y_max,dx_max):
    '''
    @description
    compute the extent of the domain to compute the steady
    state of a bubble
    '''

    # x_max of the computational domain as a integer
    # number of dx_max
    domain_x_max = int(x_max/dx_max)*dx_max
    if(domain_x_max<x_max):
        domain_x_max+=dx_max
    
    # y_max of the computational domain as a integer
    # number of dx_max
    domain_y_max = int(y_max/dx_max)*dx_max
    if(domain_y_max<y_max):
        domain_y_max+=dx_max

    # extent of the computational domain
    domain_extent = [[0 for x in range(2)] for x in range(2)]

    domain_extent[0][0] = 0.0
    domain_extent[1][0] = domain_x_max
    domain_extent[0][1] = 0.0
    domain_extent[1][1] = domain_y_max    

    return domain_extent
