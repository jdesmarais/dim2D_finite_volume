#!/usr/bin/python

'''
@description: useful functions to compute quantities
from the Van der Waals equation of state
'''

from math import sqrt


def compute_latentHeat(temperature):
    '''
    @description: approximate the latent heat
    close to the critical temperature
    '''

    return 15.08*sqrt(1.0-temperature)
