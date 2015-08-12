#!/usr/bin/python

import sys


# convert grayscale value into RGB color
def grayscale_to_RGB(grayscale_value):
    '''
    @description
    determine the RGB code for a color from
    its grayscale values
    '''
    
    if( (0<=grayscale_value) and (grayscale_value<=1) ):

        R = 1.0 - grayscale_value

    else:

        print 'error for grayscale value'
        print grayscale_value
        sys.exit(2)            

    return (R,R,R)
