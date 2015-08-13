#!/usr/bin/python

import sys
import os
import inspect

'''
@description
automatization constants used to determine the bubble
contours
'''

cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],"../../../sm_lg_domain_automatization")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from create_sm_lg_inputs import get_parameter

from library_sm_lg_inputs import (get_we,
                                  get_cv_r)

# extract length_c, dim2d_a, dim2d_b, dim2d_M, dim2d_cv, dim2d_R
# and dim2d_K from the dim2d_parameters.f fortran file
dim2dParamPath = os.path.join(os.getenv('augeanstables'),
                              'src',
                              'physical_models',
                              'dim2d',
                              'dim2d_parameters.f')

length_c  = float(get_parameter('length_c', dim2dParamPath))
dim2d_a   = float(get_parameter( 'dim2d_a', dim2dParamPath))
dim2d_b   = float(get_parameter( 'dim2d_b', dim2dParamPath))
dim2d_M   = float(get_parameter( 'dim2d_M', dim2dParamPath))
dim2d_cv  = float(get_parameter('dim2d_cv', dim2dParamPath))
dim2d_R   = float(get_parameter( 'dim2d_R', dim2dParamPath))
dim2d_K   = float(get_parameter( 'dim2d_K', dim2dParamPath))

# compute the reduced cv_r
cv_r = get_cv_r(dim2d_M,dim2d_cv,dim2d_R)

# compute the Weber number
we = get_we(length_c, dim2d_a, dim2d_b, dim2d_M, dim2d_K)
