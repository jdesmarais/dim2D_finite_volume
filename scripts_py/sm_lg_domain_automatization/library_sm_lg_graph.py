#!/usr/bin/python

'''
@description
useful functions to generate graphs for the comparison
of the small and large domains simulations
'''

import sys
import os
import subprocess
from library_sm_lg_error import generate_fortran_exe

import numpy as np
import matplotlib.pyplot as plt


cmd_subfolder = os.path.realpath(os.path.abspath(
        os.path.join(os.getenv('augeanstables'),
                     'scripts_py',
                     'scripts_erymanthianboar',
                     'fortranfile-0.2.1')
        ))
#
#
#            os.path.split(inspect.getfile( inspect.currentframe() ))[0],
#            "fortranfile-0.2.1")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from fortranfile import *


augeanstables=os.getenv('augeanstables')
extract_1D_var=os.path.join(augeanstables,
                            'scripts_fortran',
                            'error_computation',
                            'extract_1D_variable')



# generate fortran executable to extract 1D variables
# from netcdf file
def generate_extract_1D_real_exe():
    '''
    @description
    generate fortran executable to extract 1D variables
    from netcdf file
    '''
    
    exePath = extract_1D_var
    
    print 'generating the fortran executable to extract 1D variables: ...'

    generate_fortran_exe(exePath)

    print 'generating the fortran executable to extract 1D variables: done'
    print ''

    return exePath    


# extract the 1D variable from netcdf file
def extract_var(exePath, ncFile, varName):
    '''
    @description:
    extract the 1D variable from netcdf file
    '''

    #1) verify that the exePath exists
    if(not os.path.isfile(exePath)):
        print 'library_sm_lg_graph'
        print 'extract_var'
        sys.exit('***exePath '+exePath+' not found***')


    #2) verify that the ncFile exists
    if(not os.path.isfile(ncFile)):
        print 'library_sm_lg_graph'
        print 'extract_var'
        sys.exit('***ncFile '+ncFile+' not found***')


    #3) generate a binary output file from
    #   the netcdf file containing the variable
    outFile = os.path.join(os.path.dirname(ncFile),'var.out')
    if(os.path.isfile(outFile)):
        os.remove(outFile)

    cmd = exePath
    cmd+= ' -i '+ncFile
    cmd+= ' -o '+outFile
    cmd+= ' -v '+varName
    subprocess.call(cmd, shell=True)


    #4) verify that the binary output file
    #   was generated
    if(not os.path.isfile(outFile)):
        print 'library_sm_lg_error'
        print 'extract_var'
        sys.exit('***output file '+outFile+' not generated***')


    #5) extract the content of the output
    #   file using FortranFile fcts
    f = FortranFile(outFile)
    var = f.readReals('d')
    f.close()


    #6) remove the binary output file
    if(os.path.isfile(outFile)):
        os.remove(outFile)

    return var


# create a graph for the maximum error in time
def generate_graph_max_error_in_time(errorPath,
                                     var_name='max_error_mass',
                                     width=3,
                                     figPath=''):
    '''
    @description:
    create a graph for the maximum error in time
    '''
    
    #1) verify that the error_file exists
    if(not os.path.isfile(errorPath)):
        print 'library_sm_lg_graph'
        print 'generate_graph_max_error_in_time'
        sys.exit('***errorFile '+errorPath+' was not found***')


    #2) generate the executable to extract the 
    #   variables from the file
    exePath = generate_extract_1D_real_exe()


    #3) extract the time variable
    time = extract_var(exePath, errorPath, 'time')


    #4) extract the maximum of the error over
    #   space ('max_error_mass' by default)
    var = extract_var(exePath, errorPath, var_name)


    #5) rescale the time variable to be b/w [0,1]
    time_rescaled = time/time[len(time)-1]


    #6) plot the time as x-coordinate
    #   and max_error as y-coordinate
    #plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.figure(figsize=(8,4))
    plt.subplot(111)

    plt.plot(time_rescaled[:],
             var[:],
             '-',
             linewidth=width,
             color='black')
    #
    plt.xlabel(r"$ t $")
    plt.ylabel(r"Maximum error")
    #plt.xlabel(r"$\displaystyle{t}$")
    #plt.ylabel(r"$\displaystyle Maximum error$")


    #8) show plot
    plt.show()

    #9) save figure
    if(not figPath==''):
        plt.savefig(figPath)

if __name__=="__main__":
    
    errorDir = os.path.join(os.getenv('HOME'),
                             'projects',
                             'dim2d_0.99_0.1',
                             'error')

    errorPath = os.path.join(errorDir,
                             'error_max.nc')

    figPath = os.path.join(errorDir,
                           'error_max_fig.eps')

    generate_graph_max_error_in_time(errorPath,
                                     figPath=figPath)

    
