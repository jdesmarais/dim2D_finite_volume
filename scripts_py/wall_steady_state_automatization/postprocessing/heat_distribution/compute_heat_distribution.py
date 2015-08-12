#!/usr/bin/python


'''
@description: compute the heat distribution, i.e.
the amount of heat which is transported across the
domain borders, the heat used for the phase transition
and the heat used to increase the internal energy of
the fluid
'''


# library imports
#------------------------------------------------------------
import sys
import os
import getopt

from compute_energyTr import compute_energyTr
from compute_energyPh import compute_energyPh
from compute_energyRa import compute_energyRa

from library_messages import print_mg_error
from library_colors import grayscale_to_RGB

import numpy as np
import matplotlib.pyplot as plt


# functions
#------------------------------------------------------------

# display the help for the program
def display_help():
    '''
    @description:
    display program help
    '''
    print ''
    print 'compute the heat distribution from the netcdf files'
    print ''
    print 'options:'
    print '--------'
    print '-h (--help)    : display this help'
    print '-i (--input=)  : main dir with the netcdf files data*.nc'
    print ''


#parse the input options
def parse_argv(argv):
    '''
    @description:
    parse the program options
    '''

    #default values
    mainDir = 'None'


    # store the options and the arguments
    # in opts, args
    try:
        opts, args = getopt.getopt(argv,
                                   "hi:",
                                   ["help",
                                    "input="])
    except getopt.GetoptError:
        display_help()
        sys.exit(2)

    if(len(opts)==0):
        display_help()
        sys.exit(2)


    mainDir = 'None'
    outputFile = 'None'

    # options
    for opt, arg in opts:

        if opt in ("-h", "--help"):
            display_help()
            sys.exit(2)

        elif opt in ("-i", "--input"):
            mainDir = arg


    # check for directory with the netcdf files
    if(mainDir=='None'):
        print_mg_error('directory for netcdf file not provided')
        display_help()
        sys.exit(2)

    else:
        if( not os.path.isdir(mainDir)):
            print_mg_error(mainDir+' does not exist')
            display_help()
            sys.exit(2)


    return mainDir


if __name__=='__main__':


    # extract the options from the command line
    mainDir = parse_argv(sys.argv[1:])


    # determine the paths for the files containing
    # the energy transported across the domain borders
    # and the energy for the phase transition
    energyTr_filename = os.path.join(mainDir,'contours','energyTr.txt')
    energyPh_filename = os.path.join(mainDir,'contours','energyPh.txt')
    energyRa_filename = os.path.join(mainDir,'contours','energyRa.txt')

    # check whether the energy transported across the
    # domain borders has already been computed, otherwise
    # compute it
    if(not os.path.isfile(energyTr_filename)):
        
        compute_energyTr(mainDir,
                         energyTr_filename)


    # check whether the energy used by the phase transition
    # has already been computed, otherwise compute it
    if(not os.path.isfile(energyPh_filename)):
        
        temperature_filename = os.path.join(mainDir,'contours','temperature.txt')
        mass_filename        = os.path.join(mainDir,'contours','mass.txt')

        compute_energyPh(temperature_filename,
                         mass_filename,
                         energyPh_filename)


    # check whether the energy distribution has already
    # been computed, otherwise compute it
    if(not os.path.isfile(energyRa_filename)):

        ncFiles = [name for name in os.listdir(mainDir)\
                   if (os.path.isfile(os.path.join(mainDir, name)) and\
                           name.endswith('.nc') and \
                           name.startswith('data') )]

        compute_energyRa(os.path.join(mainDir,ncFiles[0]),
                         energyTr_filename,
                         energyPh_filename,
                         energyRa_filename)


    # plot the energy distribution as a function of time
    # energyRa[:,0] : file ID
    # energyRa[:,1] : time
    # energyRa[:,2] : energy flowing across the domain borders
    # energyRa[:,3] : energy used for phase transition
    # energyRa[:,4] : energy for the internal energy
        
    energyRa = np.loadtxt(energyRa_filename)

    plt.close("all")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))

    ax = fig.add_subplot(111)

    style = ['s-','-','o-']
    width = 3

    for i in range(0,len(style)):

        grayscale_value = 0.1 + 0.9*(float(i)/float(len(style)-1))

        plt.plot(
            energyRa[:,1],
            energyRa[:,2+i],
            style[i],
            linewidth=width,
            color=grayscale_to_RGB(grayscale_value))

    plt.show()
