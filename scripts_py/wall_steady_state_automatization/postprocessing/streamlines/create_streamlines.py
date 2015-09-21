#!/usr/bin/python

'''
@description: extract the streamlines from a given
netcdf file
'''

import sys
import os
import inspect


# python path
#------------------------------------------------------------
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],"../")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)


# library imports
#------------------------------------------------------------
import getopt
import subprocess
import shlex
import shutil

from library_messages import\
    print_mg_error,\
    print_mg_progress,\
    print_mg_final

import numpy as np
from matplotlib.pyplot import cm
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
    print 'extract the streamline from a given'
    print 'netcdf file'
    print ''
    print 'options:'
    print '--------'
    print '-h (--help)    : display this help'
    print '-i (--input=)  : netcdf file *.nc'
    print ''


#parse the input options
def parse_argv(argv):
    '''
    @description:
    parse the program options
    '''

    #default values
    inputFile = 'None'


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


    inputFile = 'None'

    # options
    for opt, arg in opts:

        if opt in ("-h", "--help"):
            display_help()
            sys.exit(2)

        elif opt in ("-i", "--input"):
            inputFile = arg


    # check for directory with the netcdf files
    if(inputFile=='None'):
        print_mg_error('netcdf file not provided')
        display_help()
        sys.exit(2)

    else:
        if( not os.path.isfile(inputFile)):
            print_mg_error(inputFile+' does not exist')
            display_help()
            sys.exit(2)
    
    return inputFile,


# convert to float
def convert_to_float(x):
    '''
    @description: if missing data, put 1e10
    '''

    try:
        y = float(x)

    except ValueError:
        y = 1e10

    return y


# extract netcdf variable from a netcdf file
def extract_netcdf_var(fileName,var):
    '''
    @description: extract the variable var from 
    the netcdf file fileName
    '''

    cmd='ncdump -v '+var+' '+fileName
    args = shlex.split(cmd)

    output = subprocess.Popen(args,stdout=subprocess.PIPE).communicate()[0]

    output = output.split('data:')[1].split(var+' =')[1].replace(';','').replace('}','').replace('/n','')
    output = output.split(',')
    output = [convert_to_float(i) for i in output if i !='']

    output = np.array(output)

    return output


# extract governing variables from a netcdf file
def extract_govvar_data(fileName):
    '''
    @description: extract the x-coordinates, y-coordinates
    and (mass,momentum_x,momentum_y,energy)
    '''

    x      = extract_netcdf_var(fileName,'x')
    y      = extract_netcdf_var(fileName,'y')
    mass   = extract_netcdf_var(fileName,'mass')
    qx     = extract_netcdf_var(fileName,'momentum_x')
    qy     = extract_netcdf_var(fileName,'momentum_y')
    energy = extract_netcdf_var(fileName,'energy')

    mass       = np.reshape(mass,(len(y),len(x)))
    momentum_x = np.reshape(momentum_x,(len(y),len(x)))
    momentum_y = np.reshape(momentum_y,(len(y),len(x)))
    energy     = np.reshape(energy,(len(y),len(x)))

    return x,y,mass,momentum_x,momentum_y,energy


# extract the x-coordinates, y-coordinates and velocity vector fields
# from a netcdf file
def extract_streamline_data(fileName):
    '''
    @description: extract the x-coordinates, y-coordinates
    and velocity vector fields from a netcdf file
    '''

    x,y,mass,momentum_x,momentum_y,energy =\
        extract_govvar_data(fileName)

    vx    = qx/mass
    vy    = qy/mass

    return x,y,vx,vy


if __name__=="__main__":

    # extract the options from the command line
    inputFile, = parse_argv(sys.argv[1:])


    # extract the streamlines data from the netcdf file
    x,y,vx,vy = extract_streamline_data(inputFile)

    v = np.sqrt(vx**2 + vy**2)

    # plot the streamlines
    plt.close("all")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))

    ax = fig.add_subplot(111)
    
    ax.plot(y[2:],vx[2:,0])

    plt.show()

    #plt.quiver(x, y[2:], vx[2:,:], vy[2:,:],        # data
    #           v[2:,:],                   # colour the arrows based on this array
    #           cmap=cm.Blues)        # length of the arrows
    #
    #plt.colorbar()                  # adds the colour bar
    #
    #plt.title('Quive Plot, Dynamic Colours')
    #plt.show()
