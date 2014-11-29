#!/usr/bin/python


import commands
import getopt

import os, sys, inspect
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *


#cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
#    os.path.split(
#    inspect.getfile( inspect.currentframe() ))[0],"fortranfile-0.2.1")))
#if cmd_subfolder not in sys.path:
#    sys.path.insert(0, cmd_subfolder)
#
#from fortranfile import *

#manage options
def manage_options():

    folder_path = '.'
    size_x = 0
    size_y = 0
    
    opts, extraparams = getopt.getopt(sys.argv[1:],
                                      "i:x:y:",
                                      ["file=","size_x=","size_y="])
    for o,p in opts:
        
        if o in ['-i','--file']:
            folder_path = p
            if (not os.path.isfile(p)):
                print 'file does not exist'
                sys.exit(0)

        elif o in ['-x','--size_x']:
            size_x = int(p)

        elif o in ['-y','--size_y']:
            size_y = int(p)

    return [folder_path,size_x,size_y]


#extract data
def extract_grdpts_id_data(filename,sizes):
    
    grdpts_id = np.empty([sizes[1], sizes[0]])
    grdpts_id.fill(0)

    f = open(filename)

    for i in range(0,20):
        f.readline()

    for j in range(0,sizes[1]):
        for i in range(0,sizes[0]):

            line = f.readline()
            line = line[8:13]

            grdpts_id[j,i] = int(line)
    
#    f = FortranFile(grdptid_filename)
#    grdpts_id = f.readInts()
    f.close()

    grdpts_id.resize(sizes[1],sizes[0])

    grdpts_id = grdpts_id[::-1,::]

    return [grdpts_id]


#plot data
def plot_grdpts_id_data(grdpts_id):

    #create the main figure
    fig=plt.figure(figsize=(3,6))

    #plot the gridpoint ID
    ax = fig.add_subplot(1,1,1)
    res = ax.imshow(grdpts_id[:,:], cmap=cm.spectral, interpolation='nearest', vmin=-1, vmax=4)
    fig.colorbar(res)

    return fig,ax


if __name__ == "__main__":

    #get filename and grdpts_id extensions
    [filename,size_x,size_y] = manage_options()

    #get grdpts_id
    [grdpts_id] = extract_grdpts_id_data(filename,[size_x,size_y])

    print grdpts_id

    ##plot grdpts_id
    fig,ax = plot_grdpts_id_data(grdpts_id)
    #
    ##show figure
    plt.show()
    
    
