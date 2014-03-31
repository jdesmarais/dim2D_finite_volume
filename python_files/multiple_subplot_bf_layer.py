#!/usr/bin/python

import commands
import getopt

import os, sys, inspect
import numpy as np
import matplotlib.figure as Figure
import matplotlib.pyplot as plt
import matplotlib.cm as cm


cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(inspect.getfile( inspect.currentframe() ))[0],"fortranfile-0.2.1")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from fortranfile import *
from pylab import *


def manage_options():

    folder_path = ''
    
    opts, extraparams = getopt.getopt(sys.argv[1:],
                                      "f:",
                                      ["folder="])
    for o,p in opts:
        
        if o in ['-f','--folder']:
            folder_path = p
            if (not os.path.isdir(p)):
                print 'directory path does not exist'
                sys.exit(0)

    return [folder_path]


def plot_grdpt_id(folder_path,bf_layer_location):

    #load buffer layer sizes file
    sizes_filename = folder_path+'/'+bf_layer_location+'_sizes.dat'
    f = FortranFile(sizes_filename)
    sizes = f.readInts()
    f.close()

    #load buffer layer nodes file
    nodes_filename = folder_path+'/'+bf_layer_location+'_nodes.dat'
    f = FortranFile(nodes_filename)
    nodes    = f.readReals('d')
    f.close()

    #load buffer layer gridpt id file
    grdpt_id_filename = folder_path+'/'+bf_layer_location+'_grdpt_id.dat'
    f = FortranFile(grdpt_id_filename)
    grdpt_id = f.readInts()
    f.close()

    #plot the gridpt id data
    size_x  = sizes[1]
    size_y  = sizes[0]
    size_ne = sizes[2]
    
    grdpt_id.resize(size_x,size_y)
    grdpt_id = grdpt_id[::-1,::]


    options = {'N'  : [6,3],
               'S'  : [6,3],
               'E'  : [3,6],
               'W'  : [3,6],
               'NE' : [3,3],
               'NW' : [3,3],
               'SE' : [3,3],
               'SW' : [3,3],
               }    

    import gtk
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_gtkagg import FigureCanvasGTKAgg as FigureCanvas

    #win = gtk.Window()
    #win.connect("destroy", lambda x: gtk.main_quit())
    #win.set_default_size(400,300)
    #win.set_title("Some Window")

    fig=plt.figure(figsize=(options[bf_layer_location][0],options[bf_layer_location][1])) #figsize=(12,6))
    ax = fig.add_subplot(1,1,1)
    res = ax.imshow(grdpt_id, cmap=cm.summer, interpolation='nearest')

    xticks = np.arange(1,size_y+1)
    yticks = np.arange(size_x,0,-1)
    
    plt.xticks(np.arange(size_y)  , xticks[:size_y])
    plt.yticks(np.arange(size_x+1), yticks[:size_x+1])

    cb = fig.colorbar(res)

    fig.suptitle(bf_layer_location+' : gridpoints id')

    #thismanager = get_current_fig_manager()
    #thismanager.window.SetPosition((500, 0))
    #thismanager.window.SetPosition((500, 0))
#    plt.title(bf_layer_location+' : gridpoints id')
#    plt.show()
    #fig.canvas.manager.window.Move(100,400)
    #thismanager.window.setGeometry(50,100,640, 545)


if __name__ == "__main__":

    #manage the options
    [folder_path] = manage_options()


    #define the buffer layer tested
    bf_layer_loc_table = ['N','S','E','W','NE','NW','SE','SW']
    

    #loop over the buffer layer tested
    #and visualize the gridpoints id
    for bf_layer_loc in bf_layer_loc_table:

        plot_grdpt_id(folder_path,bf_layer_loc)


    plt.show()
