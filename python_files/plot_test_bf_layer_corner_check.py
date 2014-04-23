#!/usr/bin/python


import commands
import getopt

import os, sys, inspect
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *


cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(inspect.getfile( inspect.currentframe() ))[0],"fortranfile-0.2.1")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from fortranfile import *
from pylab import *

from plot_interface import (manage_options,
                            get_nb_sublayers,
                            make_matrix_for_all_bf_layers,
                            plot_nodes_and_grdptid_with_all_bf_layers)
    
if __name__ == "__main__":
    

    #manage the options
    [folder_path] = manage_options()

    #find the maximum number of sublayers per main layer
    nb_sublayers = get_nb_sublayers(folder_path)


    for i in range(0,2):
        
         #test add_sublayer
         #=================================================================
         #combine data from several sublayers in one large matrix
         #-----------------------------------------------------------------
         interior_size_filename  = folder_path+'/interior_sizes.dat'
         interior_nodes_filename = folder_path+'/interior_nodes.dat'
         
         suffix_size    = '_sizes_'+str(i)+'.dat'
         suffix_nodes   = '_nodes_'+str(i)+'.dat'
         suffix_grdptid = '_grdpt_id_'+str(i)+'.dat'
         
         [lm_nodes,lm_grdptid] = make_matrix_for_all_bf_layers(interior_size_filename,
                                                               interior_nodes_filename,
                                                               folder_path,
                                                               nb_sublayers,
                                                               suffix_size,
                                                               suffix_nodes,
                                                               suffix_grdptid)
         #display
         #-----------------------------------------------------------------
         fig, ax = plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes,
                                                             lm_grdptid)
         fig.canvas.set_window_title("test bf_layer_corner_check"+str(i))

    
    #show all
    plt.show()
    
