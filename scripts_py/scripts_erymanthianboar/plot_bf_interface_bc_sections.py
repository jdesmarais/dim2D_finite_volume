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

from library_plot_bf_layer import (manage_options,
                                   get_nb_sublayers,
                                   make_matrix_for_all_bf_layers,
                                   plot_nodes_and_grdptid_with_all_bf_layers)
    
if __name__ == "__main__":
    

    #manage the options
    [folder_path] = manage_options()
    print folder_path


    test_windows_title = ["1: Allocation: N at corner N_E",
                          "2: Allocation: W at corner N_W",
                          "3: Allocation: E at corner N_E",
                          "4: Allocation: N at corner N_W",
                          "5: Allocation: E edge",
                          "6: Allocation: W edge",
                          "7: Allocation: S at corner S_E",
                          "8: Allocation: S at corner S_W",
                          "9: Merge: E edge",
                          "10: Reallocation: W edge",
                          "11: Additional layers in N and S"]
    

    nb_sublayers = 6
    

    #plot first tests
    for i in range(0,11):

        #extract data for the interior points and the buffer layers
    	#-----------------------------------------------------------------
        test_index = str(i+1)

        interior_size_filename      = folder_path+'/interior_sizes.dat'
        interior_grdptsid_filename  = folder_path+'/interior_grdpts_id.dat'
        interior_nodes_filename     = folder_path+'/interior_nodes.dat'

    	suffix_size    = '_sizes'+test_index+'.dat'
    	suffix_nodes   = '_nodes'+test_index+'.dat'
    	suffix_grdptid = '_grdpt_id'+test_index+'.dat'
    	
    	[lm_nodes,lm_grdptid, margin] = make_matrix_for_all_bf_layers(interior_size_filename,
                                                              interior_grdptsid_filename,
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
        fig.canvas.set_window_title(test_windows_title[i])


    #show all
    plt.show()
    
    
