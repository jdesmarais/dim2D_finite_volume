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


    test_windows_title = ["Initial detector position",
                          "Bc_interior_pt modified by interior",
                          "Bc_interior_pt modified by bf_sublayer",
                          "Bc_interior_pt modified by bf_sublayer2",
                          "Test bf_detector_list bf recombination",
                          "Test bf_detector_list recombination"]

    nb_sublayers = 6


    #plot the tests
    for i in range(0,6):
        
    	#extract data for the interior points and the buffer layers
    	#-----------------------------------------------------------------
        test_index = str(i+1)

        #if(i==3):
        #    interior_size_filename      = folder_path+'/interior_sizes'+test_index+'.dat'
        #    interior_grdptsid_filename  = folder_path+'/interior_grdpts_id'+test_index+'.dat'
        #    interior_nodes_filename     = folder_path+'/interior_nodes'+test_index+'.dat'
        #else:
        interior_size_filename      = folder_path+'/interior_sizes'+test_index+'.dat'
        interior_grdptsid_filename  = folder_path+'/interior_grdpts_id'+test_index+'.dat'
        interior_nodes_filename     = folder_path+'/interior_nodes'+test_index+'.dat'

    	suffix_size    = '_sizes'+test_index+'.dat'
    	suffix_nodes   = '_nodes'+test_index+'.dat'
    	suffix_grdptid = '_grdpt_id'+test_index+'.dat'
    	
    	[lm_nodes,lm_grdptid] = make_matrix_for_all_bf_layers(interior_size_filename,
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


    ##plot the nbc_templates
    ##------------------------
    ##create the main figure
    #fig=plt.figure(figsize=(12,6))
    #for j in range(1,8):
    #    for i in range(1,8):
    #
    #        if(i!=4 and j!=4):
    #            nbc_filename = folder_path+'/nbc_template'+str(i)+str(j)+'.dat'
    #        
    #            f = FortranFile(nbc_filename)
    #            nbc_template = f.readInts()
    #            f.close()
    #            nbc_template.resize(3,3)
    #
    #            nbc_template = nbc_template[::-1,::]
    #
    #            ax = fig.add_subplot(7,7,i+(7-j)*7)
    #            res = ax.imshow(nbc_template, cmap=cm.spectral, interpolation='nearest', vmin=-1, vmax=4)
    #            ax.get_xaxis().set_visible(False)
    #            ax.get_yaxis().set_visible(False)
    #        #fig.colorbar(res)
    #
    #fig.canvas.set_window_title('test create_template')

    
    #show all
    plt.show()
    
    
