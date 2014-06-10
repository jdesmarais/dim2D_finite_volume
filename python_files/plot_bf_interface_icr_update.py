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
                                   plot_nodes_and_grdptid_with_all_bf_layers,
                                   extract_detectors_data,
                                   add_detectors,
                                   plot_all_bf_layers)
    
if __name__ == "__main__":
    

    #manage the options
    [folder_path] = manage_options()
    print folder_path


    nb_sublayers = 6


    #plot the tests
    for i in range(0,10):
        
    	#extract data for the interior points and the buffer layers
    	#-----------------------------------------------------------------
        test_index = str(i)

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
    	
    	[lm_nodes,lm_grdptid,margin] = make_matrix_for_all_bf_layers(interior_size_filename,
                                                                     interior_grdptsid_filename,
                                                                     interior_nodes_filename,
                                                                     folder_path,
                                                                     nb_sublayers,
                                                                     suffix_size,
                                                                     suffix_nodes,
                                                                     suffix_grdptid,
                                                                     continuous=True)

        #extract detector informations
        N_detector_filename = folder_path+'/N_detectors'+test_index+'.dat'
        S_detector_filename = folder_path+'/S_detectors'+test_index+'.dat'
        E_detector_filename = folder_path+'/E_detectors'+test_index+'.dat'
        W_detector_filename = folder_path+'/W_detectors'+test_index+'.dat'

        [N_detector, S_detector, E_detector, W_detector] = extract_detectors_data(N_detector_filename,
                                                                                  S_detector_filename,
                                                                                  E_detector_filename,
                                                                                  W_detector_filename)

        [lm_grdptid] = add_detectors(N_detector,
                                     S_detector,
                                     E_detector,
                                     W_detector,
                                     lm_grdptid,
                                     margin)    
        #display
        #-----------------------------------------------------------------
        fig, ax = plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes,
                                                            lm_grdptid)
        
        #fig, ax = plot_all_bf_layers(lm_nodes,
        #                             lm_grdptid)
        
        fig.canvas.set_window_title('time step: '+str(i))

    
    #show all
    plt.show()
    
    
