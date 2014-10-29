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


    test_windows_title = ["1: Allocation: N_E and W at corner",
                          "2: Allocation: E at N_E corner and N_W",
                          "3: Allocation: E and W layers",
                          "4: Allocation: S_W and S_E corners",
                          "5: Reallocation W and merge E",
                          "6: Before neighbor dependence tests"]
    

    nb_sublayers = 6
    

    #plot first tests
    for i in range(0,6):

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


    #plot the test_get_nbf_layers_sharing_grdpts_with
    fig = plt.figure(figsize=(16,9))
    
    for i in range(6,19):
        
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
        ax = fig.add_subplot(4,4,i-5)
        res = ax.imshow(lm_nodes[0,:,:],
                        cmap=cm.spectral,
                        interpolation='nearest',
                        vmin=0.0,
                        vmax=1.0)
        
    fig.canvas.set_window_title('7: test_get_nbf_layers_sharing_grdpts_with')


    #plot the test of removing a buffer layer
    for i in range(19,20):

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
        fig.canvas.set_window_title('8: test_remove_sublayer')


    #plot the test of determining the boundary layers
    for i in range(20,21):

        #extract data for the interior points and the buffer layers
    	#-----------------------------------------------------------------
        test_index = str(i+1)

        interior_size_filename      = folder_path+'/interior_sizes_bc_sec.dat'
        interior_grdptsid_filename  = folder_path+'/interior_grdpts_id_bc_sec.dat'
        interior_nodes_filename     = folder_path+'/interior_nodes_bc_sec.dat'

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
        fig.canvas.set_window_title('9: test_determine_interior_bc_layers')


    #plot the test of determining the boundary procedures
    for i in range(21,22):

        #extract data for the interior points and the buffer layers
    	#-----------------------------------------------------------------
        test_index = str(i+1)

        interior_size_filename      = folder_path+'/interior_sizes_bc_proc.dat'
        interior_grdptsid_filename  = folder_path+'/interior_grdpts_id_bc_proc.dat'
        interior_nodes_filename     = folder_path+'/interior_nodes_bc_proc.dat'

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
        fig.canvas.set_window_title('10: test_determine_interior_bc_procedures')
    

    #plot the test of synchronizing the nodes of the buffer layers
    #with the interior domain
    for i in range(22,23):

        #extract data for the interior points and the buffer layers
    	#-----------------------------------------------------------------
        test_index = str(i+1)

        interior_size_filename      = folder_path+'/interior_sizes_sync_int.dat'
        interior_grdptsid_filename  = folder_path+'/interior_grdpts_id_sync_int.dat'
        interior_nodes_filename     = folder_path+'/interior_nodes_sync_int.dat'

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
        fig.canvas.set_window_title('11: test_sync_nodes_with_interior')
    

    #plot the test of synchronizing the nodes at the interface
    #between buffer main layers
    for i in range(23,24):

        #extract data for the interior points and the buffer layers
    	#-----------------------------------------------------------------
        test_index = str(i+1)

        interior_size_filename      = folder_path+'/interior_sizes_sync_bml.dat'
        interior_grdptsid_filename  = folder_path+'/interior_grdpts_id_sync_bml.dat'
        interior_nodes_filename     = folder_path+'/interior_nodes_sync_bml.dat'

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
        fig.canvas.set_window_title('12: test_sync_nodes_at_mainlayer_interfaces')


    #show all
    plt.show()
    
    
