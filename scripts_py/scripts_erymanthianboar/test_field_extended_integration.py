#!/usr/bin/python


import commands
import getopt

import os, sys, inspect
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from pylab import *


cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(inspect.getfile( inspect.currentframe() ))[0],
    "fortranfile-0.2.1")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from fortranfile import *
from pylab import *

from library_plot_bf_layer import (manage_options,
                                   get_nb_sublayers,
                                   make_matrix_for_all_bf_layers,
                                   plot_nodes_and_grdptid_with_all_bf_layers)
    

def extract_nodes_and_maps(
    sizes_filename,
    x_map_filename,
    y_map_filename,
    nodes_filename):

    #load the interior point sizes file
    #sizes_filename = folder_path+'/interior_sizes.dat'
    f = FortranFile(sizes_filename)
    sizes = f.readInts()
    f.close()

    #load the interior point x_map file
    #x_map_filename = folder_path+'/interior_x_map.dat'
    f = FortranFile(x_map_filename)
    x_map = f.readReals('d')
    f.close()

    #load the interior point y_map file
    #y_map_filename = folder_path+'/interior_y_map.dat'
    f = FortranFile(y_map_filename)
    y_map = f.readReals('d')
    f.close()

    #load the interior point nodes file
    #nodes_filename = folder_path+'/interior_nodes.dat'
    f = FortranFile(nodes_filename)
    nodes = f.readReals('d')
    f.close()

    #reorganize the nodes data
    size_x  = sizes[0]
    size_y  = sizes[1]
    size_ne = sizes[2]

    nodes.resize(size_ne,size_y,size_x)

    #return the interior point data
    return [sizes, x_map, y_map, nodes]


def compare_nodes_and_maps(
    x_map, y_map, nodes,
    x_map_ext, y_map_ext, nodes_ext):

    #check whether the x_map are the same
    for i in range(0,len(x_map)):
        
        [test_validated] = is_test_validated(
            x_map[i],
            x_map_ext[i])

        if(not(test_validated)):
            print 'x_map[', i, ']', x_map[i], x_map_ext[i]

            



if __name__ == "__main__":
    

    #manage the options
    [folder_path] = manage_options()
    print folder_path

    #window titles
    test_windows_title = ["0: Initial nodes",
                          "1: After one integration"]

    nb_sublayers = 6
    

    #plot the field extended: four pieces
    for i in range(0,2):

        #extract data for the interior points and the buffer layers
    	#-----------------------------------------------------------------
        test_index = str(i)

        interior_size_filename      = folder_path+'/sizes_ext_int_'+str(i)+'.dat'
        interior_grdptsid_filename  = folder_path+'/grdpts_id_ext_int_'+str(i)+'.dat'
        interior_nodes_filename     = folder_path+'/nodes_ext_int_'+str(i)+'.dat'

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
                                                                      suffix_grdptid,
                                                                      continuous=False)
    
        #display
        #-----------------------------------------------------------------
        fig, ax = plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes,
                                                            lm_grdptid)
        fig.canvas.set_window_title(test_windows_title[i])

        
    #plot the field extended: one piece
    for i in range(0,2):

        #extract data for the interior points and the buffer layers
    	#-----------------------------------------------------------------
        test_index = str(i)

        interior_size_filename      = folder_path+'/sizes_'+str(i)+'.dat'
        interior_grdptsid_filename  = folder_path+'/grdpts_id_ext_int_'+str(i)+'.dat'
        interior_nodes_filename     = folder_path+'/nodes_'+str(i)+'.dat'

    	suffix_size    = '_no_sizes'+test_index+'.dat'
    	suffix_nodes   = '_no_nodes'+test_index+'.dat'
    	suffix_grdptid = '_no_grdpt_id'+test_index+'.dat'
    	
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


    #extract the nodes from the one piece domain
    sizes_filename  = folder_path+'/sizes_0.dat'
    x_map_filename  = folder_path+'/x_map_0.dat'
    y_map_filename  = folder_path+'/y_map_0.dat'
    nodes_filename  = folder_path+'/nodes_0.dat'

    [sizes_0,
     x_map_0,
     y_map_0,
     nodes_0] = extract_nodes_and_maps(
        sizes_filename,
        x_map_filename,
        y_map_filename,
        nodes_filename)

    sizes_filename  = folder_path+'/sizes_1.dat'
    x_map_filename  = folder_path+'/x_map_1.dat'
    y_map_filename  = folder_path+'/y_map_1.dat'
    nodes_filename  = folder_path+'/nodes_1.dat'

    [sizes_1,
     x_map_1,
     y_map_1,
     nodes_1] = extract_nodes_and_maps(
        sizes_filename,
        x_map_filename,
        y_map_filename,
        nodes_filename)

    
    #extract the nodes from the four pieces domain
    sizes_filename  = folder_path+'/sizes_ext_0.dat'
    x_map_filename  = folder_path+'/x_map_ext_0.dat'
    y_map_filename  = folder_path+'/y_map_ext_0.dat'
    nodes_filename  = folder_path+'/nodes_ext_0.dat'

    [sizes_ext_0,
     x_map_ext_0,
     y_map_ext_0,
     nodes_ext_0] = extract_nodes_and_maps(
        sizes_filename,
        x_map_filename,
        y_map_filename,
        nodes_filename)

    sizes_filename  = folder_path+'/sizes_ext_1.dat'
    x_map_filename  = folder_path+'/x_map_ext_1.dat'
    y_map_filename  = folder_path+'/y_map_ext_1.dat'
    nodes_filename  = folder_path+'/nodes_ext_1.dat'

    [sizes_ext_1,
     x_map_ext_1,
     y_map_ext_1,
     nodes_ext_1] = extract_nodes_and_maps(
        sizes_filename,
        x_map_filename,
        y_map_filename,
        nodes_filename)


    #compare the data
    compare_nodes_and_maps(
        x_map_0,
        y_map_0,
        nodes_0,
        x_map_ext_0,
        y_map_ext_0,
        nodes_ext_0)

    compare_nodes_and_maps(
        x_map_1,
        y_map_1,
        nodes_1,
        x_map_ext_1,
        y_map_ext_1,
        nodes_ext_1)
