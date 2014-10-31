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

    nodes = nodes[::,::-1,::]

    #return the interior point data
    return [sizes, x_map, y_map, nodes]


def extract_one_piece(index):
    
    sizes_filename  = folder_path+'/sizes_'+str(index)+'.dat'
    x_map_filename  = folder_path+'/x_map_'+str(index)+'.dat'
    y_map_filename  = folder_path+'/y_map_'+str(index)+'.dat'
    nodes_filename  = folder_path+'/nodes_'+str(index)+'.dat'
    
    [sizes,x_map,y_map,nodes] = extract_nodes_and_maps(
        sizes_filename,
        x_map_filename,
        y_map_filename,
        nodes_filename)

    return [sizes, x_map, y_map, nodes]


def extract_four_pieces(index):
    
    sizes_filename  = folder_path+'/sizes_ext_'+str(index)+'.dat'
    x_map_filename  = folder_path+'/x_map_ext_'+str(index)+'.dat'
    y_map_filename  = folder_path+'/y_map_ext_'+str(index)+'.dat'
    nodes_filename  = folder_path+'/nodes_ext_'+str(index)+'.dat'

    [sizes_ext, x_map_ext, y_map_ext, nodes_ext] = extract_nodes_and_maps(
        sizes_filename,
        x_map_filename,
        y_map_filename,
        nodes_filename)
    
    return [sizes_ext, x_map_ext, y_map_ext, nodes_ext]


def extract_four_pieces_interface(test_index):

    interior_size_filename      = folder_path+'/sizes_ext_int_'+str(test_index)+'.dat'
    interior_grdptsid_filename  = folder_path+'/grdptsid_ext_int_'+str(test_index)+'.dat'
    interior_nodes_filename     = folder_path+'/nodes_ext_int_'+str(test_index)+'.dat'

    suffix_size    = '_sizes'+str(test_index)+'.dat'
    suffix_nodes   = '_nodes'+str(test_index)+'.dat'
    suffix_grdptid = '_grdpt_id'+str(test_index)+'.dat'
    
    [lm_nodes,lm_grdptid, margin] = make_matrix_for_all_bf_layers(interior_size_filename,
                                                                  interior_grdptsid_filename,
                                                                  interior_nodes_filename,
                                                                  folder_path,
                                                                  nb_sublayers,
                                                                  suffix_size,
                                                                  suffix_nodes,
                                                                  suffix_grdptid,
                                                                  continuous=False)
    return [lm_nodes,lm_grdptid]


def compare_nodes_and_maps(
    x_map, y_map, nodes,
    x_map_ext, y_map_ext, nodes_ext):


    test_validated = True

    #check whether the x_map are the same
    for i in range(0,len(x_map)):
        
        test_loc = is_test_validated(
            x_map[i],
            x_map_ext[i])

        if(not(test_loc)):
            print 'x_map[', i, ']', x_map[i], x_map_ext[i]

        test_validated = test_loc and test_validated

    print 'test_x_map: ',  test_validated


    test_validated = True

    #check whether the y_map are the same
    for i in range(0,len(y_map)):
        
        test_loc = is_test_validated(
            y_map[i],
            y_map_ext[i])

        if(not(test_loc)):
            print 'y_map[', i, ']', y_map[i], y_map_ext[i]

        test_validated = test_loc and test_validated

    print 'test_y_map: ', test_validated


def is_test_validated(x1,x2):
    
    test = int(x1*1000)==int(x2*1000)

    return test

            
def plot_nodes_and_maps(nodes_0,nodes_1,var=0):

    #create the main figure
    fig=plt.figure(figsize=(12,6))

    #plot the gridpoint ID
    ax = fig.add_subplot(1,2,1)
    #res = ax.imshow(nodes_0[0,:,:], cmap=cm.spectral, interpolation='nearest')#, vmin=0.0, vmax=1.0)
    res = ax.imshow(nodes_0[var,:,:], cmap=cm.spectral, interpolation='nearest', vmin=-0.9, vmax=0.9)
    fig.colorbar(res)
    
    #plot the nodes
    ax = fig.add_subplot(1,2,2)
    res = ax.imshow(nodes_1[var,:,:], cmap=cm.spectral, interpolation='nearest', vmin=-0.9, vmax=0.9)
    fig.colorbar(res)

    return fig,ax


if __name__ == "__main__":
    

    #manage the options
    [folder_path] = manage_options()
    print folder_path

    #window titles
    test_windows_title = ["0: Initial nodes",
                          "1: Time derivatives"]

    nb_sublayers = 6
    var = 3


    #extract the nodes from the one piece domain at t=0
    [sizes_0, x_map_0, y_map_0, nodes_0] = extract_one_piece(0)

    
    #extract the time dev from the one piece domain at t=0
    [sizes_1, x_map_1, y_map_1, nodes_1] = extract_one_piece(1)

    #extract the nodes from the one piece domain after first step
    #in integration
    [sizes_2, x_map_2, y_map_2, nodes_2] = extract_one_piece(2)

    #extract the nodes from the one piece domain after
    #synchronization cycle
    [sizes_3, x_map_3, y_map_3, nodes_3] = extract_one_piece(3)

    #extract the nodes from the one piece domain after
    #integration cycle
    [sizes_4, x_map_4, y_map_4, nodes_4] = extract_one_piece(4)



    #extract the nodes from the four pieces domain at t=0
    [sizes_ext_0, x_map_ext_0, y_map_ext_0, nodes_ext_0] = extract_four_pieces(0)
    [lm_nodes_0,lm_grdptid_0] = extract_four_pieces_interface(0)

    #extract the time dev from the four pieces domain at t=0
    [sizes_ext_1, x_map_ext_1, y_map_ext_1, nodes_ext_1] = extract_four_pieces(1)
    [lm_nodes_1,lm_grdptid_1] = extract_four_pieces_interface(1)

    #extract the nodes from the four pieces domain after first step
    #in integration
    [sizes_ext_2, x_map_ext_2, y_map_ext_2, nodes_ext_2] = extract_four_pieces(2)

    #extract the nodes from the four pieces domain after
    #synchronization cycle
    [sizes_ext_3, x_map_ext_3, y_map_ext_3, nodes_ext_3] = extract_four_pieces(3)
    [lm_nodes_3,lm_grdptid_3] = extract_four_pieces_interface(3)

    #extract the nodes from the four pieces domain after
    #integration cycle
    [sizes_ext_4, x_map_ext_4, y_map_ext_4, nodes_ext_4] = extract_four_pieces(4)


    #print the x_map and y_map
    print ''
    for i in range(0,len(x_map_0)):
        print 'x_map(',i,')', x_map_0[i], x_map_ext_0[i]
    print ''

    print ''
    for i in range(0,len(y_map_0)):
        print 'y_map(',i,')', y_map_0[i], y_map_ext_0[i]
    print ''


    #compare the data
    compare_nodes_and_maps(
        x_map_0, y_map_0, nodes_0,
        x_map_ext_0, y_map_ext_0, nodes_ext_0)

 
    #plot the data: time derivatives
    fig, ax = plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes_0,lm_grdptid_0,var=var,vmin=-0.9,vmax=0.9)
    fig.canvas.set_window_title('interface: initialization')

    fig, ax = plot_nodes_and_maps(nodes_0, nodes_1, var=var)
    fig.canvas.set_window_title('one piece domain')

    fig, ax = plot_nodes_and_maps(nodes_ext_0, nodes_ext_1, var=var)
    fig.canvas.set_window_title('four pieces domain')

    fig, ax = plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes_1,lm_grdptid_1,var=var,vmin=-0.9,vmax=0.9)
    fig.canvas.set_window_title('interface: time derivatives')

    
    ##plot the data: after one step
    #fig, ax = plot_nodes_and_maps(nodes_0, nodes_2, var=var)
    #fig.canvas.set_window_title('one piece domain: after one step')
    #
    #fig, ax = plot_nodes_and_maps(nodes_ext_0, nodes_ext_2, var=var)
    #fig.canvas.set_window_title('four pieces domain: after one step')
    #
    #
    ##plot the data: after nodes synchronization
    #fig, ax = plot_nodes_and_maps(nodes_0, nodes_3, var=var)
    #fig.canvas.set_window_title('one piece domain: after nodes synchronization')
    #
    #fig, ax = plot_nodes_and_maps(nodes_ext_0, nodes_ext_3, var=var)
    #fig.canvas.set_window_title('four pieces domain: after nodes synchronization')
    #
    #fig, ax = plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes_3,lm_grdptid_3, var=var, vmin=-0.9, vmax=0.9)
    #fig.canvas.set_window_title('interface: after nodes synchronization')
    #
    #
    ##plot the data: after integration cycle
    #fig, ax = plot_nodes_and_maps(nodes_0, nodes_4, var=var)
    #fig.canvas.set_window_title('one piece domain: after integration cycle')
    #
    #fig, ax = plot_nodes_and_maps(nodes_ext_0, nodes_ext_4, var=var)
    #fig.canvas.set_window_title('four pieces domain: after integration cycle')


    #show all
    plt.show()

    
