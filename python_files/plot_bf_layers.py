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


backgrd_pt = -10000
none_pt = 0
interior_pt = 1
bc_pt = 2

bc_size = 2

nodes_type = 1
grdptid_type = 2


def manage_options():

    folder_path = './'
    
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


def extract_interior_data(folder_path):

    #load the interior point sizes file
    sizes_filename = folder_path+'/interior_sizes.dat'
    f = FortranFile(sizes_filename)
    sizes = f.readInts()
    f.close()

    #load the interior point nodes file
    nodes_filename = folder_path+'/interior_nodes.dat'
    f = FortranFile(nodes_filename)
    nodes = f.readReals('d')
    f.close()

    #reorganize the nodes data
    size_x = sizes[0]
    size_y = sizes[1]
    size_ne = sizes[2]

    nodes.resize(size_ne,size_y,size_x)
    nodes = nodes[::,::,::]

    #return the interior point data
    return [sizes, nodes]


def extract_bf_data(folder_path,suffix_files,bf_layer_location):


    #load buffer layer sizes file
    sizes_filename = folder_path+'/'+bf_layer_location+suffix_files['sizes']
    f = FortranFile(sizes_filename)
    sizes = f.readInts()
    f.close()

    #load buffer layer nodes file
    nodes_filename = folder_path+'/'+bf_layer_location+suffix_files['nodes']
    f = FortranFile(nodes_filename)
    nodes    = f.readReals('d')
    f.close()

    #load buffer layer gridpt id file
    grdpt_id_filename = folder_path+'/'+bf_layer_location+suffix_files['grdptid']
    f = FortranFile(grdpt_id_filename)
    grdpt_id = f.readInts()
    f.close()

    #reorganize the gridpt id data
    size_x  = sizes[0]
    size_y  = sizes[1]
    size_ne = sizes[2]

    nodes.resize(size_ne,size_y,size_x)
    nodes = nodes[::,::,::]
    
    grdpt_id.resize(size_y,size_x)
    grdpt_id = grdpt_id[::,::]


    #return the buffer layer data
    return [sizes, nodes, grdpt_id]


def create_matrix_with_all_bf_layers(data,data_type):

    plot_buffer = {'N_'  : True,
                   'S_'  : True,
                   'E_'  : True,
                   'W_'  : True,
                   'NE' : True,
                   'NW' : True,
                   'SE' : True,
                   'SW' : True,}

    #the data are organized as follows:
    #each dictionnary entry is accessed by the id
    #of the buffer layer {N,S,E,W,NE,NW,SE,SW}
    #then each entry contains the size of the buffer
    #layer, the nodes and the gridpoint ids

    #the data type refers either to the nodes or
    #the gridpoints id

    #1) determination of size of the large matrix
    #containing all the buffer layers
    #determination of the largest size of the buffer layers
    #in the x-direction: largest number of columns
    bf_tmp_size_x = max(data['NE'][0][0],
                        data['NW'][0][0],
                        data['E_'][0][0],
                        data['W_'][0][0],
                        data['SE'][0][0],
                        data['SW'][0][0])

    bf_tmp_size_y = max(data['NE'][0][1],
                        data['NW'][0][1],
                        data['N_'][0][1],
                        data['S_'][0][1],
                        data['SE'][0][1],
                        data['SW'][0][1])
    interspace = 2

    nodes_size_x = data['interior'][0][1]
    nodes_size_y = data['interior'][0][0]
    
    lm_size_x = 2*bf_tmp_size_x + nodes_size_x + 4*interspace
    lm_size_y = 2*bf_tmp_size_y + nodes_size_y + 4*interspace
    

    #2) allocation of the large matrix and initialization
    lm = np.empty([lm_size_y, lm_size_x])
    lm.fill(backgrd_pt)

    
    #3) filling the large matrix with the buffer layers data
    #3.1) SW buffer layer
    if(plot_buffer['SW']):
    	bf_layer_size_x = data['SW'][0][0]
    	bf_layer_size_y = data['SW'][0][1]
    	i_match = interspace + bf_tmp_size_x - bf_layer_size_x
    	j_match = interspace + bf_tmp_size_y - bf_layer_size_y
    	if(data_type==nodes_type):
    	    lm[j_match:j_match+bf_layer_size_y,
    	       i_match:i_match+bf_layer_size_x] = data['SW'][data_type][0,:,:]
    	if(data_type==grdptid_type):
    	    lm[j_match:j_match+bf_layer_size_y,
    	       i_match:i_match+bf_layer_size_x] = data['SW'][data_type][:,:]

    #3.2) S  buffer layer
    if(plot_buffer['S_']):
    	bf_layer_size_x = data['S_'][0][0]
    	bf_layer_size_y = data['S_'][0][1]
    	i_match = interspace + bf_tmp_size_x + interspace + data['S_'][0][3] - 2*bc_size +1
    	j_match = interspace + bf_tmp_size_y - bf_layer_size_y
    	if(data_type==nodes_type):
    	    lm[j_match:j_match+bf_layer_size_y,
    	       i_match:i_match+bf_layer_size_x] = data['S_'][data_type][0,:,:]
    	if(data_type==grdptid_type):
    	    lm[j_match:j_match+bf_layer_size_y,
    	       i_match:i_match+bf_layer_size_x] = data['S_'][data_type][:,:]
    
    #3.3) SE buffer layer
    if(plot_buffer['SE']):
        bf_layer_size_x = data['SE'][0][0]
        bf_layer_size_y = data['SE'][0][1]
        i_match = interspace + bf_tmp_size_x + interspace + nodes_size_x + interspace
        j_match = interspace + bf_tmp_size_y - bf_layer_size_y
        if(data_type==nodes_type):
            lm[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = data['SE'][data_type][0,:,:]
        if(data_type==grdptid_type):
            lm[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = data['SE'][data_type][:,:]
        
    #3.4) E  buffer layer
    if(plot_buffer['E_']):
        bf_layer_size_x = data['E_'][0][0]
        bf_layer_size_y = data['E_'][0][1]
        i_match = interspace + bf_tmp_size_x + interspace + nodes_size_x + interspace
        j_match = interspace + bf_tmp_size_y + interspace + data['E_'][0][5] - 2*bc_size +1
        if(data_type==nodes_type):
            lm[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = data['E_'][data_type][0,:,:]
        if(data_type==grdptid_type):
            lm[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = data['E_'][data_type][:,:]
    
    #3.5) nodes
    if(data_type==nodes_type):
        i_match = interspace + bf_tmp_size_x + interspace
	j_match = interspace + bf_tmp_size_y + interspace
        lm[j_match:j_match+nodes_size_y,
           i_match:i_match+nodes_size_x] = data['interior'][1][0,:,:]

    if(data_type==grdptid_type):
	i_match = interspace + bf_tmp_size_x + interspace
	j_match = interspace + bf_tmp_size_y + interspace
	lm[j_match:j_match+nodes_size_y,
	   i_match:i_match+bc_size] = bc_pt
	lm[j_match:j_match+nodes_size_y,
	   i_match+nodes_size_x-bc_size:i_match+nodes_size_x] = bc_pt
	lm[j_match:j_match+bc_size,
	   i_match:i_match+nodes_size_x] = bc_pt
	lm[j_match+nodes_size_y-bc_size:j_match+nodes_size_y,
	   i_match:i_match+nodes_size_x] = bc_pt
	lm[j_match+bc_size:j_match+nodes_size_y-bc_size,
	   i_match+bc_size:i_match+nodes_size_x-bc_size]=interior_pt    
    
    #3.6) W  buffer layer
    if(plot_buffer['W_']):
        bf_layer_size_x = data['W_'][0][0]
        bf_layer_size_y = data['W_'][0][1]
        i_match = interspace + bf_tmp_size_x - bf_layer_size_x
        j_match = interspace + bf_tmp_size_y + interspace + data['W_'][0][5] - 2*bc_size +1
        if(data_type==nodes_type):
            lm[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = data['W_'][data_type][0,:,:]
        if(data_type==grdptid_type):
            lm[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = data['W_'][data_type][:,:]
    
    #3.7) NE buffer layer
    if(plot_buffer['NE']):
        bf_layer_size_x = data['NE'][0][0]
        bf_layer_size_y = data['NE'][0][1]
        i_match = interspace + bf_tmp_size_x + interspace + nodes_size_x + interspace
        j_match = interspace + bf_tmp_size_y + interspace + nodes_size_y + interspace
        if(data_type==nodes_type):
            lm[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = data['NE'][data_type][0,:,:]
        if(data_type==grdptid_type):
            lm[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = data['NE'][data_type][:,:]
    
    #3.8) N  buffer layer
    if(plot_buffer['N_']):
        bf_layer_size_x = data['N_'][0][0]
        bf_layer_size_y = data['N_'][0][1]

        i_match = interspace + bf_tmp_size_x + interspace + data['N_'][0][3] - bc_size -1
        j_match = interspace + bf_tmp_size_y + interspace + nodes_size_y + interspace
        if(data_type==nodes_type):    
            lm[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = data['N_'][data_type][0,:,:]
        if(data_type==grdptid_type):
            lm[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = data['N_'][data_type][:,:]
        
    #3.9) NW buffer layer
    if(plot_buffer['NW']):
        bf_layer_size_x = data['NW'][0][0]
        bf_layer_size_y = data['NW'][0][1]
        i_match = interspace + bf_tmp_size_x - bf_layer_size_x
        j_match = interspace + bf_tmp_size_y + interspace + nodes_size_y + interspace
        if(data_type==nodes_type):
            lm[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = data['NW'][data_type][0,:,:]
        if(data_type==grdptid_type):
            lm[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = data['NW'][data_type][:,:]
        

    #rearrange the large matrix for plotting
    lm = lm[::-1,::]

    return lm

    
def plot_matrix_with_all_buffer_layers(lm):
    
    #plot the matrix with all the buffer layers
    fig=plt.figure(figsize=(12,6))
    ax = fig.add_subplot(1,1,1)
    res = ax.imshow(lm, cmap=cm.summer, interpolation='nearest')

    return fig,ax

def plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes, lm_grdptid):

    #create the main figure
    fig=plt.figure(figsize=(12,6))

    #plot the gridpoint ID
    ax = fig.add_subplot(1,2,1)
    res = ax.imshow(lm_grdptid, cmap=cm.spectral, interpolation='nearest', vmin=-1, vmax=4)
    fig.colorbar(res)
    
    #plot the nodes
    ax = fig.add_subplot(1,2,2)
    res = ax.imshow(lm_nodes, cmap=cm.spectral, interpolation='nearest', vmin=0.0, vmax=1.0)
    fig.colorbar(res)    

    return fig,ax
    
if __name__ == "__main__":
    

    #manage the options
    [folder_path] = manage_options()


    #test the allocation procedure
    #=================================================================
    
    #extract data for the interior points and the buffer layers
    #-----------------------------------------------------------------
    #data container
    data = {}

    #extract the data for the interior points
    data['interior'] = extract_interior_data(folder_path)

    #extract the data of the buffer layers
    bf_layer_loc_table = ['N_','S_','E_','W_','NE','NW','SE','SW']
    
    suffix_files={}
    suffix_files['sizes']  ='_sizes.dat'
    suffix_files['nodes']  ='_nodes.dat'
    suffix_files['grdptid']='_grdpt_id.dat'

    for bf_layer_loc in bf_layer_loc_table:
        data[bf_layer_loc] = extract_bf_data(folder_path,
                                             suffix_files,
                                             bf_layer_loc)

    #create the large matrix containing the data for the gridpoint id
    lm_grdptid = create_matrix_with_all_bf_layers(data,
                                                  grdptid_type)

    #create the large matrix containing the data for the nodes
    lm_nodes = create_matrix_with_all_bf_layers(data,
                                                nodes_type)
    
    #display
    #-----------------------------------------------------------------
    fig, ax = plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes,
                                                        lm_grdptid)
    fig.canvas.set_window_title("Allocation test")


    #test the reallocation procedure
    #=================================================================
    
    #extract data for the interior points and the buffer layers
    #-----------------------------------------------------------------
    #data container
    data = {}

    #extract the data for the interior points
    data['interior'] = extract_interior_data(folder_path)

    #extract the data of the buffer layers
    bf_layer_loc_table = ['N_','S_','E_','W_','NE','NW','SE','SW']
    suffix_files={}
    suffix_files['sizes']  ='_sizes2.dat'
    suffix_files['nodes']  ='_nodes2.dat'
    suffix_files['grdptid']='_grdpt_id2.dat' 

    for bf_layer_loc in bf_layer_loc_table:
        data[bf_layer_loc] = extract_bf_data(folder_path,
                                             suffix_files,
                                             bf_layer_loc)

    #create the large matrix containing the data for the gridpoint id
    lm_grdptid = create_matrix_with_all_bf_layers(data,
                                                  grdptid_type)

    #create the large matrix containing the data for the nodes
    lm_nodes = create_matrix_with_all_bf_layers(data,
                                                nodes_type)
    
    #display
    #-----------------------------------------------------------------
    fig2, ax2 = plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes,
                                                          lm_grdptid)
    fig2.canvas.set_window_title("Reallocation test")


    #test the new gridpoints procedure
    #=================================================================
    
    #extract data for the interior points and the buffer layers
    #-----------------------------------------------------------------
    #data container
    data = {}
    
    #extract the data for the interior points
    data['interior'] = extract_interior_data(folder_path)
    
    #extract the data of the buffer layers
    bf_layer_loc_table = ['N_','S_','E_','W_','NE','NW','SE','SW']
    suffix_files={}
    suffix_files['sizes']  ='_sizes3.dat'
    suffix_files['nodes']  ='_nodes3.dat'
    suffix_files['grdptid']='_grdpt_id3.dat'
    
    for bf_layer_loc in bf_layer_loc_table:
        data[bf_layer_loc] = extract_bf_data(folder_path,
                                             suffix_files,
                                             bf_layer_loc)
    
    #create the large matrix containing the data for the gridpoint id
    lm_grdptid = create_matrix_with_all_bf_layers(data,
                                                  grdptid_type)
    
    #create the large matrix containing the data for the nodes
    lm_nodes = create_matrix_with_all_bf_layers(data,
                                                nodes_type)
    
    #display
    #-----------------------------------------------------------------
    fig3, ax3 = plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes,
                                                          lm_grdptid)
    fig3.canvas.set_window_title("New gridpoints test")
    
    
    #show all
    plt.show()
    
