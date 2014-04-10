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


def get_nb_sublayers(folder_path):

    #load the interior point sizes file
    nb_sublayers_filename = folder_path+'/sublayers_nb.dat'
    f = FortranFile(nb_sublayers_filename)
    nb_sublayers = f.readInts()
    f.close()
    return nb_sublayers


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
    size_x  = sizes[0]
    size_y  = sizes[1]
    size_ne = sizes[2]

    nodes.resize(size_ne,size_y,size_x)
    nodes = nodes[::,::,::]

    #return the interior point data
    return [sizes, nodes]


def extract_bf_layer_data(
    sizes_filename,
    nodes_filename,
    grdptid_filename):


    #load buffer layer sizes file
    f = FortranFile(sizes_filename)
    sizes = f.readInts()
    f.close()

    #load buffer layer nodes file
    f = FortranFile(nodes_filename)
    nodes    = f.readReals('d')
    f.close()

    #load buffer layer gridpt id file
    f = FortranFile(grdptid_filename)
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


def make_matrix_for_all_bf_layers(
    folder_path,
    nb_sublayers,
    suffix_size,
    suffix_nodes,
    suffix_grdptid):

    #possible mainlayers
    mainlayers_char = ['N_','S_','E_','W_','NE','NW','SE','SW']
    
    #extract the size of the sublayers
    size_x = 0
    size_y = 0
    for mainlayer in mainlayers_char:

        #test for possible filenames
        for i in range(1,nb_sublayers+1):

            #if the file exists, extract the size of the sublayer
            filename = folder_path+'/'+mainlayer+str(i)+suffix_size
            if(os.path.isfile(filename)):
                
                #extract the size of the sublayer
                f = FortranFile(filename)
                sizes = f.readInts()
                f.close()

                #compare the size with the actual extent of the
                #large matrix
                compare_procedure = {
                    'N_' : update_size_y,
                    'S_' : update_size_y,
                    'E_' : update_size_x,
                    'W_' : update_size_x,
                    'NE' : update_size_x_and_y,
                    'NW' : update_size_x_and_y,
                    'SE' : update_size_x_and_y,
                    'SW' : update_size_x_and_y}

                [size_x,size_y] = compare_procedure[mainlayer](
                    size_x,
                    size_y,
                    sizes[0],
                    sizes[1])

    #the sizes of the buffer layers are now determined
    bf_tmp_size_x = size_x
    bf_tmp_size_y = size_y

    #extract the data for the interior nodes
    [nodes_sizes,nodes] = extract_interior_data(folder_path)
    nodes_size_x = nodes_sizes[0]
    nodes_size_y = nodes_sizes[1]
    
    #interspace
    interspace = 2

    lm_size_x = 2*bf_tmp_size_x + nodes_size_x + 4*interspace
    lm_size_y = 2*bf_tmp_size_y + nodes_size_y + 4*interspace

    
    #allocation of the large matrices for the nodse and the grdptid
    lm_nodes = np.empty([lm_size_y, lm_size_x])
    lm_nodes.fill(backgrd_pt)

    lm_grdptid = np.empty([lm_size_y, lm_size_x])
    lm_grdptid.fill(backgrd_pt)
    

    #initialization of the large matrix with the nodes data
    #and the grdpts ID
    i_match = interspace + bf_tmp_size_x + interspace
    j_match = interspace + bf_tmp_size_y + interspace
    
    lm_nodes[
        j_match:j_match+nodes_size_y,
        i_match:i_match+nodes_size_x] = nodes[0,:,:]

    lm_grdptid[j_match:j_match+nodes_size_y,
       i_match:i_match+bc_size] = bc_pt
    lm_grdptid[j_match:j_match+nodes_size_y,
       i_match+nodes_size_x-bc_size:i_match+nodes_size_x] = bc_pt
    lm_grdptid[j_match:j_match+bc_size,
       i_match:i_match+nodes_size_x] = bc_pt
    lm_grdptid[j_match+nodes_size_y-bc_size:j_match+nodes_size_y,
       i_match:i_match+nodes_size_x] = bc_pt
    lm_grdptid[j_match+bc_size:j_match+nodes_size_y-bc_size,
       i_match+bc_size:i_match+nodes_size_x-bc_size]=interior_pt  

    #fill the large matrix with possible sublayers
    for mainlayer in mainlayers_char:

        #test for possible filenames
        for i in range(1,nb_sublayers+1):

            #if the file exists, extract the size, the data and
            #the grdptID of the sublayer
            filename = folder_path+'/'+mainlayer+str(i)+suffix_size
            
            if(os.path.isfile(filename)):

                print filename

                [lm_nodes,lm_grdptid] = fill_bf_layer_data(
                    folder_path,
                    mainlayer,
                    i,
                    suffix_size,
                    suffix_nodes,
                    suffix_grdptid,
                    interspace,
                    bf_tmp_size_x,
                    bf_tmp_size_y,
                    nodes_size_x,
                    nodes_size_y,
                    lm_nodes,
                    lm_grdptid)

    lm_nodes   = lm_nodes[::-1,::]
    lm_grdptid = lm_grdptid[::-1,::]
    
    return [lm_nodes,lm_grdptid]

 
def update_size_x(size_x,size_y,new_size_x,new_size_y):
    size_x = max(size_x,new_size_x)
    return [size_x,size_y]
 
def update_size_y(size_x,size_y,new_size_x,new_size_y):
    size_y = max(size_y,new_size_y)
    return [size_x,size_y]    

def update_size_x_and_y(size_x,size_y,new_size_x,new_size_y):
    size_x = max(size_x,new_size_x)
    size_y = max(size_y,new_size_y)
    return [size_x,size_y]

def fill_bf_layer_data(
    folder_path,
    mainlayer,
    sublayer,
    suffix_size,
    suffix_nodes,
    suffix_grpdtid,
    interspace,
    bf_tmp_size_x,
    bf_tmp_size_y,
    nodes_size_x,
    nodes_size_y,
    lm_nodes,
    lm_grdptid):


    #determine the filenames
    size_filename    = folder_path+'/'+mainlayer+str(sublayer)+suffix_size
    nodes_filename   = folder_path+'/'+mainlayer+str(sublayer)+suffix_nodes
    grdptid_filename = folder_path+'/'+mainlayer+str(sublayer)+suffix_grdptid


    #extract the size, nodes and grdptid of the sublayer
    [sizes,nodes,grdptid] = extract_bf_layer_data(
            size_filename,
            nodes_filename,
            grdptid_filename)
    
    bf_layer_size_x = sizes[0]
    bf_layer_size_y = sizes[1]
    

    #determine the matched indices between the large matrix
    #and the sublayers
    if(mainlayer=='SW'):        
    	i_match = interspace + bf_tmp_size_x - bf_layer_size_x
    	j_match = interspace + bf_tmp_size_y - bf_layer_size_y

    if(mainlayer=='S_'):
        i_match = interspace + bf_tmp_size_x + interspace + sizes[3] - 2*bc_size +1
    	j_match = interspace + bf_tmp_size_y - bf_layer_size_y

    if(mainlayer=='SE'):
        i_match = interspace + bf_tmp_size_x + interspace + nodes_size_x + interspace
        j_match = interspace + bf_tmp_size_y - bf_layer_size_y

    if(mainlayer=='E_'):
        i_match = interspace + bf_tmp_size_x + interspace + nodes_size_x + interspace
        j_match = interspace + bf_tmp_size_y + interspace + sizes[5] - 2*bc_size +1

    if(mainlayer=='W_'):
        i_match = interspace + bf_tmp_size_x - bf_layer_size_x
        j_match = interspace + bf_tmp_size_y + interspace + sizes[5] - 2*bc_size +1

    if(mainlayer=='NE'):
        i_match = interspace + bf_tmp_size_x + interspace + nodes_size_x + interspace
        j_match = interspace + bf_tmp_size_y + interspace + nodes_size_y + interspace

    if(mainlayer=='N_'):
        i_match = interspace + bf_tmp_size_x + interspace + sizes[3] - bc_size -1
        j_match = interspace + bf_tmp_size_y + interspace + nodes_size_y + interspace

    if(mainlayer=='NW'):
        i_match = interspace + bf_tmp_size_x - bf_layer_size_x
        j_match = interspace + bf_tmp_size_y + interspace + nodes_size_y + interspace
 
    lm_nodes[j_match:j_match+bf_layer_size_y,
             i_match:i_match+bf_layer_size_x] = nodes[0,:,:]
    
    lm_grdptid[j_match:j_match+bf_layer_size_y,
               i_match:i_match+bf_layer_size_x] = grdptid[:,:]

    return [lm_nodes,lm_grdptid]

    
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

    #find the maximum number of sublayers per main layer
    nb_sublayers = get_nb_sublayers(folder_path)


    #test the allocation procedure
    #=================================================================
    #combine data from several sublayers in one large matrix
    #-----------------------------------------------------------------
    suffix_size    = '_sizes.dat'
    suffix_nodes   = '_nodes.dat'
    suffix_grdptid = '_grdpt_id.dat'
    
    [lm_nodes,lm_grdptid] = make_matrix_for_all_bf_layers(folder_path,
                                                          nb_sublayers,
                                                          suffix_size,
                                                          suffix_nodes,
                                                          suffix_grdptid)
    
    #display
    #-----------------------------------------------------------------
    fig, ax = plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes,
                                                        lm_grdptid)
    fig.canvas.set_window_title("Allocation test")


    ##test the reallocation procedure
    ##=================================================================
    #
    ##extract data for the interior points and the buffer layers
    ##-----------------------------------------------------------------
    ##data container
    #data = {}
    #
    ##extract the data for the interior points
    #data['interior'] = extract_interior_data(folder_path)
    #
    ##extract the data of the buffer layers
    #bf_layer_loc_table = ['N_','S_','E_','W_','NE','NW','SE','SW']
    #suffix_files={}
    #suffix_files['sizes']  ='_sizes2.dat'
    #suffix_files['nodes']  ='_nodes2.dat'
    #suffix_files['grdptid']='_grdpt_id2.dat' 
    #
    #for bf_layer_loc in bf_layer_loc_table:
    #    data[bf_layer_loc] = extract_bf_data(folder_path,
    #                                         suffix_files,
    #                                         bf_layer_loc)
    #
    ##create the large matrix containing the data for the gridpoint id
    #lm_grdptid = create_matrix_with_all_bf_layers(data,
    #                                              grdptid_type)
    #
    ##create the large matrix containing the data for the nodes
    #lm_nodes = create_matrix_with_all_bf_layers(data,
    #                                            nodes_type)
    #
    ##display
    ##-----------------------------------------------------------------
    #fig2, ax2 = plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes,
    #                                                      lm_grdptid)
    #fig2.canvas.set_window_title("Reallocation test")
    #
    #
    ##test the new gridpoints procedure
    ##=================================================================
    #
    ##extract data for the interior points and the buffer layers
    ##-----------------------------------------------------------------
    ##data container
    #data = {}
    #
    ##extract the data for the interior points
    #data['interior'] = extract_interior_data(folder_path)
    #
    ##extract the data of the buffer layers
    #bf_layer_loc_table = ['N_','S_','E_','W_','NE','NW','SE','SW']
    #suffix_files={}
    #suffix_files['sizes']  ='_sizes3.dat'
    #suffix_files['nodes']  ='_nodes3.dat'
    #suffix_files['grdptid']='_grdpt_id3.dat'
    #
    #for bf_layer_loc in bf_layer_loc_table:
    #    data[bf_layer_loc] = extract_bf_data(folder_path,
    #                                         suffix_files,
    #                                         bf_layer_loc)
    #
    ##create the large matrix containing the data for the gridpoint id
    #lm_grdptid = create_matrix_with_all_bf_layers(data,
    #                                              grdptid_type)
    #
    ##create the large matrix containing the data for the nodes
    #lm_nodes = create_matrix_with_all_bf_layers(data,
    #                                            nodes_type)
    #
    ##display
    ##-----------------------------------------------------------------
    #fig3, ax3 = plot_nodes_and_grdptid_with_all_bf_layers(lm_nodes,
    #                                                      lm_grdptid)
    #fig3.canvas.set_window_title("New gridpoints test")
    
    
    #show all
    plt.show()
    
