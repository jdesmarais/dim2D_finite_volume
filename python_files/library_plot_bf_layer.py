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


#local variables defining the code for the grid points
backgrd_pt = -10000
no_pt = 0
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


def get_nb_sublayers(nb_sublayers_filename):

    #load the interior point sizes file
    #nb_sublayers_filename = folder_path+'/sublayers_nb.dat'
    f = FortranFile(nb_sublayers_filename)
    nb_sublayers = f.readInts()
    f.close()
    return nb_sublayers


def extract_interior_data(sizes_filename, grdptid_filename, nodes_filename):

    #load the interior point sizes file
    #sizes_filename = folder_path+'/interior_sizes.dat'
    f = FortranFile(sizes_filename)
    sizes = f.readInts()
    f.close()

    #load the interior point grdpt_id file
    #grdptid_filename = folder_path+'/interior_nodes.dat'
    f = FortranFile(grdptid_filename)
    grdpts_id = f.readInts()
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

    grdpts_id.resize(size_ne,size_y,size_x)
    nodes.resize(size_ne,size_y,size_x)

    #return the interior point data
    return [sizes, grdpts_id, nodes]


def extract_detectors_data(
    N_detector_filename,
    S_detector_filename,
    E_detector_filename,
    W_detector_filename):

    #N detectors
    f = FortranFile(N_detector_filename)
    N_detector = f.readInts()
    f.close()

    #S detectors
    f = FortranFile(S_detector_filename)
    S_detector = f.readInts()
    f.close()

    #E detectors
    f = FortranFile(E_detector_filename)
    E_detector = f.readInts()
    f.close()

    #W detectors
    f = FortranFile(W_detector_filename)
    W_detector = f.readInts()
    f.close()

    N_detector.resize(len(N_detector)/2,2)
    S_detector.resize(len(S_detector)/2,2)
    E_detector.resize(len(E_detector)/2,2)
    W_detector.resize(len(W_detector)/2,2)
    
    return [N_detector, S_detector, E_detector, W_detector]

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
    interior_size_filename,
    interior_grdptsid_filename,
    interior_nodes_filename,
    folder_path,
    nb_sublayers,
    suffix_size,
    suffix_nodes,
    suffix_grdptid,
    continuous=False):

    #possible mainlayers
    mainlayers_char = ['N_','S_','E_','W_','NE','NW','SE','SW']

    
    #extract the data for the interior nodes
    [nodes_sizes,grdpts_id,nodes] = extract_interior_data(
        interior_size_filename,
        interior_grdptsid_filename,
        interior_nodes_filename)
    nodes_size_x = nodes_sizes[0]
    nodes_size_y = nodes_sizes[1]
    
    
    #extract the size of the sublayers
    size_x = 0
    size_y = 0
    for mainlayer in mainlayers_char:

        #test for possible filenames
        for i in range(1,nb_sublayers+1):

            #if the file exists, extract the size of the sublayer
            filename = folder_path+'/'+mainlayer+str(i)+suffix_size
            
            if(os.path.isfile(filename)):

                #print filename used
                print filename+': analyzed'
                            
                
                #extract the size of the sublayer
                f = FortranFile(filename)
                sizes = f.readInts()
                f.close()

                #compare the size with the actual extent of the
                #large matrix
                compare_procedure = {
                    'N_' : compute_sizes_N,
                    'S_' : compute_sizes_S,
                    'E_' : compute_sizes_E,
                    'W_' : compute_sizes_W}
                
                [size_x,size_y] = compare_procedure[mainlayer](
                    size_x,
                    size_y,
                    nodes_size_x,
                    nodes_size_y,
                    sizes[3],
                    sizes[4],
                    sizes[5],
                    sizes[6])

            #else:
            #    #print filename not recognized
            #    print filename+': does not exist'

    #the sizes of the buffer layers are now determined
    bf_tmp_size_x = 25 #size_x
    bf_tmp_size_y = 25 #size_y

    #interspace
    interspace = 2
    margin = [interspace + bf_tmp_size_x+interspace-1,
              interspace + bf_tmp_size_y+interspace-1]
    

    lm_size_x = 2*bf_tmp_size_x + nodes_size_x + 4*interspace
    lm_size_y = 2*bf_tmp_size_y + nodes_size_y + 4*interspace

    
    #allocation of the large matrices for the nodes and the grdptid
    lm_nodes = np.empty([4,lm_size_y, lm_size_x])
    lm_nodes[0,:,:]=backgrd_pt #.fill(backgrd_pt)
    lm_nodes[1,:,:]=0 #.fill(backgrd_pt)
    lm_nodes[2,:,:]=0 #.fill(backgrd_pt)
    lm_nodes[3,:,:]=0 #.fill(backgrd_pt)

    lm_grdptid = np.empty([lm_size_y, lm_size_x])
    lm_grdptid.fill(backgrd_pt)
    

    #initialization of the large matrix with the nodes data
    #and the grdpts ID
    i_match = interspace + bf_tmp_size_x + interspace
    j_match = interspace + bf_tmp_size_y + interspace

    lm_grdptid[
        j_match:j_match+nodes_size_y,
        i_match:i_match+nodes_size_x] = grdpts_id[0,:,:]
    
    lm_nodes[
        :,
        j_match:j_match+nodes_size_y,
        i_match:i_match+nodes_size_x] = nodes[:,:,:]


    #fill the large matrix with possible sublayers
    for mainlayer in mainlayers_char:

        #test for possible filenames
        for i in range(1,nb_sublayers+1):

            #if the file exists, extract the size, the data and
            #the grdptID of the sublayer
            filename = folder_path+'/'+mainlayer+str(i)+suffix_size
            
            if(os.path.isfile(filename)):

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
                    lm_grdptid,
                    continuous)

    lm_nodes   = lm_nodes[::,::-1,::]
    lm_grdptid = lm_grdptid[::-1,::]

    return [lm_nodes,lm_grdptid,margin]


def add_detectors(
    N_detector, S_detector,
    E_detector, W_detector,
    lm_grdptid, margin):

    lm_grdptid = lm_grdptid[::-1,::]

    N_detector_color = 3.0
    S_detector_color = 0.0
    E_detector_color = 2.5
    W_detector_color = 1.75

    #N detectors
    for i in range(0,len(N_detector)):
        x = margin + N_detector[i,0]
        y = margin + N_detector[i,1]
        lm_grdptid[y,x] = N_detector_color+0.5*i/(len(N_detector))

    #S detectors
    for i in range(0,len(S_detector)):
        x = margin + S_detector[i,0]
        y = margin + S_detector[i,1]
        lm_grdptid[y,x] = S_detector_color+0.5*i/(len(S_detector))

    #E detectors
    for i in range(0,len(E_detector)):
        x = margin + E_detector[i,0]
        y = margin + E_detector[i,1]
        lm_grdptid[y,x] = E_detector_color+0.5*i/(len(E_detector))

    #W detectors
    for i in range(0,len(W_detector)):
        x = margin + W_detector[i,0]
        y = margin + W_detector[i,1]
        lm_grdptid[y,x] = W_detector_color+0.5*i/(len(W_detector))

    lm_grdptid = lm_grdptid[::-1,::]

    return [lm_grdptid]
 
def compute_sizes_N(size_x, size_y, nx, ny, align11, align12, align21, align22):
    new_size_x = max(0,bc_size-align11,align12-(nx-bc_size+1))
    new_size_y = max(0,align22-(ny-bc_size+1))

    size_x = max(size_x, new_size_x+2*bc_size+1)
    size_y = max(size_y, new_size_y+2*bc_size+1)
    
    return [size_x,size_y]

def compute_sizes_S(size_x, size_y, nx, ny, align11, align12, align21, align22):
    new_size_x = max(0,bc_size-align11,align12-(nx-bc_size+1))
    new_size_y = max(0,-align21+bc_size)

    size_x = max(size_x, new_size_x+2*bc_size+1)
    size_y = max(size_y, new_size_y+2*bc_size+1)
    
    return [size_x,size_y]

def compute_sizes_E(size_x, size_y, nx, ny, align11, align12, align21, align22):
    new_size_x = max(0,align12-(nx-bc_size+1))
    new_size_y = max(0,bc_size-align21,align22-(ny-bc_size+1))

    size_x = max(size_x, new_size_x+2*bc_size+1)
    size_y = max(size_y, new_size_y+2*bc_size+1)
    
    return [size_x,size_y]

def compute_sizes_W(size_x, size_y, nx, ny, align11, align12, align21, align22):
    new_size_x = max(0,bc_size-align11)
    new_size_y = max(0,bc_size-align21,align22-(ny-bc_size+1))

    size_x = max(size_x, new_size_x+2*bc_size+1)
    size_y = max(size_y, new_size_y+2*bc_size+1)
    
    return [size_x,size_y]

def fill_bf_layer_data(
    folder_path,
    mainlayer,
    sublayer,
    suffix_size,
    suffix_nodes,
    suffix_grdptid,
    interspace,
    bf_tmp_size_x,
    bf_tmp_size_y,
    nodes_size_x,
    nodes_size_y,
    lm_nodes,
    lm_grdptid,
    continuous):


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
    if(continuous):
        extra_alignment = interspace + 2*bc_size
    else:
        extra_alignment = 0

    if(mainlayer=='SW'):        
    	i_match = interspace + bf_tmp_size_x - bf_layer_size_x + extra_alignment
    	j_match = interspace + bf_tmp_size_y - bf_layer_size_y + extra_alignment

    if(mainlayer=='S_'):
        i_match = interspace + bf_tmp_size_x + interspace + sizes[3] - 2*bc_size +1
    	j_match = interspace + bf_tmp_size_y - bf_layer_size_y + extra_alignment

    if(mainlayer=='SE'):
        i_match = interspace + bf_tmp_size_x + interspace + nodes_size_x + interspace - extra_alignment
        j_match = interspace + bf_tmp_size_y - bf_layer_size_y + extra_alignment

    if(mainlayer=='E_'):
        i_match = interspace + bf_tmp_size_x + interspace + nodes_size_x + interspace - extra_alignment
        j_match = interspace + bf_tmp_size_y + interspace + sizes[5] - 2*bc_size +1 

    if(mainlayer=='W_'):
        i_match = interspace + bf_tmp_size_x - bf_layer_size_x + extra_alignment
        j_match = interspace + bf_tmp_size_y + interspace + sizes[5] - 2*bc_size +1

    if(mainlayer=='NE'):
        i_match = interspace + bf_tmp_size_x + interspace + nodes_size_x + interspace - extra_alignment
        j_match = interspace + bf_tmp_size_y + interspace + nodes_size_y + interspace - extra_alignment

    if(mainlayer=='N_'):
        i_match = interspace + bf_tmp_size_x + interspace + sizes[3] - 2*bc_size +1
        j_match = interspace + bf_tmp_size_y + interspace + nodes_size_y + interspace - extra_alignment

    if(mainlayer=='NW'):
        i_match = interspace + bf_tmp_size_x - bf_layer_size_x + extra_alignment
        j_match = interspace + bf_tmp_size_y + interspace + nodes_size_y + interspace - extra_alignment
 
    if(continuous):
        selective_copy(mainlayer,
                   lm_nodes, lm_grdptid,
                   nodes, grdptid,
                   i_match, j_match,
                   bf_layer_size_x, bf_layer_size_y)
    else:
        #print mainlayer
        #
        #print shape(nodes)
        #print bf_layer_size_y
        #print bf_layer_size_x

        #print i_match, j_match, i_match+bf_layer_size_x, j_match+bf_layer_size_y, shape(lm_nodes)

        lm_nodes[:,j_match:j_match+bf_layer_size_y,
                 i_match:i_match+bf_layer_size_x] = nodes[:,:,:]
    
        #print shape(grdptid)
        #print bf_layer_size_y
        #print bf_layer_size_x

        lm_grdptid[j_match:j_match+bf_layer_size_y,
                   i_match:i_match+bf_layer_size_x] = grdptid[:,:]

    return [lm_nodes,lm_grdptid]


def selective_copy(mainlayer,
                   lm_nodes, lm_grdptid,
                   nodes, grdptid,
                   i_match, j_match,
                   bf_layer_size_x, bf_layer_size_y):

    #copy what is sure to be copied
    if(mainlayer=='S_'):
        lm_j_min = j_match
        lm_j_max = j_match+bf_layer_size_y
        lm_i_min = i_match
        lm_i_max = i_match+bf_layer_size_x
        
        j_min = 0
        j_max = nodes.shape[1]
        i_min = 0
        i_max = nodes.shape[2]

        conditional_copy(lm_nodes, lm_grdptid,
                         nodes, grdptid,
                         lm_i_min, lm_i_max,
                         lm_j_min, lm_j_max,
                         i_min, i_max,
                         j_min, j_max)

        #lm_nodes[j_match:j_match+bf_layer_size_y-bc_size,
        #         i_match:i_match+bf_layer_size_x] = nodes[0,0:-bc_size,:]
        #
        #lm_grdptid[j_match:j_match+bf_layer_size_y-bc_size,
        #           i_match:i_match+bf_layer_size_x] = grdptid[0:-bc_size,:]

    if(mainlayer=='N_'):
        lm_j_min = j_match
        lm_j_max = j_match+bf_layer_size_y
        lm_i_min = i_match
        lm_i_max = i_match+bf_layer_size_x
        
        j_min = 0
        j_max = nodes.shape[1]
        i_min = 0
        i_max = nodes.shape[2]

        conditional_copy(lm_nodes, lm_grdptid,
                         nodes, grdptid,
                         lm_i_min, lm_i_max,
                         lm_j_min, lm_j_max,
                         i_min, i_max,
                         j_min, j_max)

        #lm_nodes[j_match+bc_size:j_match+bf_layer_size_y,
        #         i_match:i_match+bf_layer_size_x] = nodes[0,bc_size:,:]
        #
        #lm_grdptid[j_match+bc_size:j_match+bf_layer_size_y,
        #           i_match:i_match+bf_layer_size_x] = grdptid[bc_size:,:]

    if(mainlayer=='E_'):
        lm_j_min = j_match
        lm_j_max = j_match+bf_layer_size_y
        lm_i_min = i_match
        lm_i_max = i_match+bf_layer_size_x
        
        j_min = 0
        j_max = nodes.shape[1]
        i_min = 0
        i_max = nodes.shape[2]

        conditional_copy(lm_nodes, lm_grdptid,
                         nodes, grdptid,
                         lm_i_min, lm_i_max,
                         lm_j_min, lm_j_max,
                         i_min, i_max,
                         j_min, j_max)

        #lm_nodes[j_match:j_match+bf_layer_size_y,
        #         i_match+bc_size:i_match+bf_layer_size_x] = nodes[0,:,bc_size:]
        #
        #lm_grdptid[j_match:j_match+bf_layer_size_y,
        #           i_match+bc_size:i_match+bf_layer_size_x] = grdptid[:,bc_size:]

    if(mainlayer=='W_'):
        lm_j_min = j_match
        lm_j_max = j_match+bf_layer_size_y
        lm_i_min = i_match
        lm_i_max = i_match+bf_layer_size_x
        
        j_min = 0
        j_max = nodes.shape[1]
        i_min = 0
        i_max = nodes.shape[2]

        conditional_copy(lm_nodes, lm_grdptid,
                         nodes, grdptid,
                         lm_i_min, lm_i_max,
                         lm_j_min, lm_j_max,
                         i_min, i_max,
                         j_min, j_max)

        #lm_nodes[j_match:j_match+bf_layer_size_y,
        #         i_match:i_match+bf_layer_size_x-bc_size] = nodes[0,:,:-bc_size]
        #
        #lm_grdptid[j_match:j_match+bf_layer_size_y,
        #           i_match:i_match+bf_layer_size_x-bc_size] = grdptid[:,:-bc_size]


def conditional_copy(
    lm_nodes, lm_grdptid,
    nodes, grdptid,
    lm_i_min, lm_i_max,
    lm_j_min, lm_j_max,
    i_min, i_max,
    j_min, j_max):

    for (lm_j,j) in zip(range(lm_j_min,lm_j_max), range(j_min,j_max)):
        for (lm_i,i) in zip(range(lm_i_min,lm_i_max), range(i_min,i_max)):
            if(grdptid[j,i]!=no_pt):
                
                lm_grdptid[lm_j,lm_i] = grdptid[j,i]
                lm_nodes[:,lm_j,lm_i] = nodes[:,j,i]

            
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
    res = ax.imshow(lm_nodes[0,:,:], cmap=cm.spectral, interpolation='nearest', vmin=0.0, vmax=1.0)
    fig.colorbar(res)    

    return fig,ax


def plot_all_bf_layers(lm_nodes, lm_grdptid):

    #create the main figure
    fig=plt.figure(figsize=(12,6))

    #plot the gridpoint ID
    ax = fig.add_subplot(1,3,1)
    res = ax.imshow(lm_grdptid, cmap=cm.spectral, interpolation='nearest', vmin=-1, vmax=4)
    fig.colorbar(res)
    
    #plot the nodes
    ax = fig.add_subplot(1,3,2)
    res = ax.imshow(lm_nodes[0,:,:], cmap=cm.spectral, interpolation='nearest', vmin=0.0, vmax=1.0)
    fig.colorbar(res)

    #plot the velocity
    ax = fig.add_subplot(1,3,3)
    Q = quiver(lm_nodes[1,::,::]/lm_nodes[0,::,::],
               lm_nodes[2,::,::]/lm_nodes[0,::,::],
               scale=0.5)
    #qk = quiverkey(Q, 0., 1.0, 1, r'$2 \frac{m}{s}$',
    #               labelpos='E',
    #               coordinates='figure',
    #               fontproperties={'weight': 'bold'})

    return fig,ax

