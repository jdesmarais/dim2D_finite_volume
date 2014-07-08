#!/usr/bin/python


import commands
import getopt

import os, sys, inspect
import shutil
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
                                   get_filenames,
                                   make_picture_bf_layers)
    
if __name__ == "__main__":

    #manage the options
    [folder_path] = manage_options()
    print folder_path

    
    #picture_folder
    picture_folder_path = './movie_pictures'

    
    #set the maximum number of sublayers
    #per mainlayer
    nb_sublayers = 6


    #create the folder where all pictures
    #for the movie are saved
    if os.path.exists(picture_folder_path):
        shutil.rmtree(picture_folder_path)
    if not os.path.exists(picture_folder_path):
        os.makedirs(picture_folder_path)


    #extract the file names for the
    #first timestep
    i=0
    input_filenames = get_filenames(folder_path,i)


    #create a picture for each timestep written
    #in the folder_path
    while (i<100 and os.path.isfile(input_filenames[0])):

        #creation of the picture
        make_picture_bf_layers(
            picture_folder_path,
            input_filenames,
            folder_path,
            nb_sublayers,
            i,
            index_max=100,
            continuous=True)

        #get the input files for the next picture
        i+=1
        input_filenames = get_filenames(folder_path,i)
    

    
