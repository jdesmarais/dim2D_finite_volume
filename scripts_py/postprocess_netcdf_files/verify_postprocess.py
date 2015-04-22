#!/usr/bin/python

import sys
import os
import getopt
import subprocess
import shutil

from merge_netcdf_files import\
    merge_netcdf_file

from truncate_netcdf_files import\
    truncate_netcdf_file,\
    create_mg_progress,\
    create_mg_final


#parse the input options
def parse_argv(argv):
    '''
    @description:
    parse the program options
    '''

    #default values
    mainDir = 'None'
    x_min   = 'None'
    x_max   = 'None'
    y_min   = 'None'
    y_max   = 'None'


    # store the options and the arguments
    # in opts, args
    try:
        opts, args = getopt.getopt(argv,
                                   "hi:t:",
                                   ["help",
                                    "input=",
                                    "x_min=",
                                    "x_max=",
                                    "y_min=",
                                    "y_max=",
                                    "nb_timesteps="])
    except getopt.GetoptError:
        display_help()
        sys.exit(2)

    if(len(opts)==0):
        display_help()
        sys.exit(2)


    for opt, arg in opts:

        if opt in ("-h", "--help"):
            display_help()
            sys.exit(2)

        elif opt in ("-i", "--input"):
            mainDir = arg

        elif opt in ("--x_min"):
            x_min = float(arg)

        elif opt in ("--x_max"):
            x_max = float(arg)

        elif opt in ("--y_min"):
            y_min = float(arg)

        elif opt in ("--y_max"):
            y_max = float(arg)

        elif opt in ("-t","--nb_timesteps"):
            nb_timesteps = int(arg)
    
    return [mainDir,x_min,x_max,y_min,y_max,nb_timesteps]


def create_mg_progress(mg_progress):
    sys.stdout.write('%s\r' % mg_progress)
    sys.stdout.flush()


def create_mg_final(mg_progress):
    sys.stdout.write('%s' % mg_progress)
    sys.stdout.flush()
    print '\n'


def get_truncation_filenames(t,mainDir):
    '''
    @description
    get the filenames for truncating files
    '''

    large_file = os.path.join(mainDir,'data'+str(t)+'.nc')
    small_file = os.path.join(mainDir,'tr_files','data'+str(t)+'.nc')

    return [large_file,small_file]


def get_merge_filenames(t,mainDir,parallelDir):
    '''
    @description
    get the filenames for merging files
    '''

    parallel_file = os.path.join(parallelDir,'data'+str(t)+'_0.nc')
    merge_file    = os.path.join(mainDir,'data'+str(t)+'.nc')
    new_merge_file= os.path.join(parallelDir,'data'+str(t)+'.nc')

    return [parallel_file,merge_file,new_merge_file]



if __name__=="__main__":

    # main directory where the simulation files are saved
    [mainDir,x_min,x_max,y_min,y_max,nb_timesteps] = parse_argv(sys.argv[1:])

    parallelDir = os.path.join(mainDir,'parallel_files')


    # loop over the files in the truncation folder
    mg_progress = 'check_netcdf_files:...'
    create_mg_progress(mg_progress)

    for t in range(0,nb_timesteps):

        [large_file,small_file] = get_truncation_filenames(t,mainDir)
        [parallel_file,merge_file,new_merge_file] = get_merge_filenames(t,mainDir,parallelDir)

        # check whether the truncation file exists:
        # if not, it means that the merging process
        # failed and it should be restarted
        if( not os.path.isfile(small_file) ):

            # remove the merged file
            if( os.path.isfile(merge_file) ):
                os.remove(merge_file)

            # ask the program to restart the
            # merge of this file
            merge_netcdf_file(t,parallelDir)

            # ask the program to restart the
            # truncation of the merged file
            truncate_netcdf_file(mainDir,t,x_min,x_max,y_min,y_max)

        # update the output in command line
        # for truncating files
        mg_progress = 'check_netcdf_files: '+str(t)+' '
        create_mg_progress(mg_progress)

    mg_progress = 'check_netcdf_files: done   '
    create_mg_final(mg_progress)
                
                
