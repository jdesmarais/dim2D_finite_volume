#!/usr/bin/python

import sys
import os
import getopt
import subprocess
import shutil

fortranExe='truncate_netcdf_file'

fortranExePath=os.path.join(os.getenv(
        'augeanstables'),
        'scripts_fortran',
        'truncate_netcdf_files',
        fortranExe)


# display the help for the program
def display_help():
    '''
    @description:
    display program help
    '''
    print ''
    print 'truncate the netcdf files from a parallel run'
    print ''
    print 'options:'
    print '--------'
    print '-h (--help)             : display this help'
    print '-i (--input=)           : main dir of the simulation'
    print ''


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
                                   "hi:",
                                   ["help",
                                    "input=",
                                    "x_min=",
                                    "x_max=",
                                    "y_min=",
                                    "y_max="])
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
    
    return [mainDir,x_min,x_max,y_min,y_max]


def create_mg_progress(mg_progress):
    sys.stdout.write('%s\r' % mg_progress)
    sys.stdout.flush()


def create_mg_final(mg_progress):
    sys.stdout.write('%s' % mg_progress)
    sys.stdout.flush()
    print '\n'


def get_filenames(t,mainDir):
    '''
    @description
    get the filenames for truncating files
    '''

    large_file = os.path.join(mainDir,'data'+str(t)+'.nc')
    small_file = os.path.join(mainDir,'tr_files','data'+str(t)+'.nc')

    return [large_file,small_file]


def truncate_netcdf_file(mainDir,t,x_min,x_max,y_min,y_max,ne=4):
    '''
    @description
    use a fortran exe to merge the parallel files
    '''

    # get the current path
    currentpath = os.getcwd()

    # change to the main dir
    os.chdir(mainDir)
    
    # check whether the fortran exe exists
    # other copy the executable from the exe folder
    if(not os.path.isfile(fortranExe)):
        if(not os.path.isfile(fortranExePath)):
            print 'the fortran executable should be generated'
            sys.exit(2)

        shutil.copy(fortranExePath, mainDir)

    # run the fortran exe to create the merge file
    cmd = fortranExePath
    cmd+=' -xmin '+str(x_min)
    cmd+=' -xmax '+str(x_max)
    cmd+=' -ymin '+str(y_min)
    cmd+=' -ymax '+str(y_max)
    cmd+=' -e '+str(ne)
    cmd+=' -t '+str(t)
    cmd+=' > /dev/null'

    subprocess.call(cmd, shell=True)

    # change back to the original dir
    os.chdir(currentpath)


if __name__=="__main__":

    # main directory where the simulation files are saved
    [mainDir,x_min,x_max,y_min,y_max] = parse_argv(sys.argv[1:])


    #if the folder where the truncated files exists, remove it
    truncationDir = os.path.join(mainDir,'tr_files')
    if(os.path.isdir(truncationDir)):
        shutil.rmtree(truncationDir)
    os.mkdir(truncationDir)

    # loop over the parallel files, if the corresponding
    # merged file exists, do nothing, otherwise use
    # exe to merge the netcdf files
    t=0
    [large_file,small_file] = get_filenames(t,mainDir)
    
    mg_progress = 'truncate_netcdf_files:...'
    create_mg_progress(mg_progress)

    while (os.path.isfile(large_file)):
    
        truncate_netcdf_file(mainDir,t,x_min,x_max,y_min,y_max)

        # update the output in command line
        # for truncating files
        mg_progress = 'truncate_netcdf_files: '+str(t)+' '
        create_mg_progress(mg_progress)
        
        # next file
        t+=1
        [large_file,small_file] = get_filenames(t,mainDir)

    mg_progress = 'truncate_netcdf_files: done   '
    create_mg_final(mg_progress)

    sys.stdout.write(str(t)+' files truncated \n')

