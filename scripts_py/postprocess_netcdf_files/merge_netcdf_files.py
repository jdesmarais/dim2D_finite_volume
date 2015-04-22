#!/usr/bin/python

import sys
import os
import getopt
import subprocess
import shutil

fortranMergeExe='merge_netcdf_file'

fortranMergeExePath=os.path.join(os.getenv(
        'augeanstables'),
        'scripts_fortran',
        'merge_netcdf_files',
        fortranMergeExe)


# display the help for the program
def display_help():
    '''
    @description:
    display program help
    '''
    print ''
    print 'merge the netcdf files from a parallel run'
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


    # store the options and the arguments
    # in opts, args
    try:
        opts, args = getopt.getopt(argv,
                                   "hi:",
                                   ["help",
                                    "input="])
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

    
    return mainDir


def create_mg_progress(mg_progress):
    sys.stdout.write('%s\r' % mg_progress)
    sys.stdout.flush()


def create_mg_final(mg_progress):
    sys.stdout.write('%s' % mg_progress)
    sys.stdout.flush()
    print '\n'


def get_filenames(t,mainDir,parallelDir):
    '''
    @description
    get the filenames for merging files
    '''

    parallel_file = os.path.join(parallelDir,'data'+str(t)+'_0.nc')
    merge_file    = os.path.join(mainDir,'data'+str(t)+'.nc')
    new_merge_file= os.path.join(parallelDir,'data'+str(t)+'.nc')

    return [parallel_file,merge_file,new_merge_file]


def merge_netcdf_file(t,parallel_folder,nb_tiles=[8,8],ne=4,bc_size=2):
    '''
    @description
    use a fortran exe to merge the parallel files
    '''
    
    currentPath = os.getcwd()

    # go to the folder where the parallel
    # files are saved
    os.chdir(parallel_folder)

    # check whether the fortran exe exists
    # other copy the executable from the exe folder
    if(not os.path.isfile(fortranMergeExe)):
        if(not os.path.isfile(fortranMergeExePath)):
            print 'the fortran executable should be generated'
            sys.exit(2)

        shutil.copy(fortranMergeExePath, parallel_folder)

    # run the fortran exe to create the merge file
    cmd = fortranMergeExePath
    cmd+=' -x '+str(nb_tiles[0])
    cmd+=' -y '+str(nb_tiles[1])
    cmd+=' -e '+str(ne)
    cmd+=' -b '+str(bc_size)
    cmd+=' -t '+str(t)

    subprocess.call(cmd, shell=True)

    # go back the previous directory
    os.chdir(currentPath)


if __name__=="__main__":

    # main directory where the simulation files are saved
    mainDir = parse_argv(sys.argv[1:])

    # directory where the parallel files are saved
    parallelDir = os.path.join(mainDir,'parallel_files')

    # loop over the parallel files, if the corresponding
    # merged file exists, do nothing, otherwise use
    # exe to merge the netcdf files
    t=0
    [parallel_file,merge_file,new_merge_file] = get_filenames(t,mainDir,parallelDir)
    
    mg_progress = 'merge_netcdf_files:...'
    create_mg_progress(mg_progress)

    while (os.path.isfile(parallel_file)):
    
        # check whether the corresponding merged file exist
        if( not os.path.isfile(merge_file) ):

            merge_netcdf_file(t,parallelDir)

            if(os.path.isfile(new_merge_file)):
                os.rename(new_merge_file,merge_file)

        # update the output in command line
        # for merging files
        mg_progress = 'merge_netcdf_files: '+str(t)+' '
        create_mg_progress(mg_progress)
        
        # next file
        t+=1
        [parallel_file,merge_file,new_merge_file] = get_filenames(t,mainDir,parallelDir)

    mg_progress = 'merge_netcdf_files: done   '
    create_mg_final(mg_progress)

    sys.stdout.write(str(t)+' files merged \n')

