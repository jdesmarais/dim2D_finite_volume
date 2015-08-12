#!/usr/bin/python

import sys
import os
import getopt
import subprocess
import shlex
import shutil

fortranEnergyTrExe='extract_energyTr'

fortranEnergyTrExePath=os.path.join(os.getenv(
        'augeanstables'),
        'scripts_fortran',
        'energyTr_computation',
        fortranEnergyTrExe)


# error message print
def print_mg_error(error_mg):
    sys.stdout.write('****'+error_mg+'**** \n')


# progress message print
def print_mg_progress(mg_progress):
    '''
    @description:
    print a message which is overwritten
    '''
    sys.stdout.write('%s\r' % mg_progress)
    sys.stdout.flush()


# final progress message print
def print_mg_final(mg_progress):
    '''
    @description:
    print a message which is not overwritten
    '''
    sys.stdout.write('%s' % mg_progress+'\n')
    sys.stdout.flush()


# display the help for the program
def display_help():
    '''
    @description:
    display program help
    '''
    print ''
    print 'extract the energy flowing across the domain'
    print 'borders from netcdf files'
    print ''
    print 'options:'
    print '--------'
    print '-h (--help)    : display this help'
    print '-i (--input=)  : main dir with the netcdf files data*.nc'
    print '-o (--output=) : text file for the output'
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
                                   "hi:o:",
                                   ["help",
                                    "input=",
                                    "output="])
    except getopt.GetoptError:
        display_help()
        sys.exit(2)

    if(len(opts)==0):
        display_help()
        sys.exit(2)


    mainDir = 'None'
    outputFile = 'None'

    # options
    for opt, arg in opts:

        if opt in ("-h", "--help"):
            display_help()
            sys.exit(2)

        elif opt in ("-i", "--input"):
            mainDir = arg

        elif opt in ("-o", "--output"):
            outputFile = arg


    # check for directory with the netcdf files
    if(mainDir=='None'):
        print_mg_error('directory for netcdf file not provided')
        display_help()
        sys.exit(2)

    else:
        if( not os.path.isdir(mainDir)):
            print_mg_error(mainDir+' does not exist')
            display_help()
            sys.exit(2)


    # check for directory where the output file should be saved
    if(outputFile=='None'):
        print_mg_error('output path not provided')
        display_help()
        sys.exit(2)

    else:
        if( not os.path.isdir(os.path.dirname(outputFile))):
            print_mg_error('output dir does not exist')
            display_help()
            sys.exit(2)
    
    return mainDir, outputFile


# determine the netcdf files in the directory
def find_netcdf_data_files_in_dir(dirPath):
    '''
    @description:
    count the netcdf files in the directory
    '''

    # list with the netcdf files
    ncFiles = [name for name in os.listdir(dirPath)\
                   if (os.path.isfile(os.path.join(dirPath, name)) and\
                           name.endswith('.nc') and \
                           name.startswith('data') )]

    # first and last timestep of the netcdf files
    i_min = int(ncFiles[0].replace('data','').replace('.nc',''))
    i_max = i_min+len(ncFiles)-1

    return i_min,i_max


# generate the fortran executable to extract the
# energy flowing across the domain borders
def generate_fortran_exe():
    '''
    @description:
    generate the fortran exe to extract the energy
    flowing across the domain borders
    '''

    currentPath = os.getcwd()

    # if the exe does not exist, it should be generated
    if(not os.path.isfile(fortranEnergyTrExePath)):

        # change dir to where the fortran Exe is generated
        os.chdir(os.path.dirname(fortranEnergyTrExePath))

        # generate the executable
        exeOpt = os.getenv('AUGEANSTABLES_PROFILE')
        os.environ['AUGEANSTABLES_PROFILE']= 'false'
        cmd='make cleanall && make '+fortranEnergyTrExe+' make clean'
        subprocess.call(cmd, shell=True)
        os.environ['AUGEANSTABLES_PROFILE']= exeOpt

        # check whether the executable was generated
        if( not os.path.isfile(fortranEnergyTrExePath)):
            print_mg_error('fortran exe not generated: '+fortranEnergyTrExePath)
            sys.exit(2)

        # change back to the current path
        os.chdir(currentPath)

    # copy the fortran exe to the current dir
    shutil.copy(fortranEnergyTrExePath, currentPath)


# extract the energy flowing across the domain borders
# from a netcdf file using fortran exe
def extract_energyTr(ncFile):
    '''
    @description:
    extract the energy flowing across the domain borders
    from a netcdf file using fortran exe
    '''

    # extract the energyTr from output of fortran exe
    cmd='./'+fortranEnergyTrExe
    cmd+=" -i "+ncFile
    args = shlex.split(cmd)

    output = subprocess.Popen(args,stdout=subprocess.PIPE).communicate()[0]
    
    # convert to float
    try :
        energyTr= float(output.replace('energyTr:',''))

    except ValueError:

        print_mg_error('error when extracting energyTr')
        print_mg_error('fortranEnergyTrExe: '+fortranEnergyTrExe)
        print_mg_error('ncFile: '+ncFile)
        print_mg_error('output: '+output)
        sys.exit(2)

    return energyTr


# extract the time from a netcdf file using ncdump
def extract_time(ncFile):
    '''
    @description:
    extract the time from netcdf file using
    ncdump command
    '''

    # extract the energyTr from output of fortran exe
    cmd='ncdump -v time '+ncFile
    args = shlex.split(cmd)
    output = subprocess.Popen(args,stdout=subprocess.PIPE).communicate()[0]
    output = output.split('\n')[-3]

    # convert to float
    try :
        time= float(output.replace('time = ','').replace(';',''))
    
    except ValueError:
    
        print_mg_error('error when extracting time')
        print_mg_error('ncFile: '+ncFile)
        print_mg_error('output: '+output)
        sys.exit(2)

    return time


if __name__=="__main__":

    # extract the options from the command line
    mainDir, outputFile = parse_argv(sys.argv[1:])


    # determine the number of netcdf files in dir
    tmin,tmax = find_netcdf_data_files_in_dir(mainDir)


    # initialize the output file
    if(os.path.isfile(outputFile)):
        os.remove(outputFile)


    # generate the fortran exe to extract the
    # energy flowing across the domain borders
    generate_fortran_exe()


    # start the extraction of the energy flowing
    # across the domain borders
    mg_progress = 'extract energyTr: ...'
    print_mg_progress(mg_progress)


    # extract the energy flowing though the domain
    # borders for all netcdf files in the directory
    out = open(outputFile, 'w')

    for t in range(tmin,tmax+1):

        ncFile = os.path.join(mainDir,'data'+str(t)+'.nc')

        time     = extract_time(ncFile)
        energyTr = extract_energyTr(ncFile)

        out.write("%i %f %f\n" % (t,time,energyTr))

        mg_progress = 'extract energyTr: '+str(t-tmin+1)+' / '+str(tmax-tmin+1)
        print_mg_progress(mg_progress)
        
    mg_progress = 'extract energyTr: done       '
    print_mg_final(mg_progress)

    out.close()


    # remove the fortran exe
    os.remove(fortranEnergyTrExe)
    
