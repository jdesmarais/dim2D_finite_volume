#!/usr/bin/python


'''
@description: compute the heat distribution, i.e.
the amount of heat which is transported across the
domain borders, the heat used for the phase transition
and the heat used to increase the internal energy of
the fluid
'''


# library imports
#------------------------------------------------------------
import sys
import os
import getopt
import subprocess
import shlex
import numpy as np

from library_messages import print_mg_error
from compute_energyPh import check_inputFile


# functions
#------------------------------------------------------------

# display the help for the program
def display_help():
    '''
    @description:
    display program help
    '''
    print ''
    print 'compute the energy ratio from the energy supplied by'
    print 'the wall, the energy transported across the domain borders'
    print 'and the energy used for the phase transition'
    print ''
    print 'options:'
    print '--------'
    print '-h (--help)       : display this help'
    print '-w (--wall=)      : netcdf file whose header contains the wall energy'
    print '-t (--transport=) : text file with the energy transported'
    print '-p (--phase=)     : text file with the energy for phase transition'
    print '-o (--output=)    : text file for the output'
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
                                   "hi:w:t:p:o:",
                                   ["help",
                                    "wall=",
                                    "transport=",
                                    "phase=",
                                    "output="])
    except getopt.GetoptError:
        display_help()
        sys.exit(2)

    if(len(opts)==0):
        display_help()
        sys.exit(2)


    wallFile = 'None'
    transportFile = 'None'
    phaseFile = 'None'
    outputFile = 'None'


    # options
    for opt, arg in opts:

        if opt in ("-h", "--help"):
            display_help()
            sys.exit(2)

        elif opt in ("-w", "--wall"):
            wallFile = arg

        elif opt in ("-t", "--transport"):
            transportFile = arg

        elif opt in ("-p", "--phase"):
            phaseFile = arg

        elif opt in ("-o", "--output"):
            outputFile = arg


    # check the wall file
    check_inputFile(wallFile,errorMg='wall file not provided')


    # check the transport file
    check_inputFile(transportFile,errorMg='transport file not provided')


    # check the phase file
    check_inputFile(phaseFile,errorMg='phase file not provided')


    # check for directory where the output file should be saved
    if(outputFile=='None'):
        print_mg_error('output path not provided')
        display_help()
        sys.exit(2)

    else:
        if( not os.path.isdir(os.path.dirname(outputFile)) and (not os.path.dirname(outputFile)=='')):
            print_mg_error('output dir does not exist')
            display_help()
            sys.exit(2)
    
    return wallFile, transportFile, phaseFile, outputFile


# extract the energy rate supplied by the wall from the
# header of the netcdf file
def extract_energyWa(wallFile):
    '''
    @description: extract the energy rate supplied by
    the wall from the header of the netcdf file
    '''

    cmd='ncdump -h '+wallFile
    args = shlex.split(cmd)
    output = subprocess.Popen(args,stdout=subprocess.PIPE).communicate()[0]

    output = output.split('wall_maximum_heat_source')[1].split(';')[0].replace('=','')

    try:
        energyWa = float(output)

    except ValueError:
        print_mg_error('extract_energyWa')
        print_mg_error('error when converting wall energy to float')
        print_mg_error('output from ncdump: ',output)
        sys.exit(2)

    return energyWa


# compute the energy distribution between transport,
# phase transition, and internal energy
def compute_energyRa(wallFile, transportFile, phaseFile, outputFile):
    '''
    @description: compute the energy distribution
    between transport, phase transition, and
    internal energy
    '''

    # extract the heat flux provided by the wall from
    # the netcdf file, wallFile
    energyWa = extract_energyWa(wallFile)


    # extract the data for the energy flowing
    # across the domain borders
    energyTr = np.loadtxt(transportFile)


    # extract the data for the energy
    # used by the phase transition process
    energyPh = np.loadtxt(phaseFile)


    # find the file index matching the 
    # energyTr and energyPh data
    tmin = int(energyTr[0,0])
    while tmin<int(energyPh[0,0]):
        tmin+=1
    if(tmin>int(energyTr[-1,0])):
        print_mg_error('compute_energyRa')
        print_mg_error('cannot find the file index')
        print_mg_error('to match energyTr and energyPh')
        sys.exit(2)
    tmax = int(energyPh[-1,0])


    # compute the energy ratios
    if(os.path.isfile(outputFile)):
        os.remove(outputFile)

    out = open(outputFile,'w')

    for t in range(tmin,tmax+1):

        if(int(energyTr[t,1])!=int(energyPh[t-tmin,1])):
            print_mg_error('compute_energyRa')
            print_mg_error('mismatch between energyTr and energyPh data')
            sys.exit(2)

        energyTr_r = energyTr[t,2]/energyWa*100
        energyPh_r = energyPh[t-tmin,2]/energyWa*100
        energyIn_r = 100-(energyTr_r+energyPh_r)

        out.write("%i %f %f %f %f\n" % (t,energyTr[t,1],energyTr_r,energyPh_r,energyIn_r))

    out.close()


if __name__=='__main__':


    # extract the options from the command line
    wallFile, transportFile, phaseFile, outputFile = parse_argv(sys.argv[1:])

    
    # compute the energy distribution between
    # transport, phase transition and internal energy
    compute_energyRa(wallFile,
                     transportFile,
                     phaseFile,
                     outputFile)

    

    
