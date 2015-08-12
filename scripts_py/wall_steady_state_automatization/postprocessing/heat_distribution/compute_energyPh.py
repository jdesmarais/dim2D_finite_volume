#!/usr/bin/python

'''
@description: compute the energy used for the phase transition
of the vapor mass in the computational domain
'''

import sys
import os
import inspect


# python path
#------------------------------------------------------------
cmd_subfolder = os.path.realpath(os.path.abspath(os.path.join(\
    os.path.split(
    inspect.getfile( inspect.currentframe() ))[0],"../")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)


# library imports
#------------------------------------------------------------
import getopt
import numpy as np

from library_messages import\
    print_mg_error,\
    print_mg_progress,\
    print_mg_final

from library_van_der_waals import \
    compute_latentHeat


# functions
#------------------------------------------------------------

# display the help for the program
def display_help():
    '''
    @description:
    display program help
    '''
    print ''
    print 'extract the energy used by the phase transition'
    print 'from files giving the temperature at the interface'
    print 'and the total vapor mass in the system as functions'
    print 'of time'
    print ''
    print 'options:'
    print '--------'
    print '-h (--help)          : display this help'
    print '-t (--temperature=)  : file with the interface temperature as a function of time'
    print '-m (--mass=)         : file with the total vapor mass as a function of time'
    print '-o (--output=)       : text file for the output'
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
                                   "hi:t:m:o:",
                                   ["help",
                                    "temperature=",
                                    "mass=",
                                    "output="])
    except getopt.GetoptError:
        display_help()
        sys.exit(2)

    if(len(opts)==0):
        display_help()
        sys.exit(2)


    temperatureFile = 'None'
    massFile = 'None'
    outputFile = 'None'

    # options
    for opt, arg in opts:

        if opt in ("-h", "--help"):
            display_help()
            sys.exit(2)

        elif opt in ("-t", "--temperature"):
            temperatureFile = arg

        elif opt in ("-m", "--mass"):
            massFile = arg

        elif opt in ("-o", "--output"):
            outputFile = arg


    # check the temperature file
    check_inputFile(temperatureFile,errorMg='temperature file not provided')


    # check the mass file
    check_inputFile(massFile,errorMg='mass file not provided')


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
    
    return temperatureFile, massFile, outputFile


# check whether the file is provided and if it exists
def check_inputFile(filePath,errorMg='file not provided'):
    '''
    @description: check whether the filePath is provided
    and if it exists
    '''
    
    if(filePath=='None'):
        print_mg_error('file not provided')
        display_help()
        sys.exit(2)

    else:
        if( not os.path.isfile(filePath)):
            print_mg_error(filePath+' does not exist')
            display_help()
            sys.exit(2)


# compute the mass rate as a function of time
def compute_massRate(mass):
    '''
    @description: compute the mass rate as a function
    of time
    '''

    # total number of timesteps
    nt = len(mass[:,0])-2


    # initialization of the array containing
    # the time derivative of the mass
    massRate = np.empty([nt,3])

    
    # computation of the time derivative
    # of the vapor mass in the domain
    for i in range(0,nt):
        
        massRate[i,0] = mass[i+1,0]
        massRate[i,1] = mass[i+1,1]
        massRate[i,2] = (mass[i+2,2]-mass[i,2])/(mass[i+2,1]-mass[i,1])


    return massRate


# compute the latent heat rate for each time
def compute_latentHeatRate(massRate,temperature,outputFile):
    '''
    @description: compute the heat power for the
    phase transition for each time
    '''   

    # find the first fileId matching the
    # temperature and the massRate data
    tmin = int(temperature[0,0])
    while tmin<massRate[0,0]:
        tmin+=1

    if(tmin>temperature[-1,0]):
        print_mg_error('compute_latentHeatRate')
        print_mg_error('matching fileId not found')
        sys.exit(2)

    tmax = int(massRate[-1,0]+1)


    # compute the latent heat rate and
    # save it in the output file
    out = open(outputFile, 'w')

    for t in range(tmin,tmax):
        
        if((temperature[t,0]-massRate[t-tmin,0])<1.0):

            latentHeatRate = massRate[t-tmin,2]*\
                compute_latentHeat(temperature[t,2])

            out.write("%i %f %f\n" % (t,temperature[t,1],latentHeatRate))

        else:
            print_mg_error('compute_latentHeatRate')
            print_mg_error('fileID do not match')
            print_mg_error('temperatureID: '+str(temperature[t,0]))
            print_mg_error('massRateID: '+str(massRate[t-tmin,0]))
            out.close()
            sys.exit(2)

    out.close()

    return


# compute the latent heat rate for each time
def compute_energyPh(temperatureFile,massFile,outputFile):
    '''
    @description: compute the latent heat rate for each
    time
    '''

    # extract the data for the temperature as a function of time
    # temperature[:][0] : file id
    # temperature[:][1] : time
    # temperature[:][2] : temperature
    temperature = np.loadtxt(temperatureFile)

    # extract the data for the vapor mass as a function of time
    # mass[:][0] : file id
    # mass[:][1] : time
    # mass[:][2] : total vapor mass in the domain
    mass = np.loadtxt(massFile)

    
    # compute the mass rate as a function of time
    massRate = compute_massRate(mass)


    # compute the corresponding energy used for
    # phase transition as a function of time
    compute_latentHeatRate(massRate,temperature,outputFile)

    
    print_mg_final('extract_energyPh: done')


if __name__=='__main__':

    # extract the options from the command line
    temperatureFile, massFile, outputFile = parse_argv(sys.argv[1:])


    # compute the latent heat rate for each time
    compute_energyPh(temperatureFile,massFile,outputFile)    
