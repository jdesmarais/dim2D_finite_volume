#!/usr/bin/python

'''
@description
extract the bubble contours over time and provide insight on the mass
contained in the vapor phase...
'''

import sys
import os
import getopt

sys.path.append(os.environ['VISIT_PYTHON_LIB'])
import visit

from library_nc_to_vtklines import (generate_vtklines,
                                    generate_time_contour_data)

from library_contours_graph import (create_graph,
                                    create_st_graph)


import numpy as np
import matplotlib.pyplot as plt


def usage():
    '''
    @description:
    describe the usage of the python script
    '''
    
    print ''
    print '-h : display this help'
    print '-i : input folder where the data*.nc files are saved'
    print '-c : contact angle for the spherical cap approximation'
    print '-t : data files analyzed [i_min,i_max,i_step]'
    print '-p : phase checked: needed for bubble nucleation postprocessing'
    print '     otherwise, the contours do not match the interface between'
    print '     the liquid and vapor phases'
    print '-g : generate the contour files, otherwise, it is assumed that'
    print '     the contours files have already been generated'
    print '-s : show the graphs'
    print ''
    print 'ex: ./extract_bubble_contour.py -i <dir> -c <90.0> -t [0,100,10]'
    print ''
    print '    this will create a new directory in <dir> named contours'
    print '    and create files where the contour coordinates are saved'
    print '    as well as contact_length and volume files and pictures'
    print '    of the contact line at different time steps with the'
    print '    spherical cap approximation contour'

    return

def parse_opts(argv):
    '''
    @description:
    parse the options of the script
    '''

    try:
        opts, args = getopt.getopt(argv,
                                   "hi:c:t:pgs",
                                   ["help",
                                    "inputDir",
                                    "contactAngle=",
                                    "timeFrame=",
                                    "phaseCheck",
                                    "genContours",
                                    "show"])
    except getopt.GetoptError:
        usage()
        sys.exit(2)
        
    inputDirProvided      = False
    contactAngleProvided  = False
    timeFrameProvided     = False

    options = {}
    options['phaseCheck']  = False
    options['genContours'] = False
    options['show']        = False
    

    for opt, arg in opts:

        if opt in ("-h","--help"):
            usage()
            sys.exit(2)
            
        elif opt in ("-i", "--inputDir"):
            if(os.path.isdir(arg)):
                inputDirProvided = True
                options['inputDir'] = arg
            else:
                print '*** '+arg+' does not exist ***'
            
        elif opt in ("-c", "--contactAngle"):
            try:
                contactAngle = float(arg)
                contactAngleProvided = True

                options['contactAngle'] = contactAngle

            except ValueError:
                contactAngleProvided = False
                print '*** '+arg+' is not a valid contact angle ***'
            
        elif opt in ("-t", "--timeFrame"):
            try:
                timeFrame = arg.split('[')[1].split(']')[0].split(',')
                for i in range(0,len(timeFrame)):
                    timeFrame[i] = int(timeFrame[i])
                timeFrameProvided = True

                options['timeFrame'] = timeFrame

            except ValueError:
                timeFrameProvided = False
                print '*** '+arg+' is not a valid list time frame [i_min,i_max,i_step] ***'

        elif opt in ("-p","--phaseCheck"):
            options['phaseCheck'] = True

        elif opt in ("-g","--genContours"):
            options['genContours'] = True

        elif opt in ("-s","--Show"):
            options['show'] = True

    if(not(inputDirProvided and contactAngleProvided and timeFrameProvided)):
        print 'the options are not correctly provided'
        usage()
        sys.exit(2)

    print 'input_dir     : ', options['inputDir']
    print 'contact_angle : ', options['contactAngle']
    print 'time_frame    : ', options['timeFrame']
    print 'phase_check   : ', options['phaseCheck']
    print 'gen_contours  : ', options['genContours']
    print 'show          : ', options['show']
    print ''
    
    return options


def generate_st_graphs(ncFolder,
                       timeRange,
                       contactAngle,
                       contourType,
                       phase_check=False,
                       contourPer=0.1,
                       genContours=False,
                       maxNbBubbleContours='None',
                       show=True):
    '''
    @description: generate the contours of the bubble
    at different timesteps, extract the contact length
    of the bubble at the wall, the volume of the bubble
    in time and plot the contours at different timesteps
    as well as the spherical cap approximation    
    '''

    # determine the paths to the folders
    ncRootPath  = os.path.join(ncFolder,'data')
    contoursDir = os.path.join(ncFolder,'contours')
    contoursRootPath = os.path.join(contoursDir,'contours')

    
    # if there is no existing contour folder
    # create one    
    if(not os.path.isdir(contoursDir)):
        os.makedirs(contoursDir)


    ## extract the contact length and the volume
    ## as functions of time
    if(genContours):
        generate_time_contour_data(
            ncRootPath,
            contoursRootPath,
            timeRange=timeRange,
            var='mass',
            contourPer=contourPer,
            contourType=contourType,
            phase_check=phase_check,
            reflection=True)
    

    # paths for saving the contact angle and volume figures
    contact_lgh_path = os.path.join(contoursDir,'contact_lgh.txt')
    volume_path      = os.path.join(contoursDir,'volume.txt')
    dataRootPath     = contoursDir

    contactLghFigPath = os.path.join(contoursDir,'contact_lgh.eps')
    volumeFigPath     = os.path.join(contoursDir,'volume.eps')
    contoursFigPath   = os.path.join(contoursDir,'contours.eps')
    contoursStFigPath = os.path.join(contoursDir,'contours_st.eps')


    # plot the contact length as funtion of time
    create_graph(contact_lgh_path,
                 contactAngle=contactAngle,
                 xlabel='$t$',
                 ylabel='contact length',
                 figPath=contactLghFigPath,
                 width=3,
                 logScale=False,
                 show=show,
                 plot_length_eq=(not phase_check))
    
    # plot the volume as function of time
    create_graph(volume_path,
                 xlabel='$t$',
                 ylabel='volume',
                 figPath=volumeFigPath,
                 width=3,
                 logScale=False,
                 show=show)

    # plot the contour at different time steps
    #------------------------------------------------------------

    # get the first timestep with a bubble
    volumePath = dataRootPath+'/volume.txt'
    volume = np.loadtxt(volumePath)
    for i in range(0,len(volume[:,0])):
        if(volume[i,1]>0):
            start_i = i
            break

    times = []
    if(maxNbBubbleContours=='None'):
        step = timeRange[2]
    else:
        step = int(float(timeRange[1]-start_i)/float(maxNbBubbleContours))
    step = max(1,step)

    for time in range(start_i,timeRange[1],step):
        times.append(time)
    times.append(timeRange[1])

    create_st_graph(dataRootPath,
                    times,
                    contactAngle,
                    xlabel='$x$',
                    ylabel='',
                    figPath=contoursFigPath,
                    width=3,
                    logScale=False,
                    show=show)

    if(not phase_check):
        create_st_graph(dataRootPath,
                        times,
                        contactAngle,
                        xlabel='$x$',
                        ylabel='',
                        figPath=contoursStFigPath,
                        width=3,
                        logScale=False,
                        show=show,
                        sphericalCap=(not phase_check))


def find_initial_bubble(volumePath):
    '''
    @description: find the time at which the bubble appears
    and the initial volume of the bubble
    '''
    
    volume = np.loadtxt(volumePath)
    
    for (t,v) in zip(volume[:,0],volume[:,1]):
        if(v>0):
            t_i = t
            v_i = v
            break
    
    return [t_i,v_i]


if __name__=='__main__':

    #visit.ExportDBAttributes()
    #SuppressMessages(2)
    #print help(visit)
    #GetGlobalAttributes()

    options = parse_opts(sys.argv[1:])
    
    visit.Launch()
    visit.SuppressMessages(0)

    contourType = 'wall_max_gradient' #'mass', 'gradient'

    # generate the contact lengthm, the volume and
    # the bubble contours graphs
    generate_st_graphs(options['inputDir'],
                       options['timeFrame'],
                       options['contactAngle'],
                       contourType,
                       phase_check=options['phaseCheck'],
                       genContours=options['genContours'],
                       contourPer=0.1,
                       maxNbBubbleContours=10,
                       show=options['show'])

    [t_i,r_i] = find_initial_bubble(os.path.join(options['inputDir'],'contours','volume.txt'))

    print 'time bubble   : ', t_i
    print 'initial_volume: ', r_i
