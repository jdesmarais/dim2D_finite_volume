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

from library_contact_lgh import curves_to_contact_lgh

from library_volume import curves_to_volume

from library_contours_graph import (create_graph,
                                    create_st_graph,
                                    create_sph_graph)


import numpy as np
import matplotlib.pyplot as plt
from math import cos, sin, sqrt, pi


def iround(x):
    """iround(number) -> integer
    Round a number to the nearest integer."""
    return int(round(x) - .5) + (x > 0)


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
    print '-x : set the limits for the x-axis when plotting the graphs'
    print '-y : set the limits for the y-axis when plotting the graphs'
    print '-w : run the visit engine without window'
    print '-l : contourType used to determine the mass contours:'
    print ''
    print '        - wall_max_gradient: determine the location of the maximum'
    print '                             gradient at the wall and use the mass'
    print '                             density at this location'
    print '        - max_gradient     : determine the location of the maximum'
    print '                             gradient and use the mass density at'
    print '                             this location'
    print '        - mass             : use the mid-mass between the liquid'
    print '                             and vapor satured phases at the'
    print '                             simulation temperature'
    print '-r : reflection activated '
    print '-e : do not select the final time'
    print '-a : scaling for the time steps: t**a'
#    print '-v : re-compute the volumes'
    print '        '
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
                                   "hi:c:t:pgx:y:wsl:rea:",
                                   ["help",
                                    "inputDir",
                                    "contactAngle=",
                                    "timeFrame=",
                                    "phaseCheck",
                                    "genContours",
                                    "show",
                                    "no_window",
                                    "reflection",
                                    "x_figsize=",
                                    "y_figsize=",
                                    "maxNbBubbleContours=",
                                    "select_detach_time"])
    except getopt.GetoptError:
        print argv
        usage()
        sys.exit(2)
        
    print argv

    inputDirProvided      = False
    contactAngleProvided  = False
    timeFrameProvided     = False

    options = {}
    options['phaseCheck']          = False
    options['genContours']         = False
    options['show']                = False
    options['x_limits']            = 'None'
    options['y_limits']            = 'None'
    options['x_figsize']           = 12
    options['y_figsize']           = 6
    options['no_window']           = False
    options['contourType']         = 'wall_max_gradient' #'mass', 'gradient'
    options['reflection']          = False
    options['select_end_time']     = True
    options['select_detach_time']  = False
    options['scaling_timesteps']   = 'None'
    options['maxNbBubbleContours'] = 7


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

        elif opt in ("-s","--show"):
            options['show'] = True

        elif opt in ("-x"):
            try:
                x_limits = arg.split('[')[1].split(']')[0].split(',')
                for i in range(0,len(x_limits)):
                    x_limits[i] = float(x_limits[i])

                options['x_limits'] = x_limits

            except ValueError:
                print '*** '+arg+' is not a valid list x-limits [x_min,x_max] ***'

        elif opt in ("-y"):
            try:
                y_limits = arg.split('[')[1].split(']')[0].split(',')
                for i in range(0,len(y_limits)):
                    y_limits[i] = float(y_limits[i])

                options['y_limits'] = y_limits

            except ValueError:
                print '*** '+arg+' is not a valid list x-limits [x_min,x_max] ***'

        elif opt in ("-w","--no_window"):
            options['no_window'] = True

        elif opt in ("-l"):
            
            if arg in ['wall_max_gradient','max_gradient','mass']:
                options["contourType"] = arg

            else:
                print 'contourType not recognized'
                usage()
                sys.exit(2)

        elif opt in ("-r"):

            options['reflection'] = True

        elif opt in ("-e"):
            
            options['select_end_time'] = False

        elif opt in ("--select_detach_time"):
            
            options['select_detach_time'] = True

        elif opt in ("-a"):

            options['scaling_timesteps'] = int(arg)

        elif opt in ("--x_figsize"):
            options['x_figsize'] = int(arg)

        elif opt in ("--y_figsize"):
            options['y_figsize'] = int(arg)

        elif opt in ("--maxNbBubbleContours"):
            options['maxNbBubbleContours'] = int(arg)


    if(not(inputDirProvided and contactAngleProvided and timeFrameProvided)):
        print 'the options are not correctly provided'
        usage()
        sys.exit(2)

    print 'input_dir           : ', options['inputDir']
    print 'contact_angle       : ', options['contactAngle']
    print 'time_frame          : ', options['timeFrame']
    print 'phase_check         : ', options['phaseCheck']
    print 'gen_contours        : ', options['genContours']
    print 'show                : ', options['show']
    print 'x_limits            : ', options['x_limits']
    print 'y_limits            : ', options['y_limits']
    print 'x_figsize           : ', options['x_figsize']
    print 'y_figsize           : ', options['y_figsize']
    print 'maxNbBubbleContours : ', options['maxNbBubbleContours']
    print 'no_window           : ', options['no_window']
    print 'select_end_time     : ', options['select_end_time']
    print 'select_detach_time  : ', options['select_detach_time']
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
                       show=True,
                       x_limits='None',
                       y_limits='None',
                       x_figsize=12,
                       y_figsize=6,
                       reflection=False,
                       select_end_time=True,
                       select_detach_time=False,
                       scalingTimesteps='None'):
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

    
    # choose whether to create the graph with the spherical cap
    # approximation
    add_spherical_cap_approx = not phase_check

    
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
            reflection=reflection)
    

    # paths for saving the contact angle and volume figures
    contact_lgh_path = os.path.join(contoursDir,'contact_lgh.txt')
    volume_path      = os.path.join(contoursDir,'volume.txt')
    mass_path        = os.path.join(contoursDir,'mass.txt')
    contour_path     = os.path.join(contoursDir,'contour.txt')
    temperature_path = os.path.join(contoursDir,'temperature.txt')
    dataRootPath     = contoursDir

    contactLghFigPath  = os.path.join(contoursDir,'contact_lgh.eps')
    volumeFigPath      = os.path.join(contoursDir,'volume.eps')
    massFigPath        = os.path.join(contoursDir,'mass.eps')
    contourFigPath     = os.path.join(contoursDir,'mass_contour.eps')
    temperatureFigPath = os.path.join(contoursDir,'temperature.eps')

    contoursFigPath    = os.path.join(contoursDir,'contours.eps')
    contoursStFigPath  = os.path.join(contoursDir,'contours_st.eps')

    contoursFigPath = os.path.join(
        '/home/jdesmarais/Code/augeanstables/scripts_py/wall_steady_state_automatization/postprocessing',
        'figs',
        'fig'+os.path.basename(ncFolder).replace('dim2d','').replace('hca0.0_','')+'_approx.eps')
    
    draw_other_than_contour = False


    # plot the contact length as funtion of time
    if(draw_other_than_contour):

        # plot the contact length
        curves_to_contact_lgh(ncFolder)
        curves_to_volume(ncFolder)

        create_graph(contact_lgh_path,
                     contactAngle=contactAngle,
                     xlabel='$t$',
                     ylabel='contact length',
                     figPath=contactLghFigPath,
                     width=3,
                     logScale=False,
                     show=show,
                     plotLengthEq=add_spherical_cap_approx,
                     volumePath=volume_path)
    
        # plot the volume as function of time
        create_graph(volume_path,
                     xlabel='$t$',
                     ylabel='volume',
                     figPath=volumeFigPath,
                     width=3,
                     logScale=False,
                     show=False)

        # plot the mass as function of time
        create_graph(mass_path,
                     xlabel='$t$',
                     ylabel='mass',
                     figPath=massFigPath,
                     width=3,
                     logScale=False,
                     show=False)

        # plot the mass chosen to draw the
        # contours as function of time
        create_graph(contour_path,
                     xlabel='$t$',
                     ylabel='mass',
                     figPath=contourFigPath,
                     width=3,
                     logScale=False,
                     show=show)

        # plot the temperature chosen to draw the
        # contours as function of time
        if(os.path.isfile(temperature_path)):
            create_graph(temperature_path,
                         xlabel='$t$',
                         ylabel='$T$',
                         figPath=temperatureFigPath,
                         width=3,
                         logScale=False,
                         show=show)


    # plot the contour at different time steps:
    # choose the timesteps to have only maxNbBubbleContours
    #------------------------------------------------------------

    # get the first timestep with a bubble
    volumePath = dataRootPath+'/volume.txt'
    volume = np.loadtxt(volumePath)
    for i in range(0,len(volume[:,0])):
        if(volume[i,2]>0):
            start_i = i
            break
    start_i = max(start_i,timeRange[0])


    # select the last timestep before the volume is zero
    end_i = len(volume[:,0])-1
    if(select_end_time):
        for i in range(len(volume[:,0])-1,start_i,-1):
            if(volume[i,2]>0):
                end_i = volume[i,0]
                break
    end_i = min(end_i,timeRange[1])
    end_i_time = end_i
    
    nt = len(volume[:,0])


    # select the last timestep before the contact length
    # is zero (detachment)
    if(select_detach_time):
        contact_lgh = np.loadtxt(contact_lgh_path)
        end_i = len(contact_lgh[:,0])-1

        for i in range(len(contact_lgh[:,0])-1,start_i,-1):
            if(contact_lgh[i,2]>0):
                end_i = contact_lgh[i,0]
                break

        end_i*=1.3
        end_i = int(end_i)

        end_i = min(end_i,end_i_time)

    print 'end_i: ', end_i


    # select the timesteps
    if(maxNbBubbleContours=='None'):
        nbContours = end_i-start_i+1
    else:
        nbContours = maxNbBubbleContours

    if(scalingTimesteps=='None'):
        scaling=1.0
    else:
        scaling=scalingTimesteps


    # extraction of the timesteps to
    # display the contours
    times   = []
    times_t = []

    #times.append(start_i)
    #times_t.append(volume[start_i,1])

    for i in range(0,nbContours):

        step = float(end_i-start_i)*(float(i)/float(nbContours-1))**scaling
        if(step>0):
            step = max(1,iround(step))
            step-= step%2
        else:
            step = 0

        timestep = start_i + step

        times.append(timestep)
        times_t.append(volume[timestep/timeRange[2],1])

    #times = [0,82,166,248,300,400,498] #to see the break-up for 135.0, ux=0.4

    times_p = np.array(times_t)
    np.set_printoptions(precision=5)
    print 'Timesteps for contours: '
    print times
    print 'Time extracted for contours: '
    print times_p

    # create the graph with only the contours at different
    # relevant times
    create_st_graph(dataRootPath,
                    times,
                    contactAngle,
                    xlabel='$x$',
                    ylabel='',
                    figPath=contoursFigPath,
                    width=3,
                    show=show,
                    x_limits=x_limits,
                    y_limits=y_limits,
                    x_figsize=x_figsize,
                    y_figsize=y_figsize)

    # create the graph with only the last contours and the
    # spherical cap approximation
    if(add_spherical_cap_approx):
        create_sph_graph(dataRootPath,
                         times[-1],
                         contactAngle,
                         figPath=contoursStFigPath,
                         show=show,
                         x_limits=x_limits,
                         y_limits=y_limits)


def find_initial_bubble(volumePath):
    '''
    @description: find the time at which the bubble appears
    and the initial volume of the bubble
    '''
    
    volume = np.loadtxt(volumePath)
    
    for (t,v) in zip(volume[:,1],volume[:,2]):
        if(v>0):
            t_i = t
            v_i = v
            break
    
    return [t_i,v_i]


def compute_contact_length_variation(contactLghPath,
                                     volumePath,
                                     contact_angle):
    '''
    @description: find the normalized difference between the
    contact length and the equilibrium contact length computed
    from the volume and the contact angle as function of time
    '''

    print 'contact_angle: ', contact_angle

    theta  = pi*float(180-contact_angle)/180.0
    theta1 = pi-theta
        
    contactLgh = np.loadtxt(contactLghPath)
    volume     = np.loadtxt(volumePath)

    nt = len(contactLgh[:,0])

    contactLghDiff = np.empty([nt])

    for i in range(0,nt):

        volumet = volume[i,2]

        eq_R   = sqrt(volumet)*1.0/sqrt(pi-theta1+cos(theta1)*sin(theta1))
        eq_lgh = 2.0*eq_R*sin(theta1)

        contactLghEq = eq_lgh #2.0*sqrt(volumet/((pi-theta)+cos(theta)*sin(theta)))*sin(theta)

        contactLghDiff[i] = (contactLghEq-contactLgh[i,2])/contactLghEq

        #print 'eq_lgh: ', contactLghEq, contactLgh[i,2], contactLghDiff[i]

    out = open(os.path.join(os.path.dirname(volumePath),'contact_lgh_n.txt'), 'w')
    for (i,t,l) in zip(volume[:,0],volume[:,1],contactLghDiff):
        out.write("%f %f %f\n" % (i,t,l))
    out.close()


def generate_contact_length_variation(contoursDir,
                                      contact_angle,
                                      show=False):
    '''
    @description: compute the normalized difference of
    the contact length as function of time and create
    a graph
    '''

    contactLghPath = os.path.join(contoursDir,'contact_lgh.txt')
    volumePath     = os.path.join(contoursDir,'volume.txt')

    compute_contact_length_variation(contactLghPath,
                                     volumePath,
                                     contact_angle)

    contoursNPath    = os.path.join(contoursDir,'contact_lgh_n.txt')
    contoursNFigPath = os.path.join(contoursDir,'contact_lgh_n.eps')

    # extract the contact length:
    contactLgh = np.loadtxt(contoursNPath)
    print 'contact_lgh: ', contactLgh[-1]

    create_graph(contoursNPath,
                 xlabel='$t$',
                 ylabel='contact length difference',
                 figPath=contoursNFigPath,
                 width=3,
                 logScale=False,
                 show=show)


if __name__=='__main__':

    #visit.ExportDBAttributes()
    #SuppressMessages(2)
    #print help(visit)
    #GetGlobalAttributes()

    options = parse_opts(sys.argv[1:])
    
    if(options['genContours']):
        if(options['no_window']):
            visit.LaunchNowin()
        else:
            visit.Launch()
        visit.SuppressMessages(1)

    contoursPath = os.path.join(options['inputDir'],'contours')


    ## generate the contact lengthm, the volume and
    ## the bubble contours graphs
    generate_st_graphs(options['inputDir'],
                       options['timeFrame'],
                       options['contactAngle'],
                       options['contourType'],
                       phase_check=options['phaseCheck'],
                       genContours=options['genContours'],
                       contourPer=0.1,
                       maxNbBubbleContours=options['maxNbBubbleContours'],
                       show=options['show'],
                       x_limits=options['x_limits'],
                       y_limits=options['y_limits'],
                       x_figsize=options['x_figsize'],
                       y_figsize=options['y_figsize'],
                       reflection=options['reflection'],
                       select_end_time=options['select_end_time'],
                       select_detach_time=options['select_detach_time'],
                       scalingTimesteps=options['scaling_timesteps'])

    [t_i,r_i] = find_initial_bubble(os.path.join(contoursPath,'volume.txt'))

    print 'time bubble   : ', t_i
    print 'initial_volume: ', r_i

    options['contactLghVar'] = not options['phaseCheck']

    if(options['contactLghVar']):

        generate_contact_length_variation(
            contoursPath,
            options['contactAngle'],
            show=options['show'])
