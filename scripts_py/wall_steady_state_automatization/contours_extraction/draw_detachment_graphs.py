#!/usr/bin/python

'''
@description
draw graphs characteristic for the study of detachment of
a spherical cap approximation of a vapor bubble by the flow

time at which the vapor bubble is detached from the wall
maximum extent of the bubble when the bubble is detached
mass inside the bubble at the detachment time
'''

import os
import sys
import numpy as np
import matplotlib.pyplot as plt
from library_contours_graph import grayscale_to_RGB
from library_contours_detachment import (find_detachment_time,
                                         find_bubble_extent)
    

def generate_detachment_data(
    simDir,
    legend='None',
    width=3,
    xlabel='',
    ylabel='',
    figPath=''):
    '''
    @description:
    determine the detachment time, draw the contour
    of the vapor bubble at the detachment time and
    compute the maximum extent of the bubble when
    it detaches from the wall
    '''


    #============================================================
    # determine the time at which the bubble detaches from the
    # wall using the contact length as function of time
    #============================================================
    contactLghPath = os.path.join(simDir,'contact_lgh.txt')

    if( not os.path.isfile(contactLghPath)):
        print 'draw_detachment_graphs'
        print 'generate_detachment_data'
        print 'contact length file not found:'
        print contactLghPath
        sys.exit(2)

    else:
        (det,i_det,t_det) = find_detachment_time(contactLghPath)
        
    if(not det):
        print 'no detachment'
        return (0.0,0.0,0.0,0.0)

    #============================================================
    # extract the vapour data in the bubble at the detachment
    #============================================================
    massPath = os.path.join(simDir,'mass.txt')

    if( not os.path.isfile(massPath)):
        print 'draw_detachment_graphs'
        print 'generate_detachment_data'
        print 'mass file not found:'
        print massPath
        sys.exit(2)

    else:
        mass = np.loadtxt(massPath)
        for i in range(0,len(mass[:,0])):
            if(mass[i,0]==i_det):
                mass_det = mass[i,2]
                break


    #============================================================
    # extract the data describing the contour of the bubble
    # at the detachment time
    #============================================================
    contourDataPath = os.path.join(simDir,'contours'+str(i_det)+'.curve')

    if( not os.path.isfile(contourDataPath)):
        print 'draw_detachment_graphs'
        print 'generate_detachment_data'
        print 'contour data file not found:'
        print contourDataPath
        sys.exit(2)

    else:
        contourData = np.loadtxt(contourDataPath)


    #============================================================
    # extract the x-coordinate where the bubble detach
    #============================================================
    y = np.array(contourData[:,1])
    i_min, = np.unravel_index(y.argmin(),y.shape)
    x_det = contourData[i_min,0]

    
    #============================================================
    # extract the maximum extent of the bubble at detachment
    # time using the contour at this time
    #============================================================
    (max_length,segment_pts) = find_bubble_extent(contourData)
    

    #============================================================
    # plot the contour of the bubble at detachment time as well
    # as the segment for the maximum extent of the bubble
    #============================================================
    plt.close("all")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))

    ax = fig.add_subplot(111)

    plt.plot(contourData[:,0],
             contourData[:,1],
             '+-',
             linewidth=width,
             color='black')

    plt.plot(segment_pts[:,0],
             segment_pts[:,1],
             '--',
             linewidth=width,
             color='black')

    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ylabel)

    if(not figPath==''):
        plt.savefig(figPath)

    return (t_det,x_det,max_length,mass_det)


if __name__=='__main__':

    mainDir = os.path.join(os.getenv('HOME'),'projects')

    #=============================================================
    # Detachment study at different contact angle and flow
    # velocities
    #=============================================================

    # directories for the detachment study
    # with different contact angles
    contactAngleArray = [67.5]
    flowVelocityArray = [0.3]

    simDirs = []

    for contactAngle in contactAngleArray:
        for flowVelocity in flowVelocityArray:

            simDir = 'dim2d_0.95_ca'+str(contactAngle)+\
                     '_vap_v'+str(flowVelocity)+'_hca0.0_sph'

            (t_det,x_det,lgh_det,mass_det) = generate_detachment_data(os.path.join(mainDir,simDir,'contours'))

            print '=================================================='
            print 'contact_angle: ', contactAngle
            print 'flow_velocity: ', flowVelocity
            print '=================================================='
            print 'detachment data'
            print 'time   : ', t_det
            print 'x-coord: ', x_det
            print 'max_lgh: ', lgh_det
            print 'mass   : ', mass_det
            print ''
            
    plt.show()

