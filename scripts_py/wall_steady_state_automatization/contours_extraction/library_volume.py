#!/usr/bin/python


import os
import numpy as np

from library_contact_lgh import (find_transition_indices,
                                 interpolate_x_lim)


def integrate_with_trapezoidal_rule(x_data,y_data):
    '''
    description:
    integrate the volume under the curve using trapezoidal
    rule
    '''

    volume = 0
    
    for i in range(0,len(x_data)-1):

        volume+= (y_data[i] + 0.5*(y_data[i+1]-y_data[i]))*(x_data[i+1]-x_data[i])

    return volume


def compute_volume(x_data,y_data,y_lim):
    '''
    @description:
    find the volume corresponding to the contour
    '''

    # find the i corresponding to the data
    # next to the interaction with the wall
    i_transition = find_transition_indices(x_data,y_data,y_lim)

    if(len(i_transition)==0):
        volume = 0

    elif(len(i_transition)%2==0):
        
        volume = 0

        for i in range(0,len(i_transition)/2):

            # transition indices
            i1 = i_transition[0+2*i]
            i2 = i_transition[1+2*i]


            # create the (xData,yData) needed to compute the volume
            n = i2-i1+2

            xData = np.empty([n])
            yData = np.empty([n])
            
            xData[1:n-1] = x_data[i1+1:i2+1]
            yData[1:n-1] = y_data[i1+1:i2+1]
            
            xData[0]   = interpolate_x_lim([x_data[i1],x_data[i1+1]],[y_data[i1],y_data[i1+1]],y_lim)
            yData[0]   = y_lim
            
            xData[n-1] = interpolate_x_lim([x_data[i2],x_data[i2+1]],[y_data[i2],y_data[i2+1]],y_lim)
            yData[n-1] = y_lim
            

            # compute the integral using trapezoidal rule
            volume+= integrate_with_trapezoidal_rule(xData,yData)


            # remove the volume belonging to the wall
            volume-= y_lim*(xData[n-1]-xData[0])
            
    else:
        print 'problem wen computing volume, i_transition: ', i_transition
        volume = -1

    return volume


def curves_to_volume(dirPath,yWall=0.0):
    '''
    @description:
    find the contact length from the contours
    '''
    
    contourRootPath = os.path.join(dirPath,'contours')

    contactLghPath = os.path.join(contourRootPath,'contact_lgh.txt')
    contactLghData = np.loadtxt(contactLghPath)

    nt = len(contactLghData[:,0])

    volumeData = np.empty([nt])

    
    # compute the contact length for all the timesteps
    for i in range(0,nt):

        i_t = contactLghData[i,0]
        t   = contactLghData[i,1]
        
        contourFile = 'contours'+str(int(i_t))+'.curve'

        curvePath = os.path.join(contourRootPath,contourFile)

        if(os.path.isfile(curvePath)):
            curveData = np.loadtxt(curvePath)
            volumeData[i] = compute_volume(curveData[:,0],
                                           curveData[:,1],
                                           yWall)

        else:
            volumeData[i] = 0.0


    # write the volume(t) in an output file
    out = open(os.path.join(contourRootPath,'volume_contact.txt'), 'w')
    for (i,t,l) in zip(contactLghData[:,0],contactLghData[:,1],volumeData):
        out.write("%f %f %f\n" % (i,t,l))
    out.close()

