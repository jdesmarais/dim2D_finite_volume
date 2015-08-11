#!/usr/bin/python

'''
@description:
extract the contact length from the contours
'''

import os
import numpy as np


def find_transition_indices(x_data,y_data,y_lim):
    '''
    @description:
    find the indices where the y_data cross the y-lim
    '''

    transition_i = []

    i_prev = -99
    for i in range(0,len(y_data)-1):
        
        if( ((y_data[i]-y_lim)*(y_data[i+1]-y_lim)) <= 0.0 ):

            if(i!=i_prev+1 and i!=i_prev+2):
                transition_i.append(i)
                i_prev = i

    return transition_i


def interpolate_x_lim(x_data,y_data,y_lim):
    '''
    @description:
    it is possible to linearly interpolate
    x_lim such that y_lim = f(x_lim).
    The linear interpolation is done using
    (x[i],y[i]) and (x[i+1],y[i+1])
    x_data[0:1] : x[i] and x[i+1]
    y_data[0:1] : y[i] and y[i+1]
    '''

    x_lim = x_data[0] + (y_lim - y_data[0])/(y_data[1] - y_data[0])*(x_data[1] - x_data[0])
    
    return x_lim


def compute_contact_lgh(x_data,y_data,y_lim):
    '''
    @description:
    find the x-coordinates corresponding to the intersections
    between the contour and the wall and deduce the contact
    length
    '''

    # find the i corresponding to the data
    #  next to the interaction with the wall
    i_transition = find_transition_indices(x_data,y_data,y_lim)
    
    # interpolate the intersections with the wall
    if(len(i_transition)==0):
        contact_lgh = 0.0

    elif(len(i_transition)%2==0):

        contact_lgh = 0

        for i in range(0,len(i_transition)/2):
            i1 = i_transition[0+2*i]
            x1 = interpolate_x_lim(x_data[i1:i1+2],y_data[i1:i1+2],y_lim)
        
            i2 = i_transition[1+2*i]
            x2 = interpolate_x_lim(x_data[i2:i2+2],y_data[i2:i2+2],y_lim)

            contact_lgh+= abs(x2-x1)

    else:
        print 'problem when extracting the contact lgh: ', i_transition
        contact_lgh = -1

    return contact_lgh


def curves_to_contact_lgh(dirPath,yWall=0.0):
    '''
    @description:
    find the contact length from the contours
    '''
    
    contourRootPath = os.path.join(dirPath,'contours')

    volumePath = os.path.join(contourRootPath,'volume.txt')
    volumeData = np.loadtxt(volumePath)

    nt = len(volumeData[:,0])

    contactLghData = np.empty([nt])

    
    # compute the contact length for all the timesteps
    for i in range(0,nt):

        i_t = volumeData[i,0]
        t   = volumeData[i,1]
        
        contourFile = 'contours'+str(int(i_t))+'.curve'

        curvePath = os.path.join(contourRootPath,contourFile)

        if(os.path.isfile(curvePath)):
            curveData = np.loadtxt(curvePath)
            contactLghData[i] = compute_contact_lgh(curveData[:,0],
                                                    curveData[:,1],
                                                    yWall)
        else:
            contactLghData[i] = 0.0


    # write the length(t) in an output file
    out = open(os.path.join(contourRootPath,'volume2.txt'), 'w')
    for (i,t,l) in zip(volumeData[:,0],volumeData[:,1],contactLghData):
        out.write("%f %f %f\n" % (i,t,l))
    out.close()


if __name__=='__main__':
    
    dirPath = os.path.join('/home/jdesmarais/projects','dim2d_0.95_ca22.5_vap')

    curves_to_contact_lgh(dirPath)
