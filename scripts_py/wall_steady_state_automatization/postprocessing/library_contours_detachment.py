#!/usr/bin/python


import numpy as np
from math import sqrt

def find_detachment_time(contactLghPath):
    '''
    @description: find the time at which the 
    contact length becomes zero in the table of
    contact length as function of time
    '''

    # the file contains the data organized as :
    # fileIndex, time, contact-length
    #
    # the fileIndex allows to associate a 
    # mass contour and the corresponding netcdf
    # file to a specific time
    #
    # e.g. data40.nc <-> contours40.txt
    # 
    contactLgh = np.loadtxt(contactLghPath)

    print contactLghPath

    i_detachment = 0
    t_detachment = 0.
    detachment   = False

    for i in range(0,len(contactLgh[:,0])):

        if(contactLgh[i,2]==0.0):
            i_detachment = int(contactLgh[i,0])
            t_detachment = contactLgh[i,1]
            detachment = True

            if(contactLgh[i,1]==0.0):
                detachment = False
            break

    

    return (detachment, i_detachment,t_detachment)


def find_bubble_extent(contourData):
    '''
    @description: find the couple of points
    in the polyline describing the bubble mass
    contour that will lead to the longest segment
    i.e. what is the length of the largest segment
    that can fit inside the bubble
    '''

    # the data in contourData are organized
    # as:
    # - x : contourData[0,:]
    # - y : contourData[1,:]

    max_length = 0.0

    # the coordinates of the segment points
    # are organized as follow:
    # 1st point: (x,y) = ([0,0],[0,1])
    # 2nd point: (x,y) = ([1,0],[1,1])
    segment_pts = np.empty([2,2])

    # loop over the points in the polyline
    npts = len(contourData[:,0])
    for i in range(0,npts):

        x1 = contourData[i,0]
        y1 = contourData[i,1]

        # create a couple with another point in the polyline
        for j in range(i+1,npts):

            x2 = contourData[j,0]
            y2 = contourData[j,1]

            # compute the distance between the two points
            distance = (x2-x1)**2 + (y2-y1)**2

            # if the distance is larger than the maximum
            # distance save the coordinates of the segment
            # points
            if(distance > max_length):

                segment_pts[0,0] = x1
                segment_pts[1,0] = x2

                segment_pts[0,1] = y1
                segment_pts[1,1] = y2

                max_length = distance

    return (sqrt(max_length),segment_pts)
