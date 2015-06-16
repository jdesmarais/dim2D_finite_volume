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

from library_wall_graph import (create_graph,
                                create_st_graph)


import numpy as np
import matplotlib.pyplot as plt


def generate_st_graphs(ncFolder,
                       timeRange,
                       contactAngle,
                       contourPer=0.1,
                       contourType):

    # determine the paths to the folders
    ncRootPath  = os.path.join(ncFolder,'data')
    contoursDir = os.path.join(ncRoorPath,'contours')
    contoursRootPath = os.path.join(contoursDir,'contours')

    
    # if there is no existing contour folder
    # create one    
    if(!os.path.isdir(contoursDir)):
        os.makedirs(contoursDir)


    # extract the contact length and the volume
    # as functions of time    
    generate_time_contour_data(
        ncRootPath,
        contoursRootPath,
        timeRange=timeRange,
        var='mass',
        contourPer=contourPer,
        contourType=contourType,
        reflection=True)
    

    # paths for saving the contact angle and volume figures
    contact_lgh_path = os.path.join(contoursDir,'contact_lgh.txt')
    volume_path      = os.path.join(contoursDir,'volume.txt')
    dataRootPath     = contoursDir

    contactLghFigPath = os.path.join(contoursDir,'contact_lgh.eps')
    volumeFigPath     = os.path.join(contoursDir,'volume.eps')
    contoursFigPath   = os.path.join(contoursDir,'contours.eps')


    # plot the contact length as funtion of time
    create_graph(contact_lgh_path,
                 contactAngle=contactAngle,
                 xlabel='$t$',
                 ylabel='contact length($t)',
                 figPath=contactAnglePath,
                 width=3,
                 logScale=False,
                 show=True)

    # plot the volume as function of time
    create_graph(volume_path,
                 xlabel='$t$',
                 ylabel='volume($t)',
                 figPath=volumePath,
                 width=3,
                 logScale=False,
                 show=True)

    # plot the contour at different time steps
    times = []
    for time in range(timeRange[0],timeRange[1],timeRange[2]):
        times.add(time)

    create_st_graph(dataRootPath,
                    times,
                    contactAngle,
                    xlabel='$x$',
                    ylabel='',
                    figPath=contoursPath,
                    width=3,
                    logScale=False,
                    show=True)


if __name__=='__main__':

    #visit.ExportDBAttributes()
    #SuppressMessages(2)
    #print help(visit)
    #GetGlobalAttributes()

    visit.Launch()
    visit.SuppressMessages(0)

    ncFolder = os.path.join(os.getenv('projects'),'dim2d_0.999_ca112.5_vap')
    timeRange = [,,]
    contactAngle = 112.5
    contourType = 'wall_max_gradient' #'mass', 'gradient'
    
    generate_st_graphs(ncFolder,
                       timeRange,
                       contactAngle,
                       contourPer=0.1,
                       contourType)


    #[graph_data,volume] = generate_vtklines(
    #    ncPath,
    #    vtkPath,
    #    var='mass',
    #    contourType='wall_max_gradient',
    #    reflection=reflection)
    #
    #print 'volume: ', volume
    #
    ## plot the curve
    #width = 3
    #
    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')
    #if(reflection):
    #    plt.figure(figsize=(12,6))
    #else:
    #    plt.figure(figsize=(6,6))
    #ax = plt.subplot(111)
    #
    #plt.plot(
    #    graph_data[0],
    #    graph_data[1],
    #    '-',
    #    linewidth=width,
    #    color='black')
    #
    #plt.show()


    # extract the contact length and the volume as functions of time
    ncRootPath  = '/home/jdesmarais/projects/dim2d_0.999_ca112.5_vap/data'
    contourRootPath = './contours/contours'
    #
    #generate_time_contour_data(
    #    ncRootPath,
    #    contourRootPath,
    #    timeRange=[0,1002,10],
    #    var='mass',
    #    contourPer=0.1,
    #    contourType='wall_max_gradient',
    #    reflection=True)
    #
    contact_lgh_path = os.path.dirname(contourRootPath)+'/contact_lgh.txt'
    volume_path      = os.path.dirname(contourRootPath)+'/volume.txt'
    dataRootPath     = os.path.dirname(contourRootPath)
    #
    #
    # plot the contact length as funtion of time
    create_graph(contact_lgh_path,
                 contactAngle=112.5,
                 xlabel='$t$',
                 ylabel='contact length',
                 figPath='',
                 width=3,
                 logScale=True,
                 show=True)
    #
    #
    ## plot the volume as function of time
    #create_graph(volume_path,
    #             xlabel='$t$',
    #             ylabel='volume',
    #             figPath='',
    #             width=3,
    #             logScale=False,
    #             show=True)


    # plot the contour at different time steps
    create_st_graph(dataRootPath,
                    [0,250,500,750,1000],
                    112.5,
                    xlabel='$x$',
                    ylabel='',
                    figPath='',
                    width=3,
                    logScale=False,
                    show=True)
