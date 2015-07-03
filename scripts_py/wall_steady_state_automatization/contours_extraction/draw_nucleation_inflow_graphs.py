#!/usr/bin/python

from draw_nucleation_graphs import get_nucleation_qties
from library_contours_graph import create_st_graph

import os


def get_simulation_dir(T,ca,sh,v,hca,linear=False):
    '''
    @description: get the name of the directory for the simulation
    '''

    simDir = 'dim2d_'+str(T)
    simDir+= '_ca'+str(ca)+'_vap'
    simDir+= '_sh'+str(sh)
    if(linear):
        simDir+= '_vl'+str(v)
    else:
        simDir+= '_v'+str(v)
    simDir+= '_hca'+str(hca)

    return simDir


if __name__=='__main__':
    
    mainDir = '/home/jdesmarais/projects'

    T   = 0.95
    ca  = 22.5
    sh  = -0.04
    hca = 0.0

    velocity_array = [0.1,0.2,0.3,0.4,0.5]

    times = []
    times.append([0,499])
    times.append([0,499])
    times.append([0,250,500])
    times.append([0,250,416])
    times.append([0,250,316])

    x_limits = [-0.1, 0.225]
    y_limits = [ 0.0, 0.120]

    linear=True

    # get the directories for the simulations
    simDirs = []

    for i in range(0,len(velocity_array)):

        v = velocity_array[i]

        simDir = get_simulation_dir(T,ca,sh,v,hca,linear=linear)
        dirPath = os.path.join(mainDir,
                               simDir,
                               'contours')
        simDirs.append(dirPath)


    # get the nucleation time and mass
    get_nucleation_qties(simDirs)


    # draw the characteristic steps for the nucleation process
    for i in range(0,len(simDirs)):

        simDir = simDirs[i]
        t = times[i]

        create_st_graph(simDir,
                        t,
                        ca,
                        xlabel='$x$',
                        width=3,
                        show=True,
                        x_limits=x_limits,
                        y_limits=y_limits)
