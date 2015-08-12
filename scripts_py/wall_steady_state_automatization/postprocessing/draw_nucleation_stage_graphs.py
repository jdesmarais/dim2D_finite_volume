#!/usr/bin/python

'''
@description
draw graphs characteristic for the study of early stages
of the bubble nucleation

mass   = f(t)
volume = g(t)


initial time when the bubble appears, at different contact angles
initial mass when the bubble appears, at different contact angles
mass = f(t), at different contact angles
'''

import os
import numpy as np
import matplotlib.pyplot as plt
from library_contours_graph import grayscale_to_RGB


if __name__=='__main__':


    # path for the simulation data
    simDir = os.path.join(os.getenv('HOME'),
                          'projects',
                          'dim2d_0.95_ca135.0_vap_fh0.02',
                          'contours')
    
    step = 30
    


    # mass and volume extraction from files
    mass = np.loadtxt(os.path.join(simDir,'mass.txt'))
    volume = np.loadtxt(os.path.join(simDir,'volume.txt'))

    
    # plot the mass and the volume on the same graph
    # with two different axis
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    ax.set_xlabel(r'$t$')
    ax.set_ylabel(r'mass($t$)')

    ax2 = ax.twinx()
    ax2.set_ylabel(r'volume($t$)')

    ax.plot(mass[0:step,1],
            mass[0:step,2],
            'o-',
            color=grayscale_to_RGB(0.2),
            linewidth=3)

    ax2.plot(volume[0:step,1],
             volume[0:step,2],
             's-',
             color=grayscale_to_RGB(0.8),
             linewidth=3)


    # linear interpolation of the mass as
    # function of time
    i1 = [11,17]
    i2 = [16,30]

    # constant initial stage
    xi = mass[i1[0]:i1[1],1]

    av = sum(mass[i1[0]:i1[1],2])/(i1[1]-i1[0])
    av = np.ones(len(xi))*av
    ax.plot(xi,av,'r-',linewidth=3)

    # linear growth rate
    xi = mass[i2[0]:i2[1],1]
    A  = np.array([xi, np.ones(len(xi))])
    yi = mass[i2[0]:i2[1],2]

    pi = np.linalg.lstsq(A.T,yi)[0] #least square approximation
    xi = mass[i2[0]-3:i2[1],1]
    mass_line1 = pi[0]*xi + pi[1]
    ax.plot(xi,mass_line1,'b-',linewidth=3)


    # linear interpolation of the volume as
    # function of time

    # constant initial stage
    xi = volume[i1[0]:i1[1],1]

    av = sum(volume[i1[0]:i1[1],2])/(i1[1]-i1[0])
    av = np.ones(len(xi))*av
    ax2.plot(xi,av,'r-',linewidth=3)

    # linear growth rate
    xi = volume[i2[0]:i2[1],1]
    A  = np.array([xi, np.ones(len(xi))])
    yi = volume[i2[0]:i2[1],2]

    pi = np.linalg.lstsq(A.T,yi)[0] #least square approximation
    xi = volume[i2[0]:i2[1],1]
    volume_line1 = pi[0]*xi + pi[1]
    ax2.plot(xi,volume_line1,'b-',linewidth=3)

    
    ## figure with volume as function of time
    #
    #
    #plt.rc('text', usetex=True)
    #plt.rc('font', family='serif')
    #fig = plt.figure(figsize=(8,6))
    #ax = fig.add_subplot(111)
    #
    #ax.set_xlabel(r'$t$')
    #ax.set_ylabel(r'volume($t$)')
    #
    #plt.plot(volume[0:step,1],
    #         volume[0:step,2],
    #         'r+-',
    #         linewidth=3)
    #
    plt.show()

