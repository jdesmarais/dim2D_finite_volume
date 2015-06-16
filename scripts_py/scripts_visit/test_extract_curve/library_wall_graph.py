#!/usr/bin/python

from math import sqrt,pi,sin,acos
import numpy as np
import matplotlib.pyplot as plt


def get_equilibrium_length(Ri,theta):
    '''
    @description: compute the equilibirum length corresponding
    to the contact angle and the initial radius
    '''

    eq_lgh = Ri*sqrt(pi/2.0)*sin(theta)/sqrt(theta-0.5*sin(2*theta))

    return eq_lgh


def get_spherical_cap_data(x,Ri,theta):
    '''
    @description: compute the spherical cap shape corresponding
    to an equilibrium contact angle between vapor and liquid equal
    to theta
    '''
    
    eq_lgh = get_equilibrium_length(Ri,theta)
    eq_R   = Ri*sqrt(pi/2.0)*1.0/sqrt(theta-0.5*sin(2*theta))

    y = np.empty([len(x)])

    for i in range(0,len(x)):
        
        if(x[i]<-eq_lgh or x[i]>eq_lgh):
            y[i]=0
        else:
            y[i] = eq_R*(sin(acos(x[i]/eq_R)) - sin(0.5*pi-theta))

    return y


def grayscale_to_RGB(grayscale_value):
    '''
    @description
    determine the RGB code for a color from
    its grayscale values
    '''
    R = 1.0 - grayscale_value

    return (R,R,R)


def create_graph(data_path,
                 xlabel='$t$',
                 ylabel='',
                 figPath='',
                 width=3,
                 logScale=False,
                 show=True,
                 contactAngle='None'):
    '''
    @description: print a graph of the contact length
    as function of time
    '''

    data = np.loadtxt(data_path)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.figure(figsize=(8,6))

    ax = plt.subplot(111)
    
    plt.plot(
        data[:,0],
        data[:,1],
        '-',
        linewidth=width,
        color='black')

    if(contactAngle!='None'):
        Ri = data[0,1]
        theta = (180-contactAngle)*pi/180

        data[:,1] = get_equilibrium_length(Ri,theta)

        plt.plot(
            data[:,0],
            data[:,1],
            '--',
            linewidth=width,
            color='red')

    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ ylabel)

    if(logScale):
        ax.set_xscale('log')
        ax.set_yscale('log')
    
    if(show):
        plt.show()

    if(not figPath==''):
        plt.savefig(figPath)


def create_st_graph(dataRootPath,
                    timestepExtracted,
                    contactAngle,
                    xlabel='$t$',
                    ylabel='',
                    figPath='',
                    width=3,
                    logScale=False,
                    show=True):

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.figure(figsize=(12,6))

    ax = plt.subplot(111)


    # plot the bubble shape at different times
    for timestep in timestepExtracted:

        dataPath = dataRootPath+'/contours'+str(timestep)+'.curve'
        data = np.loadtxt(dataPath)

        grayscale_value = float(timestep)/float(timestepExtracted[-1])

        plt.plot(
            data[:,0],
            data[:,1],
            '-',
            linewidth=width,
            color=grayscale_to_RGB(grayscale_value))


    # plot the spherical cap if the bubble volume remains constant
    dataPath = dataRootPath+'/volume.txt'
    volume = np.loadtxt(dataPath)
    Ri = sqrt(2.0*volume[-1,1]/pi)
    theta = (180-contactAngle)*pi/180

    sph_data = get_spherical_cap_data(data[:,0],Ri,theta)

    plt.plot(
        data[:,0],
        sph_data,
        '--',
        linewidth=width,
        color='red')

    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ylabel)

    if(show):
        plt.show()

    if(not figPath==''):
        plt.savefig(figPath)


