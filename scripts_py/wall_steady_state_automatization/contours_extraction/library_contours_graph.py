#!/usr/bin/python

from math import sqrt,pi,sin,cos,acos
import numpy as np
import matplotlib.pyplot as plt
import os


def get_equilibrium_length(Ri,theta):
    '''
    @description: compute the equilibirum length corresponding
    to the contact angle and the initial radius
    '''

    eq_lgh = Ri*sqrt(pi/2.0)*sin(theta)/sqrt(theta-0.5*sin(2*theta))

    return eq_lgh


def get_spherical_cap_data(y_min,xmin,xmax,Ri,theta,filled=False):
    '''
    @description: compute the spherical cap shape corresponding
    to an equilibrium contact angle between vapor and liquid equal
    to theta
    '''
    
    if(filled):
        theta1 = pi-theta
        eq_R   = Ri*sqrt(pi/2.0)*1.0/sqrt(pi-theta1+cos(theta1)*sin(theta1))
        eq_lgh = eq_R*sin(theta1)

        x1 = [xmin,-eq_lgh]
        y1 = [y_min,y_min]

        x2 = np.arange(-eq_lgh,-eq_R,(-eq_R+eq_lgh)/300.)
        y2 = np.empty([len(x2)])
        
        for i in range(0,len(x2)):
            y2[i] = y_min+eq_R*cos(theta1)-eq_R*sin(acos(x2[i]/eq_R))

        x3 = np.arange(-eq_R,eq_R,2*eq_R/1000.)
        y3 = np.empty([len(x3)])
        
        for i in range(0,len(x3)):
            y3[i] = y_min+eq_R*cos(theta1)+eq_R*sin(acos(x3[i]/eq_R))

        x4 = np.arange(eq_R,eq_lgh,(-eq_R+eq_lgh)/300.)
        y4 = np.empty([len(x4)])
        
        for i in range(0,len(x4)):
            y4[i] = y_min+eq_R*cos(theta1)-eq_R*sin(acos(x4[i]/eq_R))

        x5 = [eq_lgh,xmax]
        y5 = [y_min,y_min]

        x1_i = len(x1)
        x2_i = x1_i+len(x2)
        x3_i = x2_i+len(x3)
        x4_i = x3_i+len(x4)
        x5_i = x4_i+len(x5)

        x = np.empty([x5_i])
        y = np.empty([len(x)])

        x[   0:x1_i] = x1
        y[   0:x1_i] = y1
        x[x1_i:x2_i] = x2
        y[x1_i:x2_i] = y2
        x[x2_i:x3_i] = x3
        y[x2_i:x3_i] = y3
        x[x3_i:x4_i] = x4
        y[x3_i:x4_i] = y4
        x[x4_i:x5_i] = x5
        y[x4_i:x5_i] = y5

    else:
        eq_lgh = get_equilibrium_length(Ri,theta)
        eq_R   = Ri*sqrt(pi/2.0)*1.0/sqrt(theta-0.5*sin(2*theta))

        x = np.arange(xmin,xmax,(xmax-xmin)/1000.)
        y = np.empty([len(x)])
        
        for i in range(0,len(x)):
            
            if(x[i]<-eq_lgh or x[i]>eq_lgh):
                y[i] = y_min
            else:
                y[i] = y_min + eq_R*(sin(acos(x[i]/eq_R)) - sin(0.5*pi-theta))
                
    return [x,y]


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
                 contactAngle='None',
                 plot_length_eq=False):
    '''
    @description: print a graph of the contact length
    as function of time
    '''

    data = np.loadtxt(data_path)

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))

    ax = fig.add_subplot(111)
    
    plt.plot(
        data[:,0],
        data[:,1],
        '-',
        linewidth=width,
        color='black')

    if(contactAngle!='None' and plot_length_eq):
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
        plt.close()

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
                    show=True,
                    sphericalCap=False):
    '''
    description: create a graph with the bubble contours
    at several timesteps
    '''

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(12,6))

    ax = fig.add_subplot(111,aspect='equal')


    y_min = 0

    # plot the bubble shape at different times
    for i in range(0,len(timestepExtracted)):
    
        timestep = timestepExtracted[i]

        dataPath = dataRootPath+'/contours'+str(timestep)+'.curve'

        if(os.path.isfile(dataPath)):

            data = np.loadtxt(dataPath)
        
            grayscale_value = 0.1+ 0.9*float(timestep)/float(timestepExtracted[-1])
        
            # plot the first bubble shape with dashed line
            if(timestep==timestepExtracted[0]):
                plotstyle = '--'
                linewidth = 2
                color = 'black'
            
            # plot the last bubble with dashed line
            elif(timestep==timestepExtracted[-1]):
                plotstyle = '--'
                linewidth = 2
                color = 'black'

            # plot the other bubbles with continuous line
            else:
                plotstyle = '-'
                linewidth = width
                color = grayscale_to_RGB(grayscale_value)
            
            plt.plot(
                data[:,0],
                data[:,1],
                plotstyle,
                linewidth=linewidth,
                color=color)

            y_min = min(y_min,min(data[:,1]))


    # plot the spherical cap
    if(sphericalCap):
        dataPath = dataRootPath+'/volume.txt'
        volume = np.loadtxt(dataPath)
        theta = (180-contactAngle)*pi/180
        Ri = sqrt(2.0*volume[-1,1]/pi)
        
        if(contactAngle<90):
            sph_data = get_spherical_cap_data(y_min,data[0,0],data[-1,0],Ri,theta,filled=True)
            
        else:
            sph_data = get_spherical_cap_data(y_min,data[0,0],data[-1,0],Ri,theta,filled=False)
            
        plt.plot(
            sph_data[0],
            sph_data[1],
            '--',
            linewidth=width,
            color='red')
    
    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ylabel)

            
    # show the graph
    if(show):
        plt.show()
        plt.close()

                
    # save in .eps format
    if(not figPath==''):
        plt.savefig(figPath)




