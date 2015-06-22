#!/usr/bin/python

from math import sqrt,pi,sin,cos,acos
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

def extract_domain_borders(dataPath):
    '''
    @description: extract the domain borders from
    the input file
    '''

    data = np.loadtxt(dataPath)

    domain_borders = {}
    domain_borders['x_min'] = data[0]
    domain_borders['x_max'] = data[1]
    domain_borders['y_min'] = data[2]
    domain_borders['y_max'] = data[3]
    domain_borders['dx']    = data[4]
    domain_borders['dy']    = data[5]

    return domain_borders


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
        eq_R   = Ri*sqrt(pi/2.0)*1.0/sqrt(theta-cos(theta)*sin(theta))

        x = np.arange(xmin,xmax,(xmax-xmin)/1000.)
        y = np.empty([len(x)])
        
        for i in range(0,len(x)):
            
            if(x[i]<-eq_lgh or x[i]>eq_lgh):
                y[i] = y_min
            else:
                y[i] = y_min + eq_R*( sin(acos(x[i]/eq_R)) - cos(theta) )
                
    return [x,y]


def grayscale_to_RGB(grayscale_value):
    '''
    @description
    determine the RGB code for a color from
    its grayscale values
    '''
    
    if( (0<=grayscale_value) and (grayscale_value<=1) ):

        R = 1.0 - grayscale_value

    else:

        print 'error for grayscale value'
        print grayscale_value
        sys.exit(2)            

    return (R,R,R)


def create_graph(data_path,
                 xlabel='$t$',
                 ylabel='',
                 figPath='',
                 width=3,
                 logScale=False,
                 show=True,
                 contactAngle='None',
                 plotLengthEq=False,
                 volumePath='None'):
    '''
    @description: print a graph of the contact length
    as function of time
    '''

    data = np.loadtxt(data_path)

    plt.close("all")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))

    ax = fig.add_subplot(111)
    
    plt.plot(
        data[:,1],
        data[:,2],
        '-',
        linewidth=width,
        color='black')

    if(contactAngle!='None' and plotLengthEq and volumePath!='None'):
        volume = np.loadtxt(volumePath)
        Ri = sqrt(2.0*volume[-1,2]/pi)
        theta = (180-contactAngle)*pi/180

        data[:,2] = get_equilibrium_length(Ri,theta)

        plt.plot(
            data[:,1],
            data[:,2],
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
                    show=True,
                    x_limits='None',
                    y_limits='None'):
    '''
    description: create a graph with the bubble contours
    at several timesteps
    '''

    plt.close("all")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(12,6))

    ax = fig.add_subplot(111,aspect='equal')

    y_min = 0.0
    y_max = 0.0

    # plot the bubble shape at different times
    for i in range(0,len(timestepExtracted)):
    
        timestep = timestepExtracted[i]

        dataPath = dataRootPath+'/contours'+str(timestep)+'.curve'

        if(os.path.isfile(dataPath)):

            data = np.loadtxt(dataPath)
        
            if(timestep==timestepExtracted[-1]):
                ratio = 1
            else:
                ratio = (float(timestep)-float(timestepExtracted[0]))/\
                        (float(timestepExtracted[-1])-float(timestepExtracted[0]))

            if(not (0<=ratio<=1)):
                print 'timestep    ', timestep
                print 'timestep[0] ', timestepExtracted[0]
                print 'timestep[-1]', timestepExtracted[-1]

            grayscale_value = 0.1 + 0.9*ratio
        
            # plot the first bubble shape with dashed line
            if(i==0):
                plotstyle = '--'
                linewidth = 2
                color = 'black'
            
            # plot the last bubble with dashed line
            elif(i==len(timestepExtracted)-1):
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
            y_max = max(y_max,max(data[:,1]))

    y_max = 1.05*y_max

    plt.ylim([y_min,y_max])

    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ylabel)


    # set the limits
    if(x_limits!='None'):
        plt.xlim(x_limits)

    if(y_limits!='None'):
        plt.ylim(y_limits)

            
    # show the graph
    if(show):
        plt.show()
        plt.close()

                
    # save in .eps format
    if(not figPath==''):
        plt.savefig(figPath)


def create_sph_graph(dataRootPath,
                     timestepExtracted,
                     contactAngle,
                     xlabel='$x$',
                     ylabel='',
                     figPath='',
                     width=3,
                     show=True,
                     x_limits='None',
                     y_limits='None'):
    '''
    description: create a contour for the last time step
    and add the spherical cap approximation
    '''

    plt.close("all")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(12,6))

    ax = fig.add_subplot(111,aspect='equal')


    # plot the contour for the last timestep and
    # compute the y_min of the contour giving
    # the wall position
    y_min = 0.0
    y_max = 0.0

    dataPath = dataRootPath+'/contours'+str(timestepExtracted)+'.curve'

    if(os.path.isfile(dataPath)):

        data = np.loadtxt(dataPath)
        
        plt.plot(
            data[:,0],
            data[:,1],
            '-',
            linewidth=width,
            color='black')

        y_min = min(y_min,min(data[:,1]))
        y_max = max(y_max,max(data[:,1]))


    # plot the spherical cap
    dataPath = dataRootPath+'/volume.txt'
    volume = np.loadtxt(dataPath)

    theta = (180-contactAngle)*pi/180
    Ri = sqrt(2.0*volume[-1,2]/pi)
        
    dataPath = dataRootPath+'/domain_borders.txt'
    domain_borders = extract_domain_borders(dataPath)    

    x_min = data[ 0,0]*1.5
    x_max = data[-1,0]*1.5
    y_min = domain_borders['y_min']+2.0*domain_borders['dy']
    y_max = y_max*1.1

    #x_min = -0.35
    #x_max =  0.35


    if(contactAngle<90):
        sph_data = get_spherical_cap_data(y_min,x_min,x_max,Ri,theta,filled=True)
        
    else:
        sph_data = get_spherical_cap_data(y_min,x_min,x_max,Ri,theta,filled=False)
        
    plt.plot(
        sph_data[0],
        sph_data[1],
        '--',
        linewidth=width,
        color='red')
    plt.ylim([y_min,y_max])

    if(x_limits!='None'):
        plt.xlim(x_limits)

    if(y_limits!='None'):
        plt.ylim(y_limits)
    
    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ylabel)

            
    # show the graph
    if(show):
        plt.show()
        plt.close()

                
    # save in .eps format
    if(not figPath==''):
        plt.savefig(figPath)

