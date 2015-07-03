#!/usr/bin/python

'''
@description
draw graphs chcracteristic for the study of bubble nucleation:

initial time when the bubble appears, at different fluxes
initial mass when the bubble appears, at different fluxes
mass = f(t), at different fluxes

initial time when the bubble appears, at different contact angles
initial mass when the bubble appears, at different contact angles
mass = f(t), at different contact angles
'''

import os
import numpy as np
import matplotlib.pyplot as plt
from library_contours_graph import grayscale_to_RGB


def create_contact_lgh_graph(simDirs,
                             legend='None',
                             width=3,
                             xlabel='',
                             ylabel='',
                             show=True,
                             logScale=False,
                             figPath=''):
    '''
    @description: create a graph with the contact length
    as function of time for different contact angles
    '''

    plt.close("all")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))


    # plot the mass as function of time
    # as the main plot
    ax = fig.add_subplot(111)

    for i in range(0,len(simDirs)):
        
        contactLgh = np.loadtxt(os.path.join(simDirs[i],'contact_lgh_n.txt'))

        # plot the mass as function
        # of time on the main graph
        grayscale_value = 0.2+ 0.8*float(i)/float(len(simDirs))
        
        #timeDev = np.empty([len(contactLgh[:,0])-2,2])
        #
        #for j in range(1,len(contactLgh[:,0])-1):
        #
        #    timeDev[j-1,0] = contactLgh[j-1,1]
        #
        #    timeDev[j-1,1] = (contactLgh[j+1,2]-contactLgh[j-1,2])/\
        #                     (contactLgh[j+1,1]-contactLgh[j-1,1])
        #
        #plt.plot(timeDev[:,0],
        #         timeDev[:,1],
        #         '-',
        #         linewidth=width,
        #         color=grayscale_to_RGB(grayscale_value))

        plt.plot(contactLgh[:,1],
                 contactLgh[:,2],
                 '+-',
                 linewidth=width,
                 color=grayscale_to_RGB(grayscale_value))

    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ylabel)

    
    if(logScale):
        #ax.set_xscale('log')
        ax.set_yscale('log')

    if(legend!='None'):
        plt.legend(legend,loc='lower right')

    ax.set_xlim([0,30])

    if(show):
        plt.show()
        plt.close()


    if(not figPath==''):
        plt.savefig(figPath)


if __name__=='__main__':

    mainDir = os.path.join(os.getenv('HOME'),'projects')

    T=0.95


    #=============================================================
    # Steady state study at different contact angles
    #=============================================================

    # directories for the nucleation study with
    # different contact angles
    contactAngleArray = [22.5,45.0,67.5,90.0,112.5,135.0]

    simDirs = []

    for contactAngle in contactAngleArray:

        simDir = 'dim2d_'+str(T)+'_ca'+str(contactAngle)+'_vap'
        simDirs.append(os.path.join(mainDir,simDir,'contours'))


    # draw a graph with the mass = f(t) for the different simulations
    create_contact_lgh_graph(simDirs,
                             legend=contactAngleArray,
                             width=3,
                             xlabel='t',
                             ylabel='r(t)',
                             show=True,
                             logScale=False)
    
