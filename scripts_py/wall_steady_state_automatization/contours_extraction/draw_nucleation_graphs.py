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


def create_mass_graph(simDirs,
                      legend='None',
                      width=3,
                      xlabel='',
                      ylabel='',
                      show=True,
                      figPath='',
                      add_zoom_above=False):
    '''
    @description: create a graph with the vapor mass
    as function of time for different contact angles
    '''

    plt.close("all")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))


    # plot the mass as function of time
    # as the main plot
    ax = fig.add_subplot(111)

    ini_time = [1.,0.]
    ini_mass = [0.,0.]

    end_i = np.empty([len(simDirs)])

    max_mass = 0.0

    for i in range(0,len(simDirs)):
        
        mass = np.loadtxt(os.path.join(simDirs[i],'mass.txt'))

        max_mass = max(max_mass,max(mass[:,2]))

        # determine the last relevant step: timestep!=0
        start_i=0
        for j in range(0,len(mass[:,0])):
            if(mass[j,2]>0.0):
                start_i=j
                break

        end_i[i]=len(mass[:,0])-1
        for j in range(start_i,len(mass[:,0])):
            if(mass[j,0]==0.0):
                end_i[i]=j-1
                break

        # extract the initial time when
        # the bubble appears
        for j in range(0,len(mass[:,0])):

            if(mass[j,2]>0.0):
                ini_time[0] = min(ini_time[0],mass[j,1])
                ini_time[1] = max(ini_time[1],mass[j,1])
                ini_mass[0] = min(ini_mass[0],mass[j,2])
                ini_mass[1] = max(ini_mass[1],mass[j,2])
                break

        # plot the mass as function
        # of time on the main graph
        grayscale_value = 0.2+ 0.8*float(i)/float(len(simDirs))
        
        plt.plot(mass[0:end_i[i],1],
                 mass[0:end_i[i],2],
                 '-',
                 linewidth=width,
                 color=grayscale_to_RGB(grayscale_value))

    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ylabel)

    if(add_zoom_above):
        ax.set_ylim(0.0,max_mass*(1.0+0.45))

    if(legend!='None'):
        plt.legend(legend,loc='lower right')

    # add a zoom where the initial bubble appears
    # on the lower right part of the plot
    ax_x_lim = ax.get_xlim()
    ax_y_lim = ax.get_ylim()

    ini_time[0] = ini_time[0]*(1.-0.60)
    ini_time[1] = ini_time[1]*(1.+0.20)

    ini_mass[0] = ini_mass[0]*(1.-0.1)
    ini_mass[1] = ini_mass[0]+(ini_time[1]-ini_time[0])*(ax_y_lim[1]-ax_y_lim[0])/(ax_x_lim[1]-ax_x_lim[0])

    width    = 0.5
    height_p = width*((ini_mass[1]-ini_mass[0])/(ax_y_lim[1]-ax_y_lim[0]))/((ini_time[1]-ini_time[0])/(ax_x_lim[1]-ax_x_lim[0]))
    height   = 0.2

    ini_mass[1] = height/height_p*ini_mass[1]

    ax_zoom = plt.axes([.15, .65, width, height]) #, axisbg='y')
    for i in range(0,len(simDirs)):

        mass = np.loadtxt(os.path.join(simDirs[i],'mass.txt'))

        grayscale_value = 0.2+ 0.8*float(i)/float(len(simDirs))

        p = plt.plot(mass[0:end_i[i],1],
                     mass[0:end_i[i],2],
                     '+-',
                     linewidth=5*width,
                     color=grayscale_to_RGB(grayscale_value))
        plt.setp(ax_zoom, xticks=[], yticks=[])

    ax_zoom.set_xlim(ini_time[0],ini_time[1])
    ax_zoom.set_ylim(ini_mass[0],ini_mass[1])

    if(show):
        plt.show()
        plt.close()

    if(not figPath==''):
        plt.savefig(figPath)


if __name__=='__main__':

    mainDir = os.path.join(os.getenv('HOME'),'projects')

    
    ##=============================================================
    ## Nucleation study at different contact angles
    ##=============================================================
    #
    ## directories for the nucleation study with
    ## different contact angles
    #contactAngleArray = [22.5,45.0,67.5,90.0,112.5,135.0]
    #
    #simDirs = []
    #
    #for contactAngle in contactAngleArray:
    #
    #    simDir = 'dim2d_0.95_ca'+str(contactAngle)+'_vap'+'_sh-0.02'
    #    simDirs.append(os.path.join(mainDir,simDir,'contours'))
    #
    #
    ## draw a graph with the mass = f(t) for the different simulations
    #create_mass_graph(simDirs,
    #                  legend=contactAngleArray,
    #                  width=3,
    #                  xlabel='t',
    #                  ylabel='mass(t)',
    #                  show=True)


    #=============================================================
    # Nucleation study at different flux intensities
    #=============================================================

    # directories for the nucleation study with
    # different contact angles
    heatFluxArray    = [-0.02,-0.04,-0.06,-0.08,-0.1]
    heatFluxArrayLeg = [ 0.02, 0.04, 0.06, 0.08, 0.1]

    simDirs = []

    for heatFlux in heatFluxArray:

        simDir = 'dim2d_0.95_ca90.0_vap_sh'+str(heatFlux)
        simDirs.append(os.path.join(mainDir,simDir,'contours'))


    # draw a graph with the mass = f(t) for the different simulations
    create_mass_graph(simDirs,
                      legend=heatFluxArrayLeg,
                      width=3,
                      xlabel='t',
                      ylabel='mass(t)',
                      show=True,
                      add_zoom_above=True)
