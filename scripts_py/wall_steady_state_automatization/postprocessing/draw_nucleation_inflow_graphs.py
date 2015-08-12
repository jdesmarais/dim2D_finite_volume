#!/usr/bin/python

from draw_nucleation_graphs import create_mass_graph, get_nucleation_qties
from library_contours_graph import create_st_graph
import os
import numpy as np
import matplotlib.pyplot as plt
from library_contours_graph import grayscale_to_RGB


def create_nucleation_time_graph(dirs,
                                 velocities,
                                 styles,
                                 legend='None',
                                 legendLoc='lower right',
                                 width=3,
                                 xlabel='$v$',
                                 ylabel='$t_n$',
                                 bothDirs=False,
                                 show=True):
    '''
    @description: create a graph with the nucleation times
    for different velocities and fluxes
    
    - dirs   : 2D array containing the directory paths
    - styles : array with the style corresponding to dirs[i]
    - leg    : legend for dirs[i]
    - width  : width of the lines
    - xlabel : xlabel for the plot
    - ylabel : ylabel for the plot
    - show   : show the plot at the end
    '''

    plt.close("all")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    i=0

    for simDirs,simVelocities,simStyle in zip(dirs,velocities,styles):

        # determine the color for the plot
        if(bothDirs):
            grayscale_value = 0.2+ 0.8*float(i%(len(dirs)/2))/float(len(dirs)/2)
        else:
            grayscale_value = 0.2+ 0.8*float(i)/float(len(dirs))
        simColor = grayscale_to_RGB(grayscale_value)        

        # array containing the data at fixed flux
        x_data = []
        y_data = []

        for simDir,simVelocity in zip(simDirs,simVelocities):

            # get the nucleation time and mass
            mass = np.loadtxt(os.path.join(simDir,'contours','mass.txt'))
            for j in range(0,len(mass[:,0])):
                if(mass[j,2]>0):
                    nucleation_time = mass[j,1]
                    nucleation_mass = mass[j,2]
                    break

            # extract the velocity corresponding to the directory
            # it is simVelocity

            # create the corresponding database for plotting
            x_data.append(simVelocity)
            y_data.append(nucleation_time)

        # plot the data
        plt.plot(x_data,
                 y_data,
                 simStyle,
                 linewidth=width,
                 color=simColor)

        i+=1

    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ylabel)

    if(legend!='None'):
        plt.legend(legend,loc=legendLoc)

    if(show):
        plt.show()
        plt.close()


def create_nucleation_time_graph2(dirsPar,
                                  dirsLin,
                                  velocities,
                                  styles,
                                  legend='None',
                                  legendLoc='lower right',
                                  width=3,
                                  xlabel='$v$',
                                  ylabel='$t_n$',
                                  show=True):
    '''
    @description: create a graph with the nucleation times
    for different velocities and fluxes
    
    - dirs   : 2D array containing the directory paths
    - styles : array with the style corresponding to dirs[i]
    - leg    : legend for dirs[i]
    - width  : width of the lines
    - xlabel : xlabel for the plot
    - ylabel : ylabel for the plot
    - show   : show the plot at the end
    '''

    plt.close("all")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)

    i=0

    for simDirsPar,simDirsLin,simVelocities,simStyle in zip(dirsPar,dirsLin,velocities,styles):

        # determine the color for the plot
        grayscale_value = 0.2+ 0.8*float(i)/float(len(dirs))
        simColor = grayscale_to_RGB(grayscale_value)        

        # array containing the data at fixed flux
        x_data = []
        y_data = []

        for simDirPar,simDirLin,simVelocity in zip(simDirsPar,simDirsLin,simVelocities):

            # get the nucleation time and mass: parabolic
            mass = np.loadtxt(os.path.join(simDirPar,'contours','mass.txt'))
            for j in range(0,len(mass[:,0])):
                if(mass[j,2]>0):
                    nucleation_time = mass[j,1]
                    nucleation_mass = mass[j,2]
                    break

            # create the corresponding database for plotting
            x_data_par.append(simVelocity)
            y_data_par.append(nucleation_time)

            
            # get the nucleation time and mass: linear
            mass = np.loadtxt(os.path.join(simDirLin,'contours','mass.txt'))
            for j in range(0,len(mass[:,0])):
                if(mass[j,2]>0):
                    nucleation_time = mass[j,1]
                    nucleation_mass = mass[j,2]
                    break

            # create the corresponding database for plotting
            x_data_par.append(simVelocity)
            y_data_par.append(nucleation_time)

        # plot the data
        plt.plot(x_data_par,
                 y_data_par,
                 simStyle+'-',
                 linewidth=width,
                 color=simColor)

        plt.plot(x_data_lin,
                 y_data_lin,
                 simStyle+'-',
                 linewidth=width,
                 color=simColor)

        i+=1

    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ylabel)

    if(legend!='None'):
        plt.legend(legend,loc=legendLoc)

    if(show):
        plt.show()
        plt.close()


if __name__=='__main__':
    
    drawNucleationInflowPar = True  #v \in \{0.1-0.5/}, fh=0.04 (parabolic)
    drawNucleationInflowLin = True  #v \in \{0.1-0.5/}, fh=0.04 (linear)
    drawNucleationFluxPar   = True  #v \in \{0.2,0.4/}, fh\in\{0.04,0.1/} (parabolic)
    drawNucleationFluxLin   = True  #v \in \{0.2,0.4/}, fh\in\{0.04,0.1/} (linear)

    drawNucleationTimePar = False
    drawNucleationTimeLin = False

    mainDir = os.path.join(os.getenv('HOME'),
                           'projects')    
    
    # directories for the nucleation study with
    # different contact angles
    contactAngle = 22.5
    heatFlux = 0.04
    velocityArray    = [0.0,0.1,0.2,0.3,0.4,0.5]

    velocityLegArrayPar = ['0.0','0.1','0.2','0.3','0.4','0.5']
    velocityLegArrayLin = velocityLegArrayPar
    
    heatFluxArray    = [0.04,0.05,0.06,0.08,0.1]
    heatFluxLegArray = [0.04,0.05,0.06,0.08,0.1]

    simDirsPar = []
    simDirsLin = []
    
    
    # directories for the simulations with the linear
    # or parabolic profiles and nucleation: varying velocity
    if(drawNucleationInflowPar or drawNucleationInflowLin):

        for velocity in velocityArray:
    
            simDir = 'dim2d_0.95_ca'+str(contactAngle)+\
                '_vap'+\
                '_fh'+str(heatFlux)

            simPathPar = simDir
            simPathLin = simDir
            
            if(velocity!=0.0):
                simPathPar+='_v'+str(velocity)            
                simPathLin+='_vl'+str(velocity)

            simPathPar+='_hca0.0'
            simPathLin+='_hca0.0'

            simDirsPar.append(os.path.join(mainDir,simPathPar,'contours'))
            simDirsLin.append(os.path.join(mainDir,simPathLin,'contours'))


    if(drawNucleationInflowPar):
        
        # borders for extracting the mass growth rate
        growthRateBorders = []
        growthRateBorders.append([90 ,400]) #v=0.0
        growthRateBorders.append([100,280]) #v=0.1
        growthRateBorders.append([100,260]) #v=0.2
        growthRateBorders.append([100,240]) #v=0.3
        growthRateBorders.append([100,220]) #v=0.4
        growthRateBorders.append([100,200]) #v=0.5

        # create the graph with mass = f(t)
        create_mass_graph(simDirsPar,
                          legend=velocityLegArrayPar,
                          width=3,
                          xlabel='t',
                          ylabel='mass(t)',
                          show=True,
                          add_zoom_above=True,
                          add_linear_interpolation=True,
                          borders_linear_interpolation=growthRateBorders,
                          legendLoc='upper right',
                          styleDashed=False)

        # extract time+mass nucleation
        get_nucleation_qties(simDirsPar)


    if(drawNucleationInflowLin):

        # borders for extracting the mass growth rate
        growthRateBorders = []
        growthRateBorders.append([90 ,400]) #v=0.0
        growthRateBorders.append([100,280]) #v=0.1
        growthRateBorders.append([110,240]) #v=0.2
        growthRateBorders.append([110,200]) #v=0.3
        growthRateBorders.append([140,200]) #v=0.4
        growthRateBorders.append([255,300]) #v=0.5
        
        # create the graph with mass = f(t)
        create_mass_graph(simDirsLin,
                          legend=velocityLegArrayLin,
                          width=3,
                          xlabel='t',
                          ylabel='mass(t)',
                          show=True,
                          add_zoom_above=True,
                          add_linear_interpolation=True,
                          borders_linear_interpolation=growthRateBorders,
                          legendLoc='upper right',
                          styleDashed=False)

        # extract time+mass nucleation
        get_nucleation_qties(simDirsLin)


    # directories for the simulations with the linear
    # or parabolic profiles and nucleation : varying fluxes 
    simDirsPar = []
    simDirsLin = []

    velocity = 0.4

    if(drawNucleationFluxLin or drawNucleationFluxPar):

        for heatFlux in heatFluxArray:
    
            simDir = 'dim2d_0.95_ca'+str(contactAngle)+\
                '_vap'+\
                '_fh'+str(heatFlux)
            
            simPathPar = simDir+\
                '_v'+str(velocity)+\
                '_hca0.0'

            simPathLin = simDir+\
                '_vl'+str(velocity)+\
                '_hca0.0'

            simDirsPar.append(os.path.join(mainDir,simPathPar,'contours'))
            simDirsLin.append(os.path.join(mainDir,simPathLin,'contours'))

    if(drawNucleationFluxPar):

        # borders for extracting the mass growth rate
        growthRateBorders = []
        growthRateBorders.append([100,180]) #fh=0.04
        growthRateBorders.append([60,190])  #fh=0.05
        growthRateBorders.append([45,160])  #fh=0.06
        growthRateBorders.append([25,140])  #fh=0.08
        growthRateBorders.append([20,45])   #fh=0.1

        # graph with the vapor mass
        create_mass_graph(simDirsPar,
                          legend=heatFluxLegArray,
                          width=3,
                          xlabel='t',
                          ylabel='mass(t)',
                          show=True,
                          add_zoom_above=True,
                          add_linear_interpolation=True,
                          borders_linear_interpolation=growthRateBorders,
                          legendLoc='upper right',
                          styleDashed=False)

        # extract time+mass nucleation
        get_nucleation_qties(simDirsPar)

    if(drawNucleationFluxLin):

        # borders for extracting the mass growth rate
        growthRateBorders = []
        growthRateBorders.append([140,200]) #fh=0.04
        growthRateBorders.append([70,150])  #fh=0.05
        growthRateBorders.append([45,140])  #fh=0.06
        growthRateBorders.append([25,120])  #fh=0.08
        growthRateBorders.append([20,70])   #fh=0.1

        # graph with the vapor mass
        create_mass_graph(simDirsLin,
                          legend=heatFluxLegArray,
                          width=3,
                          xlabel='t',
                          ylabel='mass(t)',
                          show=True,
                          add_zoom_above=True,
                          add_linear_interpolation=True,
                          borders_linear_interpolation=growthRateBorders,
                          legendLoc='upper right',
                          styleDashed=False)

        # extract time+mass nucleation
        get_nucleation_qties(simDirsLin)
        
    if(drawNucleationTimePar or drawNucleationTimeLin):

        dirs = []
        velocities = []
        leg  = [0.04,0.05,0.06,0.08,0.1]
        styles = ['o-','^-','s-','p-','.-',]


        # flux = 0.04
        dirs.append([
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.04_v0.1_hca0.0'),
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.04_v0.2_hca0.0'),
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.04_v0.3_hca0.0'),
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.04_v0.4_hca0.0'),
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.04_v0.5_hca0.0')])

        velocities.append([0.1,0.2,0.3,0.4,0.5])


        # flux = 0.05
        dirs.append([
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.05_v0.2_hca0.0'),
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.05_v0.4_hca0.0')])

        velocities.append([0.2,0.4])


        # flux = 0.06
        dirs.append([
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.06_v0.2_hca0.0'),
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.06_v0.4_hca0.0')])

        velocities.append([0.2,0.4])


        # flux = 0.08
        dirs.append([
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.08_v0.2_hca0.0'),
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.08_v0.4_hca0.0')])

        velocities.append([0.2,0.4])


        # flux = 0.1
        dirs.append([
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.1_v0.2_hca0.0'),
                os.path.join(mainDir,'dim2d_0.95_ca22.5_vap_fh0.1_v0.4_hca0.0')])

        velocities.append([0.2,0.4])


        if(drawNucleationTimePar):

            create_nucleation_time_graph(dirs,
                                         velocities,
                                         styles,
                                         legend=leg)

#        if(drawNucleationTimeLin):
#
#            for i in range(0,len(dirs)):
#
#                dirs.append([])
#                styles.append([])
#
#                styles[len(dirs)+i].append(styles[i].replace('-','--',1))
#
#                for j in range(0,len(dirs[i])):
#
#                    dirs[len(dirs)+i].append(dirs[i][j].replace('_v0','_vl0',1))
#
#            print dirs
#            print styles
#
##            create_nucleation_time_graph(dirs,
##                                         velocities,
##                                         styles,
##                                         legend=leg)
