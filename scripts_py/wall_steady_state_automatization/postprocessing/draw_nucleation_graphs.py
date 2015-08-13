#!/usr/bin/python

'''
@description
draw graphs characteristic for the study of bubble nucleation:

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
from library_colors import grayscale_to_RGB
from math import sqrt

from library_van_der_waals import compute_latentHeat


def get_nucleation_qties(simDirs):

    l = len(os.path.basename(os.path.dirname(simDirs[0])))
    for i in range(1,len(simDirs)):
        l1 = len(os.path.basename(os.path.dirname(simDirs[i])))
        if(l1>l):
            l=l1

    format  = "%"+str(l)+"s | %11s | %20s"
    format2 = "%"+str(l)+"s | %4.2f - %4.2f | %3.2e - %3.2e"

    print ''
    print(format % ('sim_dir', 'time', 'mass'))
    print(format % ('----------------', '-----------', '--------------------'))

    time = []

    for i in range(0,len(simDirs)):
        
        mass = np.loadtxt(os.path.join(simDirs[i],'mass.txt'))

        for j in range(0,len(mass[:,0])):
            if(mass[j,2]>0):
                dt = 0.5*(mass[j+1,1]-mass[j,1])
                dm = 0.5*(mass[j+1,2]-mass[j,2])
                print(format2 % (os.path.basename(os.path.dirname(simDirs[i])), mass[j,1], dt, mass[j,2], dm))

                time.append(mass[j,1])
                break

    print ''
    print ''

    return time,


def create_mass_graph(simDirs,
                      legend='None',
                      width=3,
                      xlabel='',
                      ylabel='',
                      show=True,
                      figPath='',
                      add_zoom_above=False,
                      add_linear_interpolation=True,
                      borders_linear_interpolation='None',
                      step=1,
                      legendLoc='lower right',
                      styleDashed=False,
                      xlim='None',
                      ylim='None'):
    '''
    @description: create a graph with the vapor mass
    as function of time for different contact angles
    '''

    growthRate='None'


    plt.close("all")

    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))


    # plot the mass as function of time
    # as the main plot
    ax = fig.add_subplot(111)

    ini_time = [100.,0.]
    ini_mass = [0.,0.]

    start_i = np.empty([len(simDirs)])
    end_i   = np.empty([len(simDirs)])

    max_mass = 0.0

    for i in range(0,len(simDirs)):
        
        mass = np.loadtxt(os.path.join(simDirs[i],'mass.txt'))

        max_mass = max(max_mass,max(mass[:,2]))

        # determine the last relevant step: timestep!=0
        start_i[i]=0
        for j in range(0,len(mass[:,0])):
            if(mass[j,2]>0.0):
                start_i[i]=j
                break

        end_i[i]=len(mass[:,0])-1
        for j in range(int(start_i[i]),len(mass[:,0])):
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
        if(styleDashed):
            grayscale_value = 0.2+ 0.8*float((i-i%2)/2.0)/float(len(simDirs)/2.0)
        else:
            grayscale_value = 0.2+ 0.8*float(i)/float(len(simDirs))
        
        if(styleDashed):
            if(i%2==0):
                style='-'
            else:
                style='--'
        else:
            style='-'

        plt.plot(mass[0:end_i[i]:step,1],
                 mass[0:end_i[i]:step,2],
                 style,
                 linewidth=width,
                 color=grayscale_to_RGB(grayscale_value))

    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ylabel)


    # extract the linear growth rate
    if(add_linear_interpolation):

        growthRate = []

        if(borders_linear_interpolation!='None'):
        
            for i in range(0,len(simDirs)):

                mass = np.loadtxt(os.path.join(simDirs[i],'mass.txt'))

                for j in range(0,len(borders_linear_interpolation[i])):
                    
                    i1 = borders_linear_interpolation[i][j][0]
                    i2 = borders_linear_interpolation[i][j][1]
                    
                    xi = mass[i1:i2,1]
                    A  = np.array([xi, np.ones(len(xi))])
                    yi = mass[i1:i2,2]
                    
                    pi = np.linalg.lstsq(A.T,yi)[0] #least square approximation
                    xi = mass[i1-5:i2+5,1]
                    line = pi[0]*xi + pi[1]

                    if(len(borders_linear_interpolation[i])>1 and j==0):
                        style = '--'
                    else:
                        style ='-'
                        growthRate.append(pi[0])

                    plt.plot(xi,line,'r',linestyle=style)
                    
                    print("%40s | growthrate: %3.2e" % (os.path.basename(os.path.dirname(simDirs[i])), pi[0]))

        else:

            for i in range(0,len(simDirs)):

                mass = np.loadtxt(os.path.join(simDirs[i],'mass.txt'))
                
                xi = mass[int(start_i[i]):int(end_i[i]),1]
                A  = np.array([xi, np.ones(len(xi))])
                yi = mass[int(start_i[i]):int(end_i[i]),2]
            
                pi = np.linalg.lstsq(A.T,yi)[0] #least square approximation
                line = pi[0]*xi + pi[1]
                plt.plot(xi,line,'r-')

                print("%40s | growthrate: %3.2e" % (os.path.basename(os.path.dirname(simDirs[i])), pi[0]))
                growthRate.append(pi[0])

    if(add_zoom_above):
        ax.set_ylim(0.0,max_mass*(1.0+0.2))
    else:
        ax.set_ylim(0.0,max_mass)

    if(xlim!='None'):
        ax.set_xlim(xlim[0],xlim[1])

    if(ylim!='None'):
        ax.set_ylim(ylim[0],ylim[1])

    if(legend!='None'):
        plt.legend(legend,loc=legendLoc)

    # add a zoom where the initial bubble appears
    # on the lower right part of the plot
    if(add_zoom_above):
        ax_x_lim = ax.get_xlim()
        ax_y_lim = ax.get_ylim()

        border = 0.2*(ini_time[1]-ini_time[0])

        ini_time[0]-= border
        ini_time[1]+= border
        
        ini_mass[0]  = ini_mass[0]*(1.-0.1)
        ini_mass_tmp = ini_mass[0]+(ini_time[1]-ini_time[0])*(ax_y_lim[1]-ax_y_lim[0])/(ax_x_lim[1]-ax_x_lim[0])
        
        width    = 0.5
        height_p = width*((ini_mass_tmp-ini_mass[0])/(ax_y_lim[1]-ax_y_lim[0]))/((ini_time[1]-ini_time[0])/(ax_x_lim[1]-ax_x_lim[0]))
        height   = 0.2
        
        ini_mass_tmp = height/height_p*ini_mass_tmp

        if(ini_mass_tmp>ini_mass[1]):
            ini_mass[1] = ini_mass_tmp
        
        ax_zoom = plt.axes([.15, .65, width, height]) #, axisbg='y')
        for i in range(0,len(simDirs)):
        
            mass = np.loadtxt(os.path.join(simDirs[i],'mass.txt'))
        
            grayscale_value = 0.2+ 0.8*float(i)/float(len(simDirs))
        
            p = ax_zoom.plot(mass[0:end_i[i]:step,1],
                         mass[0:end_i[i]:step,2],
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

    return growthRate,


def create_growthrate_graph(heatFlux,
                            growthRate,
                            width=3,
                            xlabel='',
                            ylabel='',
                            show=True,
                            legend='None',
                            xlim='None',
                            ylim='None'):
    '''
    @description: draw a graph with the mass growth rate
    as a function of the heat flux
    '''


    # open new figure for the mass growth rate
    plt.close("all")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)


    # plot the mass growth rate
    plt.plot(heatFlux,
             growthRate,
             'o',
             linewidth=width,
             color=grayscale_to_RGB(0.8))

    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ylabel)


    # add the linear interpolation
    xi = np.array(heatFlux)
    A  = np.array([xi, np.ones(len(xi))])
    yi = np.array(growthRate)

    pi = np.linalg.lstsq(A.T,yi)[0] #least square approximation
    xi = np.insert(xi,0,0.0)
    xi = np.append(xi,0.11)
    line = pi[0]*xi + pi[1]

    plt.plot(xi,line,'r-')

    plt.plot([-pi[1]/pi[0]],[0.0],'rs')
    plt.plot([0.0,0.11],[0.0,0.0],'k--')

    print 'equivalent_latent_heat: ', pi[0]
    print 'equivalent_temperature: ', 1.0 - (pi[0]/15.08)**2


    # graph properties
    if(xlim!='None'):
        ax.set_xlim(xlim[0],xlim[1])

    if(ylim!='None'):
        ax.set_ylim(ylim[0],ylim[1])

    if(legend!='None'):
        plt.legend(legend,loc=legendLoc)


    # show the graph
    if(show):
        plt.show()


def create_temperature_graph(simDirs,
                             legend='None',
                             width=3,
                             xlabel='',
                             ylabel='',
                             show=True,
                             figPath='',
                             borders='None',
                             legendLoc='lower right',
                             xlim='None',
                             ylim='None'):
    '''
    @description: create a graph with the temperature
    at the interface for different conditions
    '''


    # figure initialization
    plt.close("all")
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)


    for i in range(0,len(simDirs)):
        
        temperature = np.loadtxt(os.path.join(simDirs[i],'temperature.txt'))


        # data plotted
        if(borders!='None'):
            i_start = borders[i][0]
            i_end   = borders[i][1]
        else:
            
            for j in range(0,len(temperature[:,0])):
                if(temperature[j,1]>0):
                    i_start = j
                    break

            for j in range(len(temperature[:,0])-1,0,-1):
                if(temperature[j,1]>0):
                    i_end = j
                    break

        # color
        grayscale_value = 0.2+ 0.8*float(i)/float(len(simDirs))
        
        # plot
        plt.plot(temperature[i_start:i_end,1],
                 temperature[i_start:i_end,2],
                 '-',
                 linewidth=width,
                 color=grayscale_to_RGB(grayscale_value))

    ax.set_xlabel(r''+xlabel)
    ax.set_ylabel(r''+ylabel)


    # graph properties
    if(xlim!='None'):
        ax.set_xlim(xlim[0],xlim[1])

    if(ylim!='None'):
        ax.set_ylim(ylim[0],ylim[1])

    if(legend!='None'):
        plt.legend(legend,loc=legendLoc)


    # show the graph
    if(show):
        plt.show()


def compute_phase_transition_energy_supplied_rate(
    simDirs,
    heat_flux,
    mass_growthRate,
    borders):
    '''
    @description: compute the ratio of energy between the
    heat for the phase transion and the heat supplied by
    the wall
    '''
    
    print heat_flux
    print mass_growthRate

    for i in range(0,len(heat_flux)):

        temperature = np.loadtxt(os.path.join(simDirs[i],'temperature.txt'))
	
	# computation of the latent heat in the
	# phase change for the entire simulation
	latent_heat = 0.0
        T_av = 0.0

	for j in range(borders[i][0],borders[i][1]):
	    T = temperature[j,2]
	    latent_heat+= compute_latentHeat(T)
            T_av+= T
	
	# average latent heat: latent heat power per unit mass
	latent_heat /= float(borders[i][1]-borders[i][0])
        T_av        /= float(borders[i][1]-borders[i][0])
	
	# average latent heat
	latent_heat*= mass_growthRate[i]

        # equivalent temperature
        print 'borders: ', borders[i][0],borders[i][1]
        print 'av_temperature: ', T_av
        print 'eq_temperature: ', 1.0 - (latent_heat/15.08)**2
        print 'mass_growth_rate: ', mass_growthRate[i]
        print 'latent_heat: ', latent_heat/mass_growthRate[i]
        print 'Q_latent_heat: ', latent_heat

	# ratio of heat for the phase transition and
	# the heat provided by the wall
	ratio = latent_heat/heat_flux[i]
	
	print 'heat_flux: ', heat_flux[i], 'ratio phase_transition/wall: ', ratio


if __name__=='__main__':

    mainDir = os.path.join(os.getenv('HOME'),
                           'projects')


    #=============================================================
    # Options for the postprocessing graphs for the nucleation
    # on a uniform surface
    #=============================================================
    draw_angle_mass_graph = False
    draw_flux_mass_graph  = False
    draw_angle_temperature_graph = False
    draw_flux_temperature_graph = False
    compute_angle_heat_ratio = True
    compute_flux_heat_ratio = False

    
    #=============================================================
    # Nucleation study at different contact angles
    #=============================================================

    # directories for the nucleation study
    # with different contact angles
    contactAngleArray = [22.5,45.0,67.5,90.0,112.5,135.0]
    
    simDirs = []
    
    for contactAngle in contactAngleArray:
        
        simDir = 'dim2d_0.95_ca'+str(contactAngle)+'_vap'+'_fh0.02'
        simDirs.append(os.path.join(mainDir,simDir,'contours'))
    
    
    borders_linear_interpolation = []
    borders_linear_interpolation.append([[80,215]])
    borders_linear_interpolation.append([[80,215]])
    borders_linear_interpolation.append([[70,215]])
    borders_linear_interpolation.append([[50,215]])
    borders_linear_interpolation.append([[30,215]])
    borders_linear_interpolation.append([[20,170]])


    if(draw_angle_mass_graph or compute_angle_heat_ratio):
    	
    	# draw a graph with the mass = f(t) for the different simulations
    	growthRate, = create_mass_graph(simDirs,
    	                                legend=contactAngleArray,
    	                                width=3,
    	                                xlabel='$t$',
    	                                ylabel='$m_v(t)$',
    	                                show=draw_angle_mass_graph,
    	                                add_zoom_above=True,
    	                                borders_linear_interpolation=borders_linear_interpolation,
    	                                ylim=[0,0.0155])
    	
    	# extract time+mass nucleation
    	time, = get_nucleation_qties(simDirs)


    # create a graph with the temperature at the interface
    if(draw_angle_temperature_graph):

        borders = []
        for i in range(0,len(borders_linear_interpolation)):

            if(len(borders_linear_interpolation[i])>1):
                borders.append(borders_linear_interpolation[i][1])
            else:
                borders.append(borders_linear_interpolation[i][0])

        create_temperature_graph(simDirs,
                                 legend=contactAngleArray,
                                 legendLoc='upper right',
                                 width=3,
                                 xlabel='$t$',
                                 ylabel='$T(t)$',
                                 show=True,
                                 borders=borders)


    # computation of the ratio of heat used
    # for the phase transition and the total
    # heat flux provided
    if(compute_angle_heat_ratio):

        heat_flux = 0.02*np.ones(len(simDirs))
        mass_growthRate = growthRate #for fh \in [0.02,0.04]

        borders = []
        for limits in borders_linear_interpolation:
            borders.append(limits[0])
    
        compute_phase_transition_energy_supplied_rate(
            simDirs,
            heat_flux,
            mass_growthRate,
            borders)


    #=============================================================
    # Nucleation study at different flux intensities
    #=============================================================

    # directories for the nucleation study
    # with different contact angles
    heatFluxArray    = [ 0.02, 0.04, 0.06, 0.08, 0.1]
    heatFluxArrayLeg = [ 0.02, 0.04, 0.06, 0.08, 0.1]

    simDirs = []

    for heatFlux in heatFluxArray:

        simDir = 'dim2d_0.95_ca90.0_vap_fh'+str(heatFlux)
        simDirs.append(os.path.join(mainDir,simDir,'contours'))

    borders_linear_interpolation = []
    borders_linear_interpolation.append([[50,70],[60,215]])
    borders_linear_interpolation.append([[15,30],[20,110]])
    borders_linear_interpolation.append([[10,40]])
    borders_linear_interpolation.append([[5,15]])
    borders_linear_interpolation.append([[4,8]])


    if(draw_flux_mass_graph or compute_flux_heat_ratio):

        # draw a graph with the mass = f(t) for the different flux
        growthRate, = create_mass_graph(simDirs,
                                        legend=heatFluxArrayLeg,
                                        width=3,
                                        xlabel='$t$',
                                        ylabel='$m(t)$',
                                        show=draw_flux_mass_graph,
                                        add_zoom_above=True,
                                        borders_linear_interpolation=borders_linear_interpolation)

        # extract time+mass nucleation
        time, = get_nucleation_qties(simDirs)


    # create a graph with the heat flux and the growth rate
    #create_growthrate_graph(heatFluxArray,
    #                        growthRate,
    #                        width=3,
    #                        xlabel='$Q_{w,m}$',
    #                        ylabel='$\dot{m}$',
    #                        show=True,
    #                        xlim=[0,0.11])


    # create a graph with the temperature at the interface
    if(draw_flux_temperature_graph):

        borders = []
        for i in range(0,len(borders_linear_interpolation)):

            if(len(borders_linear_interpolation[i])>1):
                borders.append(borders_linear_interpolation[i][1])
            else:
                borders.append(borders_linear_interpolation[i][0])

        create_temperature_graph(simDirs,
                                 legend=heatFluxArrayLeg,
                                 legendLoc='upper right',
                                 width=3,
                                 xlabel='$t$',
                                 ylabel='$T(t)$',
                                 show=True,
                                 borders=borders)
    
    #for i in range(0,len(heatFluxArray)):
    #    
    #    print 'heat_flux: ', heatFluxArray[i], 'Q*t: ', heatFluxArray[i]*time[i]
      

    # computation of the ratio of heat used
    # for the phase transition and the total
    # heat flux provided
    if(compute_flux_heat_ratio):

        #heat_flux = heatFluxArray[0:2]    #for fh \in [0.02,0.04]
        #mass_growthRate = growthRate[0:2] #for fh \in [0.02,0.04]
        #borders = [[60,215],[20,110]]

        heat_flux = heatFluxArray    #for fh \in [0.02,0.04]
        mass_growthRate = growthRate #for fh \in [0.02,0.04]
        borders = [[60,215],[20,110],[10,40],[5,15],[4,8]]
    
        compute_phase_transition_energy_supplied_rate(
            simDirs,
            heat_flux,
            mass_growthRate,
            borders)
