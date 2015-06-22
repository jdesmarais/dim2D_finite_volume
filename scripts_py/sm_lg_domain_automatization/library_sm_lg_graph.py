#!/usr/bin/python

'''
@description
useful functions to generate graphs for the comparison
of the small and large domains simulations
'''

import sys
import os
import subprocess
from library_sm_lg_error import generate_fortran_exe

import numpy as np
import matplotlib.pyplot as plt


cmd_subfolder = os.path.realpath(os.path.abspath(
        os.path.join(os.getenv('augeanstables'),
                     'scripts_py',
                     'scripts_erymanthianboar',
                     'fortranfile-0.2.1')
        ))
#
#
#            os.path.split(inspect.getfile( inspect.currentframe() ))[0],
#            "fortranfile-0.2.1")))
if cmd_subfolder not in sys.path:
    sys.path.insert(0, cmd_subfolder)

from fortranfile import *


augeanstables=os.getenv('augeanstables')
extract_1D_var=os.path.join(augeanstables,
                            'scripts_fortran',
                            'error_computation',
                            'extract_1D_variable')


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


# generate fortran executable to extract 1D variables
# from netcdf file
def generate_extract_1D_real_exe():
    '''
    @description
    generate fortran executable to extract 1D variables
    from netcdf file
    '''
    
    exePath = extract_1D_var
    
    print 'generating the fortran executable to extract 1D variables: ...'

    generate_fortran_exe(exePath)

    print 'generating the fortran executable to extract 1D variables: done'
    print ''

    return exePath    


# extract the 1D variable from netcdf file
def extract_var(exePath, ncFile, varName):
    '''
    @description:
    extract the 1D variable from netcdf file
    '''

    #1) verify that the exePath exists
    if(not os.path.isfile(exePath)):
        print 'library_sm_lg_graph'
        print 'extract_var'
        sys.exit('***exePath '+exePath+' not found***')


    #2) verify that the ncFile exists
    if(not os.path.isfile(ncFile)):
        print 'library_sm_lg_graph'
        print 'extract_var'
        sys.exit('***ncFile '+ncFile+' not found***')


    #3) generate a binary output file from
    #   the netcdf file containing the variable
    outFile = os.path.join(os.path.dirname(ncFile),'var.out')
    if(os.path.isfile(outFile)):
        os.remove(outFile)

    cmd = exePath
    cmd+= ' -i '+ncFile
    cmd+= ' -o '+outFile
    cmd+= ' -v '+varName
    subprocess.call(cmd, shell=True)


    #4) verify that the binary output file
    #   was generated
    if(not os.path.isfile(outFile)):
        print 'library_sm_lg_error'
        print 'extract_var'
        sys.exit('***output file '+outFile+' not generated***')


    #5) extract the content of the output
    #   file using FortranFile fcts
    f = FortranFile(outFile)
    var = f.readReals('d')
    f.close()


    #6) remove the binary output file
    if(os.path.isfile(outFile)):
        os.remove(outFile)

    return var


# extract the error properties from netcdf file
def extract_max_error_in_time(errorPath, var_names=['time','max_error_mass']):
    '''
    @description:
    extract the maximum of the error over the domain as function
    in time from the path to the folder /error
    in the output [time_rescaled,error], 'time_rescaled' corresponds
    to a rescaled time from 0 to 1 and the 'error' is the error over
    space
    '''

    data = []
    
    #1) verify that the error_file exists
    if(not os.path.isfile(errorPath)):
        print 'library_sm_lg_graph'
        print 'extract_max_error_in_time'
        sys.exit('***errorFile '+errorPath+' was not found***')


    #2) generate the executable to extract the 
    #   variables from the file
    exePath = extract_1D_var
    if(not os.path.isfile(exePath)):
        exePath = generate_extract_1D_real_exe()


    #3) extract the other variables asked
    for var_name in var_names:

        var = extract_var(exePath, errorPath, var_name)
        if(var_name=='time'):
            time_rescaled = var/var[len(var)-1]
            var = time_rescaled

        data.append(var)

    return data


# create a graph with the different error in time
def create_error_graph(
    data,
    legendPosition='best',
    legendParam='None',
    graphPties='None',
    width=3,
    figPath='',
    show=True,
    logScale=True,
    plot_ylim='None'):
    '''
    @description:
    create a graph gathering the error in time for conditions

    - data       : [ [time1,error1] , [time2,error2] , ...]
                   where time1,...,error1,... are 1D-arrays

    - legendParam: [ leg1, leg2, ...]
                   text associated with each line

    - graphPties : [ [color1,lineType1], [color2,lineType2] , ... ]
                   where color1,lineType1] corresponds to the
                   set of data [time1,error1]
                   the graph properties indicates how the lines
                   should be displayed, by default, a gradient
                   of color will be used, but it can be set
                   explicitly by the user
                     - color1    : the color
                     - lineType1 : the type of line ('+', '-', '-+')
                     
    - width      : the default width for the lines

    - figPath    : the default path to save the figure
    '''
    
    
    # create the figure for the plot
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.figure(figsize=(8,4))
    ax = plt.subplot(111)


    # add the data on the plot
    if(graphPties=='None'):

        # add the data lines on the plot
        for graph_data,grah_pties in zip(data, graphPties):
        
            plt.plot(
                graph_data[0],
                graph_data[1],
                '-',
                linewidth=width,
                color='black')
    else:

        # add the data lines on the plot
        for graph_data,graph_pties in zip(data, graphPties):
        
            plt.plot(
                graph_data[0],
                graph_data[1],
                graph_pties[1],
                linewidth=width,
                color=graph_pties[0])
        
    # add the labels
    plt.xlabel(r"$ t $")
    plt.ylabel(r"Maximum error")

    # add a log scale
    if(logScale):
        ax.set_yscale('log')

    # add the plot limits
    if(not plot_ylim=='None'):
        ax.set_ylim(plot_ylim)

    #add the legend
    if(not legendParam=='None'):
        plt.legend(legendParam,loc=legendPosition)

    # show plot
    if(show):
        plt.show()

    # save figure
    if(not figPath==''):
        plt.savefig(figPath)


# create a graph with the different error in time
def create_error_graph_with_location(
    dataError,
    dataPosition,
    legendParam='None',
    graphPties='None',
    width=3,
    figPath='',
    show=True,
    logScale=False,
    plot_ylim='None'):
    '''
    @description:
    create a graph showing the maximum error as function of time
    and the location of the error in time

    - dataError   :[ time, error ]
                   where time and error are 1D-arrays

    - dataPosition:[ time, position ]
                   where time and position are 1D-arrays

    - legendParam: [ leg1, leg2, ...]
                   text associated with each line

    - graphPties : [ [color1,lineType1], [color2,lineType2] ]
                   where color1,lineType1] corresponds to the
                   set of data [time,error] and [color2,lineType2]
                   corresponds to the set of data [time,position]
                   the graph properties indicates how the lines
                   should be displayed, by default, a gradient
                   of color will be used, but it can be set
                   explicitly by the user
                     - color1    : the color
                     - lineType1 : the type of line ('+', '-', '-+')
                     
    - width      : the default width for the lines

    - figPath    : the default path to save the figure
    '''
    
    if(graphPties=='None'):
        graphPties = [['black','-'],['grey','-s']]        

    
    # create the figure for the plot
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    fig = plt.figure(figsize=(8,4))


    # add the data lines on the plot
    ax1 = plt.subplot(111)

    # plot the maximum of the error
    ax1.plot(dataError[0],
             dataError[1],
             graphPties[0][1],
             linewidth=width,
             color=graphPties[0][0])
    ax1.set_xlabel(r"$ t $")
    ax1.set_ylabel(r"Maximum error",color=graphPties[0][0])
    for tl in ax1.get_yticklabels():
        tl.set_color(graphPties[0][0])


    # plot the location of the maximum of the error
    ax2 = ax1.twinx()
    ax2.plot(dataPosition[0],
             dataPosition[1],
             graphPties[1][1],
             linewidth=width,
             color=graphPties[1][0])
    ax2.set_ylabel(r"Maximum error location ($y$-axis)",color='black')#graphPties[1][0])
    ax2.plot(dataPosition[0][0:len(dataPosition[0]-1):10],
             dataPosition[1][0:len(dataPosition[1]-1):10],
             's',
             linewidth=width,
             color=graphPties[1][0])
    for tl in ax2.get_yticklabels():
        tl.set_color(graphPties[1][0])


    # add a log scale
    if(logScale):
        ax1.set_yscale('log')

    # add the plot limits
    if(not plot_ylim=='None'):
        ax1.set_ylim(plot_ylim)

    #add the legend
    if(not legendParam=='None'):
        plt.legend(legendParam,loc='lower right')

    # show plot
    if(show):
        plt.show()

    # save figure
    if(not figPath==''):
        plt.savefig(figPath)


# create a graph with the different error in time
def create_error_graph_perturbation(
    data,
    legendPosition='best',
    legendParam='None',
    graphPties='None',
    width=3,
    figPath='',
    show=True,
    logScale=True,
    plot_ylim='None',
    xlabel='None'):
    '''
    @description:
    create a graph gathering the error max over time for different
    conditions

    - data       : [ [time1,error1] , [time2,error2] , ...]
                   where time1,...,error1,... are 1D-arrays

    - legendParam: [ leg1, leg2, ...]
                   text associated with each line

    - graphPties : [ [color1,lineType1], [color2,lineType2] , ... ]
                   where color1,lineType1] corresponds to the
                   set of data [time1,error1]
                   the graph properties indicates how the lines
                   should be displayed, by default, a gradient
                   of color will be used, but it can be set
                   explicitly by the user
                     - color1    : the color
                     - lineType1 : the type of line ('+', '-', '-+')
                     
    - width      : the default width for the lines

    - figPath    : the default path to save the figure
    '''
    
    
    # create the figure for the plot
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.figure(figsize=(8,4))
    ax = plt.subplot(111)


    # add the data on the plot
    if(graphPties=='None'):

        # add the data lines on the plot
        for graph_data,graph_pties in zip(data, graphPties):
        
            plt.plot(
                graph_data[0],
                graph_data[1],
                '-',
                linewidth=width,
                color='black')
    else:

        # add the data lines on the plot
        for graph_data,graph_pties in zip(data, graphPties):
        
            plt.plot(
                graph_data[0],
                graph_data[1],
                graph_pties[1],
                linewidth=width,
                color=graph_pties[0])
        
    # add the labels
    if(not xlabel=='None'):
        plt.xlabel(r""+xlabel)
    plt.ylabel(r"Maximum error")

    # add a log scale
    if(logScale):
        ax.set_xscale('log')

    ax.set_yscale('log')

    # add the plot limits
    if(not plot_ylim=='None'):
        ax.set_ylim(plot_ylim)

    #add the legend
    if(not legendParam=='None'):
        plt.legend(legendParam,loc=legendPosition)

    # show plot
    if(show):
        plt.show()

    # save figure
    if(not figPath==''):
        plt.savefig(figPath)


# create a graph for the maximum error in time
def generate_graph_max_error_in_time(errorPath,
                                     var_name=['time','max_error_mass'],
                                     width=3,
                                     figPath=''):
    '''
    @description:
    create a graph for the maximum error in time
    '''
    
    #1) extract the error in time
    [time_rescaled,error] = extract_max_error_in_time(
        errorPath,
        var_name)


    #2) plot the time as x-coordinate
    #   and max_error as y-coordinate
    create_error_graph(
        data=([[time_rescaled,error]]),
        legendParam=[0.99],
        graphPties=([['black','--']]),
        width=3,
        figPath=figPath,
        plot_ylim=[0.0001,0.1])


if __name__=="__main__":
    
    errorDir = os.path.join(os.getenv('HOME'),
                            'projects',
                            'dim2d_hedstrom_xy',
                            'dim2d_0.99_0.1',
                            'error')

    errorPath = os.path.join(errorDir,
                             'error_max.nc')

    figPath = os.path.join(errorDir,
                           'error_max_fig.eps')

    generate_graph_max_error_in_time(errorPath,
                                     figPath=figPath)

    
