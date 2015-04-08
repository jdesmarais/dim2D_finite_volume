#!/usr/bin/python

'''
@description
script used to create the error graphs of the
study of the impact of the DIM open b.c in 2D

the results are obtained by comparing simulations
on small domains with an equivalent simulation on
a large domain where the perturbations do not have
time to reach the edges and to come back

1) temperature study
\f$ T \in [0.95,0.99,0.995,0.999]\f$ at \f$v=0.1\f$

2) flow mean velocity study
\f$ v \in [0.05, 0.1, 0.25, 0.5]\f$ at \f$T=0.99\f$

3) threshold study
threshold \in [0.0001, 0.001, 0.01, 0.1, 0.2, 0.3]
for \f$ T \in [0.95,0.99,0.995,0.999] \f$
and \f$ v \in [0.05, 0.1, 0.25, 0.5]  \f$

'''

import os

from library_sm_lg_results import\
    get_simulation_dir

from library_sm_lg_graph import\
    extract_max_error_in_time,\
    create_error_graph,\
    grayscale_to_RGB


if __name__=="__main__":


    # main directory where the simulations are saved
    #------------------------------------------------------------
    mainDir = os.path.join(os.getenv('HOME'),
                           'projects',
                           '20150401_dim2d_bubble_transported')


    # dictionnaries to associate a temperature parameter with
    # an index
    #------------------------------------------------------------
    temperature_dict   = {'0.95':0,'0.99':1,'0.995':2,'0.999':3}
    flow_velocity_dict = {'0.05':0, '0.1':1, '0.25':2,  '0.5':3}
    md_threshold_dict  = {   '0':0,'0.05':1,  '0.1':2,  '0.2':3, '0.3':4}



    # options for the creation of the error graphs
    #------------------------------------------------------------
    # - temperatureStudy          : create error graph for the
    #                               temperature study
    # - velocityStudy             : create error graph for the
    #                               velocity study
    # - thresholdTemperatureStudy : create error graph for the
    #                               threshold study on temperature
    # - thresholdVelocityStudy    : create error graph for the
    #                               threshold study on velocity
    #------------------------------------------------------------
    temperatureStudy          = False
    velocityStudy             = False
    thresholdTemperatureStudy = True
    thresholdVelocityStudy    = False


    
    # options for the line drawing
    #------------------------------------------------------------
    # - color_temperature   : grayscale colors for the temperture
    # - color_flow_velocity : grayscale colors for the velocity
    # - style_threshold     : line style for the md_threshold study
    # - width               : width of the lines
    # - show                : show the graph before saving it
    # - logScale            : use log scale to display the error
    #------------------------------------------------------------
    color_temperature = [
        grayscale_to_RGB(0.30),
        grayscale_to_RGB(0.50),
        grayscale_to_RGB(0.75),
        grayscale_to_RGB(0.90)]

    color_flow_velocity = [
        grayscale_to_RGB(0.30),
        grayscale_to_RGB(0.50),
        grayscale_to_RGB(0.75),
        grayscale_to_RGB(0.90)]

    style_threshold = ['+','.','--','-','^']

    width = 3

    show = True

    logScale = True

    plot_ylim_T           = [0.00001,0.1]
    plot_ylim_v           = 'None'
    plot_ylim_T_threshold = [0.000001,0.1]
    plot_ylim_v_threshold = 'None'

    fig_T           = 'fig_error_temperature.eps'
    fig_v           = 'fig_error_velocity.eps'
    fig_T_threshold = 'fig_error_temperature_threshold.eps'
    fig_v_threshold = 'fig_error_velocity_threshold.eps'


    # paths for the error files
    main_err_dirs = os.path.join(mainDir)
    err_dirs = []

    for temperature in [0.95,0.99,0.995,0.999]:

        # add a new array to store the flow velocity entries
        # at this temperature
        err_dirs.append([])
        T_i = temperature_dict[str(temperature)]

        for flow_velocity in [0.05,0.1,0.25,0.5]:

            # add a new array to store the md_threshold entries
            # at this flow velocity
            err_dirs[T_i].append([])

            v_i = flow_velocity_dict[str(flow_velocity)]

            for md_threshold in [0,0.05,0.1,0.2,0.3]:

                dir_name = get_simulation_dir(temperature,flow_velocity,md_threshold)
                err_dir  = os.path.join(main_err_dirs,dir_name,'error','error_max.nc')

                err_dirs[T_i][v_i].append(err_dir)


    #============================================================
    #1) temperature study
    #============================================================
    if(temperatureStudy):

        temperature_array = [0.99,0.995,0.999]
        flow_velocity     = 0.1
        md_threshold      = 0
        
        # indices for the extraction of the error path
        v_i  = flow_velocity_dict[str(flow_velocity)]
        md_i = md_threshold_dict[str(md_threshold)]
    
        # extract the error data for the temperature study
        data = []
        legendParam = []
        graphPties = []
        for temperature in temperature_array:
            
            # get the path to the error data
            T_i  = temperature_dict[str(temperature)]
            errorPath = err_dirs[T_i][v_i][md_i]

            # extract the error data
            [time_rescaled,error] = extract_max_error_in_time(errorPath)
            
            # stored the error data for the graph
            data.append([time_rescaled,error])

            # legend
            legendParam.append(temperature)

            # color
            graphPties.append([color_temperature[T_i],'-'])


        # create the graph
        create_error_graph(
            data,
            legendParam=legendParam,
            graphPties=graphPties,
            width=width,
            figPath=os.path.join(mainDir,fig_T),
            show=show,
            logScale=logScale,
            plot_ylim=plot_ylim_T)


    #============================================================
    #2) velocity study
    #============================================================
    if(velocityStudy):

        temperature         = 0.99
        flow_velocity_array = [0.05,0.1,0.25,0.5]
        md_threshold        = 0
        
        # indices for the extraction of the error path
        T_i  = temperature_dict[str(temperature)]
        md_i = md_threshold_dict[str(md_threshold)]
    
        # extract the error data for the velocity study
        data = []
        legendParam = []
        graphPties = []
        for flow_velocity in flow_velocity_array:
            
            # get the path to the error data
            v_i  = flow_velocity_dict[str(flow_velocity)]
            errorPath = err_dirs[T_i][v_i][md_i]

            # extract the error data
            [time_rescaled,error] = extract_max_error_in_time(errorPath)
            
            # stored the error data for the graph
            data.append([time_rescaled,error])

            # legend
            legendParam.append(flow_velocity)

            # graph properties
            graphPties.append([color_flow_velocity[v_i],'-'])


        # create the graph
        create_error_graph(
            data,
            legendParam=legendParam,
            graphPties=graphPties,
            width=width,
            figPath=os.path.join(mainDir,fig_v),
            show=show,
            logScale=logScale,
            plot_ylim=plot_ylim_v)
             

    #============================================================    
    #3) threshold study
    #============================================================
    if(thresholdTemperatureStudy):

        temperature_array  = [0.995]#[0.95,0.99,0.995,0.999]
        flow_velocity      = 0.1
        md_threshold_array = [0.05, 0.1, 0.2, 0.3]
        
        v_i  = flow_velocity_dict[str(flow_velocity)]

        # extract the error data for the temperature threshold study
        data = []
        legendParam = []
        graphPties = []
        for md_threshold in md_threshold_array:

            # indices for the extraction of the error path
            md_i = md_threshold_dict[str(md_threshold)]

            for temperature in temperature_array:

                # indices for the extraction of the error path
                T_i  = temperature_dict[str(temperature)]
                errorPath = err_dirs[T_i][v_i][md_i]
                
                # extract the error data
                [time_rescaled,error] = extract_max_error_in_time(errorPath)
            
                # stored the error data for the graph
                data.append([time_rescaled,error])

                #legend
                legendParam.append(str(temperature)+'\_'+str(md_threshold))

                # create the graph properties
                graphPties.append([color_temperature[T_i],
                                   style_threshold[md_i]])

        # create the graph
        create_error_graph(
            data,
            legendParam=legendParam,
            graphPties=graphPties,
            width=width,
            figPath=os.path.join(mainDir,fig_T_threshold),
            show=show,
            logScale=logScale,
            plot_ylim=plot_ylim_T_threshold)


    if(thresholdVelocityStudy):

        temperature_array  = 0.99
        flow_velocity      = [0.05,0.1,0.25,0.5]
        md_threshold_array = [0.05, 0.1, 0.2, 0.3]
        
        T_i  = temperature_dict[str(temperature)]

        # extract the error data for the temperature threshold study
        data = []
        legendParam = []
        graphPties = []
        for md_threshold in md_threshold_array:

            # indices for the extraction of the error path
            md_i = md_threshold_dict[str(md_threshold)]

            for flow_velocity in flow_velocity_array:

                # indices for the extraction of the error path
                v_i  = flow_velocity_dict[str(flow_velocity)]
                
                # extract the error data
                [time_rescaled,error] = extract_max_error_in_time(errorPath)
            
                # stored the error data for the graph
                data.append([time_rescaled,error])

                #legend
                legendParam.append(str(flow_velocity)+'_'+str(md_threshold))

                # create the graph properties
                graphPties.append([color_flow_velocity[v_i],
                                   style_threshold[md_i]])

        # create the graph
        create_error_graph(
            data,
            legendParam=legendParam,
            graphPties=graphPties,
            width=width,
            figPath=os.path.join(mainDir,fig_v_threshold),
            show=show,
            logScale=logScale,
            plot_ylim=plot_ylim_v_threshold)
    

    
        
