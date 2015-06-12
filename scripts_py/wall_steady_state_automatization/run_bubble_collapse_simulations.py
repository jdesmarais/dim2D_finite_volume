#!/usr/bin/python

'''
@description
script used to generate the main results of the
bubble collapse closed to a wall

the results are obtained by varying the ratio between
the initial radius of the bubble and its interface width
and the contact angle at the wall between the phases

1) ratio: 1.0, 1.5, 2.0 at T=0.95 and contact angle=90.0

1) contact angle: 22.5, 45.0, 67.5, 90.0 at T=0.95 and ratio=1.0

'''

import os

from library_wall_bubblecollapse_results import generate_wall_bubblecollapse_results


if __name__=="__main__":


    # main directory where the simulations
    # are saved
    mainDir='/home/jdesmarais/projects'


    # input used as template
    #------------------------------------------------------------
    maindir_input = os.path.join(os.getenv('augeanstables'),
                                 'src','config',
                                 'default_inputs','dim2d')

    model_input = os.path.join(
        maindir_input,
        'dim2d_bubble_collapse.txt')

    ratioStudy        = True
    contactangleStudy = False


    #1) ratio study
    if(ratioStudy):

        steady_state_ac = 0
        temperature     = 0.95
        ratio_array     = [1.5,2.0]
        phase_at_center = 'vapor'
        gravity_ac      = 0
        gravity         = 0.000
        contact_angle   = 90.0
        total_nb_files  = 1200
    
        for ratio in ratio_array:
            
            [destDir, nameRun] =\
                \
                generate_wall_bubblecollapse_results(
                mainDir,
                steady_state_ac,
                temperature,
                contact_angle,
                phase_at_center,
                ratio,
                model_input,
                gravity_ac=gravity_ac,
                gravity_amp=gravity,
                total_nb_files=total_nb_files)

    #2) contact_angle study
    if(contactangleStudy):

        steady_state_ac     = 0
        temperature         = 0.95
        ratio               = 1.0
        contact_angle_array = [22.5,45.0,67.5,90.0,112.5,135.0]
        phase_at_center     = 'vapor'
        gravity_ac          = 0
        gravity             = 0.000
        contact_angle       = 90.0
        total_nb_files      = 1200
    
        for contact_angle in contact_angle_array:
            
            [destDir, nameRun] =\
                \
                generate_wall_bubblecollapse_results(
                mainDir,
                steady_state_ac,
                temperature,
                contact_angle,
                phase_at_center,
                ratio,
                model_input,
                gravity_ac=gravity_ac,
                gravity_amp=gravity,
                total_nb_files=total_nb_files)
             
    
