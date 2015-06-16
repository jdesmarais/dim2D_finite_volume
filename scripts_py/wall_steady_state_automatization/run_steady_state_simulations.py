#!/usr/bin/python

'''
@description
script used to generate the main results of the
steady state droplet/bubble close to a wall

the results are obtained by varying the temperature,
the contact angle and the gravity in the simulations

1) contact angle: 22.5, 45.0, 67.5, 90.0 at g=0.003
   at T=0.999

2) gravity: 0.003, 0.005, 0.015 at contact angle 45
   at T=0.999

'''

import os

from library_wall_st_results import generate_wall_st_results


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
        'dim2d_bubble_next_to_wall.txt')

    contactAngleStudy         = True
    gravityStudy              = False


    #1) contact angle study
    if(contactAngleStudy):

        steady_state_ac     = 1
        temperature         = 0.999
        contact_angle_array = [135.0] #[22.5,45.0,67.5,90.0,112.5,135.0]
        phase_at_center     = 'vapor'
        gravity             = 0.000
    
        for contact_angle in contact_angle_array:
            
            [destDir, nameRun] =\
                \
                generate_wall_st_results(mainDir,
                                         steady_state_ac,
                                         temperature,
                                         contact_angle,
                                         phase_at_center,
                                         model_input,
                                         gravity_ac=1,
                                         gravity_amp=gravity)

    #2) gravity study
    if(gravityStudy):

        steady_state_ac     = 1
        temperature         = 0.999
        contact_angle       = 45.0
        phase_at_center     = 'vapor'
        gravity_array       = [0.005,0.015]
    
        for gravity in gravity_array:
            
            [destDir, nameRun] =\
                \
                generate_wall_st_results(mainDir,
                                       steady_state_ac,
                                       temperature,
                                       contact_angle,
                                       phase_at_center,
                                       model_input,
                                       gravity_ac=1,
                                       gravity_amp=gravity)
             
    
