#!/usr/bin/python

import sys

bc_type_code = ['bc_nodes_choice',
                'bc_fluxes_choice',
                'bc_timedev_choice',
                'bc_flux_and_node_choice']


def get_flow_config(flow_direction,inputs_computed):
    '''
    @description
    determine the configuration of the flow
    '''

    #in order to set the correct direction of the flow
    #(N,S,E,W,NE,NW,SE,SW) from the inputs.txt, three
    #parameters are initialized in parameters_input.f
    #[1]: the direction of the flow is either horizontal
    #(x_direction), vertical (y_direction) or diagonal
    #(xy_direction)
    #[2]: whether the flow is right or left (1.0/-1.0)
    #[3]: whether the flow is upward or downwards (1.0/-1.0)
    flow_direction_code = {
        'N':  [ 'y_direction', 1.0, 1.0],
        'S':  [ 'y_direction', 1.0,-1.0],
        'E':  [ 'x_direction', 1.0, 1.0],
        'W':  [ 'x_direction',-1.0, 1.0],
        'NE': ['xy_direction', 1.0, 1.0],
        'NW': ['xy_direction',-1.0, 1.0],
        'SE': ['xy_direction', 1.0,-1.0],
        'SW': ['xy_direction',-1.0,-1.0]}

    inputs_computed['flow_direction'] = flow_direction_code[flow_direction][0]
    inputs_computed['flow_x_side']    = flow_direction_code[flow_direction][1]
    inputs_computed['flow_y_side']    = flow_direction_code[flow_direction][2]



def get_bc_config(bc_choice,inputs_computed):
    '''
    @description
    determine the configuration of the boundary conditions
       - bc_{N,S,E,W,NW,NE,SW,SE}_choice,
       - bc_{N,S,E,W,NW,NE,SW,SE}_type
       - bc_order{1,2,3,4,5,6,7,8}
    '''
    
    if(bc_choice=='reflection_xy_choice'):

        inputs_computed['bc_N_choice'] = 'reflection_y_choice'
        inputs_computed['bc_S_choice'] = 'reflection_y_choice'
        inputs_computed['bc_E_choice'] = 'reflection_x_choice'
        inputs_computed['bc_W_choice'] = 'reflection_x_choice'

        inputs_computed['bc_NE_choice'] = 'reflection_y_choice'
        inputs_computed['bc_NW_choice'] = 'reflection_y_choice'
        inputs_computed['bc_SE_choice'] = 'reflection_y_choice'
        inputs_computed['bc_SW_choice'] = 'reflection_y_choice'

        set_bc_order(['W','E','SW','S','SE','NW','N','NE'],
                     inputs_computed)


    elif(bc_choice=='periodic_xy_choice'):

        print 'library_config'
        print 'get_bc_config'
        print 'periodic not implemented'
        sys.exit(2)

        inputs_computed['bc_N_choice'] = 'periodic_y_choice'
        inputs_computed['bc_S_choice'] = 'periodic_y_choice'
        inputs_computed['bc_E_choice'] = 'periodic_x_choice'
        inputs_computed['bc_W_choice'] = 'periodic_x_choice'

        inputs_computed['bc_NE_choice'] = 'periodic_x_choice'
        inputs_computed['bc_NW_choice'] = 'periodic_x_choice'
        inputs_computed['bc_SE_choice'] = 'periodic_x_choice'
        inputs_computed['bc_SW_choice'] = 'periodic_x_choice'

        set_bc_order(['S','N','SW','SE','W','E','NW','NE'],
                     inputs_computed)


    elif(bc_choice=='wall_xy_choice'):

        print 'library_config'
        print 'get_bc_config'
        print 'wall_xy_choice not implemented'
        print 'the corners are not computed correctly'
        sys.exit(2)

       
    elif(bc_choice=='wall_S_reflection_choice'):

        inputs_computed['bc_N_choice'] = 'reflection_y_choice'
        inputs_computed['bc_S_choice'] = 'wall_choice'
        inputs_computed['bc_E_choice'] = 'reflection_x_choice'
        inputs_computed['bc_W_choice'] = 'reflection_x_choice'

        inputs_computed['bc_NW_choice'] = 'reflection_x_choice'
        inputs_computed['bc_NE_choice'] = 'reflection_x_choice'
        inputs_computed['bc_SW_choice'] = 'reflection_x_choice'
        inputs_computed['bc_SE_choice'] = 'reflection_x_choice'

        set_bc_order(['W','E','S','N','SW','SE','NW','NE'],
                     inputs_computed)


    elif(bc_choice=='wall_S_open_choice'):

        inputs_computed['bc_N_choice'] = 'hedstrom_choice'
        inputs_computed['bc_S_choice'] = 'wall_choice'
        inputs_computed['bc_E_choice'] = 'hedstrom_choice'
        inputs_computed['bc_W_choice'] = 'hedstrom_choice'

        inputs_computed['bc_NW_choice'] = 'hedstrom_choice'
        inputs_computed['bc_NE_choice'] = 'hedstrom_choice'
        inputs_computed['bc_SW_choice'] = 'wall_choice'
        inputs_computed['bc_SE_choice'] = 'wall_choice'

        set_bc_order(['SW','S','SE','W','E','NW','N','NE'],
                     inputs_computed)


    elif(bc_choice=='hedstrom_xy_choice'):
        
        inputs_computed['bc_N_choice'] = 'hedstrom_choice'
        inputs_computed['bc_S_choice'] = 'hedstrom_choice'
        inputs_computed['bc_E_choice'] = 'hedstrom_choice'
        inputs_computed['bc_W_choice'] = 'hedstrom_choice'

        inputs_computed['bc_NW_choice'] = 'hedstrom_choice'
        inputs_computed['bc_NE_choice'] = 'hedstrom_choice'
        inputs_computed['bc_SW_choice'] = 'hedstrom_choice'
        inputs_computed['bc_SE_choice'] = 'hedstrom_choice'

        set_bc_order(['SW','S','SE','W','E','NW','N','NE'],
                     inputs_computed)


    elif(bc_choice=='hedstrom_xy_corners_choice'):
        
        print 'library_config'
        print 'get_bc_config'
        print 'hedstrom_xy_corners_choice not implemented'
        print 'the corners are not implemented yet'
        sys.exit(2)

    
    elif(bc_choice=='poinsot_xy_choice'):

        print 'library_config'
        print 'get_bc_config'
        print 'poinsot_xy_choice'
        print 'the bc_operators have not been adapted yet'
        print 'for the generalized boundary operator'
        sys.exit(2)

        inputs_computed['bc_N_choice'] = 'poinsot_choice'
        inputs_computed['bc_S_choice'] = 'poinsot_choice'
        inputs_computed['bc_E_choice'] = 'poinsot_choice'
        inputs_computed['bc_W_choice'] = 'poinsot_choice'

        inputs_computed['bc_NW_choice'] = 'poinsot_choice'
        inputs_computed['bc_NE_choice'] = 'poinsot_choice'
        inputs_computed['bc_SW_choice'] = 'poinsot_choice'
        inputs_computed['bc_SE_choice'] = 'poinsot_choice'


    elif(bc_choice=='yoolodato_xy_choice'):

        print 'library_config'
        print 'get_bc_config'
        print 'yoolodato_xy_choice'
        print 'the bc_operators have not been adapted yet'
        print 'for the generalized boundary operator'
        sys.exit(2)

        inputs_computed['bc_N_choice'] = 'yoolodato_choice'
        inputs_computed['bc_S_choice'] = 'yoolodato_choice'
        inputs_computed['bc_E_choice'] = 'yoolodato_choice'
        inputs_computed['bc_W_choice'] = 'yoolodato_choice'

        inputs_computed['bc_NW_choice'] = 'yoolodato_choice'
        inputs_computed['bc_NE_choice'] = 'yoolodato_choice'
        inputs_computed['bc_SW_choice'] = 'yoolodato_choice'
        inputs_computed['bc_SE_choice'] = 'yoolodato_choice'


    else:

        print 'library_config'
        print 'get_bc_config'
        print 'bc_choice not recognized: '+str(bc_choice)
        sys.exit(2)

    set_bc_type(inputs_computed)


def set_bc_type(inputs_computed):
    '''
    @description
    initialize the type of all the boundary procedures
    '''
    
    cardinal_points = ['N','S','E','W','SW','SE','NW','NE']

    for cardinal_pt in cardinal_points:
        
        key_bc      = 'bc_'+cardinal_pt+'_choice'
        key_bc_type = 'bc_'+cardinal_pt+'_type_choice'
        inputs_computed[key_bc_type] = get_bc_type(inputs_computed[key_bc])


def set_bc_order(order,inputs_computed):
    '''
    @description
    turn the order how the b.c. are applied into
    a boundary type order understandable by the
    fortran program
    '''

    for i in range(0,8):
        key = 'bc_order'+str(i+1)
        inputs_computed[key] = get_bc_order_type(order[i])


def get_bc_type(bc_choice):
    '''
    @description
    associate the type of boundary procedure
    with its bc_choice
    '''

    if(bc_choice=='reflection_x_choice' or
       bc_choice=='reflection_y_choice' or
       bc_choice=='periodic_x_choice' or
       bc_choice=='periodic_y_choice'):
        return bc_type_code[0]

    elif(bc_choice=='wall_choice'):
        return bc_type_code[3]

    elif(bc_choice=='hedstrom_choice' or
         bc_choice=='poinsot_choice' or
         bc_choice=='yoolodato_choice'):
        return bc_type_code[2]

    else:
        print 'library_config'
        print 'get_bc_type'
        print 'bc_choice not recognized: '+str(bc_choice)
        sys.exit(2)


def get_bc_order_type(bc_type):
    '''
    @description
    turn the cardinal coordinates into a boundary type
    '''

    if(bc_type=='N'):
        return 'N_edge_type'

    elif(bc_type=='S'):
        return 'S_edge_type'

    elif(bc_type=='E'):
        return 'E_edge_type'

    elif(bc_type=='W'):
        return 'W_edge_type'

    if(bc_type=='NW'):
        return 'NW_corner_type'

    elif(bc_type=='NE'):
        return 'NE_corner_type'

    elif(bc_type=='SW'):
        return 'SW_corner_type'

    elif(bc_type=='SE'):
        return 'SE_edge_type'

    else:
        print 'library_config'
        print 'get_bc_order_type'
        print 'bc_type not recognized: '+str(bc_type)
        sys.exit(2)
    
