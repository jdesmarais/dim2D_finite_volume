      program test_bf_layer_time

        use bc_operators_class, only :
     $       bc_operators

        use bf_layer_time_class, only :
     $       bf_layer_time

        use check_data_module, only :
     $       is_real_vector_validated,
     $       is_real_matrix3D_validated,
     $       is_int_vector_validated,
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       no_pt,
     $       bc_pt,
     $       bc_interior_pt,
     $       interior_pt,
     $       
     $       align_N,align_S,
     $       align_E,align_W,
     $       
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $     
     $       no_overlap,
     $       N_overlap,
     $       S_overlap,
     $       E_overlap,
     $       W_overlap,
     $       NE_overlap,
     $       NW_overlap,
     $       SE_overlap,
     $       SW_overlap,
     $       NS_overlap,
     $       EW_overlap,
     $       
     $       BF_SUCCESS,
     $       
     $       cpt1normal_and_cpt4normal,   
     $       cpt1normal_and_cpt4not,  
     $       cpt1normal_and_cpt4overlap,
     $       cpt1not_and_cpt4normal,
     $       cpt1not_and_cpt4not,
     $       cpt1not_and_cpt4overlap,
     $       cpt1overlap_and_cpt4normal,
     $       cpt1overlap_and_cpt4not,
     $       cpt1overlap_and_cpt4overlap,
     $       
     $       cpt2normal_and_cpt3normal,
     $       cpt2normal_and_cpt3not,  
     $       cpt2normal_and_cpt3overlap,
     $       cpt2not_and_cpt3normal, 
     $       cpt2not_and_cpt3not,     
     $       cpt2not_and_cpt3overlap,
     $       cpt2overlap_and_cpt3normal,
     $       cpt2overlap_and_cpt3not,
     $       cpt2overlap_and_cpt3overlap

        use parameters_constant, only :
     $       bc_timedev_choice,
     $       hedstrom_xy_choice,
     $       N,S,E,W
        
        use parameters_input, only :
     $       x_min, x_max,
     $       y_min, y_max,
     $       nx,ny,ne,
     $       bc_size,
     $       bc_choice,
     $       bc_N_type_choice,
     $       bc_S_type_choice,
     $       bc_E_type_choice,
     $       bc_W_type_choice
        
        use parameters_kind, only :
     $       ikind,
     $       rkind
        
        use pmodel_eq_class, only :
     $       pmodel_eq

        use rk3tvd_steps_module, only :
     $       compute_1st_step,
     $       compute_1st_step_nopt
        
        use sd_operators_class, only :
     $       sd_operators
        
        use td_operators_class, only :
     $       td_operators

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled      = .true.
        test_validated = .true.


        test_loc = test_update_bc_sections(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_update_bc_sections: '',L1)', test_loc
        print '()'

        
        test_loc = test_update_integration_borders(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_update_integration_borders: '',L1)', test_loc
        print '()'


        test_loc = test_compute_time_dev(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_time_dev: '',L1)', test_loc
        print '()'


        test_loc = test_compute_integration_step(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_integration_step: '',L1)', test_loc
        print '()'


        test_loc = test_set_x_borders(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_set_x_borders: '',L1)', test_loc
        print '()'


        test_loc = test_get_x_borders(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_x_borders: '',L1)', test_loc
        print '()'


        test_loc = test_set_y_borders(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_set_y_borders: '',L1)', test_loc
        print '()'


        test_loc = test_get_y_borders(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_y_borders: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_update_bc_sections(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_time)             :: bf_layer_used
          integer(ikind), dimension(2)    :: x_borders
          integer(ikind), dimension(2)    :: y_borders
          integer       , dimension(5,60) :: test_bc_sections
          
          logical :: test_loc


          test_validated = .true.


          allocate(bf_layer_used%grdpts_id(30,26))

          bf_layer_used%localization = N

          x_borders = [5,30]
          y_borders = [3,24]

          call bf_layer_used%set_x_borders(x_borders)
          call bf_layer_used%set_y_borders(y_borders)

          
          !          ________________________________________________________
          !    |     |                                                      |    |
          !    |     |                                                      |    |
          !24- |     |                      3 3 3 3 3 3                     |    |
          !    |     |              3 3 3 3 3 2 2 2 2 3 3 3 3 3             |    |
          !22- |     |              3 2 2 2 2 2     2 2 2 2 2 3             |    |
          !    |     |      3 3 3 3 3 2                     2 3 3 3 3 3     |    |
          !20- |     |    3 3 2 2 2 2 2                     2 2 2 2 2 3 3   |    |
          !19- |     |  3 3 2 2                                     2 2 3 3 |    |
          !    |     |  3 2 2                                         2 2 3 |    |
          !17- |     |  3 2                                             2 3 |    |
          !    |     |  3 2                                             2 3 |    |
          !15- |     |  3 2 2                                         2 2 3 |    |
          !    |     |  3 3 2 2                                     2 2 3 3 |    |
          !13- |     |    3 3 2 2 2 2 2                     2 2 2 2 2 3 3   |    |
          !12- |     |      3 3 3 3 3 2                     2 3 3 3 3 3     |    |
          !    |     |              3 2                     2 3             |    |
          !10- |     |              3 2 2 2 2 2     2 2 2 2 2 3             |    |
          !    |     |              3 3 3 3 3 2     2 3 3 3 3 3             |    |
          !    |     |                      3 2     2 3                     |    |
          !    |     |          3 3 3 3 3 3 3 2     2 3 3 3 3 3 3 3         |    |
          ! 6- |     |          3 2 2 2 2 2 2 2     2 2 2 2 2 2 2 3         |    |
          !    |     |          3 2                             2 3         |    |
          ! 4- |     |          3 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 3         |    |
          !    |     |          3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3 3         |    |
          !    |     |                                                      |    |
          !    |     |                                                      |    |
          !     ---- |------------------------------------------------------|-----
          !             | | |   |   |       |       |         | |     | | |
          !             3 4 5   7   9       13      17        2223    262728
          ! --------------------------------------------------
          bf_layer_used%grdpts_id = reshape(
     $         (/
     $         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 3, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 3, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 3, 2, 1, 1, 2, 3, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 2, 1, 1, 2, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 3, 3, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 3, 0, 0, 0, 0,
     $         0, 0, 0, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 0, 0, 0,
     $         0, 0, 3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 0, 0,
     $         0, 0, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 0, 0,
     $         0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0,
     $         0, 0, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 0, 0,
     $         0, 0, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 0, 0,
     $         0, 0, 3, 3, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 3, 3, 0, 0,
     $         0, 0, 0, 3, 3, 2, 2, 2, 2, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 0, 0, 0,
     $         0, 0, 0, 0, 3, 3, 3, 3, 3, 2, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 3, 3, 3, 3, 3, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 0, 0, 3, 2, 2, 2, 2, 2, 1, 1, 2, 2, 2, 2, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 2, 2, 2, 2, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 3, 3, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
     $         0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0/),
     $         (/30,26/))

          test_bc_sections = reshape((/
     $         SW_corner_type, 7 ,  3, no_overlap                , no_overlap,
     $         S_edge_type   , 9 ,  3, 22                        , no_overlap,
     $         SE_corner_type, 23,  3, no_overlap                , no_overlap,
     $         W_edge_type   , 7 ,  5, 5                         , no_overlap,
     $         E_edge_type   , 23,  5, 5                         , no_overlap,
     $         NW_corner_type, 7 ,  6, no_overlap                , no_overlap,
     $         N_edge_type   , 9 ,  6, 12                        , no_overlap,
     $         NW_edge_type  , 13,  6, no_overlap                , no_overlap,
     $         NE_edge_type  , 17,  6, no_overlap                , no_overlap,
     $         N_edge_type   , 19,  6, 22                        , no_overlap,
     $         NE_corner_type, 23,  6, no_overlap                , no_overlap,
     $         W_edge_type   , 13,  8, 8                         , no_overlap,
     $         E_edge_type   , 17,  8, 8                         , no_overlap,
     $         SW_corner_type, 9 ,  9, no_overlap                , no_overlap,
     $         S_edge_type   , 11,  9, 12                        , no_overlap,
     $         SW_edge_type  , 13,  9, no_overlap                , no_overlap,
     $         SE_edge_type  , 17,  9, no_overlap                , no_overlap,
     $         S_edge_type   , 19,  9, 20                        , no_overlap,
     $         SE_corner_type, 21,  9, no_overlap                , no_overlap,
     $         W_edge_type   , 9 , 11, 11                        , no_overlap,
     $         E_edge_type   , 21, 11, 11                        , no_overlap,
     $         SW_corner_type, 5 , 12, cpt2normal_and_cpt3not    , no_overlap,
     $         S_edge_type   , 7 , 12, 8                         , no_overlap,
     $         SW_edge_type  , 9 , 12, no_overlap                , no_overlap,
     $         SE_edge_type  , 21, 12, no_overlap                , no_overlap,
     $         S_edge_type   , 23, 12, 24                        , no_overlap,
     $         SE_corner_type, 25, 12, cpt1normal_and_cpt4not    , no_overlap,
     $         SW_corner_type, 4 , 13, cpt2overlap_and_cpt3not   , W_overlap,
     $         SW_edge_type  , 5 , 13, cpt2normal_and_cpt3not    , SW_overlap,
     $         SE_edge_type  , 25, 13, cpt1normal_and_cpt4not    , SE_overlap,
     $         SE_corner_type, 26, 13, cpt1overlap_and_cpt4not   , no_overlap,
     $         SW_corner_type, 3 , 14, cpt2overlap_and_cpt3normal, EW_overlap,
     $         SW_edge_type  , 4 , 14, cpt2overlap_and_cpt3normal, SW_overlap,
     $         SE_edge_type  , 26, 14, cpt1overlap_and_cpt4normal, SE_overlap,
     $         SE_corner_type, 27, 14, cpt1overlap_and_cpt4normal, no_overlap,
     $         W_edge_type   , 3 , 16, 17                        , EW_overlap,
     $         E_edge_type   , 27, 16, 17                        , no_overlap,
     $         NW_corner_type, 3 , 18, cpt1normal_and_cpt4not    , EW_overlap,
     $         NW_edge_type  , 4,  18, cpt1normal_and_cpt4not    , NW_overlap,
     $         NE_edge_type  , 26, 18, cpt2normal_and_cpt3not    , NE_overlap,
     $         NE_corner_type, 27, 18, cpt2normal_and_cpt3not    , no_overlap,
     $         NW_corner_type, 4 , 19, cpt1overlap_and_cpt4not   , W_overlap,
     $         NW_edge_type  , 5 , 19, cpt1overlap_and_cpt4normal, NW_overlap,
     $         NE_edge_type  , 25, 19, cpt2overlap_and_cpt3normal, NE_overlap,
     $         NE_corner_type, 26, 19, cpt2overlap_and_cpt3not   , no_overlap,
     $         NW_corner_type,  5, 20, cpt1overlap_and_cpt4normal, no_overlap,
     $         N_edge_type   ,  7, 20, 8                         , no_overlap,
     $         NW_edge_type  ,  9, 20, no_overlap                , no_overlap,
     $         NE_edge_type  , 21, 20, no_overlap                , no_overlap,
     $         N_edge_type   , 23, 20, 24                        , no_overlap,
     $         NE_corner_type, 25, 20, cpt2overlap_and_cpt3normal, no_overlap,
     $         NW_corner_type, 9 , 22, no_overlap                , no_overlap,
     $         N_edge_type   , 11, 22, 12                        , no_overlap,
     $         NW_edge_type  , 13, 22, no_overlap                , N_overlap,
     $         NE_edge_type  , 17, 22, no_overlap                , N_overlap,
     $         N_edge_type   , 19, 22, 20                        , no_overlap,
     $         NE_corner_type, 21, 22, no_overlap                , no_overlap,
     $         NW_corner_type, 13, 23, no_overlap                , no_overlap,
     $         N_edge_type   , 15, 23, 16                        , no_overlap,
     $         NE_corner_type, 17, 23, no_overlap                , no_overlap
     $         /),
     $         (/5,60/))

          !North
          bf_layer_used%localization = N
          call bf_layer_used%update_bc_sections()

          test_loc = is_int_matrix_validated(
     $         bf_layer_used%bc_sections,
     $         test_bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc


          !South
          bf_layer_used%localization = S
          call bf_layer_used%update_bc_sections()

          test_loc = is_int_matrix_validated(
     $         bf_layer_used%bc_sections,
     $         test_bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc


          !East
          bf_layer_used%localization = E
          call bf_layer_used%update_bc_sections()

          test_loc = is_int_matrix_validated(
     $         bf_layer_used%bc_sections,
     $         test_bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc

          
          !West
          bf_layer_used%localization = W
          call bf_layer_used%update_bc_sections()

          test_loc = is_int_matrix_validated(
     $         bf_layer_used%bc_sections,
     $         test_bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc

        end function test_update_bc_sections


        function test_update_integration_borders(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_time)        :: bf_layer_used
          integer, dimension(4)      :: test_localization
          logical, dimension(2,4)    :: test_neighbors
          integer, dimension(2,2,10) :: test_borders
          integer                    :: i,k,l
          logical                    :: test_loc


          test_validated = .true.


          test_localization = [N,S,E,W]
          test_neighbors    = reshape((/
     $         .false.,.false.,
     $         .true. ,.false.,
     $         .false.,.true.,
     $         .true. ,.true./),
     $         (/2,4/))
          test_borders = reshape((/
     $         1,3,7,6,
     $         1,1,7,4,
     $         
     $         3,1,7,6,
     $         3,3,7,6,
     $         3,1,7,4,
     $         3,3,7,4,
     $         
     $         1,1,5,6,
     $         1,3,5,6,
     $         1,1,5,4,
     $         1,3,5,4/),
     $         (/2,2,10/))

          allocate(bf_layer_used%nodes(7,6,ne))

          !N,S
          do k=1,2

             !input
             call bf_layer_used%ini(test_localization(k))

             !output
             call bf_layer_used%update_integration_borders()

             !validation
             test_loc = is_int_vector_validated(
     $            bf_layer_used%x_borders,
     $            test_borders(1,:,k),
     $            detailled)
             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'',1) failed'')',k
             end if

             test_loc = is_int_vector_validated(
     $            bf_layer_used%y_borders,
     $            test_borders(2,:,k),
     $            detailled)
             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'',2) failed'')',k
             end if
             
          end do


          !E,W
          do k=3,4

             !input
             call bf_layer_used%ini(test_localization(k))

             do l=1, size(test_neighbors,2)

                call bf_layer_used%set_neighbor1_share(test_neighbors(1,l))
                call bf_layer_used%set_neighbor2_share(test_neighbors(2,l))

                !output
                call bf_layer_used%update_integration_borders()

                !validation
                if(k.eq.3) then
                   i=3+(l-1)
                else
                   i=7+(l-1)
                end if
                test_loc = is_int_vector_validated(
     $               bf_layer_used%x_borders,
     $               test_borders(1,:,i),
     $               detailled)
                test_validated = test_validated.and.test_loc
                
                if(detailled.and.(.not.test_loc)) then
                   print '(''test('',2I2,'',1) failed'')',k,l
                end if

                test_loc = is_int_vector_validated(
     $               bf_layer_used%y_borders,
     $               test_borders(2,:,i),
     $               detailled)
                test_validated = test_validated.and.test_loc

                if(detailled.and.(.not.test_loc)) then
                   print '(''test('',2I2,'',2) failed'')',k,l
                end if

             end do

          end do

        end function test_update_integration_borders



        function test_compute_time_dev(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_time)                 :: bf_layer_used
          type(td_operators)                  :: td_operators_used
          real(rkind)                         :: t
          type(sd_operators)                  :: s
          type(pmodel_eq)                     :: p_model
          type(bc_operators)                  :: bc_used
          real(rkind)   , dimension(nx)       :: x_map
          real(rkind)   , dimension(ny)       :: y_map
          real(rkind)   , dimension(nx,ny,ne) :: nodes
          real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
          real(rkind)   , dimension(nx,ny,ne) :: timedev_test

          integer(ikind) :: i,j
          integer        :: k

          !input
          t = 0.0d0          

          !alignment
          bf_layer_used%alignment = reshape((/
     $         align_E,
     $         align_S+1,
     $         align_E+nx-2*bc_size,
     $         align_S+1+ny-2*bc_size/),
     $         (/2,2/))

          !grdpts_id
          allocate(bf_layer_used%grdpts_id(nx,ny))
          bf_layer_used%grdpts_id = reshape(
     $         (/ ((interior_pt, i=1,nx),j=1,ny) /),
     $         (/nx,ny/))

          !x_map
          allocate(bf_layer_used%x_map(nx))
          x_map = (/ (x_min+(i-1)*0.1d0,i=1,nx) /)
          bf_layer_used%x_map = x_map

          !y_map
          allocate(bf_layer_used%y_map(ny))
          y_map = (/ (y_min+(j-1)*0.2d0,j=1,ny) /)
          bf_layer_used%y_map = (/ (y_min+(j-1)*0.2d0,j=1,ny) /)

          !nodes
          nodes = reshape((/
     $         1.48d0, 1.30d0, 1.35d0, 1.31d0, 1.43d0, 1.31d0,
     $         1.26d0, 1.45d0,  1.4d0, 1.29d0, 1.37d0, 1.41d0,
     $         1.46d0, 1.27d0, 1.47d0, 1.28d0, 1.25d0, 1.43d0,
     $         1.48d0, 1.26d0, 1.41d0, 1.34d0, 1.31d0, 1.39d0,
     $         1.42d0, 1.46d0, 1.38d0, 1.26d0, 1.37d0, 1.33d0,
     $         1.41d0, 1.22d0, 1.42d0, 1.23d0, 1.21d0, 1.40d0,
     $         
     $         0.128d0, 0.127d0, 0.142d0, 0.129d0, 0.136d0, 0.124d0,
     $         1.138d0, 0.148d0, 0.132d0, 0.125d0, 0.175d0, 0.123d0,
     $         0.146d0, 0.143d0, 0.145d0, 0.182d0, 0.135d0, 0.154d0,
     $         0.123d0, 0.129d0, 0.124d0, 0.162d0, 0.152d0, 0.142d0,
     $         0.168d0, 0.198d0, 0.186d0, 0.163d0, 0.126d0, 0.168d0,
     $         0.164d0, 0.134d0, 0.154d0, 0.128d0, 0.153d0, 0.145d0,
     $         
     $         0.0050d0, 0.020d0, 0.060d0, 0.056d0, 0.062d0, 0.062d0,
     $         0.0025d0, 0.001d0, 0.015d0,  0.07d0, 0.085d0, 0.011d0,
     $         0.0100d0, 0.002d0, 0.050d0,  0.08d0, 0.015d0, 0.057d0,
     $           0.08d0, 0.015d0,  0.09d0, 0.065d0, 0.042d0, 0.067d0,
     $          0.026d0,  0.03d0, 0.045d0, 0.052d0, 0.023d0, 0.051d0,
     $           0.02d0, 0.012d0, 0.098d0, 0.056d0, 0.024d0, 0.090d0,
     $         
     $         4.88d0, 4.870d0,	4.855d0, 4.834d0, 4.592d0, 4.834d0,
     $         4.85d0, 4.865d0, 4.845d0, 4.875d0, 4.815d0, 4.875d0,
     $         4.89d0, 4.870d0, 4.860d0, 4.826d0, 4.723d0, 4.826d0,
     $         4.83d0,  4.95d0,  4.62d0, 4.952d0, 4.852d0, 4.952d0,
     $         4.81d0, 4.758d0, 4.762d0,  4.95d0, 4.703d0, 4.95d0,
     $         4.98d0, 4.780d0, 4.608d0, 4.628d0, 4.237d0, 4.862d0
     $         /),
     $         (/6,6,ne/))
          allocate(bf_layer_used%nodes(nx,ny,ne))
          bf_layer_used%nodes = nodes

          interior_nodes = reshape(
     $         (/ (((100*(k-1) + 10*(j-1) + (i-1),i=1,nx),j=1,ny),k=1,ne) /),
     $         (/nx,ny,ne/))

          !bc_sections
          allocate(bf_layer_used%bc_sections(5,8))
          bf_layer_used%bc_sections = reshape((/
     $         SW_corner_type, 1           , 1           , no_overlap, no_overlap,
     $         S_edge_type   , bc_size+1   , 1           , nx-bc_size, no_overlap,
     $         SE_corner_type, nx-bc_size+1, 1           , no_overlap, no_overlap,
     $         W_edge_type   , 1           , bc_size+1   , ny-bc_size, no_overlap,
     $         E_edge_type   , nx-bc_size+1, bc_size+1   , ny-bc_size, no_overlap,
     $         NW_corner_type, 1           , ny-bc_size+1, no_overlap, no_overlap,
     $         N_edge_type   , bc_size+1   , ny-bc_size+1, nx-bc_size, no_overlap,
     $         NE_corner_type, nx-bc_size+1, ny-bc_size+1, no_overlap, no_overlap/),
     $         (/5,8/))

          !integration borders
          bf_layer_used%x_borders = [bc_size+1,nx-bc_size]
          bf_layer_used%y_borders = [bc_size+1,ny-bc_size]

          call bf_layer_used%allocate_before_timeInt()

          !output
          call bf_layer_used%compute_time_dev(
     $         td_operators_used,
     $         t, s, p_model, bc_used,
     $         interior_nodes)

          !validation
          timedev_test = td_operators_used%compute_time_dev(
     $         t, nodes, x_map, y_map,
     $         s,p_model,bc_used)

          test_validated = is_real_matrix3D_validated(
     $         bf_layer_used%bf_compute_used%time_dev,
     $         timedev_test,
     $         detailled)

        end function test_compute_time_dev


        function test_compute_integration_step(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_time)                              :: bf_layer_used
          type(td_operators)                               :: td_operators_used
          real(rkind)                                      :: t
          real(rkind)   , dimension(nx,ny,ne)              :: nodes
          real(rkind)   , dimension(nx)                    :: x_map
          real(rkind)   , dimension(ny)                    :: y_map
          type(sd_operators)                               :: s
          type(pmodel_eq)                                  :: p_model
          type(bc_operators)                               :: bc_used
          real(rkind)   , dimension(nx,ny,ne)              :: interior_nodes
          real(rkind)   , dimension(nx,ny,ne)              :: interior_nodes_tmp

          real(rkind)                         :: dt
          real(rkind)   , dimension(nx,ny,ne) :: interior_timedev

          integer(ikind) :: i,j
          logical        :: test_loc

          test_validated = .true.



          !input
          t  = 0.0d0
          dt = 0.1d0

          !alignment
          bf_layer_used%alignment = reshape((/
     $         align_E,
     $         align_S+1,
     $         align_E+nx-2*bc_size,
     $         align_S+1+ny-2*bc_size/),
     $         (/2,2/))

          !grdpts_id
          allocate(bf_layer_used%grdpts_id(nx,ny))
          bf_layer_used%grdpts_id = reshape(
     $         (/ ((interior_pt, i=1,nx),j=1,ny) /),
     $         (/nx,ny/))

          !x_map
          allocate(bf_layer_used%x_map(nx))
          x_map = (/ (x_min+(i-1)*0.1d0,i=1,nx) /)
          bf_layer_used%x_map = x_map

          !y_map
          allocate(bf_layer_used%y_map(ny))
          y_map = (/ (y_min+(j-1)*0.2d0,j=1,ny) /)
          bf_layer_used%y_map = (/ (y_min+(j-1)*0.2d0,j=1,ny) /)

          !nodes
          nodes = reshape((/
     $         1.48d0, 1.30d0, 1.35d0, 1.31d0, 1.43d0, 1.31d0,
     $         1.26d0, 1.45d0,  1.4d0, 1.29d0, 1.37d0, 1.41d0,
     $         1.46d0, 1.27d0, 1.47d0, 1.28d0, 1.25d0, 1.43d0,
     $         1.48d0, 1.26d0, 1.41d0, 1.34d0, 1.31d0, 1.39d0,
     $         1.42d0, 1.46d0, 1.38d0, 1.26d0, 1.37d0, 1.33d0,
     $         1.41d0, 1.22d0, 1.42d0, 1.23d0, 1.21d0, 1.40d0,
     $         
     $         0.128d0, 0.127d0, 0.142d0, 0.129d0, 0.136d0, 0.124d0,
     $         1.138d0, 0.148d0, 0.132d0, 0.125d0, 0.175d0, 0.123d0,
     $         0.146d0, 0.143d0, 0.145d0, 0.182d0, 0.135d0, 0.154d0,
     $         0.123d0, 0.129d0, 0.124d0, 0.162d0, 0.152d0, 0.142d0,
     $         0.168d0, 0.198d0, 0.186d0, 0.163d0, 0.126d0, 0.168d0,
     $         0.164d0, 0.134d0, 0.154d0, 0.128d0, 0.153d0, 0.145d0,
     $         
     $         0.0050d0, 0.020d0, 0.060d0, 0.056d0, 0.062d0, 0.062d0,
     $         0.0025d0, 0.001d0, 0.015d0,  0.07d0, 0.085d0, 0.011d0,
     $         0.0100d0, 0.002d0, 0.050d0,  0.08d0, 0.015d0, 0.057d0,
     $           0.08d0, 0.015d0,  0.09d0, 0.065d0, 0.042d0, 0.067d0,
     $          0.026d0,  0.03d0, 0.045d0, 0.052d0, 0.023d0, 0.051d0,
     $           0.02d0, 0.012d0, 0.098d0, 0.056d0, 0.024d0, 0.090d0,
     $         
     $         4.88d0, 4.870d0,	4.855d0, 4.834d0, 4.592d0, 4.834d0,
     $         4.85d0, 4.865d0, 4.845d0, 4.875d0, 4.815d0, 4.875d0,
     $         4.89d0, 4.870d0, 4.860d0, 4.826d0, 4.723d0, 4.826d0,
     $         4.83d0,  4.95d0,  4.62d0, 4.952d0, 4.852d0, 4.952d0,
     $         4.81d0, 4.758d0, 4.762d0,  4.95d0, 4.703d0, 4.95d0,
     $         4.98d0, 4.780d0, 4.608d0, 4.628d0, 4.237d0, 4.862d0
     $         /),
     $         (/6,6,ne/))
          allocate(bf_layer_used%nodes(nx,ny,ne))
          bf_layer_used%nodes = nodes

          interior_nodes = nodes

          !bc_sections
          allocate(bf_layer_used%bc_sections(5,8))
          bf_layer_used%bc_sections = reshape((/
     $         SW_corner_type, 1           , 1           , no_overlap, no_overlap,
     $         S_edge_type   , bc_size+1   , 1           , nx-bc_size, no_overlap,
     $         SE_corner_type, nx-bc_size+1, 1           , no_overlap, no_overlap,
     $         W_edge_type   , 1           , bc_size+1   , ny-bc_size, no_overlap,
     $         E_edge_type   , nx-bc_size+1, bc_size+1   , ny-bc_size, no_overlap,
     $         NW_corner_type, 1           , ny-bc_size+1, no_overlap, no_overlap,
     $         N_edge_type   , bc_size+1   , ny-bc_size+1, nx-bc_size, no_overlap,
     $         NE_corner_type, nx-bc_size+1, ny-bc_size+1, no_overlap, no_overlap/),
     $         (/5,8/))

          !integration borders
          bf_layer_used%x_borders = [bc_size+1,nx-bc_size]
          bf_layer_used%y_borders = [bc_size+1,ny-bc_size]

          call bf_layer_used%allocate_before_timeInt()

          call bf_layer_used%compute_time_dev(
     $         td_operators_used,
     $         t,s, p_model,bc_used,
     $         interior_nodes)

          !output
          call bf_layer_used%compute_integration_step(
     $         dt,
     $         compute_1st_step_nopt)

          !validation
          interior_timedev = td_operators_used%compute_time_dev(
     $         t, interior_nodes, x_map, y_map,
     $         s,p_model,bc_used)

          call compute_1st_step(
     $         interior_nodes,dt,
     $         interior_nodes_tmp,
     $         interior_timedev,
     $         [1,nx],
     $         [1,ny])

c$$$          print '(''timedev'')'
c$$$          do k=1,ne
c$$$             do j=1,6
c$$$                print '(6F20.14)', interior_timedev(1:6,7-j,k)
c$$$             end do
c$$$             print '()'
c$$$          end do
c$$$
c$$$          print '(''interior_nodes'')'
c$$$          do k=1,ne
c$$$             do j=1,6
c$$$                print '(6F20.14)', interior_nodes(1:6,7-j,k)
c$$$             end do
c$$$             print '()'
c$$$          end do

          test_loc = is_real_matrix3D_validated(
     $         bf_layer_used%bf_compute_used%nodes_tmp,
     $         interior_nodes_tmp,
     $         detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_real_matrix3D_validated(
     $         bf_layer_used%nodes,
     $         interior_nodes,
     $         detailled)
          test_validated = test_validated.and.test_loc          

        end function test_compute_integration_step


        function test_set_x_borders(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          

          type(bf_layer_time)          :: bf_layer_used
          integer(ikind), dimension(2) :: x_borders

          x_borders = [1,2]
          call bf_layer_used%set_x_borders(x_borders)

          test_validated = is_int_vector_validated(
     $         bf_layer_used%x_borders,
     $         x_borders,
     $         detailled)

        end function test_set_x_borders


        function test_get_x_borders(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          

          type(bf_layer_time)          :: bf_layer_used
          integer(ikind), dimension(2) :: x_borders
          integer(ikind), dimension(2) :: x_borders_test
          
          x_borders_test = [1,2]
          bf_layer_used%x_borders = x_borders_test
          x_borders = bf_layer_used%get_x_borders()

          test_validated = is_int_vector_validated(
     $         x_borders,
     $         x_borders_test,
     $         detailled)

        end function test_get_x_borders


        function test_set_y_borders(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          

          type(bf_layer_time)          :: bf_layer_used
          integer(ikind), dimension(2) :: y_borders

          y_borders = [1,2]
          call bf_layer_used%set_y_borders(y_borders)

          test_validated = is_int_vector_validated(
     $         bf_layer_used%y_borders,
     $         y_borders,
     $         detailled)

        end function test_set_y_borders


        function test_get_y_borders(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          

          type(bf_layer_time)          :: bf_layer_used
          integer(ikind), dimension(2) :: y_borders
          integer(ikind), dimension(2) :: y_borders_test
          
          y_borders_test = [1,2]
          bf_layer_used%y_borders = y_borders_test
          y_borders = bf_layer_used%get_y_borders()

          test_validated = is_int_vector_validated(
     $         y_borders,
     $         y_borders_test,
     $         detailled)

        end function test_get_y_borders


        subroutine check_inputs()

          implicit none

          logical :: test_parameters

          test_parameters = bc_choice.eq.hedstrom_xy_choice

          test_parameters = test_parameters.and.(bc_N_type_choice.eq.bc_timedev_choice)
          test_parameters = test_parameters.and.(bc_S_type_choice.eq.bc_timedev_choice)
          test_parameters = test_parameters.and.(bc_E_type_choice.eq.bc_timedev_choice)
          test_parameters = test_parameters.and.(bc_W_type_choice.eq.bc_timedev_choice)
          test_parameters = test_parameters.and.(nx.eq.6)
          test_parameters = test_parameters.and.(ny.eq.6)
          
        end subroutine check_inputs


      end program test_bf_layer_time
