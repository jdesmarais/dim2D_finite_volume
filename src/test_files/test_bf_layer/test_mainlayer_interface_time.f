      program test_mainlayer_interface_time

        use bf_layer_class, only :
     $       bf_layer

        use check_data_module, only :
     $       is_int_matrix_validated

        use mainlayer_interface_time_class, only :
     $       mainlayer_interface_time        

        use parameters_bf_layer, only :
     $       align_N,
     $       align_W,
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       no_overlap,
     $       S_overlap

        use parameters_constant, only :
     $       N

        use parameters_kind, only :
     $       ikind

        implicit none

        
        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        test_loc = test_extract_grdpts_id_for_merge(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_extract_grdpts_id_for_merge: '',L1)', test_loc
        print '()'


        test_loc = test_update_bc_sections(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_update_bc_sections: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_loc


        contains


        function test_update_bc_sections(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mainlayer_interface_time)              :: mainlayer_interface_used
          type(bf_layer)                              :: bf_layer_used
          integer(ikind), dimension(:,:), allocatable :: test_bc_sections


          test_validated = .true.


          ! first test:
          ! we test whether the simple merge works
          !
          !
          !3 3 3 3 3 3     3 3 3 3 3 3    3 3 3 3 3 3     3 3 3 3 3 3  
          !3 2 2 2_2_3_____3_2_2 2 2 3    3 2 2 2 2_3_____3_2 2 2 2 3
          !3 2    |2 3|3 3|3 2|    2 3    3 2    |2 3 3 3 3 2|    2 3
          !3 2    |2_2|2_2|2_2|    2 3 -> 3 2    |2_2_2_2_2_2|    2 3
          !3 2                     2 3    3 2                     2 3
          !3 2                     2 3    3 2                     2 3
          !2 2                     2 2    2 2                     2 2
          !------------------------------------------------------------

          ! input
          bf_layer_used%alignment = reshape((/
     $         align_W+5, align_N, align_W+14, align_N+4/),
     $         (/2,2/))

          bf_layer_used%localization = N

          allocate(bf_layer_used%grdpts_id(14,9))

          bf_layer_used%grdpts_id = reshape((/
     $         1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     $         1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     $         2,2,1,1,1,1,1,1,1,1,1,1,2,2,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,2,2,2,2,2,2,1,1,2,3,
     $         3,2,1,1,2,3,3,3,3,2,1,1,2,3,
     $         3,2,2,2,2,3,0,0,3,2,2,2,2,3,
     $         3,3,3,3,3,3,0,0,3,3,3,3,3,3/),
     $         (/14,9/))

          bf_layer_used%x_borders = [1,14]
          bf_layer_used%y_borders = [3,9]

          call bf_layer_used%update_bc_sections()

          allocate(test_bc_sections(5,11))
          test_bc_sections = reshape((/
     $         NW_edge_type  , 1, 3, no_overlap, no_overlap,
     $         NE_edge_type  , 13,3, no_overlap, no_overlap,
     $         W_edge_type   , 1, 5, 7         , no_overlap,
     $         E_edge_type   , 13,5, 7         , no_overlap,
     $         N_edge_type   , 5, 6, 10        , no_overlap,
     $         NW_corner_type, 1, 8, no_overlap, no_overlap,
     $         N_edge_type   , 3, 8, 4         , no_overlap,
     $         NE_corner_type, 5, 8, no_overlap, no_overlap,
     $         NW_corner_type, 9, 8, no_overlap, no_overlap,
     $         N_edge_type   , 11,8, 12        , no_overlap,
     $         NE_corner_type, 13,8, no_overlap, no_overlap/),
     $         (/5,11/))


          ! output
          call mainlayer_interface_used%update_bc_sections(bf_layer_used)


          ! validation
          test_loc = is_int_matrix_validated(
     $         bf_layer_used%bc_sections,
     $         test_bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test 1 failed'')'
          end if


          ! second test:
          ! we test whether the simple merge works
          !
          !
          !3 3 3 3 3 3____ 3 3 3 3 3 3    3 3 3 3_3_3_____3_3_3 3 3 3  
          !3 2 2 2_2_3 3 3_3_2_2 2 2 3    3 2 2 2|2 3 3 3 3 2|2 2 2 3
          !3 2    |2_2_2_2_2_2|    2 3    3 2    |2_2_2_2_2_2|    2 3
          !3 2                     2 3 -> 3 2                     2 3
          !3 2                     2 3    3 2                     2 3
          !3 2                     2 3    3 2                     2 3
          !2 2                     2 2    2 2                     2 2
          !------------------------------------------------------------
          bf_layer_used%grdpts_id = reshape((/
     $         1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     $         1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     $         2,2,1,1,1,1,1,1,1,1,1,1,2,2,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,2,2,2,2,2,2,1,1,2,3,
     $         3,2,2,2,2,3,3,3,3,2,2,2,2,3,
     $         3,3,3,3,3,3,0,0,3,3,3,3,3,3/),
     $         (/14,9/))

          bf_layer_used%x_borders = [1,14]
          bf_layer_used%y_borders = [3,9]

          call bf_layer_used%update_bc_sections()

          test_bc_sections = reshape((/
     $         NW_edge_type  , 1, 3, no_overlap, no_overlap,
     $         NE_edge_type  , 13,3, no_overlap, no_overlap,
     $         W_edge_type   , 1, 5, 7         , no_overlap,
     $         E_edge_type   , 13,5, 7         , no_overlap,
     $         N_edge_type   , 5, 7, 10        , no_overlap,
     $         NW_corner_type, 1, 8, no_overlap, no_overlap,
     $         N_edge_type   , 3, 8, 4         , no_overlap,
     $         NE_corner_type, 5, 8, no_overlap, S_overlap,
     $         NW_corner_type, 9, 8, no_overlap, S_overlap,
     $         N_edge_type   , 11,8, 12        , no_overlap,
     $         NE_corner_type, 13,8, no_overlap, no_overlap/),
     $         (/5,11/))


          ! output
          call mainlayer_interface_used%update_bc_sections(bf_layer_used)


          ! validation
          test_loc = is_int_matrix_validated(
     $         bf_layer_used%bc_sections,
     $         test_bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test 2 failed'')'
          end if

        end function test_update_bc_sections


        function test_extract_grdpts_id_for_merge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mainlayer_interface_time) :: mainlayer_interface_used
          type(bf_layer)                 :: bf_layer_used

          integer(ikind), dimension(2,2) :: gen_coords
          integer       , dimension(3,4) :: grdpts_id
          integer       , dimension(3,4) :: test_grdpts_id


          ! input
          bf_layer_used%alignment = reshape((/
     $         align_W+5,align_N,align_W+6,align_N+1/),
     $         (/2,2/))

          allocate(bf_layer_used%grdpts_id(6,6))

          bf_layer_used%grdpts_id = reshape((/
     $         1,1,1,1,1,1,
     $         1,1,1,1,1,1,
     $         2,2,1,1,2,2,
     $         3,2,1,1,2,3,
     $         3,2,2,2,2,3,
     $         3,3,3,3,3,3/),
     $         (/6,6/))

          gen_coords = reshape((/
     $         align_W+2,align_N,align_W+4,align_N+3/),(/2,2/))

          test_grdpts_id = reshape((/
     $         2,2,2,
     $         3,3,2,
     $         0,3,2,
     $         0,3,3/),
     $         (/3,4/))

          ! output
          call mainlayer_interface_used%extract_grdpts_id_for_merge(
     $         bf_layer_used,
     $         gen_coords,
     $         grdpts_id)


          ! validation
          test_validated = is_int_matrix_validated(
     $         grdpts_id,
     $         test_grdpts_id,
     $         detailled)

        end function test_extract_grdpts_id_for_merge

      end program test_mainlayer_interface_time
