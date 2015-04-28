      program test_mainlayer_interface_icr

        use bf_layer_class, only :
     $       bf_layer

        use mainlayer_interface_icr_class, only :
     $       mainlayer_interface_icr

        use parameters_bf_layer, only :
     $       align_N,
     $       align_W,
     $       N_edge_type,
     $       no_overlap

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


        test_loc = test_analyze_bc_section_edge(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_bc_section_edge: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated


        contains


        function test_analyze_bc_section_edge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mainlayer_interface_icr) :: mainlayer_interface_used
          type(bf_layer)                :: bf_layer_used
          integer(ikind), dimension(5)  :: bc_section

          logical :: test_loc


          test_validated = .true.


          ! first test:
          ! we test whether a long edge is not analyzed
          !
          !    3 3 3 3 3 3 3 3 3 3 3 3 3 3
          !    3 2 2 2 2 2 2 2 2 2 2 2 2 3
          !    3 2                     2 3
          !    3 2                     2 3
          !    3 2                     2 3
          !    3 2                     2 3
          !    2 2                     2 2
          !--------------------------------------------
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
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,2,2,2,2,2,2,2,2,2,2,2,3,
     $         3,3,3,3,3,3,3,3,3,3,3,3,3,3/),
     $         (/14,9/))

          bf_layer_used%x_borders = [1,14]
          bf_layer_used%y_borders = [3,9]

          bc_section = [N_edge_type,3,6,12,no_overlap]

          test_loc = .not.(
     $         mainlayer_interface_used%analyze_bc_section_edge(
     $         bc_section,
     $         bf_layer_used))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test 1 failed'')'
          end if

          
          ! second test:
          ! we test whether an edge with no anti-
          ! corners of type (merge+overlap) is
          ! present at the sides of an edge-like
          ! bc_section
          !
          !  3 3 3 3 3 3     3 3 3 3 3 3  
          !  3 2 2 2_2_3_____3_2_2 2 2 3
          !  3 2    |2 3|3 3|3 2|    2 3
          !  3 2    |2_2|2_2|2_2|    2 3
          !  3 2                     2 3
          !  3 2                     2 3
          !  2 2                     2 2
          !--------------------------------------------
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

          call bf_layer_used%update_bc_sections()

          bc_section = [N_edge_type,7,6,8,no_overlap]
          
          test_loc = mainlayer_interface_used%analyze_bc_section_edge(
     $         bc_section,
     $         bf_layer_used)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test 2 failed'')'
          end if


          ! third test:
          ! we test whether an edge with only one
          ! anti-corner of type (merge+overlap) is
          ! present at the sides of an edge-like
          ! bc_section is adapted
          !
          !  3 3 3 3 3 3     
          !  3 2 2 2_2_3_____3_3_3 3 3 3
          !  3 2    |2 3|3 3|3 2|2 2 2 3
          !  3 2    |2_2|2_2|2_2|    2 3
          !  3 2                     2 3
          !  3 2                     2 3
          !  2 2                     2 2
          !--------------------------------------------
          bf_layer_used%grdpts_id = reshape((/
     $         1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     $         1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     $         2,2,1,1,1,1,1,1,1,1,1,1,2,2,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,2,2,2,2,2,2,1,1,2,3,
     $         3,2,1,1,2,3,3,3,3,2,2,2,2,3,
     $         3,2,2,2,2,3,0,0,3,3,3,3,3,3,
     $         3,3,3,3,3,3,0,0,0,0,0,0,0,0/),
     $         (/14,9/))

          call bf_layer_used%update_bc_sections()
          
          test_loc = mainlayer_interface_used%analyze_bc_section_edge(
     $         bc_section,
     $         bf_layer_used)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test 3 failed'')'
          end if


          ! fourth test:
          ! we test whether an edge with only one
          ! anti-corner of type (merge+overlap) is
          ! present at the sides of an edge-like
          ! bc_section is adapted
          !
          !                  3 3 3 3 3 3
          !  3 3 3 3 3_3_____3_2_2 2 2 3
          !  3 2 2 2|2 3|3 3|3 2|    2 3
          !  3 2    |2_2|2_2|2_2|    2 3
          !  3 2                     2 3
          !  3 2                     2 3
          !  2 2                     2 2
          !--------------------------------------------
          bf_layer_used%grdpts_id = reshape((/
     $         1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     $         1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     $         2,2,1,1,1,1,1,1,1,1,1,1,2,2,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,2,2,2,2,2,2,1,1,2,3,
     $         3,2,2,2,2,3,3,3,3,2,1,1,2,3,
     $         3,3,3,3,3,3,0,0,3,2,2,2,2,3,
     $         0,0,0,0,0,0,0,0,3,3,3,3,3,3/),
     $         (/14,9/))

          call bf_layer_used%update_bc_sections()
          
          test_loc =  mainlayer_interface_used%analyze_bc_section_edge(
     $         bc_section,
     $         bf_layer_used)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test 4 failed'')'
          end if
          

          ! fifth test:
          ! we test whether an edge with two
          ! anti-cornersof type (merge+overlap) is
          ! present at the sides of an edge-like
          ! bc_section is adapted
          !
          !
          !  3 3 3 3 3_3_____3_3 3 3 3 3
          !  3 2 2 2|2 3|3 3|3 2|2 2 2 3
          !  3 2    |2_2|2_2|2_2|    2 3
          !  3 2                     2 3
          !  3 2                     2 3
          !  2 2                     2 2
          !--------------------------------------------
          bf_layer_used%grdpts_id = reshape((/
     $         1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     $         1,1,1,1,1,1,1,1,1,1,1,1,1,1,
     $         2,2,1,1,1,1,1,1,1,1,1,1,2,2,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,2,2,2,2,2,2,1,1,2,3,
     $         3,2,2,2,2,3,3,3,3,2,2,2,2,3,
     $         3,3,3,3,3,3,0,0,3,3,3,3,3,3,
     $         0,0,0,0,0,0,0,0,0,0,0,0,0,0/),
     $         (/14,9/))

          call bf_layer_used%update_bc_sections()
          
          test_loc = mainlayer_interface_used%analyze_bc_section_edge(
     $         bc_section,
     $         bf_layer_used)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test 5 failed'')'
          end if

        end function test_analyze_bc_section_edge

      end program test_mainlayer_interface_icr
