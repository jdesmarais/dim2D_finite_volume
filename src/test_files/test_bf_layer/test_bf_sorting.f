      program test_bf_sorting

        use bf_sorting_module, only :
     $       bubble_sort_grdpts

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_kind, only :
     $       ikind

        implicit none

        
        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        test_loc = test_bubble_sort_grdpts(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_bubble_sort_grdpts: '',L1)', test_loc


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_bubble_sort_grdpts(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(2,5) :: grdpts_in
          integer, dimension(2,5) :: grdpts_test


          !input
          grdpts_in = reshape((/
     $         2,1,
     $         4,2,
     $         2,-1,
     $         1,-4,
     $         0,2/),
     $         (/2,5/))

          grdpts_test = reshape((/
     $         1,-4,
     $         2,-1,
     $         2,1,
     $         0,2,
     $         4,2/),
     $         (/2,5/))

          !output
          call bubble_sort_grdpts(grdpts_in)

          !validation
          test_validated = is_int_matrix_validated(
     $         grdpts_in,
     $         grdpts_test,
     $         detailled)

        end function test_bubble_sort_grdpts

      end program test_bf_sorting
