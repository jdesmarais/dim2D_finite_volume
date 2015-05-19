      program test_file_merge_param

        use check_data_module, only :
     $     is_int_matrix_validated

        use file_merge_param_class, only :
     $     file_merge_param

        use parameters_kind, only :
     $       rkind


        implicit none

        
        logical :: test_loc
        logical :: test_validated
        logical :: detailled


        detailled = .true.
        test_validated = .true.


        test_loc = test_set_ranks(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_set_ranks: '',L1)', test_loc

        test_loc = test_determine_indices_for_extracting_and_saving(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_determine_indices_for_extracting_and_saving: '',L1)', test_loc

        contains


        function test_determine_indices_for_extracting_and_saving(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(file_merge_param)      :: fm_used
          real(rkind), dimension(2,2) :: borders_lg_domain
          real(rkind), dimension(2,2) :: borders_extraction
          real(rkind), dimension(2)   :: grid_spacings
          integer    , dimension(2)   :: nb_pts_per_tile
          integer                     :: bc_size

          integer    , dimension(2,2) :: extracted_from_file_test
          integer    , dimension(2,2) :: saved_in_lg_domain_test

          borders_lg_domain(1,1) = 1.0
          borders_lg_domain(1,2) = 31.0
          borders_lg_domain(2,1) = 0.0
          borders_lg_domain(2,2) = 39.0

          borders_extraction(1,1) = 19.0
          borders_extraction(1,2) = 25.0
          borders_extraction(2,1) = 24.0
          borders_extraction(2,2) = 27.0

          grid_spacings(1) = 2.0
          grid_spacings(2) = 3.0

          nb_pts_per_tile(1) = 8
          nb_pts_per_tile(2) = 9

          bc_size = 2

          extracted_from_file_test(1,1) = 2
          extracted_from_file_test(1,2) = 5
          extracted_from_file_test(2,1) = 4
          extracted_from_file_test(2,2) = 5

          saved_in_lg_domain_test(1,1) = 1
          saved_in_lg_domain_test(1,2) = 4
          saved_in_lg_domain_test(2,1) = 1
          saved_in_lg_domain_test(2,2) = 2


          call fm_used%set_ranks(2,1,2)

          call fm_used%determine_indices_for_extracting_and_saving(
     $         borders_lg_domain,
     $         borders_extraction,
     $         grid_spacings,
     $         nb_pts_per_tile,
     $         bc_size)


          test_validated = is_int_matrix_validated(
     $         fm_used%extracted_from_file,
     $         extracted_from_file_test,
     $         detailled)

          test_validated = test_validated.and.is_int_matrix_validated(
     $         fm_used%saved_in_lg_domain,
     $         saved_in_lg_domain_test,
     $         detailled)

        end function test_determine_indices_for_extracting_and_saving


        function test_set_ranks(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(file_merge_param) :: fm_used
          integer, dimension(3)  :: ranks_x
          integer, dimension(3)  :: ranks_y
          integer, dimension(3)  :: ranks_test
          integer                :: nb_tiles_y

          integer :: k
          logical :: test_loc


          test_validated = .true.


          ranks_x = [0,3,6]
          ranks_y = [2,3,5]
          nb_tiles_y = 4

          ranks_test = [2,15,29]


          do k=1,size(ranks_x,1)
             
             call fm_used%set_ranks(ranks_x(k),ranks_y(k),nb_tiles_y)
             test_loc = fm_used%rank.eq.ranks_test(k)
             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test '',I1,'' failed'')', k
                print '(''rank: '',I2,'' -> '', I2)',
     $               fm_used%rank,
     $               ranks_test(k)
             end if

          end do

        end function test_set_ranks

      end program test_file_merge_param
