      program test_nf90_operators_mt

        use check_data_module, only :
     $     is_real_vector_validated,
     $     is_real_matrix_validated,
     $     is_int_vector_validated,
     $     is_int_matrix_validated

        use nf90_operators_mt_module, only :
     $       nf90_get_coordinate_borders,
     $       nf90_get_lg_domain_param,
     $       get_rank_tr_corner_files

        use parameters_kind, only :
     $       rkind

        implicit none

        
        logical :: test_loc
        logical :: test_validated
        logical :: detailled

        detailled = .true.
        test_validated = .true.


        test_loc = test_nf90_get_coordinate_borders(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_nf90_get_coordinate_borders: '',L1)', test_loc


        test_loc = test_nf90_get_lg_domain_param(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_nf90_get_lg_domain_param: '',L1)', test_loc


        test_loc = test_get_rank_tr_corner_files(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_rank_tr_corner_files: '',L1)', test_loc


        contains


        function test_get_rank_tr_corner_files(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(2,2) :: borders_lg_domain
          real(rkind), dimension(2,2) :: borders_tr_domain
          real(rkind), dimension(2)   :: map_spacing
          integer    , dimension(2)   :: nb_tiles
          integer    , dimension(2)   :: nb_pts_per_tile
          integer                     :: bc_size
          integer    , dimension(2,2) :: corner_ranks
          integer    , dimension(2,2) :: corner_ranks_test


          borders_lg_domain = reshape((/1.0d0,0.0d0,31.0d0,39.0d0/),(/2,2/))
          !borders_tr_domain = reshape((/19.0d0,24.0d0,25.0d0,27.0d0/),(/2,2/))
          !corner_ranks_test = reshape((/1,1,2,1/),(/2,2/))

          !borders_tr_domain = reshape((/5.0d0,6.0d0,11.0d0,18.0d0/),(/2,2/))
          !corner_ranks_test = reshape((/0,0,0,0/),(/2,2/))

          borders_tr_domain = reshape((/19.0d0,6.0d0,21.0d0,21.0d0/),(/2,2/))
          corner_ranks_test = reshape((/1,0,2,1/),(/2,2/))

          map_spacing       = [2.0d0,3.0d0]
          nb_tiles          = [3,2]
          nb_pts_per_tile   = [8,9]
          bc_size           = 2

          

          call get_rank_tr_corner_files(
     $         borders_lg_domain,
     $         borders_tr_domain,
     $         map_spacing,
     $         nb_tiles,
     $         nb_pts_per_tile,
     $         bc_size,
     $         corner_ranks)

          test_validated = is_int_matrix_validated(
     $         corner_ranks,
     $         corner_ranks_test,
     $         detailled)

        end function test_get_rank_tr_corner_files


        function test_nf90_get_lg_domain_param(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          character*(*), parameter :: SW_corner_filename = 'data0_0.nc'
          character*(*), parameter :: NE_corner_filename = 'data0_63.nc'

          real(rkind), dimension(2,2) :: borders_lg_domain
          real(rkind), dimension(2)   :: map_spacing
          integer    , dimension(2)   :: nb_pts_per_tile

          real(rkind), dimension(2,2) :: borders_lg_domain_test
          real(rkind), dimension(2)   :: map_spacing_test
          integer    , dimension(2)   :: nb_pts_per_tile_test


          test_validated = .true.

          borders_lg_domain_test = reshape((/
     $         -5.53455d0,-5.53455d0,
     $          5.53455d0, 5.53455d0/),
     $         (/2,2/))

          map_spacing_test(1) = -5.52965d0+5.53455d0
          map_spacing_test(2) = -5.52965d0+5.53455d0             

          nb_pts_per_tile_test = [286,286]


          call nf90_get_lg_domain_param(
     $         SW_corner_filename,
     $         NE_corner_filename,
     $         borders_lg_domain,
     $         map_spacing,
     $         nb_pts_per_tile)


          test_loc = is_real_matrix_validated(
     $         borders_lg_domain,
     $         borders_lg_domain_test,
     $         detailled)
          test_validated = test_validated.and.test_loc

          if(detailled.and.(.not.test_loc)) then
             print '(''borders_lg_domain failed'')'
          end if

          test_loc = is_real_vector_validated(
     $         map_spacing,
     $         map_spacing_test,
     $         detailled)
          test_validated = test_validated.and.test_loc

          if(detailled.and.(.not.test_loc)) then
             print '(''map_spacing failed'')'
          end if

          test_loc = is_int_vector_validated(
     $         nb_pts_per_tile,
     $         nb_pts_per_tile_test,
     $         detailled)
          test_validated = test_validated.and.test_loc

          if(detailled.and.(.not.test_loc)) then
             print '(''nb_pts_per_tile failed'')'
          end if

        end function test_nf90_get_lg_domain_param


        function test_nf90_get_coordinate_borders(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          character*(*), parameter :: filename = 'data0_1.nc'

          real(rkind), dimension(2,2) :: map_borders
          real(rkind), dimension(2)   :: map_spacing

          real(rkind), dimension(2,2) :: map_borders_test
          real(rkind), dimension(2)   :: map_spacing_test

          map_borders_test(1,:) = [-5.53455d0,-4.13805d0]
          map_borders_test(2,:) = [-4.15275d0,-2.75625d0]

          map_spacing_test(1) = -4.14785d0+4.15275d0
          map_spacing_test(2) = -5.52965d0+5.53455d0             


          call nf90_get_coordinate_borders(
     $         filename,
     $         map_borders, map_spacing)

          test_validated = is_real_matrix_validated(
     $            map_borders,
     $            map_borders_test,
     $            detailled)

          test_validated = test_validated.and.
     $         is_real_vector_validated(
     $         map_spacing,
     $         map_spacing_test,
     $         detailled)

       end function test_nf90_get_coordinate_borders

      end program test_nf90_operators_mt
