      program test_nf90_operators_truncate

        use check_data_module, only :
     $     is_int_matrix_validated

        use netcdf

        use nf90_operators_truncate_module, only :
     $       nf90_truncate_files,
     $       nf90_put_var_truncate,
     $       nf90_def_var_truncate,
     $       nf90_get_extraction_indices,
     $       nf90_def_header_truncate,
     $       nf90_handle_err

        use parameters_cst, only :
     $       NOT_SUCCESS,
     $       SUCCESS

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        
        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        test_loc = test_nf90_def_header_truncate(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_nf90_def_header_truncate: '',L1)', test_loc
        print '()'


        test_loc = test_nf90_get_extraction_indices(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_nf90_get_extraction_indices: '',L1)', test_loc
        print '()'


        test_loc = test_nf90_def_var_truncate(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_nf90_def_var_truncate: '',L1)', test_loc
        print '()'


        test_loc = test_nf90_put_var_truncate(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_nf90_put_var_truncate: '',L1)', test_loc
        print '()'


        call test_nf90_truncate_files()


        print '(''test_validated: '',L1)', test_validated


        contains


        subroutine test_nf90_truncate_files()

          implicit none


          integer    , parameter      :: ne=4
          integer    , parameter      :: nb_timesteps=2
          real(rkind), dimension(2,2) :: borders_tr


          ! input
          borders_tr(1,1) = -1.5078d0
          borders_tr(1,2) =  1.7591d0
          borders_tr(2,1) = -1.8309d0
          borders_tr(2,2) =  1.8309d0


          ! output
          call nf90_truncate_files(
     $         ne,
     $         nb_timesteps,
     $         borders_tr)


          print '(''check the files:'')'
          print '(''   - ./tr_files/data0.nc'')'
          print '(''   - ./tr_files/data1.nc'')'
          print '()'
          

        end subroutine test_nf90_truncate_files


        function test_nf90_put_var_truncate(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          character*(*), parameter :: filename_lg = 'data0.nc'
          character*(*), parameter :: filename_tr = './test_files/data1.nc'

          integer      , dimension(2,2) :: borders_indices

          integer :: ncid_lg
          integer :: ncid_tr

          integer, parameter     :: ne=4
          integer, dimension(3)  :: coordinates_id
          integer, dimension(ne) :: var_id

          real(rkind), dimension(1)        :: time
          real(rkind), dimension(92)       :: x_map
          real(rkind), dimension(103)      :: y_map
          real(rkind), dimension(92,103,4) :: nodes
          

          integer :: ierror


          ! input
          call nf90_def_header_truncate(
     $         filename_lg,
     $         filename_tr,
     $         ncid_lg,
     $         ncid_tr,
     $         ierror=ierror)

          borders_indices(1,1) = 12
          borders_indices(1,2) = 103
          borders_indices(2,1) = 3
          borders_indices(2,2) = 105


          call nf90_def_var_truncate(
     $         ncid_lg,
     $         ncid_tr,
     $         borders_indices,
     $         coordinates_id,
     $         var_id)

          ! output
          call nf90_put_var_truncate(
     $         ncid_lg,
     $         ncid_tr,
     $         borders_indices,
     $         coordinates_id,
     $         var_id,
     $         time,
     $         x_map,
     $         y_map,
     $         nodes,
     $         ierror)
          

          call nf90_handle_err(
     $         NF90_CLOSE(ncid_lg))

          call nf90_handle_err(
     $         NF90_CLOSE(ncid_tr))

          test_validated = ierror.eq.SUCCESS
          if(detailled.and.(.not.test_validated)) then
             print '(''put_var failed'')'
          end if

        end function test_nf90_put_var_truncate


        function test_nf90_def_var_truncate(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          character*(*), parameter :: filename_lg = 'data0.nc'
          character*(*), parameter :: filename_tr = './test_files/data1.nc'

          integer      , dimension(2,2) :: borders_indices

          integer :: ncid_lg
          integer :: ncid_tr

          integer, parameter     :: ne=4
          integer, dimension(3)  :: coordinates_id
          integer, dimension(ne) :: var_id

          integer :: ierror


          ! input
          call nf90_def_header_truncate(
     $         filename_lg,
     $         filename_tr,
     $         ncid_lg,
     $         ncid_tr,
     $         ierror=ierror)

          borders_indices(1,1) = 12
          borders_indices(1,2) = 103
          borders_indices(2,1) = 3
          borders_indices(2,2) = 105


          ! output
          call nf90_def_var_truncate(
     $         ncid_lg,
     $         ncid_tr,
     $         borders_indices,
     $         coordinates_id,
     $         var_id)

          call nf90_handle_err(
     $         NF90_CLOSE(ncid_lg))

          call nf90_handle_err(
     $         NF90_CLOSE(ncid_tr))

          test_validated = ierror.eq.SUCCESS
          if(detailled.and.(.not.test_validated)) then
             print '(''def_var failed'')'
          end if


        end function test_nf90_def_var_truncate


        function test_nf90_get_extraction_indices(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          character*(*), parameter      :: filename_lg = 'data0.nc'
          real(rkind)  , dimension(2,2) :: borders_tr
          integer      , dimension(2,2) :: borders_indices
          integer      , dimension(2,2) :: test_borders_indices

          integer :: ierror
          integer :: ncid


          ! input
          call nf90_handle_err(
     $         NF90_OPEN(
     $         filename_lg,
     $         NF90_NOWRITE,
     $         ncid))

          borders_tr(1,1) = -1.5078d0
          borders_tr(1,2) =  1.7591d0
          borders_tr(2,1) = -1.8309d0
          borders_tr(2,2) =  1.8309d0
          
          test_borders_indices(1,1) = 12
          test_borders_indices(1,2) = 103
          test_borders_indices(2,1) = 3
          test_borders_indices(2,2) = 105


          ! output
          call nf90_get_extraction_indices(
     $         ncid,
     $         borders_tr,
     $         borders_indices,
     $         ierror)

          call nf90_handle_err(NF90_CLOSE(ncid))


          ! validation
          test_validated = is_int_matrix_validated(
     $         borders_indices,
     $         test_borders_indices,
     $         detailled)

        end function test_nf90_get_extraction_indices


        function test_nf90_def_header_truncate(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          character*(*), parameter :: filename_lg = 'data0.nc'
          character*(*), parameter :: filename_tr = './test_files/data1.nc'
          integer :: ncid_lg
          integer :: ncid_tr
          integer :: ierror


          ! output
          call nf90_def_header_truncate(
     $         filename_lg,
     $         filename_tr,
     $         ncid_lg,
     $         ncid_tr,
     $         ierror=ierror)

          call nf90_handle_err(NF90_CLOSE(ncid_lg))
          call nf90_handle_err(NF90_CLOSE(ncid_tr))

          ! validation
          test_validated = ierror.eq.SUCCESS
          if(detailled.and.(.not.test_validated)) then
             print '(''def header failed'')'
          end if

        end function test_nf90_def_header_truncate

      end program test_nf90_operators_truncate
