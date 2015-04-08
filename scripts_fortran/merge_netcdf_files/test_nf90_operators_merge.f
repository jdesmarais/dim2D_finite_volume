      program test_nf90_operators_merge

        use netcdf

        use nf90_operators_merge_module, only :
     $       nf90_merge_files,
     $       nf90_merge_files_at_timestep,
     $       nf90_put_var_merge,
     $       nf90_def_var_merge,
     $       nf90_def_header_merge,
     $       nf90_handle_err,
     $       get_filename

        use parameters_cst, only :
     $       SUCCESS

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        logical :: test_loc
        logical :: test_validated
        logical :: detailled


        detailled = .true.
        test_validated = .true.

        
        test_loc = test_nf90_def_header_merge(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_nf90_def_header_merge: '',L1)', test_loc
        print ''


        test_loc = test_nf90_def_var_merge(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_nf90_def_var_merge: '',L1)', test_loc
        print ''


        test_loc = test_nf90_put_var_merge(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_nf90_put_var_merge: '',L1)', test_loc
        print ''

        
        test_loc = test_nf90_merge_files_at_timestep(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_nf90_merge_files_at_timestep: '',L1)', test_loc
        print ''


        test_loc = test_nf90_merge_files()
        test_validated = test_validated.and.test_loc
        print '(''test_nf90_merge_files: '',L1)', test_loc
        print ''


        print '(''test_validated: '',L1)', test_validated

        contains

        function test_nf90_def_header_merge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          character*(*), parameter :: filename_rank0 = 'def_header_rank0.nc'
          character*(*), parameter :: filename_merge = 'def_header.nc'
          integer                  :: ncid_merge
          integer                  :: ierror
          

          ! output
          call nf90_def_header_merge(
     $         filename_rank0,
     $         filename_merge,
     $         ncid_merge,
     $         ierror=ierror)
          
          call nf90_handle_err(NF90_CLOSE(ncid_merge))


          ! validation
          test_validated = ierror.eq.SUCCESS
          if(detailled.and.(.not.test_validated)) then
             print '(''nf90_def_header_merge: failed'')'
          end if

        end function test_nf90_def_header_merge


        function test_nf90_def_var_merge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          character*(*), parameter :: filename_rank0 = 'def_var_rank0.nc'
          character*(*), parameter :: filename_merge = 'def_var.nc'
          integer                  :: ncid_merge
          integer                  :: ierror

          integer, parameter    :: ne=4
          integer, parameter    :: bc_size=2
          integer, dimension(3) :: coordinates_id
          integer, dimension(ne):: var_id
                    

          ! input
          call nf90_def_header_merge(
     $         filename_rank0,
     $         filename_merge,
     $         ncid_merge,
     $         ierror=ierror)

          ! output
          call nf90_def_var_merge(
     $         filename_rank0,
     $         ncid_merge,
     $         bc_size,[8,8],
     $         coordinates_id,
     $         var_id,
     $         ierror=ierror)
          
          call nf90_handle_err(
     $         NF90_CLOSE(ncid_merge))


          ! validation
          test_validated = ierror.eq.SUCCESS
          if(detailled.and.(.not.test_validated)) then
             print '(''nf90_def_var_merge: failed'')'
          end if

        end function test_nf90_def_var_merge


        function test_nf90_put_var_merge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          character*(*), parameter :: filename_rank0 = 'data0_0.nc'
          character*(*), parameter :: filename_merge = 'put_var.nc'
          integer                  :: ncid_merge
          integer                  :: ierror

          integer, parameter    :: ne=4
          integer, parameter    :: bc_size=2
          integer, dimension(3) :: coordinates_id
          integer, dimension(ne):: var_id
                    
          integer, dimension(2), parameter :: nb_tiles = [8,8]
          integer              , parameter :: timestep=0
          character(len=16)                :: filename_rank

          real(rkind), dimension(1)                  :: time
          real(rkind), dimension(:)    , allocatable :: x_map
          real(rkind), dimension(:)    , allocatable :: y_map
          real(rkind), dimension(:,:,:), allocatable :: nodes

          integer :: size_x
          integer :: size_y
          integer :: k


          ! input
          call nf90_def_header_merge(
     $         filename_rank0,
     $         filename_merge,
     $         ncid_merge,
     $         ierror=ierror)

          call nf90_def_var_merge(
     $         filename_rank0,
     $         ncid_merge,
     $         bc_size,nb_tiles,
     $         coordinates_id,
     $         var_id,
     $         size_x=size_x,
     $         size_y=size_y,
     $         ierror=ierror)

          allocate(x_map(size_x))
          allocate(y_map(size_y))
          allocate(nodes(size_x,size_y,ne))

          ! output
          do k=0, nb_tiles(1)*nb_tiles(2)-1

             call get_filename(
     $            filename_rank,
     $            timestep,
     $            k)

             call nf90_put_var_merge(
     $            trim(filename_rank),
     $            k,
     $            ncid_merge,
     $            bc_size,
     $            nb_tiles,
     $            coordinates_id,
     $            var_id,
     $            time,
     $            x_map,
     $            y_map,
     $            nodes,
     $            ierror=ierror)

          end do

          deallocate(x_map)
          deallocate(y_map)
          deallocate(nodes)
          
          call nf90_handle_err(NF90_CLOSE(ncid_merge))


          ! validation
          test_validated = ierror.eq.SUCCESS
          if(detailled.and.(.not.test_validated)) then
             print '(''nf90_put_var_merge: failed'')'
          end if

        end function test_nf90_put_var_merge


        function test_nf90_merge_files_at_timestep(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          character*(*), parameter :: filename_rank0 = 'data0_0.nc'
          character*(*), parameter :: filename_merge = 'data0.nc'
          integer                  :: ncid_merge
          integer                  :: ierror

          integer, parameter    :: ne=4
          integer, parameter    :: bc_size=2
          integer, dimension(3) :: coordinates_id
          integer, dimension(ne):: var_id
                    
          integer, dimension(2), parameter :: nb_tiles = [8,8]
          integer              , parameter :: timestep=0

          real(rkind), dimension(1)                  :: time
          real(rkind), dimension(:)    , allocatable :: x_map
          real(rkind), dimension(:)    , allocatable :: y_map
          real(rkind), dimension(:,:,:), allocatable :: nodes

          integer :: size_x
          integer :: size_y


          ! input
          call nf90_def_header_merge(
     $         filename_rank0,
     $         filename_merge,
     $         ncid_merge,
     $         ierror=ierror)

          call nf90_def_var_merge(
     $         filename_rank0,
     $         ncid_merge,
     $         bc_size,nb_tiles,
     $         coordinates_id,
     $         var_id,
     $         size_x=size_x,
     $         size_y=size_y,
     $         ierror=ierror)

          call nf90_handle_err(NF90_CLOSE(ncid_merge))

          allocate(x_map(size_x))
          allocate(y_map(size_y))
          allocate(nodes(size_x,size_y,ne))


          ! output
          call nf90_merge_files_at_timestep(
     $         nb_tiles,
     $         ne,
     $         bc_size,
     $         timestep,
     $         time,
     $         x_map,
     $         y_map,
     $         nodes)

          deallocate(x_map)
          deallocate(y_map)
          deallocate(nodes)


          ! validation
          test_validated = ierror.eq.SUCCESS
          if(detailled.and.(.not.test_validated)) then
             print '(''nf90_merge_files_at_timestep: failed'')'
          end if

        end function test_nf90_merge_files_at_timestep


        function test_nf90_merge_files()
     $       result(test_validated)

          implicit none

          logical             :: test_validated

          integer              , parameter :: ne=4
          integer              , parameter :: bc_size=2                    
          integer, dimension(2), parameter :: nb_tiles = [8,8]
          integer              , parameter :: nb_timesteps=2


          call nf90_merge_files(
     $         nb_tiles,
     $         ne,
     $         bc_size,
     $         nb_timesteps)

          test_validated = .true.

        end function test_nf90_merge_files

      end program test_nf90_operators_merge
