      program test_bf_restart

        use bf_restart_module, only :
     $     get_restart_alignment,
     $     get_nb_detectors,
     $     read_detectors_from_file

        use netcdf

        use nf90_operators_module
        
        use nf90_operators_read_module

        use pmodel_eq_class, only :
     $       pmodel_eq

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       rkind

        character*(*), parameter :: interior_filename = 'data137.nc'
        character*(*), parameter :: bf_filename = 'E_1_137.nc'
        character*(*), parameter :: dct_filename = 'detectors136.curve'
                
        real(rkind), dimension(2)  :: x_borders
        real(rkind), dimension(2)  :: y_borders
        

        !test: nf90_read_borders
        call test_nf90_read_borders(
     $       bf_filename,
     $       x_borders,
     $       y_borders)

        !test: get_restart_alignment
        call test_get_restart_alignment(
     $       interior_filename,
     $       x_borders,
     $       y_borders)

        !test: get_nb_detectors
        call test_get_nb_detectors(dct_filename)

        !test: read_detectors_from_file
        call test_read_detectors_from_file(dct_filename)
        

        contains


        subroutine test_nf90_read_borders(
     $       bf_filename,
     $       x_borders,
     $       y_borders)

          implicit none

          character*(*)            , intent(in)  :: bf_filename
          real(rkind), dimension(2), intent(out) :: x_borders
          real(rkind), dimension(2), intent(out) :: y_borders
        
          integer                    :: ncid
          type(pmodel_eq)            :: p_model
          integer    , dimension(3)  :: coordinates_id
          integer    , dimension(ne) :: data_id
          integer    , dimension(2)  :: sizes

          !open file
          call nf90_open_file_for_reading(bf_filename,ncid)
          call nf90_get_varid(ncid, p_model, coordinates_id, data_id)
          
          !read borders
          call nf90_read_borders(
     $         ncid,
     $         coordinates_id, 
     $         x_borders,
     $         y_borders,
     $         sizes)

          !close file
          call nf90_close_file(ncid)
          
          !print borders
          print '(''test_nf90_read_borders'')'
          print '(''----------------------------------------'')'
          print *, 'size_x: ', sizes(1)
          print *, 'size_y: ', sizes(2)
          print *, 'x_min:  ', x_borders(1)
          print *, 'x_max:  ', x_borders(2)
          print *, 'y_min:  ', y_borders(1)
          print *, 'y_max:  ', y_borders(2)
          print '()'

       end subroutine test_nf90_read_borders


       subroutine test_get_restart_alignment(
     $     interior_filename,
     $     x_borders,
     $     y_borders)

         implicit none

         character*(*)              , intent(in) :: interior_filename
         real(rkind)  , dimension(2), intent(in) :: x_borders
         real(rkind)  , dimension(2), intent(in) :: y_borders

         type(pmodel_eq)                  :: p_model
         integer                          :: ncid
         integer    , dimension(3)        :: coordinates_id
         integer    , dimension(ne)       :: data_id
         real(rkind)                      :: time
         real(rkind), dimension(nx,ny,ne) :: interior_nodes
         real(rkind), dimension(nx)       :: interior_x_map
         real(rkind), dimension(ny)       :: interior_y_map

         integer, dimension(2,2) :: bf_alignment

         
         !open interior_filename
         call nf90_open_file_for_reading(
     $        interior_filename,
     $        ncid)

         !read data inside interior
         call nf90_get_varid(
     $        ncid,
     $        p_model,
     $        coordinates_id,
     $        data_id)

         call nf90_get_var_model(
     $        ncid,
     $        coordinates_id,
     $        data_id,
     $        time,
     $        interior_nodes,
     $        interior_x_map,
     $        interior_y_map)

         !close file
         call nf90_close_file(ncid)


         !get the alignment compared to the interior
         bf_alignment = get_restart_alignment(
     $        interior_x_map,
     $        interior_y_map,
     $        x_borders,
     $        y_borders)


         !print the alignment
         print '()'
         print '(''test_get_restart_alignment: '')'
         print '(''----------------------------------------'')'
         print *, 'bf_alignment(1,1): ', bf_alignment(1,1)
         print *, 'bf_alignment(1,2): ', bf_alignment(1,2)
         print *, 'bf_alignment(2,1): ', bf_alignment(2,1)
         print *, 'bf_alignment(2,2): ', bf_alignment(2,2)

         print *, 'interior_x_map(bf_alignment(1,1)-bc_size)',
     $        interior_x_map(bf_alignment(1,1)-bc_size)

         print *, 'interior_y_map(bf_alignment(2,1)-bc_size)',
     $        interior_x_map(bf_alignment(2,1)-bc_size)

         print *, 'interior_y_map(bf_alignment(2,2)+bc_size)',
     $        interior_y_map(bf_alignment(2,2)+bc_size)
         print '()'

         !get the alignment compared to the interior
         bf_alignment = get_restart_alignment(
     $        interior_x_map,
     $        interior_y_map,
     $        [-0.54d0, 0.54d0],
     $        [-0.54d0, 0.54d0])


         !print the alignment
         print *, 'bf_alignment(1,1): ', bf_alignment(1,1)
         print *, 'bf_alignment(1,2): ', bf_alignment(1,2)
         print *, 'bf_alignment(2,1): ', bf_alignment(2,1)
         print *, 'bf_alignment(2,2): ', bf_alignment(2,2)
         print '()'

       end subroutine test_get_restart_alignment


       subroutine test_get_nb_detectors(dct_filename)

         implicit none

         character*(*), intent(in) :: dct_filename

         integer, dimension(4) :: nb_detectors

         nb_detectors = get_nb_detectors(dct_filename)
         print '(''test_get_nb_detectors: '')'
         print '(''----------------------------------------'')'
         print *, 'N_detectors: ', nb_detectors(1)
         print *, 'S_detectors: ', nb_detectors(2)
         print *, 'E_detectors: ', nb_detectors(3)
         print *, 'W_detectors: ', nb_detectors(4)
         print '()'

       end subroutine test_get_nb_detectors


       subroutine test_read_detectors_from_file(dct_filename)

         implicit none

         character*(*), intent(in) :: dct_filename

         real(rkind), dimension(:,:), allocatable :: N_detectors
         real(rkind), dimension(:,:), allocatable :: S_detectors
         real(rkind), dimension(:,:), allocatable :: E_detectors
         real(rkind), dimension(:,:), allocatable :: W_detectors

         call read_detectors_from_file(
     $        dct_filename,
     $        N_detectors,
     $        S_detectors,
     $        E_detectors,
     $        W_detectors)

         print '(''test_read_detectors_from_file'')'
         print '(''----------------------------------------'')'
         print '(''N_detectors:'')'
         print '(2F8.4)', N_detectors(:,1)
         print '(''      ...'')'
         print '(2F8.4)', N_detectors(:,size(N_detectors,2))
         print '()'

         print '(''S_detectors:'')'
         print '(2F8.4)', S_detectors(:,1)
         print '(''      ...'')'
         print '(2F8.4)', S_detectors(:,size(S_detectors,2))
         print '()'

         print '(''E_detectors:'')'
         print '(2F8.4)', E_detectors(:,1)
         print '(''      ...'')'
         print '(2F8.4)', E_detectors(:,size(E_detectors,2))
         print '()'

         print '(''W_detectors:'')'
         print '(2F8.4)', W_detectors(:,1)
         print '(''      ...'')'
         print '(2F8.4)', W_detectors(:,size(W_detectors,2))
         print '()'

       end subroutine test_read_detectors_from_file

      end program test_bf_restart
