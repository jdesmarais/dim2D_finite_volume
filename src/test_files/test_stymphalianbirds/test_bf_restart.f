      program test_bf_restart

        use bf_restart_module, only :
     $     get_restart_alignment

        use netcdf

        use nf90_operators_module
        
        use nf90_operators_read_module

        use pmodel_eq_class, only :
     $       pmodel_eq

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       rkind

        character*(*), parameter   :: interior_filename = 'data137.nc'
        character*(*), parameter   :: bf_filename       = 'E_1_137.nc'
                
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
          
          print *, 'test_nf90_read_borders: 1'

          !read borders
          call nf90_read_borders(
     $         ncid,
     $         coordinates_id, 
     $         x_borders,
     $         y_borders,
     $         sizes)

          print *, 'test_nf90_read_borders: 2'

          !close file
          call nf90_close_file(ncid)
          
          !print borders
          print *, 'size_x: ', sizes(1)
          print *, 'size_y: ', sizes(2)
          print *, 'x_min:  ', x_borders(1)
          print *, 'x_max:  ', x_borders(2)
          print *, 'y_min:  ', y_borders(1)
          print *, 'y_max:  ', y_borders(2)

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

       end subroutine test_get_restart_alignment

      end program test_bf_restart
