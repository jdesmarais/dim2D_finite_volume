      !> @file
      !> module containing subroutines to read existing
      !> netcdf files
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module containing subroutines to read existing
      !> netcdf files
      !
      !> @date
      !> 05_12_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module nf90_operators_read_module

        use netcdf

        use nf90_operators_module, only :
     $       nf90_handle_err

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only : 
     $       pmodel_eq

        
        implicit none

        private
        public ::
     $       nf90_open_file_for_reading,
     $       nf90_get_varid,
     $       nf90_get_var_model,
     $       nf90_get_var_model_opt,
     $       nf90_read_borders


        contains


        !open a netcdf file for reading
        subroutine nf90_open_file_for_reading(filename,ncid)

          implicit none

          character*(*), intent(in)  :: filename
          integer      , intent(out) :: ncid

          integer :: retval

          retval = NF90_OPEN(trim(filename), NF90_NOWRITE, ncid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

        end subroutine nf90_open_file_for_reading


        !get the varid of the data variables based on their name
        subroutine nf90_get_varid(ncid,p_model,coordinates_id,data_id)

          implicit none

          integer               , intent(in)  :: ncid
          type(pmodel_eq)       , intent(in)  :: p_model
          integer, dimension(3) , intent(out) :: coordinates_id
          integer, dimension(ne), intent(out) :: data_id

          integer                          :: retval
          character(len=10), dimension(ne) :: name_var
          integer                          :: k


          !1) get the varid for the 'time' variable
          retval = NF90_INQ_VARID(ncid, "time",coordinates_id(1))
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !2) get the varid for the 'x' variable
          retval = NF90_INQ_VARID(ncid, "x",coordinates_id(2))
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !3) get the varid for the 'y' variable
          retval = NF90_INQ_VARID(ncid, "y",coordinates_id(3))
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !4) get the name of the governing variables
          !   for the file
          name_var = p_model%get_var_name()
          
          do k=1,ne
             
             !get the varid for the 'name_var(k)' variable
             retval = NF90_INQ_VARID(ncid, trim(name_var(k)),data_id(k))
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

          end do

        end subroutine nf90_get_varid


        !get the coordinates and the governing variable fields
        subroutine nf90_get_var_model(
     $     ncid,
     $     coordinates_id,
     $     data_id,
     $     time, nodes, x_map, y_map)

          implicit none

          integer                            , intent(in)  :: ncid
          integer    , dimension(3)          , intent(in)  :: coordinates_id
          integer    , dimension(ne)         , intent(in)  :: data_id
          real(rkind)                        , intent(out) :: time
          real(rkind), dimension(nx,ny,ne)   , intent(out) :: nodes
          real(rkind), dimension(nx)         , intent(out) :: x_map
          real(rkind), dimension(ny)         , intent(out) :: y_map

          real(rkind), dimension(1) :: time_table
          integer                   :: k
          integer                   :: retval


          !1) get the time
          retval = NF90_GET_VAR(ncid, coordinates_id(1), time_table)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          time = time_table(1)

          print '(''restart time: ok'')'


          !2) get the x_map
          retval = NF90_GET_VAR(ncid, coordinates_id(2), x_map)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          print '(''restart x_map: ok'')'


          !3) get the y_map
          retval = NF90_GET_VAR(ncid, coordinates_id(3), y_map)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          print '(''restart y_map: ok'')'


          !4) get the governing variables
          do k=1,ne

             retval = NF90_GET_VAR(
     $            ncid,
     $            data_id(k),
     $            nodes(:,:,k),
     $            START=[1,1,1],
     $            COUNT=[1,nx,ny])
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)
             
             print '(''restart data('',I2,''): ok'')', k

          end do

        end subroutine nf90_get_var_model


        !get the coordinates and the governing variable fields
        subroutine nf90_get_var_model_nopt(
     $     ncid,
     $     coordinates_id,
     $     data_id,
     $     time, nodes, x_map, y_map,
     $     COUNT)

          implicit none

          integer                      , intent(in)    :: ncid
          integer    , dimension(3)    , intent(in)    :: coordinates_id
          integer    , dimension(ne)   , intent(in)    :: data_id
          real(rkind)                  , intent(out)   :: time
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind), dimension(:)    , intent(inout) :: x_map
          real(rkind), dimension(:)    , intent(inout) :: y_map
          integer    , dimension(3)    , intent(in)    :: COUNT

          real(rkind), dimension(1) :: time_table
          integer                   :: k
          integer                   :: retval


          !1) get the time
          retval = NF90_GET_VAR(ncid, coordinates_id(1), time_table)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          time = time_table(1)

          print '(''restart time: ok'')'


          !2) get the x_map
          retval = NF90_GET_VAR(ncid, coordinates_id(2), x_map)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          print '(''restart x_map: ok'')'


          !3) get the y_map
          retval = NF90_GET_VAR(ncid, coordinates_id(3), y_map)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          print '(''restart y_map: ok'')'


          !4) get the governing variables
          do k=1,ne

             retval = NF90_GET_VAR(
     $            ncid,
     $            data_id(k),
     $            nodes(:,:,k),
     $            START=[1,1,1],
     $            COUNT=COUNT)
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)
             
             print '(''restart data('',I2,''): ok'')', k

          end do

        end subroutine nf90_get_var_model_nopt


        !read the borders of the x- and y- coordinates
        subroutine nf90_read_borders(
     $     ncid,
     $     coordinates_id, 
     $     x_borders,
     $     y_borders,
     $     sizes)

          implicit none

          integer                  , intent(in)  :: ncid
          integer, dimension(3)    , intent(in)  :: coordinates_id
          real(rkind), dimension(2), intent(out) :: x_borders
          real(rkind), dimension(2), intent(out) :: y_borders
          integer    , dimension(2), intent(out) :: sizes


          real(rkind), dimension(1)  :: x_min
          real(rkind), dimension(1)  :: x_max
          real(rkind), dimension(1)  :: y_min
          real(rkind), dimension(1)  :: y_max

          integer :: retval


          !get the extension of the x variable
          retval = nf90_inquire_dimension(
     $         ncid,
     $         coordinates_id(2),
     $         len=sizes(1))
          call nf90_handle_err(retval)
          
          !get the extension of the y variable
          retval = nf90_inquire_dimension(
     $         ncid,
     $         coordinates_id(3),
     $         len=sizes(2))
          call nf90_handle_err(retval)
          
          !get x_min
          retval = nf90_get_var(
     $         ncid,
     $         coordinates_id(2),
     $         x_min,
     $         START=[1],
     $         COUNT=[1])
          call nf90_handle_err(retval)
          
          !get x_max
          retval = nf90_get_var(
     $         ncid,
     $         coordinates_id(2),
     $         x_max,
     $         START=[sizes(1)],
     $         COUNT=[1])
          call nf90_handle_err(retval)
          
          !get y_min
          retval = nf90_get_var(
     $         ncid,
     $         coordinates_id(3),
     $         y_min,
     $         START=[1],
     $         COUNT=[1])
          call nf90_handle_err(retval)

          !get y_max
          retval = nf90_get_var(
     $         ncid,
     $         coordinates_id(3),
     $         y_max,
     $         START=[sizes(2)],
     $         COUNT=[1])
          call nf90_handle_err(retval)


          !x-borders
          x_borders(1) = x_min(1)
          x_borders(2) = x_max(1)

          !y-borders
          y_borders(1) = y_min(1)
          y_borders(2) = y_max(1)

        end subroutine nf90_read_borders

      end module nf90_operators_read_module
