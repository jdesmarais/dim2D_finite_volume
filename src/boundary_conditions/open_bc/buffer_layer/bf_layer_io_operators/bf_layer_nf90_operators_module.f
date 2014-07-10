      !> @file
      !> module encapsulating the subroutines to write netcdf files
      !> from the buffer layer data
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the subroutines to write netcdf files
      !> from the buffer layer data
      !
      !> @date
      ! 10_07_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_nf90_operators_module

        use bf_layer_errors_module, only : error_mainlayer_id
        use dim2d_eq_class        , only : dim2d_eq
        use netcdf                
        use parameters_bf_layer   , only : no_pt
        use parameters_constant   , only : institut,
     $                                     prog_version,
     $                                     ref,
     $                                     convention,
     $                                     N,S,E,W
        use parameters_input      , only : ne
        use parameters_kind       , only : ikind, rkind

        implicit none

        private
        public :: print_bf_layer_on_netcdf
      
        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the grdpts_id and the nodes, data typical from
        !> a buffer layer on an external netcdf file
        !
        !> @date
        !> 10_07_2014 - initial version - J.L. Desmarais
        !
        !>@param filename
        !> file name for the netcdf output
        !
        !>@param grdpts_id
        !> array where the role of each grid point is stored (no_pt,
        !> bc_pt, bc_interior_pt, interior_pt)
        !
        !> @param nodes
        !> array where the governing variables are saved at each grid
        !> point
        !
        !> @param t
        !> time corresponding to the data for the grdpts and the nodes
        !--------------------------------------------------------------
        subroutine print_bf_layer_on_netcdf(
     $       filename, p_model,
     $       bf_loc, bf_order,
     $       x_start, y_start, dx, dy,
     $       grdpts_id, nodes, time)

          implicit none

          character(*)                  , intent(in)    :: filename
          type(dim2d_eq)                , intent(in)    :: p_model
          integer                       , intent(in)    :: bf_loc
          integer                       , intent(in)    :: bf_order
          real(rkind)                   , intent(in)    :: x_start
          real(rkind)                   , intent(in)    :: y_start
          real(rkind)                   , intent(in)    :: dx
          real(rkind)                   , intent(in)    :: dy
          integer     , dimension(:,:)  , intent(in)    :: grdpts_id
          real(rkind) , dimension(:,:,:), intent(inout) :: nodes
          real(rkind)                   , intent(in)    :: time

          
          real(rkind) :: NF90_FILL_MYREAL
          integer     :: ncid

          integer, dimension(3)  :: coords_id
          integer                :: grdptsid_id
          integer, dimension(ne) :: nodes_id


          !determination of the convention for the
          !missing value
          NF90_FILL_MYREAL = set_nf90_fill_myreal()
          

          !set the gridpoints corresponding to no_pt
          !with a missing_value such that they are
          !recognized by the visualization software
          call set_missing_data(grdpts_id, nodes, NF90_FILL_MYREAL)
          

          !open the netcdf file
          call nf90_open_file(filename, ncid)
          

          !define the header of the file
          call nf90_write_header(
     $         ncid,
     $         get_bf_layer_title(bf_loc, bf_order))


          !define the variables saved in the file
          call bf_layer_nf90_def_var(
     $         ncid,
     $         p_model, size(grdpts_id,1), size(grdpts_id,2),
     $         NF90_FILL_MYREAL,
     $         coords_id, grdptsid_id, nodes_id)
          

          !save the variables in the file
          call bf_layer_nf90_put_var(
     $         ncid,
     $         coords_id, grdptsid_id, nodes_id,
     $         x_start, y_start, dx, dy, 
     $         grdpts_id, nodes, time)


          !close the file
          call nf90_close_file(ncid)

        end subroutine print_bf_layer_on_netcdf


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the default fill (missing_value) for the variables of
        !> type real(rkind) written on netcdf outputs
        !
        !> @date
        !> 10_07_2014 - initial version - J.L. Desmarais
        !
        !> @result NF90_FILL_MYREAL
        !> default fill (missing value) for the variables of type
        !> real(rkind) written on netcdf outputs
        !--------------------------------------------------------------
        function set_nf90_fill_myreal() result(NF90_FILL_MYREAL)

          implicit none

          real(rkind) :: NF90_FILL_MYREAL

          select case(rkind)
            case(4)
               NF90_FILL_MYREAL = NF90_FILL_FLOAT
            case(8)
               NF90_FILL_MYREAL = NF90_FILL_DOUBLE
          end select
          
        end function set_nf90_fill_myreal


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set as missing values in the nodes array the grid points
        !> corresponding to no_pt in the grdpts_id array
        !
        !> @date
        !> 10_07_2014 - initial version - J.L. Desmarais
        !
        !>@param grdpts_id
        !> array where the role of each grid point is stored (no_pt,
        !> bc_pt, bc_interior_pt, interior_pt)
        !
        !> @param nodes
        !> array where the governing variables are saved at each grid
        !> point
        !
        !> @param missing_value
        !> default value used to identify the missing variables to type
        !> real(rkind)
        !--------------------------------------------------------------
        subroutine set_missing_data(grdpts_id, nodes, missing_value)

          implicit none

          integer     , dimension(:,:)  , intent(in)    :: grdpts_id
          real(rkind) , dimension(:,:,:), intent(inout) :: nodes
          real(rkind)                   , intent(in)    :: missing_value

          integer(ikind) :: i,j
          integer        :: k
          
          !replacing the missing values in nodes by
          !the convention
          do j=1, size(grdpts_id,2)
             do i=1, size(grdpts_id,1)
                if(grdpts_id(i,j).eq.no_pt) then
                   do k=1, ne
                      nodes(i,j,k) = missing_value
                   end do
                end if
             end do
          end do

        end subroutine set_missing_data

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> handle the netcdf errors when working with a netcdf file
        !
        !> @date
        !> 10_07_2014 - initial version - J.L. Desmarais
        !
        !>@param nf90_err
        !> integer identifying the netcdf error caught
        !--------------------------------------------------------------
        subroutine nf90_handle_err(nf90_err)

          implicit none

          integer, intent(in) ::  nf90_err
        
          if (nf90_err.ne.NF90_NOERR) then
             print *, 'Error: ', NF90_STRERROR(nf90_err)
             stop 2
          end if

        end subroutine nf90_handle_err


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> open a netcdf file
        !
        !> @date
        !> 10_07_2014 - initial version - J.L. Desmarais
        !
        !>@param filename
        !> filename for the netcdf file open
        !
        !>@param ncid
        !> integer identifying uniquely the netcdf file opened
        !--------------------------------------------------------------
        subroutine nf90_open_file(filename, ncid)

          implicit none

          character(*), intent(in)  :: filename
          integer     , intent(out) :: ncid

          integer :: retval

          retval = NF90_CREATE(trim(filename), NF90_NETCDF4, ncid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

        end subroutine nf90_open_file


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> write the default header for the netcdf file when
        !> saving the buffer layer data
        !
        !> @date
        !> 10_07_2014 - initial version - J.L. Desmarais
        !
        !>@param ncid
        !> integer identifying uniquely the netcdf file opened
        !
        !>@param title
        !> title saved in the netcdf file
        !--------------------------------------------------------------
        subroutine nf90_write_header(ncid,title)

          implicit none

          integer     , intent(in) :: ncid
          character(*), intent(in) :: title


          character(len=29)     :: history
          integer, dimension(8) :: now_data
          integer               :: retval


          !retrieve the current date and time from the system
          !--------------------------------------------------
          call date_and_time(values=now_data)
          write (history,
     $         '(''date'', 1x, i4.4, ''/'', i2.2, ''/'', i2.2, 1x,
     $           ''time'', 1x, i2.2, '':'', i2.2, '':'', i2.2)')
     $         now_data(1), now_data(2), now_data(3),
     $         now_data(5), now_data(6), now_data(7)


          !write the header
          !-----------------
          retval = nf90_put_att(ncid,nf90_global,'title',trim(title))
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = nf90_put_att(ncid,nf90_global,'history',history)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'institution',institut)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'source',prog_version)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'references',ref)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'convention',convention)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

        end subroutine nf90_write_header


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> define the type and extends of the variables saved
        !> in the netcdf file
        !
        !> @date
        !> 10_07_2014 - initial version - J.L. Desmarais
        !
        !>@param ncid
        !> integer identifying uniquely the netcdf file opened
        !
        !>@param size_x
        !> extent of the data along the x-direction
        !
        !>@param size_y
        !> extent of the data along the y-direction
        !
        !>@param coords_id
        !> integer identifying where in the netcdf file the coordinates
        !> along the x- and y-directions are to be saved
        !
        !>@param grdptsid_id
        !> integer identifying where in the netcdf file the grdpts_id
        !> array is to be saved
        !
        !>@param nodes_id
        !> integer identifying where in the netcdf file the nodes
        !> array is to be saved
        !--------------------------------------------------------------
        subroutine bf_layer_nf90_def_var(
     $     ncid, p_model, size_x, size_y, missing_data,
     $     coords_id, grdptsid_id, nodes_id)

          implicit none

          integer               , intent(in)  :: ncid
          type(dim2d_eq)        , intent(in)  :: p_model
          integer               , intent(in)  :: size_x
          integer               , intent(in)  :: size_y
          real(rkind)           , intent(in)  :: missing_data
          integer, dimension(3) , intent(out) :: coords_id
          integer               , intent(out) :: grdptsid_id
          integer, dimension(ne), intent(out) :: nodes_id

          !definition of the coordinates saved in the netcdf file
          character*(*), parameter :: T_NAME = 'time'
          character*(*), parameter :: X_NAME = 'x'
          character*(*), parameter :: Y_NAME = 'y'

          character*(*), parameter :: UNITS    = 'units'
          character*(*), parameter :: T_UNITS  = 's/s'
          character*(*), parameter :: X_UNITS  = 'm/m'
          character*(*), parameter :: Y_UNITS  = 'm/m'

          character*(*), parameter :: AXIS     = 'axis'
          character*(*), parameter :: X_AXIS   = 'X'
          character*(*), parameter :: Y_AXIS   = 'Y'
          
          character*(*), parameter :: LONG_NAME  ='long_name'
          character*(*), parameter :: T_LONG_NAME='reduced time'
          character*(*), parameter :: X_LONG_NAME='reduced x-coordinate'
          character*(*), parameter :: Y_LONG_NAME='reduced y-coordinate'

          integer, parameter       :: NDIMS = 3
          integer, dimension(NDIMS):: dimids


          !identification of the memory locations to save the
          !netcdf variables
          integer :: NF_MYREAL

          character(len=10), dimension(ne) :: name_var
          character(len=33), dimension(ne) :: longname_var
          character(len=23), dimension(ne) :: unit_var

          integer :: t_dimid
          integer :: x_dimid
          integer :: y_dimid
      
          integer :: t_varid
          integer :: x_varid
          integer :: y_varid

          integer :: k
          integer :: retval


          !define the type of variables stored
          select case(RKIND)
            case(4)
               NF_MYREAL=NF90_FLOAT
            case(8)
               NF_MYREAL=NF90_DOUBLE
            case default
               print '(''nf90_operators_wr_class :'')'
               print '(''nf90_def_var_model'')'
               stop 'NF_MYREAL'
          end select


          !define the dimensions
          retval = NF90_DEF_DIM(ncid, T_NAME, 1, t_dimid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_DEF_DIM(ncid, X_NAME, size_x, x_dimid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_DEF_DIM(ncid, Y_NAME, size_y, y_dimid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !define the coordinates variables
          retval = NF90_DEF_VAR(ncid, T_NAME, NF_MYREAL, t_dimid,t_varid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_DEF_VAR(ncid, X_NAME, NF_MYREAL, x_dimid,x_varid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = NF90_DEF_VAR(ncid, Y_NAME, NF_MYREAL, y_dimid,y_varid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          

          !assign the units to the variables
          retval = NF90_PUT_ATT(ncid, t_varid, UNITS, T_UNITS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_PUT_ATT(ncid, x_varid, UNITS, X_UNITS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = NF90_PUT_ATT(ncid, y_varid, UNITS, Y_UNITS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          

          !assign the type of axis to the coordinates
          retval = NF90_PUT_ATT(ncid, x_varid, AXIS, X_AXIS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = NF90_PUT_ATT(ncid, y_varid, AXIS, Y_AXIS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          

          !assign the long name for the description of the coordinates
          retval = NF90_PUT_ATT(ncid, t_varid, LONG_NAME, T_LONG_NAME)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_PUT_ATT(ncid, x_varid, LONG_NAME, X_LONG_NAME)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = NF90_PUT_ATT(ncid, y_varid, LONG_NAME, Y_LONG_NAME)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !assign the coordinates id
          coords_id(1) = t_varid
          coords_id(2) = x_varid
          coords_id(3) = y_varid


          !define the extension in x and y directions
          dimids(1) = t_dimid
          dimids(2) = x_dimid
          dimids(3) = y_dimid


          !define the role of the grid points
          retval = NF90_DEF_VAR(
     $         ncid,
     $         'grdpts_id',
     $         NF90_INT,
     $         dimids,
     $         grdptsid_id)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          !assign units to the role of grid points
          retval = NF90_PUT_ATT(
     $         ncid,
     $         grdptsid_id,
     $         UNITS,
     $         '-')
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          !assign long_name to the role of grid points
          retval = NF90_PUT_ATT(
     $         ncid,
     $         grdptsid_id,
     $         LONG_NAME,
     $         'role of the buffer layer grid points')
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          
          !define the main variables of the governing equations
          !DEC$ FORCEINLINE RECURSIVE
          name_var     = p_model%get_var_name()
          !DEC$ FORCEINLINE RECURSIVE
          longname_var = p_model%get_var_longname()
          !DEC$ FORCEINLINE RECURSIVE
          unit_var     = p_model%get_var_unit()

          do k=1, ne

             !define the governing variables
             retval = NF90_DEF_VAR(
     $            ncid,
     $            trim(name_var(k)),
     $            NF_MYREAL,
     $            dimids,
     $            nodes_id(k))
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             !assign units to the governing variables
             retval = NF90_PUT_ATT(
     $            ncid,
     $            nodes_id(k),
     $            UNITS,
     $            unit_var(k))
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             !assign long_name to the governing variables
             retval = NF90_PUT_ATT(
     $            ncid,
     $            nodes_id(k),
     $            LONG_NAME,
     $            longname_var(k))
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             !assign missing value to
             !the governing variables
             retval = nf90_put_att(
     $            ncid,
     $            nodes_id(k),
     $            'missing_data',
     $            missing_data)
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             !assign fill value to
             !the governing variables
             retval = nf90_put_att(
     $            ncid,
     $            nodes_id(k),
     $            '_FillValue',
     $            missing_data)
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

          end do
          
          !stop the definition of the variables saved in the file
          retval = NF90_ENDDEF(ncid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)          

        end subroutine bf_layer_nf90_def_var


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> define the type and extends of the variables saved
        !> in the netcdf file
        !
        !> @date
        !> 10_07_2014 - initial version - J.L. Desmarais
        !
        !>@param ncid
        !> integer identifying uniquely the netcdf file opened
        !
        !>@param coords_id
        !> integer identifying where in the netcdf file the coordinates
        !> along the x- and y-directions are to be saved
        !
        !>@param grdptsid_id
        !> integer identifying where in the netcdf file the grdpts_id
        !> array is to be saved
        !
        !>@param nodes_id
        !> integer identifying where in the netcdf file the nodes
        !> array is to be saved
        !
        !>@param grdpts_id
        !> array where the role of each grid point is stored (no_pt,
        !> bc_pt, bc_interior_pt, interior_pt)
        !
        !> @param nodes
        !> array where the governing variables are saved at each grid
        !> point
        !
        !> @param time
        !> time corresponding to the data for the grdpts and the nodes
        !--------------------------------------------------------------
        subroutine bf_layer_nf90_put_var(
     $     ncid, coords_id, grdptsid_id, nodes_id,
     $     x_start, y_start, dx, dy,
     $     grdpts_id, nodes, time)

          implicit none

          integer                      , intent(in) :: ncid
          integer    , dimension(3)    , intent(in) :: coords_id
          integer                      , intent(in) :: grdptsid_id
          integer    , dimension(ne)   , intent(in) :: nodes_id
          real(rkind)                  , intent(in) :: x_start
          real(rkind)                  , intent(in) :: y_start
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          integer    , dimension(:,:)  , intent(in) :: grdpts_id
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: time

          real(rkind), dimension(1)              :: time_table
          integer                                :: t_varid
          integer                                :: x_varid
          integer                                :: y_varid
          real(rkind), dimension(:), allocatable :: x_map
          real(rkind), dimension(:), allocatable :: y_map
          integer    , dimension(3)              :: start_op
          integer    , dimension(3)              :: count_op
          integer                                :: i,j,k
          integer                                :: retval


          !define time and variable id
          time_table(1) = time

          t_varid = coords_id(1)
          x_varid = coords_id(2)
          y_varid = coords_id(3)

          
          !write the time
          retval = NF90_PUT_VAR(ncid, t_varid, time_table)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !write the x_map coordinates
          allocate(x_map(size(grdpts_id,1)))
          do i=1, size(grdpts_id,1)
             x_map(i) = x_start + (i-1)*dx
          end do
          retval = NF90_PUT_VAR(ncid, x_varid, x_map)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          deallocate(x_map)


          !write the y_map coordinates
          allocate(y_map(size(grdpts_id,2)))
          do j=1, size(grdpts_id,2)
             y_map(j) = y_start + (j-1)*dy
          end do
          retval = NF90_PUT_VAR(ncid, y_varid, y_map)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          deallocate(y_map)
          

          !start of the writing
          start_op = [1,1,1]
          count_op = [1,size(grdpts_id,1),size(grdpts_id,2)]

          
          !write the grdpts_id
          retval = NF90_PUT_VAR(
     $         ncid,
     $         grdptsid_id,
     $         grdpts_id,
     $         START=start_op,
     $         COUNT=count_op)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !write the variables of the governing equations          
          do k=1, ne

             retval = NF90_PUT_VAR(
     $            ncid,
     $            nodes_id(k),
     $            nodes(:,:,k),
     $            START=start_op,
     $            COUNT=count_op)
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

          end do

        end subroutine bf_layer_nf90_put_var


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> close the netcdf file
        !
        !> @date
        !> 10_07_2014 - initial version - J.L. Desmarais
        !
        !>@param ncid
        !> integer identifying uniquely the netcdf file opened
        !--------------------------------------------------------------
        subroutine nf90_close_file(ncid)

          implicit none

          integer, intent(in) :: ncid

          integer :: retval

          retval = NF90_CLOSE(ncid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

        end subroutine nf90_close_file


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the title for the netcdf file corresponding
        !> to the buffer layer written
        !
        !> @date
        !> 10_07_2014 - initial version - J.L. Desmarais
        !
        !>@param bf_localiaztion
        !> integer identifying the cardinal coordinate of the buffer
        !> layer
        !
        !>@param bf_order
        !> integer identifying the buffer layer in the buffer main layer
        !
        !>@return title
        !> title for the netcdf file giving information about the
        !> localization of the buffer layer
        !--------------------------------------------------------------
        function get_bf_layer_title(
     $     bf_localization, bf_order) result(title)

          implicit none

          integer, intent(in) :: bf_localization
          integer, intent(in) :: bf_order
          character(len=24)   :: title


          character(len=5)    :: cardinal_coord

          select case(bf_localization)
            case(N)
               cardinal_coord = 'North'
            case(S)
               cardinal_coord = 'South'
            case(E)
               cardinal_coord = 'East '
            case(W)
               cardinal_coord = 'West '
            case default
               call error_mainlayer_id(
     $              'bf_layer_nf90_operators_module',
     $              'get_bf_layer_title',
     $              bf_localization)
          end select
          
          write(title,'(A5,'' buffer layer : '',I3)')
     $         cardinal_coord, bf_order

        end function get_bf_layer_title

      end module bf_layer_nf90_operators_module
