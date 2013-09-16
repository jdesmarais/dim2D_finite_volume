      !> @file
      !> module containing subroutines to manipulate
      !> netcdf files
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module containing subroutines to manipulate
      !> netcdf files
      !
      !> @date
      !> 14_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module nf90_operators_module

        use field_class        , only : field
        use netcdf
        use parameters_constant, only : prog_version
        use parameters_kind    , only : rkind, ikind
        use parameters_input   , only : npx,npy,nx,ny,ne,
     $                                  x_min,x_max,y_min,y_max,
     $                                  t_max,dt,detail_print
        use dim2d_eq_class     , only : dim2d_eq
        use dim2d_parameters   , only : viscous_r, Re, We, Pr, cv_r

        implicit none

        private
        public :: nf90_handle_err,
     $            nf90_write_header,
     $            nf90_def_var_model,
     $            nf90_put_var_model

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> handle the error encountered while manipulating
        !> netcdf files
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
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
        !> function giving the name of the file
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param ncid
        !> integer identifying the netcdf file
        !
        !>@param p_model
        !> physical model
        !
        !>@param t
        !> reduced simulation time
        !--------------------------------------------------------------
        subroutine nf90_write_header(ncid, p_model)

          implicit none

          integer       , intent(in) :: ncid
          type(dim2d_eq), intent(in) :: p_model


          character(len=10) :: title
          character(len=29) :: history
          character*(*)     :: institut
          character*(*)     :: source
          character*(*)     :: ref
          character*(*)     :: convention

          parameter (institut   = 'Eindhoven university of technology')
          parameter (source     = prog_version)
          parameter (ref        = 'desmaraisjulien@gmail.com')
          parameter (convention = 'cf-1.6')

          integer, dimension(8) :: now_data
          integer               :: retval


          !<retrieve the model name
          title = p_model%get_model_name()

          !<retrieve the current date and time from the system
          call date_and_time(values=now_data)
          write (history,
     $         '(''date'', 1x, i4.4, ''/'', i2.2, ''/'', i2.2, 1x,
     $           ''time'', 1x, i2.2, '':'', i2.2, '':'', i2.2)')
     $         now_data(1), now_data(2), now_data(3),
     $         now_data(5), now_data(6), now_data(7)


          !<write the header
          retval = nf90_put_att(ncid,nf90_global,'title',trim(title))
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = nf90_put_att(ncid,nf90_global,'history',history)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'institution',institut)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'source',source)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'references', ref)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'convention',convention)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !<write the characteristic parameters of the simulation
          retval = nf90_put_att(ncid,nf90_global,'x_min',x_min)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'x_max',x_max)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'y_min',y_min)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'y_max',y_max)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'t_max',t_max)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'dt',dt)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'detail_print',detail_print)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          

          !<write the characteristic parameters for the DIM
          retval = nf90_put_att(ncid,nf90_global,'viscous_r',viscous_r)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'Re',Re)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'We',We)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'Pr',Pr)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = nf90_put_att(ncid,nf90_global,'cv_r',cv_r)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

        end subroutine nf90_write_header


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine defining the variables saved in
        !> the netcdf file
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !
        !>@param p_model
        !> physical model
        !
        !>@param coordinate_id
        !> integer table identifying the location of the coordinate
        !> tables in the netcdf file
        !
        !>@param data_id
        !> integer table identifying the location of the main varibles
        !> of the governing equations in the netcdf file
        !--------------------------------------------------------------
        subroutine nf90_def_var_model(
     $     ncid,
     $     p_model,
     $     bc_size,
     $     coordinates_id,
     $     data_id)

          implicit none

          integer               , intent(in)    :: ncid
          type(dim2d_eq)        , intent(in)    :: p_model
          integer               , intent(in)    :: bc_size
          integer, dimension(3) , intent(inout) :: coordinates_id
          integer, dimension(ne), intent(inout) :: data_id


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

          integer :: NF_MYREAL

          character(len=10), dimension(ne) :: name_var
          character(len=32), dimension(ne) :: longname_var
          character(len=10), dimension(ne) :: unit_var

          integer :: t_dimid
          integer :: x_dimid
          integer :: y_dimid
      
          integer :: t_varid
          integer :: x_varid
          integer :: y_varid

          integer(ikind) :: nx_domain
          integer(ikind) :: ny_domain

          integer :: k
          integer :: retval


          !<define the type of variables stored
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


          !<define the dimensions
          nx_domain = npx*(nx-2*bc_size) + 2*bc_size
          ny_domain = npy*(ny-2*bc_size) + 2*bc_size

          retval = NF90_DEF_DIM(ncid, T_NAME, 1, t_dimid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_DEF_DIM(ncid, X_NAME, nx_domain, x_dimid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_DEF_DIM(ncid, Y_NAME, ny_domain, y_dimid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !<define the coordinates variables
          retval = NF90_DEF_VAR(ncid, T_NAME, NF_MYREAL, t_dimid,t_varid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_DEF_VAR(ncid, X_NAME, NF_MYREAL, x_dimid,x_varid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = NF90_DEF_VAR(ncid, Y_NAME, NF_MYREAL, y_dimid,y_varid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          

          !<assign the units to the variables
          retval = NF90_PUT_ATT(ncid, t_varid, UNITS, T_UNITS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_PUT_ATT(ncid, x_varid, UNITS, X_UNITS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = NF90_PUT_ATT(ncid, y_varid, UNITS, Y_UNITS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          

          !<assign the type of axis to the coordinates
          retval = NF90_PUT_ATT(ncid, x_varid, AXIS, X_AXIS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = NF90_PUT_ATT(ncid, y_varid, AXIS, Y_AXIS)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          

          !<assign the long name for the description of the coordinates
          retval = NF90_PUT_ATT(ncid, t_varid, LONG_NAME, T_LONG_NAME)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

          retval = NF90_PUT_ATT(ncid, x_varid, LONG_NAME, X_LONG_NAME)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)
          
          retval = NF90_PUT_ATT(ncid, y_varid, LONG_NAME, Y_LONG_NAME)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !<assign the coordinates id
          coordinates_id(1) = t_varid
          coordinates_id(2) = x_varid
          coordinates_id(3) = y_varid


          !<define the extension in x and y directions
          dimids(1) = t_dimid
          dimids(2) = x_dimid
          dimids(3) = y_dimid

          
          !<define the main variables of the governing equations
          name_var     = p_model%get_var_name()
          longname_var = p_model%get_var_longname()
          unit_var     = p_model%get_var_unit()

          do k=1, ne

             !<define the netcdf variables
             retval = NF90_DEF_VAR(ncid, trim(name_var(k)), NF_MYREAL,
     $            dimids, data_id(k))
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             !assign the units to the variables
             retval = NF90_PUT_ATT(ncid, data_id(k), UNITS, unit_var(k))
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

             !assign the long_name to the variables
             retval = NF90_PUT_ATT(ncid, data_id(k), LONG_NAME,
     $            longname_var(k))
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)

          end do
          
          !<stop the definition of the variables saved in the file
          retval = NF90_ENDDEF(ncid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)

        end subroutine nf90_def_var_model


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> put the main variables in the netcdf file
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !> @param ncid
        !> integer identifying the netcdf file opened
        !
        !> @param coordinates_id
        !> table with the id identifying the coordinates
        !
        !> @param data_id
        !> table with the id identifying the data
        !
        !> @param time
        !> reduced simulation time
        !
        !> @param field_used
        !> object encapsulating the main variables
        !
        !> @param start
        !> integer, dimension(2) identifying the sw point
        !> where the partial filling starts
        !
        !> @param count
        !> integer, dimension(2) identifying the profile
        !> of the gridpoints written
        !
        !> @param limit
        !> integer, dimension(4) identifying the borders
        !> of the sub-tile written out of the computational
        !> domain
        !----------------------------------------------------------------
        subroutine nf90_put_var_model(
     $     ncid,
     $     coordinates_id,
     $     data_id,
     $     time, field_used,
     $     start, count, limit)

          implicit none
          
          integer                        , intent(in) :: ncid
          integer    , dimension(3)      , intent(in) :: coordinates_id
          integer    , dimension(ne)     , intent(in) :: data_id
          real(RKIND)                    , intent(in) :: time
          class(field)                   , intent(in) :: field_used
          integer, dimension(2), optional, intent(in) :: start
          integer, dimension(2), optional, intent(in) :: count
          integer, dimension(4), optional, intent(in) :: limit


          integer    , dimension(3) :: start_op, count_op
          real(RKIND), dimension(1) :: time_table
          integer :: t_varid
          integer :: x_varid
          integer :: y_varid
          integer :: k
          integer :: retval


          time_table(1) = time

          t_varid = coordinates_id(1)
          x_varid = coordinates_id(2)
          y_varid = coordinates_id(3)

          
          !< write the time
          retval = NF90_PUT_VAR(ncid, t_varid, time_table)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !< depending if the start and count optional variables
          !> are present, we should write the whole domain or
          !> only a part identified by (start, count) with the
          !> exception of the time
          !-----------------------------------------------------
          if(.not.present(start)) then

             start_op = [1,1,1]
             count_op = [1,nx,ny]

             !< write the x_map coordinates
             retval = NF90_PUT_VAR(ncid, x_varid, field_used%x_map)
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)


             !< write the y_map coordinates
             retval = NF90_PUT_VAR(ncid, y_varid, field_used%y_map)
             !DEC$ FORCEINLINE RECURSIVE
             call nf90_handle_err(retval)
             

             !< write the variables of the governing equations
             do k=1, ne

                retval = NF90_PUT_VAR(
     $               ncid,
     $               data_id(k),
     $               field_used%nodes(:,:,k),
     $               START=start_op,
     $               COUNT=count_op)
                !DEC$ FORCEINLINE RECURSIVE
                call nf90_handle_err(retval)

             end do

         else

            !< write the x_map coordinates
            retval = NF90_PUT_VAR(
     $           ncid,
     $           x_varid,
     $           field_used%x_map(limit(1):limit(2)),
     $           START=[start(1)],
     $           COUNT=[count(1)])
            !DEC$ FORCEINLINE RECURSIVE
            call nf90_handle_err(retval)


            !< write the y_map coordinates
            retval = NF90_PUT_VAR(
     $           ncid,
     $           y_varid,
     $           field_used%y_map(limit(3):limit(4)),
     $           START=[start(2)],
     $           COUNT=[count(2)])
            !DEC$ FORCEINLINE RECURSIVE
            call nf90_handle_err(retval)

            
            !< write the variables of the governing equations
            do k=1, ne

               retval = NF90_PUT_VAR(
     $              ncid,
     $              data_id(k),
     $              field_used%nodes(limit(1):limit(2),limit(3):limit(4),k),
     $              START=[1, start(1), start(2)], 
     $              COUNT=[1, count(1), count(2)])
               !DEC$ FORCEINLINE RECURSIVE
               call nf90_handle_err(retval)

            end do

          end if

        end subroutine nf90_put_var_model

      end module nf90_operators_module
