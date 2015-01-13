      module nf90_operators_error_module

        use compare_type_module, only :
     $       compare_reals,
     $       compare_doubles

        use netcdf

        use parameters_cst, only :
     $       NOT_SUCCESS,
     $       SUCCESS

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none


        private
        public :: 
     $       nf90_put_var_error,
     $       nf90_def_var_error,
     $       nf90_def_header_error,
     $       nf90_get_time,
     $       nf90_get_gov_var,
     $       nf90_get_maps,
     $       nf90_handle_err


        contains
        
        !put the variables in the error file
        subroutine nf90_put_var_error(
     $       ncid_error,
     $       coordinates_id,
     $       error_id,
     $       time,
     $       x_map,
     $       y_map,
     $       error,
     $       error_max_id,
     $       error_max,
     $       x_error_max_id,
     $       x_error_max,
     $       y_error_max_id,
     $       y_error_max)

          implicit none

          integer                                , intent(in) :: ncid_error
          integer    , dimension(3)              , intent(in) :: coordinates_id
          integer    , dimension(:)              , intent(in) :: error_id
          real(rkind)                            , intent(in) :: time
          real(rkind), dimension(:)              , intent(in) :: x_map
          real(rkind), dimension(:)              , intent(in) :: y_map
          real(rkind), dimension(:,:,:)          , intent(in) :: error
          integer    , dimension(:)    , optional, intent(in) :: error_max_id
          real(rkind), dimension(:)    , optional, intent(in) :: error_max
          integer    , dimension(:)    , optional, intent(in) :: x_error_max_id
          real(rkind), dimension(:)    , optional, intent(in) :: x_error_max
          integer    , dimension(:)    , optional, intent(in) :: y_error_max_id
          real(rkind), dimension(:)    , optional, intent(in) :: y_error_max

          
          integer :: retval
          integer :: k


          ! put the time
          retval = NF90_PUT_VAR(ncid_error,coordinates_id(1),time)
          call nf90_handle_err(retval)


          ! put the x-coordinates
          retval = NF90_PUT_VAR(ncid_error,coordinates_id(2),x_map)
          call nf90_handle_err(retval)
          

          ! put the y-coordinates
          retval = NF90_PUT_VAR(ncid_error,coordinates_id(3),y_map)
          call nf90_handle_err(retval)

          
          ! save the relative error
          do k=1, size(error,3)

             retval = NF90_PUT_VAR(
     $            ncid_error,
     $            error_id(k),
     $            error(:,:,k),
     $            START=[1,1,1],
     $            COUNT=[1,size(error,1),size(error,2)])
             call nf90_handle_err(retval)

          end do


          ! if present, save the error_max
          if(present(error_max_id).and.present(error_max)) then

             do k=1, size(error,3)

                retval = NF90_PUT_VAR(ncid_error,error_max_id(k),error_max(k))
                call nf90_handle_err(retval)

             end do

          end if


          ! if present, save the x_error_max
          if(present(x_error_max_id).and.present(x_error_max)) then

             do k=1, size(error,3)

                retval = NF90_PUT_VAR(ncid_error,x_error_max_id(k),x_error_max(k))
                call nf90_handle_err(retval)

             end do

          end if


          ! if present, save the y_error_max
          if(present(y_error_max_id).and.present(y_error_max)) then

             do k=1, size(error,3)

                retval = NF90_PUT_VAR(ncid_error,y_error_max_id(k),y_error_max(k))
                call nf90_handle_err(retval)

             end do

          end if

        end subroutine nf90_put_var_error


        !definition of the variables saved in the error file
        subroutine nf90_def_var_error(
     $       filename_sm_domain,
     $       ncid_error,
     $       nx,ny,ne,
     $       coordinates_id,
     $       error_id,
     $       error_max_id,
     $       x_error_max_id,
     $       y_error_max_id)

          implicit none

          character*(*)                                     , intent(in)  :: filename_sm_domain
          integer                                           , intent(in)  :: ncid_error
          integer                                           , intent(in)  :: nx
          integer                                           , intent(in)  :: ny
          integer                                           , intent(in)  :: ne
          integer      , dimension(3)                       , intent(out) :: coordinates_id
          integer      , dimension(:), allocatable          , intent(out) :: error_id
          integer      , dimension(:), allocatable, optional, intent(out) :: error_max_id
          integer      , dimension(:), allocatable, optional, intent(out) :: x_error_max_id
          integer      , dimension(:), allocatable, optional, intent(out) :: y_error_max_id

          integer :: retval

          integer, dimension(3) :: sizes
          integer               :: NF_MYREAL
          integer               :: ncid_sm_domain
          integer, dimension(3) :: dims_id
          integer               :: k,l

          integer             :: lenNameVar
          character(len=60)   :: nameVar
          character(len=1024) :: nameVarModified

          integer             :: lenCharVarAtt
          integer             :: nAttsVar
          character(len=60)   :: nameAtt
          character(len=1024) :: charVarAtt
          character(len=1024) :: charVarAttModified

          character(len=30), dimension(:), allocatable :: nameGovVar


          !0) sizes
          sizes = [1,nx,ny]


          !1) define the type of variables stored
          select case(RKIND)
            case(4)
               NF_MYREAL=NF90_FLOAT
            case(8)
               NF_MYREAL=NF90_DOUBLE
            case default
               print '(''error_netcdf_file_module'')'
               print '(''nf90_def_var_error'')'
               stop 'NF_MYREAL'
          end select

          
          !2) open the small domain netcdf file
          retval = NF90_OPEN(
     $         trim(filename_sm_domain),
     $         NF90_NOWRITE,
     $         ncid_sm_domain)
          call nf90_handle_err(retval)


          !3) extract the name, units and longname
          !   of the coordinates
          do k=1,3
             
             ! extract the name of the coordinate
             ! and save it in the error_file
             retval = NF90_INQUIRE_VARIABLE(
     $            ncid_sm_domain,
     $            k,
     $            name=nameVar,
     $            nAtts=nAttsVar)
             call nf90_handle_err(retval)
             

             ! define the dimension with the name
             retval = NF90_DEF_DIM(
     $            ncid_error,
     $            trim(nameVar),
     $            sizes(k),
     $            dims_id(k))
             call nf90_handle_err(retval)


             ! define the corresponding variable
             retval = NF90_DEF_VAR(
     $            ncid_error,
     $            trim(nameVar),
     $            NF_MYREAL,
     $            dims_id(k),
     $            coordinates_id(k))
             call nf90_handle_err(retval)


             ! extract the attributes of the
             ! coordinates and save them in
             ! the error file
             do l=1,nAttsVar

                retval = NF90_INQ_ATTNAME(
     $               ncid_sm_domain,
     $               k,
     $               l,
     $               nameAtt)
                call nf90_handle_err(retval)

                retval = NF90_GET_ATT(
     $               ncid_sm_domain,
     $               k,
     $               trim(nameAtt),
     $               charVarAtt)
                call nf90_handle_err(retval)

                retval = NF90_PUT_ATT(
     $               ncid_error,
     $               coordinates_id(k),
     $               trim(nameAtt),
     $               trim(charVarAtt))
                call nf90_handle_err(retval)
                
             end do

          end do


          !4) extract and define the governing 
          !   variables
          !   [1,3]: coordinates
          !   [4,:]: governing variables
          allocate(error_id(ne))
          allocate(nameGovVar(ne))

          do k=1, ne

             ! extract the governing variable
             ! from the netcdf file of the
             ! small domain
             retval = NF90_INQUIRE_VARIABLE(
     $            ncid_sm_domain,
     $            3+k,
     $            name=nameVar,
     $            nAtts=nAttsVar)
             call nf90_handle_err(retval)


             ! save the name of the governing
             ! variable
             lenNameVar                      = len(trim(nameVar))
             nameVarModified                 = ''
             nameVarModified(1:6)            = 'error_'
             nameVarModified(7:7+lenNameVar) = trim(nameVar)
             nameGovVar(k)                   = trim(nameVarModified)


             ! define the governing variable to
             ! be saved in the netcdf file 
             ! containing the error
             retval = NF90_DEF_VAR(
     $            ncid_error,
     $            nameGovVar(k),
     $            NF_MYREAL,
     $            dims_id,
     $            error_id(k))
             call nf90_handle_err(retval)

             
             ! extract and assign the units to
             ! the governing variables
             ! extract the attributes of the
             ! coordinates and save them in
             ! the error file
             do l=1,nAttsVar

                retval = NF90_INQ_ATTNAME(
     $               ncid_sm_domain,
     $               3+k,
     $               l,
     $               nameAtt)
                call nf90_handle_err(retval)

                retval = NF90_GET_ATT(
     $               ncid_sm_domain,
     $               3+k,
     $               trim(nameAtt),
     $               charVarAtt)
                call nf90_handle_err(retval)

                select case(trim(nameAtt))
                  case('units')
                     charVarAttModified = '-'

                     retval = NF90_PUT_ATT(
     $                    ncid_error,
     $                    error_id(k),
     $                    trim(nameAtt),
     $                    trim(charVarAttModified))
                     call nf90_handle_err(retval)

                  case('long_name')
                     lenCharVarAtt            = len(trim(charVarAtt))
                     charVarAttModified       = ''
                     charVarAttModified(1:18) = 'relative error of '
                     charVarAttModified(19:19+lenCharVarAtt) = trim(charVarAtt)

                     
                     retval = NF90_PUT_ATT(
     $                    ncid_error,
     $                    error_id(k),
     $                    trim(nameAtt),
     $                    trim(charVarAttModified))
                     call nf90_handle_err(retval)
                     
                  case default

                     retval = NF90_PUT_ATT(
     $                    ncid_error,
     $                    error_id(k),
     $                    trim(nameAtt),
     $                    trim(charVarAtt))
                     call nf90_handle_err(retval)

                end select
                
             end do

          end do


          !5) if asked by the user, define also
          !   the variable for the error_max
          if(present(error_max_id)) then

             allocate(error_max_id(ne))

             do k=1, ne


                ! define the length of the governing variable
                lenNameVar = len(trim(nameGovVar(k)))

                ! define the variable name for the error max
                nameVar = ''
                nameVar(1:4) = 'max_'
                nameVar(5:5+lenNameVar) = trim(nameGovVar(k))
                
                retval = NF90_DEF_VAR(
     $               ncid_error,
     $               trim(nameVar),
     $               NF_MYREAL,
     $               [dims_id(1)],
     $               error_max_id(k))
                call nf90_handle_err(retval)

                ! assign the unit
                retval = NF90_PUT_ATT(
     $               ncid_error,
     $               error_max_id(k),
     $               'units',
     $               '-')
                call nf90_handle_err(retval)

                ! assign the description of the variable
                nameVar = ''
                nameVar = 'maximum of '
                nameVar(12:12+lenNameVar) = trim(nameGovVar(k))
                retval = NF90_PUT_ATT(
     $               ncid_error,
     $               error_max_id(k),
     $               'long_name',
     $               nameVar)
                call nf90_handle_err(retval)

             end do

          end if


          !5) if asked by the user, define also
          !   the variable for the x_error_max
          if(present(x_error_max_id)) then

             allocate(x_error_max_id(ne))

             do k=1, ne

                ! define the variable name for the x_error_max
                lenNameVar = len(trim(nameGovVar(k)))
                nameVar = ''
                nameVar(1:6) = 'x_max_'
                nameVar(7:7+lenNameVar) = trim(nameGovVar(k))
                
                retval = NF90_DEF_VAR(
     $               ncid_error,
     $               trim(nameVar),
     $               NF_MYREAL,
     $               [dims_id(1)],
     $               x_error_max_id(k))
                call nf90_handle_err(retval)

                ! assign the unit
                retval = NF90_PUT_ATT(
     $               ncid_error,
     $               x_error_max_id(k),
     $               'units',
     $               'm/m')
                call nf90_handle_err(retval)

                ! assign the description of the variable
                nameVar = ''
                nameVar = 'x-coordinate of max( '
                nameVar(22:22+lenNameVar) = trim(nameGovVar(k))
                nameVar(22+lenNameVar+1:22+lenNameVar+1) = ')'
                
                retval = NF90_PUT_ATT(
     $               ncid_error,
     $               x_error_max_id(k),
     $               'long_name',
     $               nameVar)
                call nf90_handle_err(retval)

             end do

          end if


          !6) if asked by the user, define also
          !   the variable for the y_error_max
          if(present(y_error_max_id)) then

             allocate(y_error_max_id(ne))

             do k=1, ne

                ! define the variable name for the y_error_max
                lenNameVar = len(trim(nameGovVar(k)))
                nameVar = ''
                nameVar(1:6) = 'y_max_'
                nameVar(7:7+lenNameVar) = trim(nameGovVar(k))

                retval = NF90_DEF_VAR(
     $               ncid_error,
     $               trim(nameVar),
     $               NF_MYREAL,
     $               [dims_id(1)],
     $               y_error_max_id(k))
                call nf90_handle_err(retval)

                ! assign the unit
                retval = NF90_PUT_ATT(
     $               ncid_error,
     $               y_error_max_id(k),
     $               'units',
     $               'm/m')
                call nf90_handle_err(retval)

                ! assign the description of the variable
                nameVar = ''
                nameVar = 'y-coordinate of max( '
                nameVar(22:22+lenNameVar) = trim(nameGovVar(k))
                nameVar(22+lenNameVar+1:22+lenNameVar+1) = ')'
                
                retval = NF90_PUT_ATT(
     $               ncid_error,
     $               y_error_max_id(k),
     $               'long_name',
     $               nameVar)
                call nf90_handle_err(retval)

             end do

          end if


          !7) stop the definition of the
          !   variables saved in the error
          !   file
          retval = NF90_ENDDEF(ncid_error)
          call nf90_handle_err(retval)


          !8) close small domain netcdf file
          retval = NF90_CLOSE(ncid_sm_domain)
          call nf90_handle_err(retval)

        end subroutine nf90_def_var_error


        !create header informations from the small and large domain
        !headers
        subroutine nf90_def_header_error(
     $     filename_sm_domain,
     $     filename_lg_domain,
     $     filename_error,
     $     ncid_error,
     $     ierror)
        
          implicit none

          character*(*)          , intent(in)  :: filename_sm_domain
          character*(*)          , intent(in)  :: filename_lg_domain
          character*(*)          , intent(in)  :: filename_error
          integer                , intent(out) :: ncid_error
          integer      , optional, intent(out) :: ierror

          integer               :: retval
          integer               :: ncid_sm_domain
          integer               :: ncid_lg_domain

          integer               :: len_title
          character(len=26)     :: title
          character(len=29)     :: history
          integer, dimension(8) :: now_data

          integer               :: size_header_sm_domain
          integer               :: size_header_lg_domain

          integer               :: k

          character(len=1028)   :: header_attname_sm_domain
          integer               :: header_att_type
          integer               :: header_att_len

          character(len=1028)   :: att_char_var
          integer               :: att_int_var
          real                  :: att_float_var
          real*8                :: att_double_var
          real                  :: att_float_var2
          real*8                :: att_double_var2

        
          if(present(ierror)) then
             ierror = NOT_SUCCESS
          end if


          !1) open netcdf file for saving the header
          !   as a combination of the small and large
          !   domains headers
          retval = NF90_CREATE(
     $         trim(filename_error),
     $         NF90_NETCDF4,
     $         ncid_error)
          call nf90_handle_err(retval)


          !2) get the date for the header of the
          !   error file
          call date_and_time(values=now_data)
          write (history,
     $         '(''date'', 1x, i4.4, ''/'', i2.2, ''/'', i2.2, 1x,
     $           ''time'', 1x, i2.2, '':'', i2.2, '':'', i2.2)')
     $         now_data(1), now_data(2), now_data(3),
     $         now_data(5), now_data(6), now_data(7)


          !3) open second netcdf file for reading
          !   small domain header
          retval = NF90_OPEN(
     $         trim(filename_sm_domain),
     $         NF90_NOWRITE,
     $         ncid_sm_domain)
          call nf90_handle_err(retval)

          
          !4) inquire the number of global
          !   variables (header variables)
          retval = NF90_INQUIRE(
     $         ncid_sm_domain,
     $         nAttributes=size_header_sm_domain)
          call nf90_handle_err(retval)


          !5) open second netcdf file for reading
          !   large domain header
          retval = NF90_OPEN(
     $         trim(filename_lg_domain),
     $         NF90_NOWRITE,
     $         ncid_lg_domain)
          call nf90_handle_err(retval)


          !6) inquire the number of global
          !   variables for the large domain
          !   (header variables)
          retval = NF90_INQUIRE(
     $         ncid_lg_domain,
     $         nAttributes=size_header_lg_domain)
          call nf90_handle_err(retval)


          !7) inquire the name of the global
          !   variables for the small and large
          !   domains
          do k=1, size_header_sm_domain

             ! inquire the name of the global
             ! variables in the header
             retval = NF90_INQ_ATTNAME(
     $            ncid_sm_domain,
     $            NF90_GLOBAL,
     $            k,
     $            header_attname_sm_domain)
             call nf90_handle_err(retval)


             ! inquire the type and length of
             ! the global variable
             retval = NF90_INQUIRE_ATTRIBUTE(
     $            ncid_sm_domain,
     $            NF90_GLOBAL,
     $            trim(header_attname_sm_domain),
     $            xtype=header_att_type,
     $            len=header_att_len)
             call nf90_handle_err(retval)


             ! save the global variable in the 
             ! new header
             select case(header_att_type)
               case(NF90_CHAR)
                  retval = nf90_get_att(
     $                 ncid_sm_domain,
     $                 NF90_GLOBAL,
     $                 trim(header_attname_sm_domain),
     $                 att_char_var)

               case(NF90_INT)
                  retval = nf90_get_att(
     $                 ncid_sm_domain,
     $                 NF90_GLOBAL,
     $                 trim(header_attname_sm_domain),
     $                 att_int_var)

               case(NF90_FLOAT)
                  retval = nf90_get_att(
     $                 ncid_sm_domain,
     $                 NF90_GLOBAL,
     $                 trim(header_attname_sm_domain),
     $                 att_float_var)

               case(NF90_DOUBLE)
                  retval = nf90_get_att(
     $                 ncid_sm_domain,
     $                 NF90_GLOBAL,
     $                 trim(header_attname_sm_domain),
     $                 att_double_var)
             
               case default
                  print '(''error_netcdf_file_module'')'
                  print '(''create_common_header_from_netcdf_files'')'
                  stop ''

             end select


             ! depending on the header attribute, the
             ! attribute is directly saved in the header
             ! of the error file or it is replaced or
             ! renamed
             ! (only modify the
             ! x_min,x_max,y_min,y_max attributes)
             select case(trim(header_attname_sm_domain))

               case('title')
                  
                  len_title = len(trim(att_char_var))
                  title(1:len_title) = trim(att_char_var)
                  title(len_title+1:len_title+16) = ': relative_error'
                  
                  retval = nf90_put_att(
     $                 ncid_error,
     $                 nf90_global,
     $                 'title',
     $                 trim(title))


               case('history')

                  retval = nf90_put_att(
     $                 ncid_error,
     $                 nf90_global,
     $                 'history',
     $                 history)

               case('x_min')
                  
                  call nf90_put_sm_lg_att(
     $                 ncid_error,
     $                 ncid_lg_domain,
     $                 att_float_var,
     $                 att_double_var,
     $                 header_att_type,
     $                 'x_min',
     $                 'x_min_sm_domain',
     $                 'x_min_lg_domain')

               case('x_max')

                  call nf90_put_sm_lg_att(
     $                 ncid_error,
     $                 ncid_lg_domain,
     $                 att_float_var,
     $                 att_double_var,
     $                 header_att_type,
     $                 'x_max',
     $                 'x_max_sm_domain',
     $                 'x_max_lg_domain')

               case('y_min')
                  
                  call nf90_put_sm_lg_att(
     $                 ncid_error,
     $                 ncid_lg_domain,
     $                 att_float_var,
     $                 att_double_var,
     $                 header_att_type,
     $                 'y_min',
     $                 'y_min_sm_domain',
     $                 'y_min_lg_domain')

               case('y_max')

                  call nf90_put_sm_lg_att(
     $                 ncid_error,
     $                 ncid_lg_domain,
     $                 att_float_var,
     $                 att_double_var,
     $                 header_att_type,
     $                 'y_max',
     $                 'y_max_sm_domain',
     $                 'y_max_lg_domain')

               ! for the other attributes, they should be
               ! compared to check whether the small and
               ! the large domain simulations were run
               ! with the same properties
               case default

                  select case(header_att_type)

                    case(NF90_INT)
                       retval = nf90_put_att(
     $                      ncid_error,
     $                      NF90_GLOBAL,
     $                      trim(header_attname_sm_domain),
     $                      att_int_var)
                       call nf90_handle_err(retval)


                    case(NF90_CHAR)
                       retval = nf90_put_att(
     $                      ncid_error,
     $                      NF90_GLOBAL,
     $                      trim(header_attname_sm_domain),
     $                      att_char_var)
                       call nf90_handle_err(retval)


                    case(NF90_FLOAT)
                       retval = nf90_get_att(
     $                      ncid_lg_domain,
     $                      NF90_GLOBAL,
     $                      trim(header_attname_sm_domain),
     $                      att_float_var2)

                       if(.not.compare_reals(att_float_var,att_float_var2)) then
                          print *, 'WARNING attribute ',
     $                          trim(header_attname_sm_domain),
     $                          ': do not match'
                       end if

                       retval = nf90_put_att(
     $                      ncid_error,
     $                      NF90_GLOBAL,
     $                      trim(header_attname_sm_domain),
     $                      att_float_var)
                       call nf90_handle_err(retval)
                          
                       
                    case(NF90_DOUBLE)
                       retval = nf90_get_att(
     $                      ncid_lg_domain,
     $                      NF90_GLOBAL,
     $                      trim(header_attname_sm_domain),
     $                      att_double_var2)
                       
                       if(.not.compare_doubles(att_double_var,att_double_var2)) then
                          print *, 'WARNING attribute ',
     $                          trim(header_attname_sm_domain),
     $                          ': do not match'
                       end if

                       retval = nf90_put_att(
     $                      ncid_error,
     $                      NF90_GLOBAL,
     $                      trim(header_attname_sm_domain),
     $                      att_double_var)
                       call nf90_handle_err(retval)


                  end select                 

             end select

          end do
          

          !8) close netcdf files
          retval = NF90_CLOSE(ncid_sm_domain)
          call nf90_handle_err(retval)
          
          retval = NF90_CLOSE(ncid_lg_domain)
          call nf90_handle_err(retval)


          if(present(ierror)) then
             ierror = SUCCESS
          end if


        end subroutine nf90_def_header_error
        

        !put the header attribute of the small and large domain
        !simulations in the error file
        subroutine nf90_put_sm_lg_att(
     $     ncid_error,
     $     ncid_lg_domain,
     $     att_real_var,
     $     att_double_var,
     $     type_sm_domain,
     $     attname,
     $     attname_sm_domain,
     $     attname_lg_domain)

          implicit none

          integer      , intent(in) :: ncid_error
          integer      , intent(in) :: ncid_lg_domain
          real         , intent(in) :: att_real_var
          real*8       , intent(in) :: att_double_var
          integer      , intent(in) :: type_sm_domain
          character*(*), intent(in) :: attname
          character*(*), intent(in) :: attname_sm_domain
          character*(*), intent(in) :: attname_lg_domain
        
          
          integer :: retval
          real    :: att_real_var2
          real*8  :: att_double_var2


          select case(type_sm_domain)

            case(NF90_FLOAT)
               retval = nf90_put_att(
     $              ncid_error,
     $              NF90_GLOBAL,
     $              attname_sm_domain,
     $              att_real_var)

               retval = nf90_get_att(
     $              ncid_lg_domain,
     $              NF90_GLOBAL,
     $              attname,
     $              att_real_var2)

               retval = nf90_put_att(
     $              ncid_error,
     $              NF90_GLOBAL,
     $              attname_lg_domain,
     $              att_real_var2)
               
            case(NF90_DOUBLE)
               retval = nf90_put_att(
     $              ncid_error,
     $              NF90_GLOBAL,
     $              attname_sm_domain,
     $              att_double_var)

               retval = nf90_get_att(
     $              ncid_lg_domain,
     $              NF90_GLOBAL,
     $              attname,
     $              att_double_var2)

               retval = nf90_put_att(
     $              ncid_error,
     $              NF90_GLOBAL,
     $              attname_lg_domain,
     $              att_double_var2)

          end select

        end subroutine nf90_put_sm_lg_att


        !get the time from the netcdf file
        function nf90_get_time(
     $       filename,
     $       ierror)
     $       result(time)

          implicit none

          character*(*)          , intent(in)  :: filename
          integer      , optional, intent(out) :: ierror
          real(rkind)                          :: time

          
          integer                   :: retval
          integer                   :: ncid
          integer                   :: time_varid
          real(rkind), dimension(1) :: time_table

          
          if(present(ierror)) then
             ierror = NOT_SUCCESS
          end if
          

          !1) open the netcdf file for reading
          retval = NF90_OPEN(trim(filename), NF90_NOWRITE, ncid)
          call nf90_handle_err(retval)


          !2) get the varid for time
          retval = NF90_INQ_VARID(ncid, 'time', time_varid)
          call nf90_handle_err(retval)

          
          !3) extract time
          retval = NF90_GET_VAR(ncid,time_varid,time_table(:))
          call nf90_handle_err(retval)
          time = time_table(1)


          !4) close netcdf file
          retval = NF90_CLOSE(ncid)
          call nf90_handle_err(retval)


          if(present(ierror)) then
             ierror = SUCCESS
          end if

        end function nf90_get_time


        !extract a var from the netcdf file
        subroutine nf90_get_gov_var(
     $       filename,
     $       var_name,
     $       var,
     $       start,
     $       count,
     $       ierror)

          implicit none

          
          character*(*)                               , intent(in)  :: filename
          character*(*), dimension(:)                 , intent(in)  :: var_name
          real(rkind)  , dimension(:,:,:), allocatable, intent(out) :: var
          integer      , dimension(2)    , optional   , intent(in)  :: start
          integer      , dimension(2)    , optional   , intent(in)  :: count
          integer                        , optional   , intent(out) :: ierror

          integer                                   :: ierror_op
          integer                                   :: retval
          integer                                   :: ncid
          integer       , dimension(:), allocatable :: varid
          integer(ikind)                            :: size_x_map
          integer(ikind)                            :: size_y_map
          integer       , dimension(3)              :: start_op
          integer       , dimension(3)              :: count_op
          integer                                   :: i


          ierror_op = NOT_SUCCESS


          !1) open the netcdf file for reading
          retval = NF90_OPEN(trim(filename), NF90_NOWRITE, ncid)
          call nf90_handle_err(retval)


          !2) get the varid for the variable 'var_name(:)'
          allocate(varid(size(var_name,1)))
          do i=1, size(var_name,1)

             retval = NF90_INQ_VARID(ncid, trim(var_name(i)), varid(i))
             call nf90_handle_err(retval)

          end do


          !3) get the extension for the variable 'var_name(:)'
          !   if no restriction is asked, the size of the array
          !   extracted matchs the size of the x and y-maps
          !   2 is the dimID of the 'x-coordinate'
          !   3 is the dimID of the 'y-coordinate'
          if((.not.present(start)).or.(.not.present(count))) then
             retval = NF90_INQUIRE_DIMENSION(ncid, 2, len=size_x_map)
             call nf90_handle_err(retval)
             retval = NF90_INQUIRE_DIMENSION(ncid, 3, len=size_y_map)
             call nf90_handle_err(retval)             

             start_op = [1,1,1]
             count_op = [1,size_x_map,size_y_map]

          !  otherwise, the size of the array is given by the optional
          !  argument 'count'
          else
             size_x_map = count(1)
             size_y_map = count(2)

             start_op = [1,start(1),start(2)]
             count_op = [1,count(1),count(2)]

          end if
          allocate(var(size_x_map,size_y_map,size(var_name,1)))


          !4) extract the arrays of var_name
          do i=1, size(varid,1)

             retval = NF90_GET_VAR(
     $            ncid,
     $            varid(i),
     $            var(:,:,i),
     $            START=start_op,
     $            COUNT=count_op)
             call nf90_handle_err(retval)

          end do


          !5) close netcdf file
          retval = NF90_CLOSE(ncid)
          call nf90_handle_err(retval)


          ierror_op = SUCCESS

          if(present(ierror)) then
             ierror = ierror_op
          end if

        end subroutine nf90_get_gov_var

      
        !extract the x_map and the y_map from the netcdf file
        subroutine nf90_get_maps(
     $       filename,
     $       x_map,
     $       y_map)

          implicit none

          character*(*)                         , intent(in)  :: filename
          real(rkind), dimension(:), allocatable, intent(out) :: x_map
          real(rkind), dimension(:), allocatable, intent(out) :: y_map

          integer        :: retval
          integer        :: ncid
          integer        :: x_varid
          integer        :: y_varid
          integer(ikind) :: size_x_map
          integer(ikind) :: size_y_map


          !1) open netcdf file for reading
          retval = NF90_OPEN(trim(filename), NF90_NOWRITE, ncid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !2) get the varid for the 'x' variable
          retval = NF90_INQ_VARID(ncid, "x", x_varid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !3) get the varid for the 'y' variable
          retval = NF90_INQ_VARID(ncid, "y", y_varid)
          !DEC$ FORCEINLINE RECURSIVE
          call nf90_handle_err(retval)


          !4) 2 is the dimID of the x-coordinate
          !   get the extent of the 'x_map'
          !   and the extract the 'x_map'
          retval = NF90_INQUIRE_DIMENSION(ncid, 2, len=size_x_map)
          call nf90_handle_err(retval)
          allocate(x_map(size_x_map))
          retval = NF90_GET_VAR(ncid,x_varid,x_map)
          call nf90_handle_err(retval)
                    

          !4) 3 is the dimID of the y-coordinate
          !   get the extent of the 'y_map'
          !   and the extract the 'y_map'
          retval = NF90_INQUIRE_DIMENSION(ncid, 3, len=size_y_map)
          call nf90_handle_err(retval)
          allocate(y_map(size_y_map))
          retval = NF90_GET_VAR(ncid,y_varid,y_map)
          call nf90_handle_err(retval)


          !5) close netcdf file
          retval = NF90_CLOSE(ncid)
          call nf90_handle_err(retval)

        end subroutine nf90_get_maps


        !print the netcdf error received
        subroutine nf90_handle_err(nf90_err)

          implicit none

          integer, intent(in) ::  nf90_err
        
          if (nf90_err.ne.NF90_NOERR) then
             print *, 'Error: ', NF90_STRERROR(nf90_err)
             stop 2
          end if

        end subroutine nf90_handle_err

      end module nf90_operators_error_module
