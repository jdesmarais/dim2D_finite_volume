      module nf90_operators_error_max_module
      
        use netcdf

        use nf90_operators_error_module, only :
     $       nf90_put_att_value,
     $       nf90_get_att_value,
     $       nf90_handle_err,
     $       nf90_get_myreal

        use parameters_cst, only :
     $       SUCCESS,
     $       NOT_SUCCESS

        use parameters_kind, only :
     $       rkind


        implicit none

        
        private
        public ::
     $       nf90_put_var_error_max,
     $       nf90_def_var_error_max,
     $       nf90_def_header_error_max,
     $       nf90_extract_error_max_in_time,
     $       nf90_extract_error_max,
     $       nf90_extract_error_var


        contains

        ! put the variables in the 'error_max.nc' file
        subroutine nf90_put_var_error_max(
     $       ncid,
     $       time_id,
     $       error_max_id,
     $       t_error_max_in_time_id,
     $       error_max_in_time_id,
     $       time,
     $       error_max,
     $       t_error_max_in_time,
     $       error_max_in_time,
     $       x_error_max_id,
     $       y_error_max_id,
     $       x_error_max,
     $       y_error_max)

          implicit none

          integer                              , intent(in) :: ncid
          integer                              , intent(in) :: time_id
          integer    , dimension(:)            , intent(in) :: error_max_id
          integer    , dimension(:)            , intent(in) :: t_error_max_in_time_id
          integer    , dimension(:)            , intent(in) :: error_max_in_time_id
                                                
          real(rkind), dimension(:)            , intent(in) :: time
          real(rkind), dimension(:,:)          , intent(in) :: error_max          
          real(rkind), dimension(:)            , intent(in) :: t_error_max_in_time
          real(rkind), dimension(:)            , intent(in) :: error_max_in_time
          
          integer    , dimension(:)  , optional, intent(in) :: x_error_max_id
          integer    , dimension(:)  , optional, intent(in) :: y_error_max_id
          real(rkind), dimension(:,:), optional, intent(in) :: x_error_max
          real(rkind), dimension(:,:), optional, intent(in) :: y_error_max


          integer :: retval
          integer :: k


          ! put the time
          retval = NF90_PUT_VAR(ncid,time_id,time)
          call nf90_handle_err(retval)


          ! put the maximum error of the governing
          ! variables over space
          do k=1, size(error_max_id,1)

             retval = NF90_PUT_VAR(
     $            ncid,
     $            error_max_id(k),
     $            error_max(k,:))
             call nf90_handle_err(retval)

          end do
   

          ! put the time at which the maximum error
          ! of the governing variables over space
          ! and time is reached
          do k=1, size(t_error_max_in_time,1)

             retval = NF90_PUT_VAR(
     $            ncid,
     $            t_error_max_in_time_id(k),
     $            t_error_max_in_time(k))
             call nf90_handle_err(retval)

          end do
          
          
          ! put the the maximum error of the
          ! governing variables over space
          ! and time
          do k=1, size(error_max_in_time,1)

             retval = NF90_PUT_VAR(
     $            ncid,
     $            error_max_in_time_id(k),
     $            error_max_in_time(k))
             call nf90_handle_err(retval)

          end do
          

          ! if asked by the user, put the x-coordinate
          ! where the maximum of the error of the governing
          ! variable over space is reached
          if(present(x_error_max_id).and.present(x_error_max)) then

             do k=1, size(x_error_max_id,1)
                
                retval = NF90_PUT_VAR(
     $               ncid,
     $               x_error_max_id(k),
     $               x_error_max(k,:))
                call nf90_handle_err(retval)
                
             end do

          end if


          ! if asked by the user, put the y-coordinate
          ! where the maximum of the error of the governing
          ! variable over space is reached
          if(present(y_error_max_id).and.present(y_error_max)) then

             do k=1, size(y_error_max_id,1)
                
                retval = NF90_PUT_VAR(
     $               ncid,
     $               y_error_max_id(k),
     $               y_error_max(k,:))
                call nf90_handle_err(retval)
                
             end do

          end if

        end subroutine nf90_put_var_error_max


        ! define the variables that will be stored in the
        ! 'error_max.nc' file
        subroutine nf90_def_var_error_max(
     $       ncid,
     $       nt,
     $       var_name,
     $       time_id,
     $       error_max_id,
     $       t_error_max_in_time_id,
     $       error_max_in_time_id,
     $       x_error_max_id,
     $       y_error_max_id)

          implicit none

          integer                                          , intent(in)  :: ncid
          integer                                          , intent(in)  :: nt
          character(*), dimension(:)                       , intent(in)  :: var_name
          integer                                          , intent(out) :: time_id
          integer     , dimension(:), allocatable          , intent(out) :: error_max_id
          integer     , dimension(:), allocatable          , intent(out) :: t_error_max_in_time_id
          integer     , dimension(:), allocatable          , intent(out) :: error_max_in_time_id
          integer     , dimension(:), allocatable, optional, intent(out) :: x_error_max_id
          integer     , dimension(:), allocatable, optional, intent(out) :: y_error_max_id
          

          integer :: retval

          integer               :: NF_MYREAL
          integer, dimension(2) :: dim_id


          !1) define the type of variables stored
          NF_MYREAL = nf90_get_myreal(rkind)


          !2) define the time as the main coordinate
          retval = NF90_DEF_DIM(
     $         ncid,
     $         'time',
     $         nt,
     $         dim_id(1))
          call nf90_handle_err(retval)


          !3) define the variable 'time' coorresponding
          !   to the coordinate 'time'
          retval = NF90_DEF_VAR(
     $         ncid,
     $         'time',
     $         NF_MYREAL,
     $         dim_id(1),
     $         time_id)
          call nf90_handle_err(retval)


          !4) define the attributes of 'time'
          retval = NF90_PUT_ATT(
     $         ncid,
     $         time_id,
     $         'units',
     $         's/s')
          call nf90_handle_err(retval)

          retval = NF90_PUT_ATT(
     $         ncid,
     $         time_id,
     $         'long_name',
     $         'reduced time')
          call nf90_handle_err(retval)


          !5) define 'max' as the second coordinate
          !   for the maximum of the error over time
          !   and space
          retval = NF90_DEF_DIM(
     $         ncid,
     $         'max',
     $         1,
     $         dim_id(2))
          call nf90_handle_err(retval)

          
          !5) define the maximum of the error
          !   over space of the governing variables
          !   as variables depending on time
          allocate(error_max_id(size(var_name,1)))
          call nf90_def_var_as_timeFct(
     $         ncid,
     $         var_name,
     $         NF_MYREAL,
     $         [dim_id(1)],
     $         'max_error_',
     $         '',
     $         '-',
     $         'maximum of error_',
     $         ' over space',
     $         error_max_id)


          !6) define the time where the 
          !   maximum of the error
          !   over space and over time
          allocate(t_error_max_in_time_id(size(var_name,1)))
          call nf90_def_var_as_timeFct(
     $         ncid,
     $         var_name,
     $         NF_MYREAL,
     $         [dim_id(2)],
     $         't_max_error_in_time_',
     $         '',
     $         '-',
     $         'time where the maximum over time and space of error_',
     $         ' is reached',
     $         t_error_max_in_time_id)

          
          !7) define the maximum of the error
          !   over space and over time
          allocate(error_max_in_time_id(size(var_name,1)))
          call nf90_def_var_as_timeFct(
     $         ncid,
     $         var_name,
     $         NF_MYREAL,
     $         [dim_id(2)],
     $         'max_error_in_time_',
     $         '',
     $         '-',
     $         'maximum of error_',
     $         ' over time and space',
     $         error_max_in_time_id)


          !8) define the x-coordinate where the
          !   maximum error is reached
          !   as function of time
          if(present(x_error_max_id)) then

             allocate(x_error_max_id(size(var_name,1)))
             call nf90_def_var_as_timeFct(
     $            ncid,
     $            var_name,
     $            NF_MYREAL,
     $            [dim_id(1)],
     $            'x_max_error_',
     $            '',
     $            '-',
     $            'x-coordinate where the maximum over space of error_',
     $            ' is reached',
     $            x_error_max_id)
             
          end if

          
          !9) define the y-coordinate where the
          !   maximum error is reached
          !   as function of time
          if(present(y_error_max_id)) then

             allocate(y_error_max_id(size(var_name,1)))
             call nf90_def_var_as_timeFct(
     $            ncid,
     $            var_name,
     $            NF_MYREAL,
     $            [dim_id(1)],
     $            'y_max_error_',
     $            '',
     $            '-',
     $            'y-coordinate where the maximum over space of error_',
     $            ' is reached',
     $            y_error_max_id)
             
          end if


          !10) stop the definition of the variables
          !    saved in the 'error_max.nc' file
          retval = NF90_ENDDEF(ncid)
          call nf90_handle_err(retval)

        end subroutine nf90_def_var_error_max


        subroutine nf90_def_var_as_timeFct(
     $     ncid,
     $     var_name,
     $     var_type,
     $     dim_ids,
     $     var_prefix,
     $     var_suffix,
     $     units,
     $     longname_prefix,
     $     longname_suffix,
     $     var_id)

          implicit none

          integer                    , intent(in)  :: ncid
          character*(*), dimension(:), intent(in)  :: var_name
          integer                    , intent(in)  :: var_type
          integer      , dimension(:), intent(in)  :: dim_ids
          character*(*)              , intent(in)  :: var_prefix
          character*(*)              , intent(in)  :: var_suffix
          character*(*)              , intent(in)  :: units
          character*(*)              , intent(in)  :: longname_prefix
          character*(*)              , intent(in)  :: longname_suffix
          integer      , dimension(:), intent(out) :: var_id
          
          integer             :: k
          character(len=1024) :: var_name_modified
          integer             :: retval

          do k=1, size(var_name,1)

             !define the governing variable
             call join_strings(
     $            [var_prefix,var_name(k),var_suffix],
     $            var_name_modified)

             retval = NF90_DEF_VAR(
     $            ncid,
     $            trim(var_name_modified),
     $            var_type,
     $            dim_ids,
     $            var_id(k))
             call nf90_handle_err(retval)


             !define the units
             retval = NF90_PUT_ATT(
     $            ncid,
     $            var_id(k),
     $            'units',
     $            trim(units))
             call nf90_handle_err(retval)


             !define the 'long_name'
             call join_strings(
     $            [longname_prefix,var_name(k),longname_suffix],
     $            var_name_modified)

             retval = NF90_PUT_ATT(
     $            ncid,
     $            var_id(k),
     $            'long_name',
     $            trim(var_name_modified))
             call nf90_handle_err(retval)

          end do

        end subroutine nf90_def_var_as_timeFct


        subroutine join_strings(
     $     words,
     $     sentence)

          implicit none

          character*(*), dimension(:), intent(in) :: words
          character(len=1024)                     :: sentence

          integer :: i_start
          integer :: k
          integer :: len_word

          
          sentence = ''
          i_start  = 1

          do k=1, size(words,1)

             len_word = len(trim(words(k)))
             sentence(i_start:i_start+len_word-1) = trim(words(k))
             i_start = i_start+len_word

          end do

        end subroutine join_strings


        subroutine nf90_def_header_error_max(
     $       errorDir,
     $       ncid)

          implicit none

          character*(*), intent(in) :: errorDir
          integer      , intent(in) :: ncid

          character(len=1024)   :: filename
          integer               :: retval
          integer               :: ncid_error_file
          integer               :: size_header_error_file

          character(len=29)     :: history
          integer, dimension(8) :: now_data
          character(len=60)     :: title

          integer               :: k

          character(len=1028)   :: att_name
          integer               :: att_type
          character(len=1028)   :: att_char_var
          integer               :: att_int_var
          real                  :: att_real_var
          real*8                :: att_double_var


          !1) open the 'error0.nc' netcdf file to
          !   copy most of the header except the
          !   'history' and 'title' global attributes
          call get_error_filename(
     $         errorDir,
     $         0,
     $         filename)

          retval = NF90_OPEN(
     $         trim(filename),
     $         NF90_NOWRITE,
     $         ncid_error_file)
          call nf90_handle_err(retval)


          !2) get the date for the header of the
          !   error file
          call date_and_time(values=now_data)
          write (history,
     $         '(''date'', 1x, i4.4, ''/'', i2.2, ''/'', i2.2, 1x,
     $           ''time'', 1x, i2.2, '':'', i2.2, '':'', i2.2)')
     $         now_data(1), now_data(2), now_data(3),
     $         now_data(5), now_data(6), now_data(7)


          !3) inquire the number of global
          !   variables (header variables) in
          !   'error.nc'
          retval = NF90_INQUIRE(
     $         ncid_error_file,
     $         nAttributes=size_header_error_file)
          call nf90_handle_err(retval)


          !4) inquire the name of the global
          !   variables for the error file
          !   store them in the new file except
          !   'history' and 'title' that are
          !   modified
          do k=1, size_header_error_file

             ! inquire the name of the global
             ! variables in the header
             retval = NF90_INQ_ATTNAME(
     $            ncid_error_file,
     $            NF90_GLOBAL,
     $            k,
     $            att_name)
             call nf90_handle_err(retval)


             ! inquire the type and length of
             ! the global variable
             retval = NF90_INQUIRE_ATTRIBUTE(
     $            ncid_error_file,
     $            NF90_GLOBAL,
     $            trim(att_name),
     $            xtype=att_type)
             call nf90_handle_err(retval)


             ! inquire the value of the global
             ! variable
             retval = nf90_get_att_value(
     $            ncid_error_file,
     $            NF90_GLOBAL,
     $            att_type,
     $            trim(att_name),
     $            att_char_var,
     $            att_int_var,
     $            att_real_var,
     $            att_double_var)
             call nf90_handle_err(retval)
                                    

             ! depending on the header attribute, the
             ! attribute is directly saved in the header
             ! of the error file or it is replaced or
             ! renamed
             ! (only modify the
             ! x_min,x_max,y_min,y_max attributes)
             select case(trim(att_name))

               case('title')
                  
                  title( 1:23) = 'DIM2D: maximum error : '
                  title(24:24+28) = 'small/large domain comparison'

                  retval = nf90_put_att(
     $                 ncid,
     $                 NF90_GLOBAL,
     $                 'title',
     $                 trim(title))


               case('history')

                  retval = nf90_put_att(
     $                 ncid,
     $                 NF90_GLOBAL,
     $                 'history',
     $                 history)


               ! for the other attributes, they are simply
               ! copy in the netcdf file
               case default

                  call nf90_put_att_value(
     $                 ncid,
     $                 NF90_GLOBAL,
     $                 att_type,
     $                 att_name,
     $                 att_char_var,
     $                 att_int_var,
     $                 att_real_var,
     $                 att_double_var)

             end select

          end do


          !5) close the netcdf file whose header
          !   has been copied
          retval = NF90_CLOSE(ncid_error_file)
          call nf90_handle_err(retval)

        end subroutine nf90_def_header_error_max


        subroutine nf90_extract_error_max_in_time(
     $       errorDir,
     $       var_name,
     $       nb_files,
     $       time,
     $       error_max,
     $       t_error_max_in_time,
     $       error_max_in_time,
     $       x_error_max,
     $       y_error_max,
     $       ierror)

          implicit none

          character*(*)                                       , intent(in)  :: errorDir
          character*(*), dimension(:)                         , intent(in)  :: var_name
          integer                                             , intent(in)  :: nb_Files
          real(rkind)  , dimension(:)            , allocatable, intent(out) :: time
          real(rkind)  , dimension(:,:)          , allocatable, intent(out) :: error_max
          real(rkind)  , dimension(:)            , allocatable, intent(out) :: t_error_max_in_time
          real(rkind)  , dimension(:)            , allocatable, intent(out) :: error_max_in_time
          real(rkind)  , dimension(:,:), optional, allocatable, intent(out) :: x_error_max
          real(rkind)  , dimension(:,:), optional, allocatable, intent(out) :: y_error_max
          integer                      , optional             , intent(out) :: ierror

          
          logical             :: dir_exists
          integer             :: len_destDir
          character(len=1024) :: filename
          integer             :: i


          !1) verify that the directoy exists
          INQUIRE (DIRECTORY=errorDir, EXIST=dir_exists)
          if(.not.dir_exists) then
             print '(''nf90_operators_error_max_module'')'
             print '(''nf90_extract_error_max_in_time'')'
             print '(''dir errorDir does not exists'')'
             print *, 'errorDir: ', errorDir
             if(present(ierror)) then
                ierror = NOT_SUCCESS
             end if
             stop ''
          end if


          !2) extract time, error_max, error_max_in_time, 
          !   x_error_max and y_error_max from netcdf files
          allocate(time(nb_files))

          allocate(error_max(size(var_name,1),nb_files))

          allocate(t_error_max_in_time(size(var_name,1)))
          allocate(error_max_in_time(size(var_name,1)))

          len_destDir = len(trim(errorDir))

          if(.not.present(x_error_max)) then

             if(.not.present(y_error_max)) then

                do i=0, nb_files-1
                   
                   ! determine the name of the error file
                   call get_error_filename(
     $                  errorDir,
     $                  i,
     $                  filename,
     $                  len_destDir_op=len_destDir)
                   
                   ! determine the error max, x_error_max, y_error_max
                   call nf90_extract_error_max(
     $                  trim(filename),
     $                  var_name,
     $                  time(i+1),
     $                  error_max(:,i+1))

                   ! determine the error_max_in_time
                   call determine_error_max_in_time(
     $                  time(i+1),
     $                  error_max(:,i+1),
     $                  i,
     $                  t_error_max_in_time,
     $                  error_max_in_time)

                end do

             else

                allocate(y_error_max(size(var_name,1),nb_files))

                do i=0, nb_files-1
                   
                   ! determine the name of the error file
                   call get_error_filename(
     $                  errorDir,
     $                  i,
     $                  filename,
     $                  len_destDir_op=len_destDir)
                   
                   ! determine the error max, x_error_max, y_error_max
                   call nf90_extract_error_max(
     $                  trim(filename),
     $                  var_name,
     $                  time(i+1),
     $                  error_max(:,i+1),
     $                  y_error_max=y_error_max(:,i+1))

                   ! determine the error_max_in_time
                   call determine_error_max_in_time(
     $                  time(i+1),
     $                  error_max(:,i+1),
     $                  i,
     $                  t_error_max_in_time,
     $                  error_max_in_time)

                end do

             end if

             
          else
             
             allocate(x_error_max(size(var_name,1),nb_files))

             if(.not.present(y_error_max)) then

                do i=0, nb_files-1
                   
                   ! determine the name of the error file
                   call get_error_filename(
     $                  errorDir,
     $                  i,
     $                  filename,
     $                  len_destDir_op=len_destDir)
                   
                   ! determine the error max, x_error_max, y_error_max
                   call nf90_extract_error_max(
     $                  trim(filename),
     $                  var_name,
     $                  time(i+1),
     $                  error_max(:,i+1),
     $                  x_error_max=x_error_max(:,i+1))

                   ! determine the error_max_in_time
                   call determine_error_max_in_time(
     $                  time(i+1),
     $                  error_max(:,i+1),
     $                  i,
     $                  t_error_max_in_time,
     $                  error_max_in_time)

                end do

             else

                allocate(y_error_max(size(var_name,1),nb_files))

                do i=0, nb_files-1
                   
                   ! determine the name of the error file
                   call get_error_filename(
     $                  errorDir,
     $                  i,
     $                  filename,
     $                  len_destDir_op=len_destDir)
                   
                   ! determine the error max, x_error_max, y_error_max
                   call nf90_extract_error_max(
     $                  trim(filename),
     $                  var_name,
     $                  time(i+1),
     $                  error_max(:,i+1),
     $                  x_error_max=x_error_max(:,i+1),
     $                  y_error_max=y_error_max(:,i+1))

                   ! determine the error_max_in_time
                   call determine_error_max_in_time(
     $                  time(i+1),
     $                  error_max(:,i+1),
     $                  i,
     $                  t_error_max_in_time,
     $                  error_max_in_time)

                end do

             end if

          end if

        end subroutine nf90_extract_error_max_in_time


        ! determine the error max in time
        subroutine determine_error_max_in_time(
     $     time,
     $     error_max,
     $     i,
     $     t_error_max_in_time,
     $     error_max_in_time)

          implicit none

          real(rkind)              , intent(in)    :: time
          real(rkind), dimension(:), intent(in)    :: error_max
          integer                  , intent(in)    :: i
          real(rkind), dimension(:), intent(inout) :: t_error_max_in_time
          real(rkind), dimension(:), intent(inout) :: error_max_in_time


          integer :: k


          if(i.eq.0) then

             error_max_in_time = error_max

          else
             do k=1, size(error_max_in_time,1)

                if(error_max_in_time(k).lt.error_max(k)) then

                   error_max_in_time(k)   = error_max(k)
                   t_error_max_in_time(k) = time

                end if

             end do

          end if

        end subroutine determine_error_max_in_time


        ! get the filename for the error file
        subroutine get_error_filename(
     $     destDir,
     $     timestep,
     $     filename,
     $     len_destDir_op)

          implicit none

          character*(*)      , intent(in)  :: destDir
          integer            , intent(in)  :: timestep
          character(len=1024), intent(out) :: filename
          integer, optional  , intent(in)  :: len_destDir_op


          integer           :: len_destDir
          integer           :: len_timestep
          character(len=16) :: format_string
          character(len=10) :: timestep_str
          integer           :: len_suffix

          ! determine the space for the directory path
          if(present(len_destDir_op)) then
             len_destDir = len_destDir_op
          else
             len_destDir = len(trim(destDir))
          end if


          ! create the string for the timestep
          if(timestep.eq.0) then
             len_timestep = 1
          else
             len_timestep = floor(log10(real(timestep)))+1
          end if

          write (format_string, "(A2,I1,A4)")
     $         '(I', len_timestep, ',A3)'

          write (timestep_str, trim(format_string)) timestep, '.nc'

          len_suffix = len(trim(timestep_str))


          ! create the filename
          filename = ''

          filename(1:len_destDir) = trim(destDir)
          filename(len_destDir+1:len_destDir+6) = '/error'
          filename(len_destDir+7:len_destDir+7+len_suffix-1) = trim(timestep_str)

        end subroutine get_error_filename


        subroutine nf90_extract_error_max(
     $       filename,
     $       var_name,
     $       time,
     $       error_max,
     $       x_error_max,
     $       y_error_max,
     $       ierror)

          implicit none

          character*(*)                        , intent(in)  :: filename
          character*(*), dimension(:)          , intent(in)  :: var_name
          real(rkind)                          , intent(out) :: time
          real(rkind)  , dimension(:)          , intent(out) :: error_max
          real(rkind)  , dimension(:), optional, intent(out) :: x_error_max
          real(rkind)  , dimension(:), optional, intent(out) :: y_error_max
          integer                    , optional, intent(out) :: ierror

          
          logical                   :: file_exists
          integer                   :: retval
          integer                   :: ncid
          integer                   :: time_id
          real(rkind), dimension(1) :: time_table
          

          !0) verify that the file 'filename' exists
          INQUIRE(FILE=trim(filename), EXIST=file_exists)
          if(.not.file_exists) then
             print '(''nf90_operators_error_max_module'')'
             print '(''nf90_extract_error_max'')'
             print '(''file filename does not exists'')'
             print *, 'filename: ', filename
             if(present(ierror)) then
                ierror = NOT_SUCCESS
             end if
             stop ''
          end if


          !1) open the netcdf file for reading
          retval = NF90_OPEN(trim(filename), NF90_NOWRITE, ncid)
          call nf90_handle_err(retval)


          !2) get the timeid for the time variable
          !   then, extract the time
          retval = NF90_INQ_VARID(
     $            ncid,
     $            'time',
     $            time_id)
          call nf90_handle_err(retval)

          retval = NF90_GET_VAR(
     $         ncid,
     $         time_id,
     $         time_table)
          call nf90_handle_err(retval)

          time = time_table(1)


          !3) get the varid for the error_max variables:
          !   the name of the maximum error over space for the
          !   governing variable var_name(i) is given by:
          !   'max_error_'+var_name(i)
          !   
          !   extract the error_max variable
          call nf90_extract_error_var(
     $         ncid,
     $         var_name,
     $         'max_error_',
     $         error_max)


          !4) if asked by the user, get the var_id for the
          !   x-coordinate of the maximum error over space:
          !   the name of the maximum error over space for the
          !   governing variable var_name(i) is given by:
          !   'x_max_error_'+var_name(i)
          !   
          !   extract the x_error_max variable
          if(present(x_error_max)) then

             call nf90_extract_error_var(
     $            ncid,
     $            var_name,
     $            'x_max_error_',
     $            x_error_max)
             
          end if


          !5) if asked by the user, get the var_id for the
          !   y-coordinate of the maximum error over space:
          !   the name of the maximum error over space for the
          !   governing variable var_name(i) is given by:
          !   'y_max_error_'+var_name(i)
          !   
          !   extract the y_error_max variable
          if(present(y_error_max)) then

             call nf90_extract_error_var(
     $            ncid,
     $            var_name,
     $            'y_max_error_',
     $            y_error_max)
             
          end if


          !6) close netcdf file
          retval = NF90_CLOSE(ncid)
          call nf90_handle_err(retval)


          if(present(ierror)) then
             ierror = SUCCESS
          end if

        end subroutine nf90_extract_error_max


        subroutine nf90_extract_error_var(
     $     ncid,
     $     var_name,
     $     prefix,
     $     error_var,
     $     ierror)

          implicit none

          integer                    , intent(in)  :: ncid
          character*(*), dimension(:), intent(in)  :: var_name
          character*(*)              , intent(in)  :: prefix
          real(rkind)  , dimension(:), intent(out) :: error_var
          integer      , optional    , intent(out) :: ierror

          integer           :: retval
          integer           :: i
          integer           :: len_prefix
          integer           :: len_var_name
          character(len=60) :: var_name_modified
          integer           :: var_id


          ! verify that the size of the error_var array is the
          ! same as the size of the var_name array
          if(size(var_name,1).ne.size(error_var,1)) then

             print '(''nf90_operators_error_max_module'')'
             print '(''nf90_extract_error_max'')'
             print '(''the size of error_var does not match'')'
             print '(''the size of var_name'')'
             if(present(ierror)) then
                ierror = NOT_SUCCESS
             end if             
             stop ''
          end if


          len_prefix   = len(trim(prefix))

          do i=1, size(var_name,1)

             len_var_name = len(trim(var_name(i)))

             var_name_modified = ''
             var_name_modified(1:len_prefix) = trim(prefix)
             var_name_modified(len_prefix+1:len_prefix+len_var_name) = trim(var_name(i))
             
             retval = NF90_INQ_VARID(
     $            ncid,
     $            trim(var_name_modified),
     $            var_id)
             call nf90_handle_err(retval)

             retval = NF90_GET_VAR(
     $            ncid,
     $            var_id,
     $            error_var(i))
             call nf90_handle_err(retval)
                
          end do

        end subroutine nf90_extract_error_var

      end module nf90_operators_error_max_module
