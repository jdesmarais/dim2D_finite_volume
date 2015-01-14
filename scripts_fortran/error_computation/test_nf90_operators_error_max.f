      program test_nf90_nf90_operators_error_max

        use netcdf

        use nf90_operators_error_module, only :
     $       nf90_handle_err
        
        use nf90_operators_error_max_module, only :
     $       nf90_extract_error_max,
     $       nf90_extract_error_max_in_time,
     $       nf90_def_header_error_max,
     $       nf90_def_var_error_max,
     $       nf90_put_var_error_max

        use parameters_kind, only :
     $       ikind, rkind

        implicit none

        
        character*(*), parameter  :: filename = 'error0.nc'
        character*(*), dimension(4), parameter :: var_name =
     $       ['mass','momentum_x','momentum_y','energy']
        real(rkind)               :: time
        real(rkind), dimension(4) :: error_max
        real(rkind), dimension(4) :: x_error_max
        real(rkind), dimension(4) :: y_error_max
        integer                   :: ierror
        
        
        character*(*), parameter  :: errorDir = './error_test'
        real(rkind), dimension(:)  , allocatable :: time_array
        real(rkind), dimension(:,:), allocatable :: error_max_time
        real(rkind), dimension(:)  , allocatable :: t_error_max_in_time
        real(rkind), dimension(:)  , allocatable :: error_max_in_time
        real(rkind), dimension(:,:), allocatable :: x_error_max_time
        real(rkind), dimension(:,:), allocatable :: y_error_max_time
        integer                                  :: nb_files


        integer  :: retval
        integer  :: ncid


        integer                            :: time_id
        integer, dimension(:), allocatable :: error_max_id
        integer, dimension(:), allocatable :: t_error_max_in_time_id
        integer, dimension(:), allocatable :: error_max_in_time_id
        integer, dimension(:), allocatable :: x_error_max_id
        integer, dimension(:), allocatable :: y_error_max_id


        !test: nf90_extract_error_max
        print '(''test: nf90_extract_error_max'')'
        call nf90_extract_error_max(
     $       filename,
     $       var_name,
     $       time,
     $       error_max,
     $       x_error_max,
     $       y_error_max,
     $       ierror=ierror)

        
        print '(''max_error:   '',4F8.4)', error_max
        print '(''x_max_error: '',4F8.4)', x_error_max
        print '(''y_max_error: '',4F8.4)', y_error_max
        print ''


        !test: nf90_extract_error_max_in_time
        print '(''test: nf90_extract_error_max_in_time'')'
        nb_files = 2
        call nf90_extract_error_max_in_time(
     $       errorDir,
     $       var_name,
     $       nb_files,
     $       time_array,
     $       error_max_time,
     $       t_error_max_in_time,
     $       error_max_in_time,
     $       x_error_max=x_error_max_time,
     $       y_error_max=y_error_max_time,
     $       ierror=ierror)

        print '(''time_array: '',2F8.3)', time_array
        print '(''error_max_in_time(:,1): '',2F8.3)', error_max_time(1,:)
        print ''


        !test: nf90_def_header_error_max
        print '(''test: nf90_def_header_error_max'')'

        retval = NF90_CREATE(
     $         'error_max.nc',
     $         NF90_NETCDF4,
     $         ncid)
        call nf90_handle_err(retval)

        call nf90_def_header_error_max(
     $       errorDir,
     $       ncid)


        !test: nf90_def_var_error_max
        print '(''test: nf90_def_var_error_max'')'

        call nf90_def_var_error_max(
     $       ncid,
     $       nb_files,
     $       var_name,
     $       time_id,
     $       error_max_id,
     $       t_error_max_in_time_id,
     $       error_max_in_time_id,
     $       x_error_max_id=x_error_max_id,
     $       y_error_max_id=y_error_max_id)


        !test: nf90_put_var_error_max
        print '(''test: nf90_put_var_error_max'')'

        call nf90_put_var_error_max(
     $       ncid,
     $       time_id,
     $       error_max_id,
     $       t_error_max_in_time_id,
     $       error_max_in_time_id,
     $       time_array,
     $       error_max_time,
     $       t_error_max_in_time,
     $       error_max_in_time,
     $       x_error_max_id=x_error_max_id,
     $       y_error_max_id=y_error_max_id,
     $       x_error_max=x_error_max_time,
     $       y_error_max=y_error_max_time)


        retval = NF90_CLOSE(ncid)
        call nf90_handle_err(retval)

      end program test_nf90_nf90_operators_error_max
