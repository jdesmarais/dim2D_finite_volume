      program compute_error_max_file

        use cmd_operators_error_max_class, only :
     $       cmd_operators_error_max

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

        
        type(cmd_operators_error_max) :: cmd_operators_used
        character(len=1024)           :: errorDir
        character(len=1024)           :: filename
        integer                       :: nb_files
        character(10), dimension(4)   :: var_name

        integer :: ierror

        real(rkind), dimension(:)  , allocatable :: time_array
        real(rkind), dimension(:,:), allocatable :: error_max_time
        real(rkind), dimension(:)  , allocatable :: t_error_max_over_time
        real(rkind), dimension(:)  , allocatable :: error_max_over_time
        real(rkind), dimension(:,:), allocatable :: x_error_max_time
        real(rkind), dimension(:,:), allocatable :: y_error_max_time

        integer :: retval
        integer :: ncid

        integer                            :: time_id
        integer, dimension(:), allocatable :: error_max_id
        integer, dimension(:), allocatable :: t_error_max_over_time_id
        integer, dimension(:), allocatable :: error_max_over_time_id
        integer, dimension(:), allocatable :: x_error_max_id
        integer, dimension(:), allocatable :: y_error_max_id


        var_name = ['mass', 'momentum_x', 'momentum_y', 'energy']


        !0) analyze the command line arguments to extract
        !   the path for the directory where the error
        !   files to be analyzed are located, the path for
        !   the output file where the maximum of the error
        !   is saved, and the total number of files to be
        !   analyzed
        call cmd_operators_used%analyse_cmd_line_arg()
        
        if(.not.cmd_operators_used%are_inputs_provided()) then
           print '(''***not all inputs are provided***'')'
           print ''
           call cmd_operators_used%display_help()
           stop ''
        end if

        call cmd_operators_used%get_directory_error_files(errorDir)
        call cmd_operators_used%get_filename_error_max(filename)
        call cmd_operators_used%get_nb_files(nb_files)


        !1) create array collecting the maximum of the
        !   error for the governing variables in time
        call nf90_extract_error_max_in_time(
     $       trim(errorDir),
     $       var_name,
     $       nb_files,
     $       time_array,
     $       error_max_time,
     $       t_error_max_over_time,
     $       error_max_over_time,
     $       x_error_max=x_error_max_time,
     $       y_error_max=y_error_max_time,
     $       ierror=ierror)

        
        !2) open the netcdf file where the arrays will
        !   be written
        retval = NF90_CREATE(
     $       trim(filename),
     $       NF90_NETCDF4,
     $       ncid)
        call nf90_handle_err(retval)


        !3) create the header of the netcdf file
        call nf90_def_header_error_max(
     $       errorDir,
     $       ncid)


        !4) define the variables saved in the
        !   netcdf file
        call nf90_def_var_error_max(
     $       ncid,
     $       nb_files,
     $       var_name,
     $       time_id,
     $       error_max_id,
     $       t_error_max_over_time_id,
     $       error_max_over_time_id,
     $       x_error_max_id=x_error_max_id,
     $       y_error_max_id=y_error_max_id)


        !5) put the arrays in the netcdf file
        call nf90_put_var_error_max(
     $       ncid,
     $       time_id,
     $       error_max_id,
     $       t_error_max_over_time_id,
     $       error_max_over_time_id,
     $       time_array,
     $       error_max_time,
     $       t_error_max_over_time,
     $       error_max_over_time,
     $       x_error_max_id=x_error_max_id,
     $       y_error_max_id=y_error_max_id,
     $       x_error_max=x_error_max_time,
     $       y_error_max=y_error_max_time)


        !6) close the netcdf file
        retval = NF90_CLOSE(ncid)
        call nf90_handle_err(retval)


        !7) deallocate the arrays
        deallocate(time_array)
        deallocate(error_max_time)
        deallocate(t_error_max_over_time)
        deallocate(error_max_over_time)
        deallocate(x_error_max_time)
        deallocate(y_error_max_time)

        deallocate(error_max_id)
        deallocate(t_error_max_over_time_id)
        deallocate(error_max_over_time_id)
        deallocate(x_error_max_id)
        deallocate(y_error_max_id)

        !8) error_max file created
        print '(''output file created'')'

      end program compute_error_max_file
