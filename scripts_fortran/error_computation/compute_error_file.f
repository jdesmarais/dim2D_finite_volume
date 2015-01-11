      program compute_error_file

        use cmd_operators_error_class, only :
     $       cmd_operators_error

        use compare_type_module, only :
     $       compare_rkind

        use compute_error_module, only :
     $       compute_relative_error,
     $       get_index_coord

        use netcdf

        use nf90_operators_error_module, only :
     $       nf90_put_var_error,
     $       nf90_def_var_error,
     $       nf90_def_header_error,
     $       nf90_get_time,
     $       nf90_get_gov_var,
     $       nf90_get_maps,
     $       nf90_handle_err
        
        use parameters_cst, only :
     $       NOT_SUCCESS,
     $       SUCCESS

        use parameters_kind, only :
     $       ikind,
     $       rkind


        implicit none

        integer, parameter :: bc_size = 2

        character*(*), dimension(4), parameter :: var_name =
     $       ['mass','momentum_x','momentum_y','energy']

        type(cmd_operators_error) :: cmd_operators_used

        character(len=1024) :: filename_sm_domain
        character(len=1024) :: filename_lg_domain
        character(len=1024) :: filename_error    

        integer :: ierror

        real(rkind) :: time_sm_domain
        real(rkind) :: time_lg_domain
        real(rkind) :: time

        real(rkind), dimension(:), allocatable :: x_map_sm_domain
        real(rkind), dimension(:), allocatable :: y_map_sm_domain
        real(rkind), dimension(:), allocatable :: x_map_lg_domain
        real(rkind), dimension(:), allocatable :: y_map_lg_domain
        real(rkind), dimension(:), allocatable :: x_map
        real(rkind), dimension(:), allocatable :: y_map

        real(rkind) :: x_min_sm_domain
        real(rkind) :: x_max_sm_domain
        real(rkind) :: y_min_sm_domain
        real(rkind) :: y_max_sm_domain        

        integer(ikind) :: i_min
        integer(ikind) :: i_max
        integer(ikind) :: j_min
        integer(ikind) :: j_max

        real(rkind), dimension(:,:,:), allocatable :: var_sm_domain
        real(rkind), dimension(:,:,:), allocatable :: var_lg_domain

        real(rkind), dimension(:,:,:), allocatable :: error
        real(rkind), dimension(:)    , allocatable :: error_max
        real(rkind), dimension(:)    , allocatable :: x_error_max
        real(rkind), dimension(:)    , allocatable :: y_error_max   

        integer                            :: retval
        integer                            :: ncid_error
        integer, dimension(3)              :: coordinates_id
        integer, dimension(:), allocatable :: error_id
        integer, dimension(:), allocatable :: error_max_id
        integer, dimension(:), allocatable :: x_error_max_id
        integer, dimension(:), allocatable :: y_error_max_id


        !0) analyse the command line arguments to get the
        !   small, large domains and error filenames
        call cmd_operators_used%analyse_cmd_line_arg()

        if(.not.cmd_operators_used%are_files_provided()) then
           print '(''compute_error_file'')'
           print '(''***filenames are not provided***'')'
           stop ''
        else
           call cmd_operators_used%get_filename_sm_domain(
     $          filename_sm_domain)
           call cmd_operators_used%get_filename_lg_domain(
     $          filename_lg_domain)
           call cmd_operators_used%get_filename_error(
     $          filename_error)
        end if

        !1) extract the time from the small and the large domains
        !   and compare if the two files correspond to the same
        !   time
        time_sm_domain = nf90_get_time(
     $       trim(filename_sm_domain),
     $       ierror)
        if(ierror.ne.SUCCESS) then
           print '(''compute_error_file'')'
           print '(''error in extracting time from small domain'')'
           stop ''
        end if

        time_lg_domain = nf90_get_time(
     $       trim(filename_lg_domain),
     $       ierror)
        if(ierror.ne.SUCCESS) then
           print '(''compute_error_file'')'
           print '(''error in extracting time from large domain'')'
           stop ''
        end if

        if(compare_rkind(time_sm_domain,time_lg_domain)) then
           time = time_sm_domain
        else
           print '(''compute_error_file'')'
           print '(''***the files do not match the same time***'')'
           stop ''
        end if


        !2) extract the maps from the two simulations and find the 
        !   indices for the x_map and y_map in the large domain
        !   simulation such that these smaller x_map(i_min:i_max)
        !   y_map(j_min,j_max) match the maps of the simulation
        !   on the small domain
        call nf90_get_maps(
     $       trim(filename_sm_domain),
     $       x_map_sm_domain,
     $       y_map_sm_domain)

        call nf90_get_maps(
     $       trim(filename_lg_domain),
     $       x_map_lg_domain,
     $       y_map_lg_domain)

        x_min_sm_domain = x_map_sm_domain(bc_size+1)
        x_max_sm_domain = x_map_sm_domain(size(x_map_sm_domain,1)-bc_size)
        y_min_sm_domain = y_map_sm_domain(bc_size+1)
        y_max_sm_domain = y_map_sm_domain(size(y_map_sm_domain,1)-bc_size)

        i_min = get_index_coord(
     $       x_min_sm_domain,
     $       x_map_lg_domain,
     $       side=.false.,
     $       ierror=ierror)
        if(ierror.ne.SUCCESS) then
           print '(''compute_error_file'')'
           print '(''***error when matching x_min in large domain***'')'
           stop ''
        end if

        i_max = get_index_coord(
     $       x_max_sm_domain,
     $       x_map_lg_domain,
     $       side=.true.,
     $       ierror=ierror)
        if(ierror.ne.SUCCESS) then
           print '(''compute_error_file'')'
           print '(''***error when matching x_max in large domain***'')'
           stop ''
        end if

        j_min = get_index_coord(
     $       y_min_sm_domain,
     $       y_map_lg_domain,
     $       side=.false.,
     $       ierror=ierror)
        if(ierror.ne.SUCCESS) then
           print '(''compute_error_file'')'
           print '(''***error when matching y_min in large domain***'')'
           stop ''
        end if

        j_max = get_index_coord(
     $       y_max_sm_domain,
     $       y_map_lg_domain,
     $       side=.true.,
     $       ierror=ierror)
        if(ierror.ne.SUCCESS) then
           print '(''compute_error_file'')'
           print '(''***error when matching y_max in large domain***'')'
           stop ''
        end if


        !3) only retain the coordinate maps for [x_min:x_max]x[y_min,y_max]
        allocate(x_map(i_max-i_min+1))
        x_map = x_map_sm_domain(bc_size+1:size(x_map_sm_domain,1)-bc_size)

        allocate(y_map(j_max-j_min+1))
        y_map = y_map_sm_domain(bc_size+1:size(y_map_sm_domain,1)-bc_size)

        
        !4) extract the governing variables corresponding to
        !   [x_min:x_max]x[y_min:y_max]
        call nf90_get_gov_var(
     $       trim(filename_sm_domain),
     $       var_name,
     $       var_sm_domain,
     $       start=[bc_size+1,bc_size+1],
     $       count=[size(x_map_sm_domain)-2*bc_size,
     $              size(y_map_sm_domain)-2*bc_size],
     $       ierror=ierror)
        if(ierror.ne.SUCCESS) then
           print '(''compute_error_file'')'
           print '(''***error extracting gov var on small domain***'')'
           stop ''
        end if

        call nf90_get_gov_var(
     $       trim(filename_lg_domain),
     $       var_name,
     $       var_lg_domain,
     $       start=[i_min,j_min],
     $       count=[i_max-i_min+1,j_max-j_min+1],
     $       ierror=ierror)
        if(ierror.ne.SUCCESS) then
           print '(''compute_error_file'')'
           print '(''***error extracting gov var on small domain***'')'
           stop ''
        end if


        !5) compute the relative error for each governing variable
        call compute_relative_error(
     $       var_sm_domain,
     $       var_lg_domain,
     $       error,
     $       error_max=error_max,
     $       x_error_max=x_error_max,
     $       y_error_max=y_error_max,
     $       x_map=x_map,
     $       y_map=y_map,
     $       ierror=ierror)


        !6) deallocate unnecessary variables
        deallocate(x_map_sm_domain)
        deallocate(y_map_sm_domain)
        deallocate(var_sm_domain)

        deallocate(x_map_lg_domain)
        deallocate(y_map_lg_domain)
        deallocate(var_lg_domain)


        !7) define the header for the error file by
        !   combining the header of the small and
        !   large domain files
        call nf90_def_header_error(
     $       trim(filename_sm_domain),
     $       trim(filename_lg_domain),
     $       trim(filename_error),
     $       ncid_error,
     $       ierror=ierror)
        if(ierror.ne.SUCCESS) then
           print '(''compute_error_file'')'
           print '(''***error defining the error file header***'')'
           stop ''
        end if


        !8) define the variables stored in the error
        !   file
        call nf90_def_var_error(
     $       trim(filename_sm_domain),
     $       ncid_error,
     $       size(x_map,1),
     $       size(y_map,1),
     $       size(var_name,1),
     $       coordinates_id,
     $       error_id,
     $       error_max_id,
     $       x_error_max_id,
     $       y_error_max_id)


        !9) put the variables in the error file
        call nf90_put_var_error(
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


        !10) close the error file
        retval = NF90_CLOSE(ncid_error)
        call nf90_handle_err(retval)

      end program compute_error_file
