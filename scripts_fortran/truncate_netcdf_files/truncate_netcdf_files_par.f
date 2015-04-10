      program truncate_netcdf_files_par

        use cmd_operators_truncate_class, only :
     $       cmd_operators_truncate

        use mpi

        use nf90_operators_truncate_module, only :
     $       nf90_truncate_files

        use parameters_kind, only :
     $       rkind

        implicit none

        type(cmd_operators_truncate) :: cmd_op
        real(rkind), dimension(2,2)  :: borders_tr
        integer                      :: ne
        integer                      :: nb_timesteps

        integer :: ierror
        integer :: rank
        integer :: nb_procs

        
        ! initialize the MPI processes, get the rank and
        ! the number of processors
        call MPI_INIT(ierror)
        if(ierror.ne.MPI_SUCCESS) then
           print '(''MPI_INIT failed'')'
           stop ''
        end if

        call MPI_COMM_RANK(MPI_COMM_WORLD, rank, ierror)
        if(ierror.ne.MPI_SUCCESS) then
           call MPI_FINALIZE(ierror)
           print '(''MPI_COMM_RANK failed'')'
           stop ''
        end if

        call MPI_COMM_SIZE(MPI_COMM_WORLD, nb_procs, ierror)
        if(ierror.ne.MPI_SUCCESS) then
           call MPI_FINALIZE(ierror)
           print '(''MPI_COMM_SIZE failed'')'
           stop ''
        end if


        ! analyse the command line arguments
        call cmd_op%analyse_cmd_line_arg()
        
        borders_tr   = cmd_op%get_borders()
        ne           = cmd_op%get_ne()
        nb_timesteps = cmd_op%get_nb_timesteps()


        if(rank.eq.0) then
           print *, 'ne: ', ne
           print *, 'borders_tr: ', borders_tr
           print *, 'nb_timesteps: ', nb_timesteps
           print *, 'nb_procs: ', nb_procs
        end if


        ! truncate all files corresponding to the simulation
        call nf90_truncate_files(
     $       ne,
     $       nb_timesteps,
     $       borders_tr,
     $       timestep_start=rank,
     $       timestep_increment=nb_procs)

        
        ! ask all processes to wait
        call MPI_BARRIER(MPI_COMM_WORLD, ierror)
        if(ierror.ne.MPI_SUCCESS) then
           print '(''MPI_BARRIER failed'')'
           stop ''
        end if


        ! finalize the MPI processes
        call MPI_FINALIZE(ierror)
        if(ierror.ne.MPI_SUCCESS) then
           print '(''MPI_FINALIZE failed'')'
           stop ''
        end if

      end program truncate_netcdf_files_par
