      program merge_netcdf_files_par
      
        use cmd_operators_merge_class, only :
     $     cmd_operators_merge

        use mpi

        use nf90_operators_merge_module, only :
     $     nf90_merge_files

        implicit none

        type(cmd_operators_merge) :: cmd_op
        integer, dimension(2)     :: nb_tiles
        integer                   :: ne
        integer                   :: bc_size
        integer                   :: nb_timesteps

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
        
        nb_tiles(1)  = cmd_op%get_nb_tiles_x()
        nb_tiles(2)  = cmd_op%get_nb_tiles_y()
        ne           = cmd_op%get_ne()
        bc_size      = cmd_op%get_bc_size()
        nb_timesteps = cmd_op%get_nb_timesteps()


        if(rank.eq.0) then
           print *, 'nb_tiles: ', nb_tiles
           print *, 'ne: ', ne
           print *, 'bc_size: ', bc_size
           print *, 'nb_timesteps: ', nb_timesteps
           print *, 'nb_procs: ', nb_procs
        end if


        ! merge all files corresponding to the simulation
        call nf90_merge_files(
     $       nb_tiles,
     $       ne,bc_size,
     $       nb_timesteps,
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

      end program merge_netcdf_files_par
