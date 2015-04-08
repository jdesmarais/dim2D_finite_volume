      program merge_netcdf_files
      
        use cmd_operators_merge_class, only :
     $     cmd_operators_merge

        use nf90_operators_merge_module, only :
     $     nf90_merge_files

        implicit none

        type(cmd_operators_merge) :: cmd_op
        integer, dimension(2)     :: nb_tiles
        integer                   :: ne
        integer                   :: bc_size
        integer                   :: nb_timesteps

        
        ! analyse the command line arguments
        call cmd_op%analyse_cmd_line_arg()
        
        nb_tiles(1)  = cmd_op%get_nb_tiles_x()
        nb_tiles(2)  = cmd_op%get_nb_tiles_y()
        ne           = cmd_op%get_ne()
        bc_size      = cmd_op%get_bc_size()
        nb_timesteps = cmd_op%get_nb_timesteps()


        ! merge all files corresponding to the simulation
        call nf90_merge_files(
     $       nb_tiles,
     $       ne,bc_size,
     $       nb_timesteps)

      end program merge_netcdf_files
