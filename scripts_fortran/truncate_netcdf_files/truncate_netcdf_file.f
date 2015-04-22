      program truncate_netcdf_file

        use cmd_operators_truncate_class, only :
     $       cmd_operators_truncate

        use nf90_operators_truncate_module, only :
     $       nf90_truncate_files

        use parameters_kind, only :
     $       rkind

        implicit none

        type(cmd_operators_truncate) :: cmd_op
        real(rkind), dimension(2,2)  :: borders_tr
        integer                      :: ne
        integer                      :: nb_timesteps

        
        ! analyse the command line arguments
        call cmd_op%analyse_cmd_line_arg()
        
        borders_tr   = cmd_op%get_borders()
        ne           = cmd_op%get_ne()
        nb_timesteps = cmd_op%get_nb_timesteps()


        print *, 'ne: ', ne
        print *, 'borders_tr: ', borders_tr
        print *, 'nb_timesteps: ', nb_timesteps


        ! truncate all files corresponding to the simulation
        call nf90_truncate_files(
     $       ne,
     $       nb_timesteps+1,
     $       borders_tr,
     $       timestep_start=nb_timesteps,
     $       timestep_increment=1)

      end program truncate_netcdf_file
