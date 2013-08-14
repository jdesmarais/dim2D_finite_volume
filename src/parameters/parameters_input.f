      !> @file
      !> module containing the user choices defined as constants
      !> and propagated at compilation time
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> user choices required at compilation time
      !
      !> @date
      !> 13_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module parameters_input

        use parameters_constant

        !<initial conditions choice
        integer, parameter :: ic_choice = drop_retraction

        !<boundary conditions choice
        integer, parameter :: bc_choice = periodic_xy_choice

        !<output choice
        integer, parameter :: io_choice   = netcdf_choice

      end module parameters_input
