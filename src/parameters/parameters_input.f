      !> @file
      !> module containing the constants used in the
      !> program: they will be propagated at compilation
      !> time to select the user choice
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

        !<boundary conditions choice
        integer, parameter :: bc_choice = periodic_xy_choice

      end module parameters_input
