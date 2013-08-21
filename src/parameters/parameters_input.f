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
      !> 20_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module parameters_input

        use parameters_constant
        use parameters_kind, only : ikind

        implicit none

        !<mpi choice
        integer, parameter :: npx = 2 !<number of processors along x
        integer, parameter :: npy = 2 !<number of processors along y

        !<size of the main tables
        !<careful, choose ne according to the physical model
        integer(ikind), parameter :: nx = 84
        integer(ikind), parameter :: ny = 84
        integer       , parameter :: ne = 4

        !<initial conditions choice
        integer, parameter :: ic_choice = drop_retraction

        !<boundary conditions choice
        integer, parameter :: bc_choice = periodic_xy_choice

        !<output choice
        integer, parameter :: io_choice   = netcdf_choice

      end module parameters_input
