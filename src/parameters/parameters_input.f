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
      !> 20_08_2013 - initial version   - J.L. Desmarais
      !-----------------------------------------------------------------
      module parameters_input

        use parameters_constant
        use parameters_kind, only : ikind, rkind

        implicit none

        !<computational field dimensions
        real(rkind), parameter :: x_min = -0.4000000000d0
        real(rkind), parameter :: x_max = 0.4000000000d0
        real(rkind), parameter :: y_min = -0.4000000000d0
        real(rkind), parameter :: y_max = 0.4000000000d0
        
        !<computational times
        real(rkind), parameter :: t_max = 0.0002100000d0
        real(rkind), parameter :: dt = 0.0000007000d0
        
        !<output writing
        real(rkind), parameter :: detail_print = 0.0000000000d0

        !<mpi choice
        integer, parameter :: npx = 1 !<number of processors along x
        integer, parameter :: npy = 2 !<number of processors along y

        !<size of the main tables
        !<careful, choose ne according to the physical model
        integer(ikind), parameter :: ntx = 804
        integer(ikind), parameter :: nty = 808

        integer(ikind), parameter :: nx = ntx/npx
        integer(ikind), parameter :: ny = nty/npy
        integer       , parameter :: ne = 4

        !<initial conditions choice
        integer, parameter :: ic_choice = drop_retraction

        !<boundary conditions choice
        integer, parameter :: bc_choice = periodic_xy_choice

        !<output choice
        integer, parameter :: io_choice   = netcdf_choice

      end module parameters_input
