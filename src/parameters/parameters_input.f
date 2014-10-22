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

        !< debug option allowing extra checks in the code
        logical    , parameter :: debug = .true.        

        !<computational field dimensions
        real(rkind), parameter :: x_min = -3.0000000000d0
        real(rkind), parameter :: x_max = 0.0000000000d0
        real(rkind), parameter :: y_min = 0.0000000000d0
        real(rkind), parameter :: y_max = 3.0000000000d0
        
        !<computational times
        real(rkind), parameter :: t_max = 1.5000000000d0
        real(rkind), parameter :: dt = 0.02d0
        
        !<output writing
        real(rkind), parameter :: detail_print = 0.0150000000d0

        !<mpi choice
        integer, parameter :: npx = 1 !<number of processors along x
        integer, parameter :: npy = 1 !<number of processors along y

        !<size of the main tables
        !<careful, choose ne according to the physical model
        integer(ikind), parameter :: ntx = 5
        integer(ikind), parameter :: nty = 5

        integer(ikind), parameter :: nx = ntx/npx
        integer(ikind), parameter :: ny = nty/npy
        integer       , parameter :: ne = 4
        integer       , parameter :: bc_size = 2

        !<initial conditions choice
        !--------------------------------------------
        !ic_choice:
        !--------------------------------------------
        !
        !for NS equations
        !--------------------------------------------
        !steady_state       : constant everywhere
        !peak               : peak in the center of the domain
        !vortex             : vortex in the center of the domain
        !
        !for DIM equations
        !--------------------------------------------
        !steady_state       : constant everywhere
        !drop_retraction    : ellipsoidal droplet
        !bubble_ascending   : initial bubble
        !homogeneous_liquid : constant liquid density
        !phase_separation   : unstable mass density
        !--------------------------------------------
        integer, parameter :: flow_direction = x_direction
        integer, parameter :: ic_choice = vortex

        !<body forces choice
        integer, parameter :: gravity_choice = no_gravity_choice

        !<boundary conditions choice
        integer, parameter :: bc_choice = periodic_xy_choice
        integer, parameter :: bcx_type_choice = bc_nodes_choice
        integer, parameter :: bcy_type_choice = bc_nodes_choice

        !<output choice
        integer, parameter :: io_choice   = netcdf_choice

        !< boundary conditions parameters
        !--------------------------------------------
        !search_nb_dt: for the open boundary conditions,
        !              number of timesteps checked in
        !              advance by the increasing detector
        !
        !search_dcr  : radius expressed as number of grid
        !              points to check around the line for
        !              removing a buffer layer
        !
        !sigma_P     : relaxation coefficient used when
        !              applying the non-reflecting outflow
        !              pressure b.c.
        !--------------------------------------------
        real(rkind), parameter :: search_nb_dt = 0.5 !1.0 !0.0001000000d0
        integer    , parameter :: search_dcr = 4
        real(rkind), parameter :: sigma_P =  0.278d0 !0.25d0
        integer    , parameter :: obc_type_N = always_outflow
        integer    , parameter :: obc_type_S = always_inflow
        integer    , parameter :: obc_type_E = always_outflow
        integer    , parameter :: obc_type_W = always_outflow

      end module parameters_input
