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
        real(rkind), parameter :: x_min = -0.1519000000d0
        real(rkind), parameter :: x_max = 0.1519000000d0
        real(rkind), parameter :: y_min = -0.1519000000d0
        real(rkind), parameter :: y_max = 0.1519000000d0
        
        !<computational times
        real(rkind), parameter :: t_max = 2.9831000000d0 !10.0d0
        real(rkind), parameter :: dt = 0.0004000000d0
        
        !<output writing
        real(rkind), parameter :: detail_print = 0.1341000000d0
        logical    , parameter :: write_domain_extension = .true.
        logical    , parameter :: write_detectors = .true.

        !<mpi choice
        integer, parameter :: npx = 1 !<number of processors along x
        integer, parameter :: npy = 1 !<number of processors along y

        !<size of the main tables
        !<careful, choose ne according to the physical model
        integer(ikind), parameter :: ntx = 67
        integer(ikind), parameter :: nty = 67

        integer(ikind), parameter :: nx = ntx/npx
        integer(ikind), parameter :: ny = nty/npy
        integer       , parameter :: ne = 4
        integer       , parameter :: bc_size = 2

        !<initial conditions choice
        !--------------------------------------------
        !flow_direction:
        !--------------------------------------------
        !x_direction  : from left to right
        !y_direction  : from bottom to up
        !xy_direction : from SW to NE corner
        !
        !--------------------------------------------
        !flow_x_side:
        !--------------------------------------------
        !-1.0d0 : from right to left
        !+1.0d0 : from left  to right
        !
        !--------------------------------------------
        !flow_y_side:
        !--------------------------------------------
        !-1.0d0 : from top to bottom
        !+1.0d0 : from bottom to top
        !
        !--------------------------------------------
        !flow_velocity:
        !--------------------------------------------
        !reduced velocity of the mean flow
        !
        !--------------------------------------------
        !T0:
        !--------------------------------------------
        !reduced temperature of the mean flow
        !
        !--------------------------------------------
        !ic_choice:
        !--------------------------------------------
        !
        !for wave2d equations
        !--------------------------------------------
        !peak               : peak in the center of the domain
        !negative_spot      : negative field in the center of the domain
        !
        !for NS equations
        !--------------------------------------------
        !steady_state       : constant everywhere
        !peak               : peak in the center of the domain
        !vortex             : vortex in the center of the domain
        !sym_x              : symmetry compared to the y-axis
        !sym_y              : symmetry compared to the x-axis
        !
        !for DIM equations
        !--------------------------------------------
        !steady_state       : constant everywhere
        !drop_retraction    : ellipsoidal droplet
        !bubble_ascending   : initial bubble
        !homogeneous_liquid : constant liquid density
        !phase_separation   : unstable mass density
        !--------------------------------------------
        integer    , parameter :: flow_direction = x_direction
        real(rkind), parameter :: flow_x_side = 1.0000000000d0
        real(rkind), parameter :: flow_y_side = 1.0000000000d0
        real(rkind), parameter :: flow_velocity = 0.1000000000d0
        
        real(rkind), parameter :: T0 = 0.9500000000d0

        integer    , parameter :: ic_choice = bubble_transported

        !<body forces choice
        integer, parameter :: gravity_choice = no_gravity_choice
        integer, parameter :: wave_forcing = no_wave_forcing

        !<boundary conditions choice
        integer, parameter :: bc_choice = hedstrom_xy_corners_choice

        !<output choice
        integer, parameter :: io_choice = netcdf_choice


        !< boundary conditions parameters
        !-----------------------------------------------------
        !constrained by the bc_choice
        !-----------------------------------------------------
        !
        !bc_N_type_choice : type of boundary condition applied
        !                   at the North boundary
        !bc_S_type_choice : type of boundary condition applied
        !                   at the South boundary
        !bc_E_type_choice : type of boundary condition applied
        !                   at the East boundary
        !bc_W_type_choice : type of boundary condition applied
        !                   at the West boundary
        !-----------------------------------------------------
        integer    , parameter :: bc_N_type_choice = bc_timedev_choice
        integer    , parameter :: bc_S_type_choice = bc_timedev_choice
        integer    , parameter :: bc_E_type_choice = bc_timedev_choice
        integer    , parameter :: bc_W_type_choice = bc_timedev_choice


        !-----------------------------------------------------
        !for the increase of the computational domain
        !-----------------------------------------------------
        !bf_openbc_md_threshold_ac : control whether the increase
        !                            of the computational domain
        !                            is also activated by the
        !                            value of the mass density at
        !                            the edge
        !
        !bf_openbc_md_threshold : the increase of the computational
        !                         domain is triggered if the mass
        !                         density is inside
        !                         mid = (\rho_vap+\rho_liq)/2
        !                         thr_vap = threshold*(mid-\rho_vap)
        !                         thr_liq = threshold*(\rho_liq-mid)
        !                         [\rho_vap+thr_vap, \rho_liq-thr_liq]
        !
        !-----------------------------------------------------
        logical    , parameter :: bf_openbc_md_threshold_ac = .true.
        real(rkind), parameter :: bf_openbc_md_threshold = 0.0001000000d0


        !-----------------------------------------------------
        !only effective for Yoo and Lodato b.c.
        !-----------------------------------------------------
        !
        !sigma_P    : relaxation coefficient used when
        !             applying the non-reflecting outflow
        !             pressure b.c.
        !
        !obc_type_N : type of boundary condition applied
        !             at the North boundary (always_outflow,
        !             always_inflow, ask_flow)
        !obc_type_S : type of boundary condition applied
        !             at the South boundary (always_outflow,
        !             always_inflow, ask_flow)
        !obc_type_E : type of boundary condition applied
        !             at the East boundary (always_outflow,
        !             always_inflow, ask_flow)
        !obc_type_W : type of boundary condition applied
        !             at the West boundary (always_outflow,
        !             always_inflow, ask_flow)
        !-----------------------------------------------------
        real(rkind), parameter :: sigma_P = 0.25d0 !0.278d0
        
        integer    , parameter :: obc_type_N = always_outflow
        integer    , parameter :: obc_type_S = always_outflow
        integer    , parameter :: obc_type_E = always_outflow
        integer    , parameter :: obc_type_W = always_inflow

      end module parameters_input
