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

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        !< debug option allowing extra checks in the code
        logical    , parameter :: debug = .true.        

        !<computational field dimensions
        real(rkind), parameter :: x_min = 0.0000000000d0
        real(rkind), parameter :: x_max = 0.6489000000d0
        real(rkind), parameter :: y_min = 0.0000000000d0
        real(rkind), parameter :: y_max = 0.5000000000d0
        
        !<computational times
        real(rkind), parameter :: t_max = 100.0000000000d0 !10.0d0
        real(rkind), parameter :: dt = 0.0008000000d0
        
        !<output writing
        real(rkind), parameter :: detail_print = 0.0063500000d0
        logical    , parameter :: write_domain_extension = .true.
        logical    , parameter :: write_detectors = .true.

        !<mpi choice
        integer, parameter :: npx = 1 !<number of processors along x
        integer, parameter :: npy = 1 !<number of processors along y

        !<size of the main tables
        !<careful, choose ne according to the physical model
        integer(ikind), parameter :: ntx = 57
        integer(ikind), parameter :: nty = 45

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
        !
        !--------------------------------------------
        !dim2d_lowTemperature:
        !--------------------------------------------
        !option to force the use of low temperature laws
        !for the saturated liquid and vapor mass densities
        !and the interface length
        !--------------------------------------------
        integer    , parameter :: flow_direction = y_direction
        real(rkind), parameter :: flow_x_side = 1.0000000000d0
        real(rkind), parameter :: flow_y_side = 1.0000000000d0
        real(rkind), parameter :: flow_velocity = 0.0000000000d0
        
        real(rkind), parameter :: T0 = 0.9990000000d0

        integer    , parameter :: ic_choice = bubble_next_to_wall

        integer    , parameter :: phase_at_center = liquid
        real(rkind), parameter :: wall_micro_contact_angle = 45.0d0

        logical    , parameter :: ic_perturbation_ac = .false.
        real(rkind), parameter :: ic_perturbation_amp = 0.0000000000d0

        logical    , parameter :: li_perturbation_ac = .false.
        real(rkind), parameter :: li_perturbation_amp = 0.0000000000d0

        logical    , parameter :: dim2d_lowTemperature = .false.

        !<body forces choice
        integer    , parameter :: gravity_ac  = .true.
        real(rkind), parameter :: gravity_amp = 0.003d0
        integer    , parameter :: wave_forcing = no_wave_forcing

        !<boundary conditions choice
        integer, parameter :: bc_choice = wall_x_simplified_choice

        !<output choice
        integer, parameter :: io_choice = netcdf_choice
        logical, parameter :: io_onefile_per_proc = .true.

        !< boundary conditions parameters
        !-----------------------------------------------------
        !type of boundary conditions applied at the edge
        !(constrained by the bc_choice parameter)
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
        integer    , parameter :: bc_N_type_choice = bc_flux_and_node_choice
        integer    , parameter :: bc_S_type_choice = bc_flux_and_node_choice
        integer    , parameter :: bc_E_type_choice = bc_flux_and_node_choice
        integer    , parameter :: bc_W_type_choice = bc_flux_and_node_choice


        !-----------------------------------------------------
        !for the open boundary conditions
        !-----------------------------------------------------
        !obc_eigenqties_strategy : control how the eigenquantities
        !                          are computed at the edge for the
        !                          open boundary conditions:
        ! 
        !                          1) obc_eigenqties_bc
        !                          2) obc_eigenqties_lin
        !
        !obc_edge_xy_strategy : control which strategy is used
        !                       when computing the gridpoints
        !                       at the anti-corner boundary
        !                       section
        !
        !                       1) obc_edge_xy_corner
        !                       2) obc_edge_xy_flux
        !                       3) obc_edge_xy_diag_flux
        !
        !
        !obc_edge_flux_strategy: control whether the capillarity
        !                        terms are included when computing
        !                        the one-side fluxes used for the
        !                        open boundary conditions
        !
        !                        1) obc_edge_flux_capillarity
        !                        2) obc_edge_flux_no_capillarity
        !
        !------------------------------------------------------------
        !openbc_perturbation_T0_ac : integer
        !------------------------------------------------------------
        ! 1: activate the perturbation of the temperature used to 
        !    determine the far field values
        ! 0: do not activate the perturbation of the temperature used
        !    to determine the far field values
        !
        !------------------------------------------------------------
        !openbc_perturbation_vx0_ac : integer
        !------------------------------------------------------------
        ! 1: activate the perturbation of the x-component of the
        !    velocity used to determine the far field values
        ! 0: do not activate the perturbation of the x-component of
        !    the velocity used to determine the far field values
        !
        !------------------------------------------------------------
        !openbc_perturbation_vy0_ac : integer
        !------------------------------------------------------------
        ! 1: activate the perturbation of the y-component of the
        !    velocity used to determine the far field values
        ! 0: do not activate the perturbation of the y-component of
        !    the velocity used to determine the far field values
        !
        !------------------------------------------------------------
        !openbc_perturbation_T0_amp : real
        !------------------------------------------------------------
        ! amplitude of the perturbation applied to the temperature
        ! used to compute the far field values
        !
        !------------------------------------------------------------
        !openbc_perturbation_vx0_amp : real
        !------------------------------------------------------------
        ! amplitude of the perturbation applied to the x-component of
        ! the velocity used to compute the far field values
        !
        !------------------------------------------------------------
        !openbc_perturbation_vy0_amp : real
        !------------------------------------------------------------
        ! amplitude of the perturbation applied to the y-component of
        ! the velocity used to compute the far field values
        !------------------------------------------------------------
        integer    , parameter :: obc_eigenqties_strategy = obc_eigenqties_lin
        integer    , parameter :: obc_edge_xy_strategy    = obc_edge_xy_flux
        integer    , parameter :: obc_edge_flux_strategy  = obc_edge_flux_capillarity
        logical    , parameter :: obc_edge_overlap_ac     = .true.
        logical    , parameter :: obc_crenel_removal_ac   = .true. !no_edge_limit (pb at interfaces between bf_layers)
        integer    , parameter :: obc_dct_distance = 5

        logical    , parameter :: obc_perturbation_T0_ac = .false.
        logical    , parameter :: obc_perturbation_vx0_ac = .false.
        logical    , parameter :: obc_perturbation_vy0_ac = .false.

        real(rkind), parameter :: obc_perturbation_T0_amp = 0.0000000000d0
        real(rkind), parameter :: obc_perturbation_vx0_amp = 0.0000000000d0
        real(rkind), parameter :: obc_perturbation_vy0_amp = 0.0000000000d0


        !------------------------------------------------------------
        !for the Yoo and Lodato open boundary conditions
        !------------------------------------------------------------
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
        !------------------------------------------------------------
        real(rkind), parameter :: sigma_P = 0.25d0 !0.278d0
        
        integer    , parameter :: obc_type_N = always_outflow
        integer    , parameter :: obc_type_S = always_outflow
        integer    , parameter :: obc_type_E = always_outflow
        integer    , parameter :: obc_type_W = always_inflow


        !< domain adaptation parameters
        !-----------------------------------------------------
        !determination of the directiosn in which the
        !computational domain can be extended
        !-----------------------------------------------------
        !
        !adapt_N_choice : choose whether the North boundary 
        !                 can be extended
        !
        !adapt_S_choice : choose whether the South boundary 
        !                 can be extended
        !
        !adapt_E_choice : choose whether the East boundary 
        !                 can be extended
        !
        !adapt_W_choice : choose whether the West boundary 
        !                 can be extended
        !-----------------------------------------------------
        integer, parameter :: adapt_N_choice = adapt_domain_choice
        integer, parameter :: adapt_S_choice = adapt_domain_choice
        integer, parameter :: adapt_E_choice = adapt_domain_choice
        integer, parameter :: adapt_W_choice = adapt_domain_choice


        !------------------------------------------------------------
        ! criterion to decide whether nodes are activated
        !------------------------------------------------------------
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
        !------------------------------------------------------------
        logical    , parameter :: bf_openbc_md_threshold_ac = .false.
        real(rkind), parameter :: bf_openbc_md_threshold = 0.0000000000d0


        !------------------------------------------------------------
        !debugging options
        !------------------------------------------------------------
        !debug_restart_for_geometry :
        !    control whether the restart option is only used to get
        !    the geometry of the previous computational domain
        !
        !debug_adapt_computational_domain :
        !    control whether the edges of the computational domain
        !    are adapted once the simulation starts
        !    should be set to .true. by default
        !
        !debug_geometry_update :
        !    control whether the new grid points are computed when
        !    increasing the computational domain (only use for tests,
        !    should be set to .false. by default)
        !
        !debug_initialize_nodes :
        !    the nodes are initialized with debug_real
        !
        !debug_initialize_timedev :
        !    the time derivatives are initialized with debug_real
        !------------------------------------------------------------
        logical    , parameter :: debug_restart_for_geometry = .false.
        logical    , parameter :: debug_adapt_computational_domain = .true.
        logical    , parameter :: debug_geometry_update = .false.

        logical    , parameter :: debug_initialize_nodes    = .true.
        logical    , parameter :: debug_initialize_bc_nodes = .false.
        logical    , parameter :: debug_initialize_timedev  = .true.
        real(rkind), parameter :: debug_real=1e30


        !------------------------------------------------------------
        !steady state simulation options
        !------------------------------------------------------------
        !steady_state_simulation :
        !    logical stating whether the simulation should be run as
        !    if it is a steady state computation (no time limit)
        !
        !steady_state_limit :
        !    parameter checked such that the simulation is considered
        !    steady state
        !------------------------------------------------------------
        logical    , parameter :: steady_state_simulation = .true.
        real(rkind), parameter :: steady_state_limit = 1.0e-12

      end module parameters_input
