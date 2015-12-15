      !> @file
      !> module containing the user input choices defined
      !> as constants and propagated at compilation time
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module containing the user input choices defined
      !> as constants and propagated at compilation time
      !
      !> @date
      !> 20_08_2013 - initial version   - J.L. Desmarais
      !-----------------------------------------------------------------
      module parameters_input

        use netcdf

        use parameters_constant

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        !============================================================
        !============================================================
        ! 1) main simulation parameters
        ! 2) initial conditions choice
        ! 3) boundary conditions
        ! 4) domain adaptation parameters
        ! 5) debugging options
        ! 6) steady state options
        !============================================================
        !============================================================


        !============================================================
        ! 1) main simulation parameters
        !============================================================
        ! size and resolution of the space domain, resolution of the
        ! time integration, number of processors and tiles,
        ! print granularity for the outputs
        !============================================================

        !>@brief debug
        !> debug option allowing extra checks in the code
        !> \li .true.: extra checks
        !> \li .false.: no checks
        !--------------------------------------------------
        logical    , parameter :: debug = .true.        

        ! computational field dimensions
        !--------------------------------------------------
        real(rkind), parameter :: x_min = -1.8500000000d0 !<@brief minimum x-coordinate for the x-map
        real(rkind), parameter :: x_max = 1.8500000000d0  !<@brief maximum x-coordinate for the x-map
        real(rkind), parameter :: y_min = -1.8500000000d0  !<@brief minimum y-coordinate for the y-map
        real(rkind), parameter :: y_max = 1.8500000000d0  !<@brief maximum y-coordinate for the y-map
        
        ! computational times
        real(rkind), parameter :: t_max = 0.0008000000d0 !<@brief maximum simulation time
        real(rkind), parameter :: dt = 0.0008000000d0     !<@brief time step
        
        ! output writing
        real(rkind), parameter :: detail_print = 1.0000000000d0   !<@brief percentage of time steps written in output: 0.0d0=no file writen, 0.5=write output every two time steps, 1.0=write all the time steps
        logical    , parameter :: write_domain_extension = .true. !<@brief write the buffer layers (domain extension)

        ! mpi choice
        integer, parameter :: npx = 1 !<@brief number of tiles along x when the domain is computed by several processors: domain=(npx x npy) tiles
        integer, parameter :: npy = 1 !<@brief number of tiles along y when the domain is computed by several processors: domain=(npx x npy) tiles

        ! size of the main tables
        ! careful, choose ne according to the physical model
        integer(ikind), parameter :: ntx = 108 !<@brief total number of grid-points along the x-direction for the computational domain
        integer(ikind), parameter :: nty = 108 !<@brief total number of grid-points along the y-direction for the computational domain

        integer(ikind), parameter :: nx = ntx/npx !<@brief number of grid-points along the x-direction for one tile
        integer(ikind), parameter :: ny = nty/npy !<@brief number of grid-points along the y-direction for one tile
        integer       , parameter :: ne = 4       !<@brief number of governing equations for the physical model chosen: for DIM, 4 governing equations \f$(\rho,q_x,q_y,\rho E)\f$
        integer       , parameter :: bc_size = 2  !<@brief number of ghost grid points for the boundary (depends on the space discretization operators chosen)

        !> @brief type of i/o operators
        !> - nectdf_choice : netcdf files as output
        !--------------------------------------------------
        integer, parameter :: io_choice = netcdf_choice

        !> @brief depending on the file system, it is possible
        !> that all processors access the same file at the same
        !> time (parallel file system), otherwise, to prevent
        !> the i/o overhead, every processor is writing its own
        !> file
        !> \li .true.: one file per processor
        !> \li .false.: one major file (parallel file system)
        !--------------------------------------------------        
        logical, parameter :: io_onefile_per_proc = .true.


        !============================================================
        ! 2) initial conditions choice
        !============================================================
        ! flow direction
        ! flow_x_side
        ! flow_y_side
        ! flow_velocity
        ! flow_profile
        ! T0
        ! ic_choice
        ! phase_at_center
        ! dim2d_lowTemperature
        ! gravity_ac
        ! gravity_amp
        ! wave_forcing
        !============================================================
        ! 2.1) perturbation
        ! 2.2) additional parameters for a bubble
        ! 2.3) additional parameters for bubble near nucleating bubble
        ! 2.4) additional parameters for two bubbles advected by flow
        !============================================================

        !>@brief
        !> direction of the main flow
        !> - x_direction  : from left to right
        !> - y_direction  : from bottom to up
        !> - xy_direction : from SW to NE corner
        !--------------------------------------------------
        integer    , parameter :: flow_direction = y_direction

        !> @brief
        !> direction of the main flow along the x-direction
        !> - \f$-1.0d0\f$ : from right to left
        !> - \f$+1.0d0\f$ : from left to right
        !--------------------------------------------------
        real(rkind), parameter :: flow_x_side = 1.0000000000d0

        !> @brief
        !> direction of the main flow along the y-direction
        !> - \f$-1.0d0\f$ : from top to bottom
        !> - \f$+1.0d0\f$ : from bottom to top
        !--------------------------------------------------
        real(rkind), parameter :: flow_y_side = 1.0000000000d0
        
        !> @brief
        !> velocity of the mean flow
        !--------------------------------------------------
        real(rkind), parameter :: flow_velocity = 0.1000000000d0

        !> @brief profile of the flow
        !> - linear_profile: linear velocity profile
        !> - parabolic profile: parabolic velocity profile
        !--------------------------------------------------
        integer    , parameter :: flow_profile = parabolic_profile

        !>@brief
        !> temperature of the mean flow
        !--------------------------------------------------
        real(rkind), parameter :: T0 = 0.9990000000d0

        !>@brief initial conditions
        !> - for wave2d equations
        !>      - peak          : peak in the center of the domain
        !>      - negative_spot : negative field in the center of the domain
        !> - for NS equations
        !>      - steady_state  : constant everywhere
        !>      - peak          : peak in the center of the domain
        !>      - vortex        : vortex in the center of the domain
        !>      - sym_x         : symmetry compared to the y-axis
        !>      - sym_y         : symmetry compared to the x-axis
        !> - for DIM equations
        !>      - steady_state       : constant everywhere
        !>      - drop_retraction    : ellipsoidal droplet
        !>      - bubble_ascending   : initial bubble
        !>      - homogeneous_liquid : constant liquid density
        !>      - phase_separation   : unstable mass density
        !--------------------------------------------------
        integer    , parameter :: ic_choice = bubbles_transported

        !>@brief
        !> phase present in the center of the bubble/droplet
        !> - vapor : create a bubble in saturated liquid 
        !> - liquid : create a droplet in saturated vapor 
        !--------------------------------------------------
        integer    , parameter :: phase_at_center = vapor

        !>@brief
        !> range over the interpolation for the mass densities
        !> of the vapor and liquid phases and the width of the
        !> interface is valid:
        !> \li .false.: \f$ T \in [0.995,0.999]\f$
        !> \li .true.: \f$ T \in [0.95,0.995]\f$
        !--------------------------------------------------
        logical    , parameter :: dim2d_lowTemperature = .false.

        !>@brief 
        !> add gravitational forces in the governing equations
        !--------------------------------------------------
        logical    , parameter :: gravity_ac = .false.

        !>@brief
        !> amplitude of the gravitional acceleration imposed
        !--------------------------------------------------
        real(rkind), parameter :: gravity_amp = 0.0000000000d0

        !>@brief 
        !> source term for the wave physical model to
        !> generate waves in the center of the computational
        !> domain
        !> - no_wave_forcing : no source term
        !> - oscillatory_forcing : sinusoid forcing
        !> - intermittent_oscillatory_forcing : sinusoid forcing then nothing
        !> - moving_oscillatory_forcing : sinusoid forcing and the center of the source is moving spatially
        !--------------------------------------------------
        integer    , parameter :: wave_forcing = no_wave_forcing

        
        !============================================================ 
        ! 2.1) perturbation of the initial conditions to estimate
        !      the stability properties
        !============================================================
        ! ic_perturbation_ac
        ! ic_perturbation_amp
        ! li_perturbation_amp
        !============================================================

        !> @brief
        !> activate the superposition of perturbations
        !> over the initial conditions
        !> \li .false. : no perturbation
        !> \li .true. : additional perturbations
        !--------------------------------------------------
        logical    , parameter :: ic_perturbation_ac = .false.

        !> @brief
        !> amplitude of the perturbations imposed over
        !> the initial conditions
        !--------------------------------------------------
        real(rkind), parameter :: ic_perturbation_amp = 0.0000000000d0

        !>@brief
        !> activate the superposition of perturbations
        !> over the width of the interface initially prescribed
        !> \li .false. : no perturbation
        !> \li .true. : perturbation activated
        !--------------------------------------------------
        logical    , parameter :: li_perturbation_ac = .false.

        !>@brief
        !> amplitude of the perturbations imposed over
        !> the width of the interface
        !> \f[ \tilde{L_i} = (1.0+\textrm{li\_perturbation\_ac})*L_i \f]
        !> where \f$ \tilde{L_i} \f$ is the perturbed width of the interface and 
        !> \f$ L_i \f$ is the theoretical width of the interface
        !--------------------------------------------------
        real(rkind), parameter :: li_perturbation_amp = 0.0000000000d0


        !============================================================
        ! 2.2) parameters if a bubble is initialized in the domain
        !============================================================
        ! ratio_bubble_interface
        !============================================================

        !> @brief ratio between the radius of the initial
        !> bubble and the width of the interface
        !-----------------------------------------------------
        real(rkind), parameter :: ratio_bubble_interface = 2.0000000000d0


        !==================================================
        ! 2.3) wall inflow bubble parameters
        !==================================================
        ! parameters related to the inflow bubble
        ! that can be introduced in the bubble nucleation i.c.
        !==================================================

        !>@brief add a bubble in the flow around the
        !> nucleation spot as an attempt to reproduce bubble jet
        !  - .true.  : add a bubble in the inflow for the bubble nucleation i.c.
        !  - .false. : do not add a bubble in the inflow for the bubble nucleation i.c.
        !--------------------------------------------------
        logical    , parameter :: inflow_bubble_ac = .false.

        !> @brief position of the inflow bubble along the x-axis
        !--------------------------------------------------
        real(rkind), parameter :: inflow_bubble_x_center = 0.0d0        

        !> @brief position of the inflow bubble along the y-axis
        !--------------------------------------------------
        real(rkind), parameter :: inflow_bubble_y_center = 0.125d0

        !>@brief radius of the inflow bubble
        !--------------------------------------------------
        real(rkind), parameter :: inflow_bubble_radius   = 0.05d0        


        !==================================================
        ! 2.4) two bubbles advected by the flow
        !==================================================
        ! li_separation
        !==================================================

        !>@brief
        !> Separation length between two bubbles for the
        !> I.C. with two bubbles advected. It is expressed
        !> as a fraction of the width of the interface at the
        !> initial temperature
        !--------------------------------------------------
        integer    , parameter :: li_separation = 3.0000000000d0


        !============================================================
        ! 3) boundary conditions
        !============================================================
        ! bc_choice
        ! bc_N_choice
        ! bc_S_choice
        ! bc_E_choice
        ! bc_W_choice
        ! bc_NW_choice
        ! bc_NE_choice
        ! bc_SE_choice
        ! bc_SW_choice        
        !============================================================
        ! 3.1) order in which the boundary conditions are computed
        ! 3.2) type of boundary conditions
        ! 3.3) additional parameters for wall boundary conditions
        ! 3.4) additional parameters for wall heat source
        ! 3.5) additional parameters for the open boundary conditions
        ! 3.6) additional parameters for Yoo and Lodato open boundary conditions
        !============================================================

        !> @brief
        !> configuration for the boundary conditions
        !> - periodic_xy_choice : periodic along x and y directions
        !> - reflection_xy_choice : reflection along x and y directions
        !> - hedstrom_xy_choice : open (hedstrom type) along x and y directions
        !> - poinsot_xy_choice : open (Poinsot and Lele type) along x and y directions
        !> - yoolodato_xy_choice : open (Yoo and Lodato type) along x and y directions
        !> - wall_xy_choice : wall along the x and y directions
        !> - wall_S_reflection_choice : wall on the bottom, reflection for the N,E,W borders
        !> - wall_S_open_choice : wall on the bottom, open for the N,E,W borders
        !> - half_wall_S_open_choice : wall on the bottom, open on the E and N, reflection for W 
        !
        !> @warning
        !> since only a limited number of configuration can
        !> be implemented in this way, a new implementation is
        !> modifying how boundary conditions are chosen. Now, the
        !> boundary conditions are specified for each boundary
        !> layer (N,S,E,W,SE,SW,NE,NW)
        !--------------------------------------------------
        integer, parameter :: bc_choice = hedstrom_xy_choice

        !> @brief boundary condition configuration for
        !> the North layer
        !--------------------------------------------------
        integer, parameter :: bc_N_choice = hedstrom_choice

        !> @brief boundary condition configuration for
        !> the South layer
        !--------------------------------------------------
        integer, parameter :: bc_S_choice = hedstrom_choice

        !> @brief boundary condition configuration for
        !> the East layer
        !--------------------------------------------------
        integer, parameter :: bc_E_choice = hedstrom_choice

        !> @brief boundary condition configuration for
        !> the West layer
        !--------------------------------------------------
        integer, parameter :: bc_W_choice = hedstrom_choice

        !> @brief boundary condition configuration for
        !> the North-West corner
        !--------------------------------------------------
        integer, parameter :: bc_NW_choice = hedstrom_choice

        !> @brief boundary condition configuration for
        !> the North-East corner
        !--------------------------------------------------
        integer, parameter :: bc_NE_choice = hedstrom_choice

        !> @brief boundary condition configuration for
        !> the South-East corner
        !--------------------------------------------------
        integer, parameter :: bc_SE_choice = hedstrom_choice

        !> @brief boundary condition configuration for
        !> the South-West corner
        !--------------------------------------------------
        integer, parameter :: bc_SW_choice = hedstrom_choice


        !============================================================
        ! 3.1) order in which the boundary conditions are computed
        !============================================================
        ! bc_order1
        ! bc_order2
        ! bc_order3
        ! bc_order4
        ! bc_order5
        ! bc_order6
        ! bc_order7
        ! bc_order8
        !============================================================

        !> @brief specify the order in which the boundary
        !> conditions are applied (to avoid conflicts):
        !> first boundary layer computed
        !--------------------------------------------------
        integer, parameter :: bc_order1 = SW_corner_type

        !> @brief specify the order in which the boundary
        !> conditions are applied (to avoid conflicts):
        !> second boundary layer computed
        !--------------------------------------------------
        integer, parameter :: bc_order2 = S_edge_type

        !> @brief specify the order in which the boundary
        !> conditions are applied (to avoid conflicts):
        !> third boundary layer computed
        !--------------------------------------------------
        integer, parameter :: bc_order3 = SE_corner_type

        !> @brief specify the order in which the boundary
        !> conditions are applied (to avoid conflicts):
        !> fourth boundary layer computed
        !--------------------------------------------------
        integer, parameter :: bc_order4 = W_edge_type

        !> @brief specify the order in which the boundary
        !> conditions are applied (to avoid conflicts):
        !> fifth boundary layer computed
        !--------------------------------------------------
        integer, parameter :: bc_order5 = E_edge_type

        !> @brief specify the order in which the boundary
        !> conditions are applied (to avoid conflicts):
        !> sixth boundary layer computed
        !--------------------------------------------------
        integer, parameter :: bc_order6 = NW_corner_type

        !> @brief specify the order in which the boundary
        !> conditions are applied (to avoid conflicts):
        !> seventh boundary layer computed
        !--------------------------------------------------
        integer, parameter :: bc_order7 = N_edge_type

        !> @brief specify the order in which the boundary
        !> conditions are applied (to avoid conflicts):
        !> eighth boundary layer computed
        !--------------------------------------------------
        integer, parameter :: bc_order8 = NE_corner_type


        !============================================================
        ! 3.2) type of boundary conditions
        !============================================================
        ! in order to align correctly the subroutines at
        ! compilation time, it is important to specify
        ! type of boundary conditions applied at the edge
        ! (constrained by the bc_choice parameter)
        !============================================================
        ! bc_N_type_choice
        ! bc_S_type_choice
        ! bc_E_type_choice
        ! bc_W_type_choice
        ! bc_NE_type_choice
        ! bc_NW_type_choice
        ! bc_SE_type_choice
        ! bc_SW_type_choice
        !============================================================

        !> @brief type of boundary condition applied at
        !> the North boundary (bc_nodes_choice, bc_fluxes_choice,
        !> bc_timedev_choice, bc_flux_and_node_choice)
        !--------------------------------------------------
        integer, parameter :: bc_N_type_choice = bc_timedev_choice

        !> @brief type of boundary condition applied at
        !> the South boundary (bc_nodes_choice, bc_fluxes_choice,
        !> bc_timedev_choice, bc_flux_and_node_choice)
        !--------------------------------------------------
        integer, parameter :: bc_S_type_choice = bc_timedev_choice

        !> @brief type of boundary condition applied at
        !> the East boundary (bc_nodes_choice, bc_fluxes_choice,
        !> bc_timedev_choice, bc_flux_and_node_choice)
        !--------------------------------------------------
        integer, parameter :: bc_E_type_choice = bc_timedev_choice
        
        !> @brief type of boundary condition applied at
        !> the West boundary (bc_nodes_choice, bc_fluxes_choice,
        !> bc_timedev_choice, bc_flux_and_node_choice)
        !--------------------------------------------------
        integer, parameter :: bc_W_type_choice = bc_timedev_choice
        
        !> @brief type of boundary condition applied at
        !> the North-West boundary (bc_nodes_choice, bc_fluxes_choice,
        !> bc_timedev_choice, bc_flux_and_node_choice)
        !--------------------------------------------------
        integer, parameter :: bc_NW_type_choice = bc_timedev_choice

        !> @brief type of boundary condition applied at
        !> the North-East boundary (bc_nodes_choice, bc_fluxes_choice,
        !> bc_timedev_choice, bc_flux_and_node_choice)
        !--------------------------------------------------
        integer, parameter :: bc_NE_type_choice = bc_timedev_choice

        !> @brief type of boundary condition applied at
        !> the South-West boundary (bc_nodes_choice, bc_fluxes_choice,
        !> bc_timedev_choice, bc_flux_and_node_choice)
        !--------------------------------------------------
        integer, parameter :: bc_SW_type_choice = bc_timedev_choice

        !> @brief type of boundary condition applied at
        !> the South-East boundary (bc_nodes_choice, bc_fluxes_choice,
        !> bc_timedev_choice, bc_flux_and_node_choice)
        !--------------------------------------------------
        integer, parameter :: bc_SE_type_choice = bc_timedev_choice


        !============================================================
        ! 3.3) additional parameters for wall boundary conditions
        !============================================================
        ! wall_surface_type
        ! wall_micro_contact_angle
        ! wall_heater_center
        ! wall_heater_length
        ! wall_heater_variation_angle_length
        ! wall_heater_micro_contact_angle
        !============================================================

        !> @brief type of wall surface
        !> - uniform_surface : the contact angle is the same everywhere
        !> - surface_with_heaters : the contact angle varies at the location of the heaters
        !-----------------------------------------------------
        integer    , parameter :: wall_surface_type = 0.0
        
        !> @brief wall micro contact angle (in degrees)
        !-----------------------------------------------------
        real(rkind), parameter :: wall_micro_contact_angle = 5.0000000000d0

        !> @brief location of the center of the heater, \f$ x_c \f$
        !--------------------------------------------------
        real(rkind), parameter :: wall_heater_center = 0.0000000000d0

        !> @brief extent of the wall heater, \f$ l_h\f$
        !> Therefore, the heater corresponds to:
        !> \f[ x \in [ x_c-l_h, x_c+l_h] \f]
        !--------------------------------------------------
        real(rkind), parameter :: wall_heater_length = 2.0000000000d0

        !> @brief characteristic length over which the contact
        !> length is varying at the wall
        !--------------------------------------------------
        real(rkind), parameter :: wall_heater_variation_angle_length = 0.5000000000d0
 
        !> @brief contact angle imposed on the heater (in degrees)
        !--------------------------------------------------
        real(rkind), parameter :: wall_heater_micro_contact_angle = 45.0000000000d0


        !============================================================
        ! 3.4) additional parameters for wall heat source
        !============================================================
        ! wall_heat_source_choice
        ! wall_maximum_heat_flux
        ! wall_heat_source_center
        ! wall_heat_source_variance
        ! wall_extra_heat_source_choice
        ! wall_maximum_extra_heat_flux
        ! wall_extra_heat_source_center
        ! wall_extra_heat_source_variance
        !============================================================

        !>@brief choice of heat source at the wall to impose
        !> the heat flux at the wall, \f$ Q_w\f$
        !> - no_heat_source       : no heat is provided \f$ Q_w=0 \f$
        !> - constant_heat_source : the heat sourcve is constant, \f$ Q_w= Q_{w,m} \f$
        !> - gaussian_heat_source : Gaussian profile for the heat source, \f$ Q_w(x) = \frac{Q_{w,m}}{\sigma_w \sqrt{2 \pi}} \exp\left( \frac{-(x-x_{w,c})^2}{2 \sigma_w^2} \right) \f$
        !--------------------------------------------------
        integer    , parameter :: wall_heat_source_choice = no_heat_source
        
        !>@brief maximum heat flux at the wall,
        !> \f$ Q_{w,m} \f$
        !--------------------------------------------------
        real(rkind), parameter :: wall_maximum_heat_flux = 0.0000000000d0        

        !>@brief center if gaussian heat source,
        !> \f$ x_{c,w}\f$
        !--------------------------------------------------
        real(rkind), parameter :: wall_heat_source_center   = wall_heater_center

        !>@brief variance if gaussian heat source,
        !> \f$ \sigma_{w,m}\f$
        !--------------------------------------------------
        real(rkind), parameter :: wall_heat_source_variance = 0.5d0*wall_heater_length

        !> @brief choice of extra heat source at the wall
        !--------------------------------------------------
        integer    , parameter :: wall_extra_heat_source_choice = no_heat_source

        !> @brief maximum extra heat flux at the wall
        !--------------------------------------------------
        real(rkind), parameter :: wall_maximum_extra_heat_flux = 0.0000000000d0

        !> @brief center if gaussian extra_heat source
        !--------------------------------------------------
        real(rkind), parameter :: wall_extra_heat_source_center   = wall_heater_center

        !> @brief variance if gaussian extra_heat source
        !-----------------------------------------------------
        real(rkind), parameter :: wall_extra_heat_source_variance = 0.5d0*wall_heater_length
        

        !============================================================
        ! 3.5) additional parameters for the open boundary conditions
        !============================================================
        ! obc_eigenqties_strategy
        ! obc_edge_xy_strategy
        ! obc_edge_flux_strategy
        ! obc_edge_overlap_ac
        ! obc_crenel_removal_ac
        ! obc_dct_distance
        ! obc_perturbation_T0_ac
        ! obc_perturbation_vx0_ac
        ! obc_perturbation_vy0_ac
        ! obc_perturbation_T0_amp
        ! obc_perturbation_vx0_amp
        ! obc_perturbation_vy0_amp
        !============================================================

        !>@brief control how the eigenquantities are computed
        !> at the edge for the open boundary conditions:
        !> - obc_eigenqties_bc: using the value at the boundary points
        !> - obc_eigenqties_lin: using the value in the far field
        !--------------------------------------------------
        integer    , parameter :: obc_eigenqties_strategy = obc_eigenqties_bc

        !> @brief control which strategy is used when computing
        !> the gridpoints at the anti-corner boundary section
        !> - obc_edge_xy_corner: as a corner
        !> - obc_edge_xy_flux: some points as a corner, some points as an edge
        !> - obc_edge_xy_diag_flux: computing the fluxes using diagonal fluxes
        !--------------------------------------------------
        integer    , parameter :: obc_edge_xy_strategy    = obc_edge_xy_flux

        !>@brief control whether the capillarity terms are included when computing
        !> the one-side fluxes used for the open boundary conditions
        !> - obc_edge_flux_capillarity: add the capillary terms
        !> - obc_edge_flux_no_capillarity: remove the capillary terms
        !--------------------------------------------------
        integer    , parameter :: obc_edge_flux_strategy  = obc_edge_flux_capillarity

        !>@brief control whether bc_sections are rearranged to
        !> combined edge and anti-corner types of bc_sections
        !--------------------------------------------------
        logical    , parameter :: obc_edge_overlap_ac     = .true.

        !>@brief control whether bc_sections are rearranged to
        !> remove crenel overlap
        !--------------------------------------------------
        logical    , parameter :: obc_crenel_removal_ac   = .true. !no_edge_limit (pb at interfaces between bf_layers)

        !> @brief control the distance between the detectors
        !> and the edge of the computational domain
        !--------------------------------------------------
        integer    , parameter :: obc_dct_distance = 5

        !>@brief activate or not the perturbation of the temperature
        !> imposed in the far field compared to the initial conditions
        !------------------------------------------------------------
        logical    , parameter :: obc_perturbation_T0_ac = .false.

        !> @brief activate or not the perturbation of the velocity along
        !> the x-direction imposed in the far field compared to the
        !> initial conditions
        !------------------------------------------------------------
        logical    , parameter :: obc_perturbation_vx0_ac = .false.

        !> @brief activate or not the perturbation of the velocity along
        !> the y-direction imposed in the far field compared to the
        !> initial conditions
        !------------------------------------------------------------
        logical    , parameter :: obc_perturbation_vy0_ac = .false.

        !> @brief amplitude of the perturbation for the temperature
        !> imposed in the far field
        !> \f[ \tilde{T}_f = (1.0 + \epsilon_T)*T_f \f]
        !> where \f$T_f\f$ is the temperature for the initial
        !> conditions, and \f$\epsilon_T\f$ is the amplitude of
        !> the perturbations imposed
        !------------------------------------------------------------
        real(rkind), parameter :: obc_perturbation_T0_amp = 0.0000000000d0

        !> @brief amplitude of the perturbation for the velocity
        !> in the x-direction imposed in the far field
        !> \f[ \tilde{v_x}_f = (1.0 + \epsilon_vx)*{v_x}_f \f]
        !> where \f${v_x}_f\f$ is the veloity in the x-direction
        !> for the initial conditions, and \f$\epsilon_{vx}\f$ is
        !> the amplitude of the perturbations imposed
        !------------------------------------------------------------
        real(rkind), parameter :: obc_perturbation_vx0_amp = 0.0000000000d0

        !> @brief amplitude of the perturbation for the velocity
        !> in the y-direction imposed in the far field
        !> \f[ \tilde{v_y}_f = (1.0 + \epsilon_vy)*{v_y}_f \f]
        !> where \f${v_y}_f\f$ is the veloity in the y-direction
        !> for the initial conditions, and \f$\epsilon_{vy}\f$ is
        !> the amplitude of the perturbations imposed
        !------------------------------------------------------------
        real(rkind), parameter :: obc_perturbation_vy0_amp = 0.0000000000d0

        
        !============================================================
        ! 3.6) additional parameters for Yoo and Lodato open boundary conditions
        !============================================================
        ! sigma_P
        ! obc_type_N
        ! obc_type_S
        ! obc_type_E
        ! obc_type_W
        !============================================================

        !>@brief relaxation coefficient used when applying the
        !> non-reflecting outflow pressure b.c.
        !> \f[ \sigma_P (P - P_{\infty}) \f]
        !--------------------------------------------------
        real(rkind), parameter :: sigma_P = 0.25d0 !0.278d0

        !>@brief type of boundary condition applied at the
        !> North boundary
        !> - always_outflow: imposed outflow
        !> - always_inflow: imposed inflow
        !> - ask_flow: check whether the boundary is of outflow of inflow type using the velocity at the boundary
        !--------------------------------------------------
        integer    , parameter :: obc_type_N = always_outflow

        !>@brief type of boundary condition applied at the
        !> South boundary
        !> - always_outflow: imposed outflow
        !> - always_inflow: imposed inflow
        !> - ask_flow: check whether the boundary is of outflow of inflow type using the velocity at the boundary
        !--------------------------------------------------
        integer    , parameter :: obc_type_S = always_outflow

        !>@brief type of boundary condition applied at the
        !> East boundary
        !> - always_outflow: imposed outflow
        !> - always_inflow: imposed inflow
        !> - ask_flow: check whether the boundary is of outflow of inflow type using the velocity at the boundary
        !--------------------------------------------------
        integer    , parameter :: obc_type_E = always_outflow

        !>@brief type of boundary condition applied at the
        !> West boundary
        !> - always_outflow: imposed outflow
        !> - always_inflow: imposed inflow
        !> - ask_flow: check whether the boundary is of outflow of inflow type using the velocity at the boundary
        !--------------------------------------------------
        integer    , parameter :: obc_type_W = always_inflow


        !============================================================
        ! 4) domain adaptation parameters
        !============================================================
        ! 4.1) directions in which the computational domain is extended
        ! 4.2) criterion to decide whether nodes are activated
        !============================================================

        !============================================================
        ! 4.1) directions in which the computational domain is extended
        !============================================================
        ! adapt_N_choice
        ! adapt_S_choice
        ! adapt_E_choice
        ! adapt_W_choice
        !============================================================

        !> @brief choose whether the North boundary can be
        !> extended
        !> - fixed_domain_choice: no boundary adaptation
        !> - adapt_domain_choice: adapt the boundary to the needs
        !--------------------------------------------------
        integer, parameter :: adapt_N_choice = adapt_domain_choice
        
        !> @brief choose whether the South boundary can be
        !> extended
        !> - fixed_domain_choice: no boundary adaptation
        !> - adapt_domain_choice: adapt the boundary to the needs
        !--------------------------------------------------
        integer, parameter :: adapt_S_choice = adapt_domain_choice

        !> @brief choose whether the East boundary can be
        !> extended
        !> - fixed_domain_choice: no boundary adaptation
        !> - adapt_domain_choice: adapt the boundary to the needs
        !--------------------------------------------------
        integer, parameter :: adapt_E_choice = adapt_domain_choice

        !> @brief choose whether the West boundary can be
        !> extended
        !> - fixed_domain_choice: no boundary adaptation
        !> - adapt_domain_choice: adapt the boundary to the needs
        !--------------------------------------------------
        integer, parameter :: adapt_W_choice = adapt_domain_choice


        !============================================================
        ! 4.2) criterion to decide whether nodes are activated
        ! for the extension of the computational domain
        !============================================================
        ! bf_openbc_md_threshold_ac
        ! bf_openbc_md_threshold
        !============================================================

        !>@brief control whether the increase of the computational
        !> domain is also activated by the value of the mass density
        !> at the edge and not just whether the speed of sound is
        !> negative
        !--------------------------------------------------
        logical    , parameter :: bf_openbc_md_threshold_ac = .false.

        !>@brief the increase of the computational domain is
        !> triggered if:
        !> \f[ \rho \in [\rho_{\textrm{vap}}+\epsilon_{\textrm{vap}},
        !>               \rho_{\textrm{liq}}-\epsilon_{\textrm{liq}}]
        !> \f]
        !> where
        !> \f[ \rho_{\textrm{mid}} = \frac{\rho_{\textrm{liq}}+\rho_{\textrm{vap}}}{2} \f]
        !> \f[ \epsilon_{\textrm{vap}} = \epsilon*(\rho_{\textrm{mid}}-\rho_\textrm{vap}) \f]
        !> \f[ \epsilon_{\textrm{liq}} = \epsilon*(\rho_{\textrm{liq}}-\rho_\textrm{mid}) \f]
        !> and \f$ \epsilon \f$ is bf_openbc_md_threshold
        !------------------------------------------------------------
        real(rkind), parameter :: bf_openbc_md_threshold = 0.0000000000d0


        !============================================================
        ! 5) debugging options
        !============================================================
        ! debug_restart_for_geometry
        ! debug_adapt_computational_domain
        ! debug_geometry_update
        ! debug_initialize_nodes
        ! debug_initialize_bc_nodes
        ! debug_initialize_timedev
        ! debug_real
        !============================================================
        
        !> @brief control whether the restart option is
        !> only used to get the geometry of the previous
        !> computational domain
        !--------------------------------------------------
        logical    , parameter :: debug_restart_for_geometry = .false.

        !> @brief control whether the edges of the computational domain
        !> are adapted once the simulation starts
        !> should be set to .true. by default
        !--------------------------------------------------
        logical    , parameter :: debug_adapt_computational_domain = .false.
        
        !> @brief control whether the new grid points are computed when
        !> increasing the computational domain (only use for tests,
        !> should be set to .false. by default)
        logical    , parameter :: debug_geometry_update = .false.

        !> @brief the nodes are initialized with debug_real
        !--------------------------------------------------
        logical    , parameter :: debug_initialize_nodes    = .true.

        !> @brief the boundary nodes are initialized with debug_real
        !------------------------------------------------------------
        logical    , parameter :: debug_initialize_bc_nodes = .true.

        !> @brief the time derivatives are initialized with debug_real
        !------------------------------------------------------------
        logical    , parameter :: debug_initialize_timedev  = .true.

        !> @brief default value when initializing nodes
        !---------------------------------------------------
        real(rkind), parameter :: debug_real=NF90_FILL_DOUBLE !1e30


        !============================================================
        ! 5) steady state options
        !============================================================
        ! steady state simulation options
        !============================================================

        !>@brief state whether the simulation should be run as
        !> if it is a steady state computation (no time limit)
        !--------------------------------------------------
        logical    , parameter :: steady_state_simulation = .false.

        !>@brief parameter checked such that the simulation
        !> is considered steady state
        !> \f[ \textrm{max}_{(x,y) \in D} \left(\left| \frac{\partial v}{\partial t} \right|\right) < \epsilon \f]
        !> where
        !> \f$ D \f$ is the computational field, \f$v\f$ the
        !> conservative variable, and \f$ \epsilon \f$ the steady_state_limit
        !------------------------------------------------------------
        real(rkind), parameter :: steady_state_limit = 1.0e-12

      end module parameters_input
