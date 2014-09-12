      !> @file
      !> class encapsulating subroutines to apply open boundary
      !> conditions at the edge of the computational domain using
      !> Yoo and Lodato boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to apply open boundary
      !> conditions at the edge of the computational domain using
      !> Yoo and Lodato boundary conditions
      !
      !> @date
      !> 10_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_class

        use bc_operators_default_class, only :
     $       bc_operators_default

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

        use lodi_corner_inflow_inflow_class, only :
     $       lodi_corner_inflow_inflow

        use lodi_corner_inflow_outflow_class, only :
     $       lodi_corner_inflow_outflow

        use lodi_corner_outflow_outflow_class, only :
     $       lodi_corner_outflow_outflow

        use lodi_edge_inflow_class, only :
     $       lodi_edge_inflow

        use lodi_edge_outflow_class, only :
     $       lodi_edge_outflow

        use lodi_timedev_xy_module, only :
     $       compute_timedev_x_edge,
     $       compute_timedev_y_edge,
     $       compute_timedev_corner

        use lodi_transverse_module, only :
     $       compute_edge_fluxes,
     $       compute_lodi_terms,
     $       compute_dev_from_flux_x,
     $       compute_dev_from_flux_y

        use ns2d_fluxes_module, only :
     $       flux_x_mass_density,
     $       flux_y_mass_density,
     $       
     $       flux_x_inviscid_momentum_x,
     $       flux_x_inviscid_momentum_y,
     $       flux_x_inviscid_total_energy,
     $       flux_y_inviscid_momentum_x,
     $       flux_y_inviscid_momentum_y,
     $       flux_y_inviscid_total_energy,
     $       
     $       flux_x_viscid_momentum_x,
     $       flux_x_viscid_momentum_y,
     $       flux_x_viscid_total_energy,
     $       flux_y_viscid_momentum_x,
     $       flux_y_viscid_momentum_y,
     $       flux_y_viscid_total_energy

        use ns2d_parameters, only :
     $       epsilon

        use ns2d_prim_module, only :
     $       compute_cons_lodi_matrix_x,
     $       compute_cons_lodi_matrix_y

        use pmodel_eq_class, only :
     $       pmodel_eq

        use parameters_constant, only :
     $       bc_timedev_choice,
     $       left,
     $       right,
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       obc_type_N,
     $       obc_type_S,
     $       obc_type_E,
     $       obc_type_W

        use parameters_kind, only :
     $       rkind,ikind

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_L1,
     $       gradient_x_x_oneside_R1,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_L1,
     $       gradient_y_y_oneside_R1,
     $       gradient_y_y_oneside_R0

        use sd_operators_x_oneside_L0_class, only :
     $       sd_operators_x_oneside_L0

        use sd_operators_x_oneside_L1_class, only :
     $       sd_operators_x_oneside_L1

        use sd_operators_x_oneside_R1_class, only :
     $       sd_operators_x_oneside_R1

        use sd_operators_x_oneside_R0_class, only :
     $       sd_operators_x_oneside_R0

        use sd_operators_y_oneside_L0_class, only :
     $       sd_operators_y_oneside_L0

        use sd_operators_y_oneside_L1_class, only :
     $       sd_operators_y_oneside_L1

        use sd_operators_y_oneside_R1_class, only :
     $       sd_operators_y_oneside_R1

        use sd_operators_y_oneside_R0_class, only :
     $       sd_operators_y_oneside_R0

        
        implicit none


        private
        public ::
     $       bc_operators,
     $       compute_fluxes_and_lodi_at_the_edges_2ndorder


        !> @class bc_operators
        !> class encapsulating subroutines to apply
        !> open boundary conditions in the x and
        !> y directions at the edge of the computational
        !> domain
        !
        !> @param ini
        !> initialize the bcx_type and bcy_type
        !> attributes of the boundary conditions
        !
        !> @param apply_bc_on_timedev
        !> apply the open boundary conditions for the time derivatives
        !---------------------------------------------------------------
        type, extends(bc_operators_default) :: bc_operators

          integer :: oneside_flow_N
          integer :: oneside_flow_S
          integer :: oneside_flow_E
          integer :: oneside_flow_W

          contains

          procedure, pass :: ini
          procedure, pass :: set_obc_type
          procedure, pass :: apply_bc_on_timedev => apply_bc_on_timedev_2ndorder

        end type bc_operators        
      

        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the main attributes
        !> of the boundary conditions
        !
        !> @date
        !> 10_09_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary conditions initialized
        !
        !>@param p_model
        !> physical model to know the type of the main variables
        !--------------------------------------------------------------
        subroutine ini(this,p_model)
        
          implicit none

          class(bc_operators), intent(inout) :: this
          type(pmodel_eq)    , intent(in)    :: p_model
          
          integer :: neq

          neq = p_model%get_eq_nb()

          this%bcx_type = bc_timedev_choice
          this%bcy_type = bc_timedev_choice

          this%oneside_flow_N = obc_type_N
          this%oneside_flow_S = obc_type_S
          this%oneside_flow_E = obc_type_E
          this%oneside_flow_W = obc_type_W

        end subroutine ini

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the main attributes
        !> of the boundary conditions
        !
        !> @date
        !> 10_09_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary conditions initialized
        !
        !>@param obc_type
        !> integer table containing the 
        !--------------------------------------------------------------
        subroutine set_obc_type(this,obc_type)

          implicit none

          class(bc_operators)  , intent(inout) :: this
          integer, dimension(4), intent(in)    :: obc_type

          this%oneside_flow_N = obc_type(N)
          this%oneside_flow_S = obc_type(S)
          this%oneside_flow_E = obc_type(E)
          this%oneside_flow_W = obc_type(W)

        end subroutine set_obc_type


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the time derivatives with space operators
        !> that are 2nd order accurate in the interior
        !> and 1st order accurate at the boundary
        !
        !> @date
        !> 10_09_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array of grid point data
        !
        !>@param dx
        !> space step along the x-axis
        !
        !>@param dy
        !> space step along the y-axis
        !
        !>@param p_model
        !> governing equations of the physical model
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param timedev
        !> time derivatives
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev_2ndorder(
     $     this,
     $     p_model,
     $     t,nodes,x_map,y_map,
     $     flux_x,flux_y,
     $     timedev)
        
          implicit none
        
          class(bc_operators)               , intent(in)    :: this
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev


          real(rkind) :: dx
          real(rkind) :: dy


          real(rkind), dimension(nx-2*bc_size,bc_size,ne) :: transverse_lodi_N
          real(rkind), dimension(nx-2*bc_size,bc_size,ne) :: transverse_lodi_S
          real(rkind), dimension(bc_size,ny-2*bc_size,ne) :: transverse_lodi_E
          real(rkind), dimension(bc_size,ny-2*bc_size,ne) :: transverse_lodi_W
          real(rkind), dimension(nx-2*bc_size,bc_size,ne) :: viscous_lodi_N
          real(rkind), dimension(nx-2*bc_size,bc_size,ne) :: viscous_lodi_S
          real(rkind), dimension(bc_size,ny-2*bc_size,ne) :: viscous_lodi_E
          real(rkind), dimension(bc_size,ny-2*bc_size,ne) :: viscous_lodi_W


          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)


          !compute the fluxes at the edge of the computational
          !domain. Distinction is made between the viscid and
          !the inviscid fluxes to be able to compute the lodi
          !transverse and viscous lodi vectors
          call compute_fluxes_and_lodi_at_the_edges_2ndorder(
     $         nodes,dx,dy,
     $         transverse_lodi_N, transverse_lodi_S,
     $         transverse_lodi_E, transverse_lodi_W,
     $         viscous_lodi_N, viscous_lodi_S,
     $         viscous_lodi_E, viscous_lodi_W,
     $         flux_x, flux_y)

          
          !compute the time derivatives at the boundary
          call compute_timedev_at_the_edges_2ndorder(
     $         p_model,
     $         t,nodes,x_map,y_map,
     $         transverse_lodi_N, transverse_lodi_S,
     $         transverse_lodi_E, transverse_lodi_W,
     $         viscous_lodi_N, viscous_lodi_S,
     $         viscous_lodi_E, viscous_lodi_W,
     $         flux_x, flux_y,
     $         timedev,
     $         this%oneside_flow_N,
     $         this%oneside_flow_S,
     $         this%oneside_flow_E,
     $         this%oneside_flow_W)

        end subroutine apply_bc_on_timedev_2ndorder


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the LODI transverse and
        !> viscous vectors as well as the fluxes at the
        !> edges of the computational domain
        !
        !> @date
        !> 10_09_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array of grid point data
        !
        !>@param dx
        !> space step along the x-axis
        !
        !>@param dy
        !> space step along the y-axis
        !
        !>@param transverse_lodi_W
        !> LODI transverse vector for the western edge
        !
        !>@param transverse_lodi_E
        !> LODI transverse vector for the eastern edge
        !
        !>@param transverse_lodi_S
        !> LODI transverse vector for the southern edge
        !
        !>@param transverse_lodi_N
        !> LODI transverse vector for the northern edge
        !
        !>@param viscous_lodi_W
        !> LODI viscous vector for the western edge
        !
        !>@param viscous_lodi_E
        !> LODI viscous vector for the eastern edge
        !
        !>@param viscous_lodi_S
        !> LODI viscous vector for the southern edge
        !
        !>@param viscous_lodi_N
        !> LODI viscous vector for the northern edge
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !-------------------------------------------------------------
        subroutine compute_fluxes_and_lodi_at_the_edges_2ndorder(
     $     nodes,dx,dy,
     $     transverse_lodi_N, transverse_lodi_S,
     $     transverse_lodi_E, transverse_lodi_W,
     $     viscous_lodi_N, viscous_lodi_S,
     $     viscous_lodi_E, viscous_lodi_W,
     $     flux_x, flux_y)

          implicit none

          real(rkind), dimension(nx,ny,ne)               , intent(in)    :: nodes
          real(rkind)                                    , intent(in)    :: dx
          real(rkind)                                    , intent(in)    :: dy
          real(rkind), dimension(nx-2*bc_size,bc_size,ne), intent(out)   :: transverse_lodi_N
          real(rkind), dimension(nx-2*bc_size,bc_size,ne), intent(out)   :: transverse_lodi_S
          real(rkind), dimension(bc_size,ny-2*bc_size,ne), intent(out)   :: transverse_lodi_E
          real(rkind), dimension(bc_size,ny-2*bc_size,ne), intent(out)   :: transverse_lodi_W
          real(rkind), dimension(nx-2*bc_size,bc_size,ne), intent(out)   :: viscous_lodi_N
          real(rkind), dimension(nx-2*bc_size,bc_size,ne), intent(out)   :: viscous_lodi_S
          real(rkind), dimension(bc_size,ny-2*bc_size,ne), intent(out)   :: viscous_lodi_E
          real(rkind), dimension(bc_size,ny-2*bc_size,ne), intent(out)   :: viscous_lodi_W
          real(rkind), dimension(nx+1,ny,ne)             , intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne)             , intent(inout) :: flux_y


          type(sd_operators_x_oneside_L0) :: s_x_L0
          type(sd_operators_x_oneside_L1) :: s_x_L1
          type(sd_operators_x_oneside_R1) :: s_x_R1
          type(sd_operators_x_oneside_R0) :: s_x_R0
          type(sd_operators_y_oneside_L0) :: s_y_L0
          type(sd_operators_y_oneside_L1) :: s_y_L1
          type(sd_operators_y_oneside_R1) :: s_y_R1
          type(sd_operators_y_oneside_R0) :: s_y_R0


          real(rkind), dimension(nx-2*bc_size+1,bc_size,ne) :: edge_inviscid_flux_N
          real(rkind), dimension(nx-2*bc_size+1,bc_size,ne) :: edge_inviscid_flux_S
          real(rkind), dimension(bc_size,ny-2*bc_size+1,ne) :: edge_inviscid_flux_E
          real(rkind), dimension(bc_size,ny-2*bc_size+1,ne) :: edge_inviscid_flux_W

          real(rkind), dimension(nx-2*bc_size+1,bc_size,ne) :: edge_viscid_flux_N
          real(rkind), dimension(nx-2*bc_size+1,bc_size,ne) :: edge_viscid_flux_S
          real(rkind), dimension(bc_size,ny-2*bc_size+1,ne) :: edge_viscid_flux_E
          real(rkind), dimension(bc_size,ny-2*bc_size+1,ne) :: edge_viscid_flux_W

          
          integer(ikind) :: j


          !fluxes computation
          !-------------------------------------------------------------

          !compute the S layers
          call compute_edge_fluxes(
     $         nodes,
     $         s_y_L0,
     $         dx, dy,
     $         1, nx-2*bc_size+1, bc_size,
     $         1,              1, 0,
     $         epsilon,
     $         flux_x_mass_density,
     $         flux_x_inviscid_momentum_x,
     $         flux_x_inviscid_momentum_y,
     $         flux_x_inviscid_total_energy,
     $         flux_x_viscid_momentum_x,
     $         flux_x_viscid_momentum_y,
     $         flux_x_viscid_total_energy,
     $         edge_inviscid_flux_S,
     $         edge_viscid_flux_S,
     $         flux_x)

          call compute_edge_fluxes(
     $         nodes,
     $         s_y_L1,
     $         dx, dy,
     $         1, nx-2*bc_size+1, bc_size,
     $         2,              2, 0,
     $         epsilon,
     $         flux_x_mass_density,
     $         flux_x_inviscid_momentum_x,
     $         flux_x_inviscid_momentum_y,
     $         flux_x_inviscid_total_energy,
     $         flux_x_viscid_momentum_x,
     $         flux_x_viscid_momentum_y,
     $         flux_x_viscid_total_energy,
     $         edge_inviscid_flux_S,
     $         edge_viscid_flux_S,
     $         flux_x)

          !compute the E and W layers
          do j=1, ny-2*bc_size+1

             !W layer
             call compute_edge_fluxes(
     $            nodes,
     $            s_x_L0,
     $            dx, dy,
     $            1, 1, 0,
     $            j, j, bc_size,
     $            epsilon,
     $            flux_y_mass_density,
     $            flux_y_inviscid_momentum_x,
     $            flux_y_inviscid_momentum_y,
     $            flux_y_inviscid_total_energy,
     $            flux_y_viscid_momentum_x,
     $            flux_y_viscid_momentum_y,
     $            flux_y_viscid_total_energy,
     $            edge_inviscid_flux_W,
     $            edge_viscid_flux_W,
     $            flux_y)

             call compute_edge_fluxes(
     $            nodes,
     $            s_x_L1,
     $            dx, dy,
     $            2, 2, 0,
     $            j, j, bc_size,
     $            epsilon,
     $            flux_y_mass_density,
     $            flux_y_inviscid_momentum_x,
     $            flux_y_inviscid_momentum_y,
     $            flux_y_inviscid_total_energy,
     $            flux_y_viscid_momentum_x,
     $            flux_y_viscid_momentum_y,
     $            flux_y_viscid_total_energy,
     $            edge_inviscid_flux_W,
     $            edge_viscid_flux_W,
     $            flux_y)

             !E layer
             call compute_edge_fluxes(
     $            nodes,
     $            s_x_R1,
     $            dx, dy,
     $            1, 1, nx-bc_size,
     $            j, j, bc_size,
     $            epsilon,
     $            flux_y_mass_density,
     $            flux_y_inviscid_momentum_x,
     $            flux_y_inviscid_momentum_y,
     $            flux_y_inviscid_total_energy,
     $            flux_y_viscid_momentum_x,
     $            flux_y_viscid_momentum_y,
     $            flux_y_viscid_total_energy,
     $            edge_inviscid_flux_E,
     $            edge_viscid_flux_E,
     $            flux_y)

             call compute_edge_fluxes(
     $            nodes,
     $            s_x_R0,
     $            dx, dy,
     $            2, 2, nx-bc_size,
     $            j, j, bc_size,
     $            epsilon,
     $            flux_y_mass_density,
     $            flux_y_inviscid_momentum_x,
     $            flux_y_inviscid_momentum_y,
     $            flux_y_inviscid_total_energy,
     $            flux_y_viscid_momentum_x,
     $            flux_y_viscid_momentum_y,
     $            flux_y_viscid_total_energy,
     $            edge_inviscid_flux_E,
     $            edge_viscid_flux_E,
     $            flux_y)

          end do

          !compute the N layers
          call compute_edge_fluxes(
     $         nodes,
     $         s_y_R1,
     $         dx, dy,
     $         1, nx-2*bc_size+1, bc_size,
     $         1,              1, ny-bc_size,
     $         epsilon,
     $         flux_x_mass_density,
     $         flux_x_inviscid_momentum_x,
     $         flux_x_inviscid_momentum_y,
     $         flux_x_inviscid_total_energy,
     $         flux_x_viscid_momentum_x,
     $         flux_x_viscid_momentum_y,
     $         flux_x_viscid_total_energy,
     $         edge_inviscid_flux_N,
     $         edge_viscid_flux_N,
     $         flux_x)

          call compute_edge_fluxes(
     $         nodes,
     $         s_y_R0,
     $         dx, dy,
     $         1, nx-2*bc_size+1, bc_size,
     $         2,              2, ny-bc_size,
     $         epsilon,
     $         flux_x_mass_density,
     $         flux_x_inviscid_momentum_x,
     $         flux_x_inviscid_momentum_y,
     $         flux_x_inviscid_total_energy,
     $         flux_x_viscid_momentum_x,
     $         flux_x_viscid_momentum_y,
     $         flux_x_viscid_total_energy,
     $         edge_inviscid_flux_N,
     $         edge_viscid_flux_N,
     $         flux_x)


          !LODI computation
          !-------------------------------------------------------------
          !S layer
          call compute_lodi_terms(
     $         nodes,
     $         dx,
     $         bc_size, 0,
     $         epsilon,
     $         edge_inviscid_flux_S,
     $         edge_viscid_flux_S,
     $         compute_cons_lodi_matrix_y,
     $         compute_dev_from_flux_x,
     $         transverse_lodi_S,
     $         viscous_lodi_S)

          !W layer
          call compute_lodi_terms(
     $         nodes,
     $         dy,
     $         0, bc_size,
     $         epsilon,
     $         edge_inviscid_flux_W,
     $         edge_viscid_flux_W,
     $         compute_cons_lodi_matrix_x,
     $         compute_dev_from_flux_y,
     $         transverse_lodi_W,
     $         viscous_lodi_W)

          !E layer
          call compute_lodi_terms(
     $         nodes,
     $         dy,
     $         nx-bc_size, bc_size,
     $         epsilon,
     $         edge_inviscid_flux_E,
     $         edge_viscid_flux_E,
     $         compute_cons_lodi_matrix_x,
     $         compute_dev_from_flux_y,
     $         transverse_lodi_E,
     $         viscous_lodi_E)
          
          !N layer
          call compute_lodi_terms(
     $         nodes,
     $         dx,
     $         bc_size, ny-bc_size,
     $         epsilon,
     $         edge_inviscid_flux_N,
     $         edge_viscid_flux_N,
     $         compute_cons_lodi_matrix_y,
     $         compute_dev_from_flux_x,
     $         transverse_lodi_N,
     $         viscous_lodi_N)

        end subroutine compute_fluxes_and_lodi_at_the_edges_2ndorder


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary time derivatives
        !> from the LODI transverse and viscous vectors as well
        !> as the fluxes at the edges of the computational domain
        !
        !> @date
        !> 11_09_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> physical model
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array of grid point data
        !
        !>@param x_map
        !> map of the x-coordinates
        !
        !>@param y_map
        !> map of the y-coordinates
        !
        !>@param transverse_lodi_N
        !> LODI transverse vector for the northern edge
        !
        !>@param transverse_lodi_S
        !> LODI transverse vector for the southern edge
        !
        !>@param transverse_lodi_E
        !> LODI transverse vector for the eastern edge
        !
        !>@param transverse_lodi_W
        !> LODI transverse vector for the western edge
        !
        !>@param viscous_lodi_N
        !> LODI viscous vector for the northern edge
        !
        !>@param viscous_lodi_S
        !> LODI viscous vector for the southern edge
        !
        !>@param viscous_lodi_E
        !> LODI viscous vector for the eastern edge
        !
        !>@param viscous_lodi_W
        !> LODI viscous vector for the western edge
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param flow_user_N
        !> user flow configuration choice for the northern boundary
        !
        !>@param flow_user_S
        !> user flow configuration choice for the northern boundary
        !
        !>@param flow_user_E
        !> user flow configuration choice for the northern boundary
        !
        !>@param flow_user_W
        !> user flow configuration choice for the northern boundary
        !-------------------------------------------------------------
        subroutine compute_timedev_at_the_edges_2ndorder(
     $     p_model,
     $     t,nodes,x_map,y_map,
     $     transverse_lodi_N, transverse_lodi_S,
     $     transverse_lodi_E, transverse_lodi_W,
     $     viscous_lodi_N, viscous_lodi_S,
     $     viscous_lodi_E, viscous_lodi_W,
     $     flux_x, flux_y,
     $     timedev,
     $     flow_user_N,
     $     flow_user_S,
     $     flow_user_E,
     $     flow_user_W)

           implicit none           

           type(pmodel_eq)                                , intent(in)    :: p_model
           real(rkind)                                    , intent(in)    :: t
           real(rkind), dimension(nx,ny,ne)               , intent(in)    :: nodes
           real(rkind), dimension(nx)                     , intent(in)    :: x_map
           real(rkind), dimension(ny)                     , intent(in)    :: y_map
           real(rkind), dimension(nx-2*bc_size,bc_size,ne), intent(in)    :: transverse_lodi_N
           real(rkind), dimension(nx-2*bc_size,bc_size,ne), intent(in)    :: transverse_lodi_S
           real(rkind), dimension(bc_size,ny-2*bc_size,ne), intent(in)    :: transverse_lodi_E
           real(rkind), dimension(bc_size,ny-2*bc_size,ne), intent(in)    :: transverse_lodi_W
           real(rkind), dimension(nx-2*bc_size,bc_size,ne), intent(in)    :: viscous_lodi_N
           real(rkind), dimension(nx-2*bc_size,bc_size,ne), intent(in)    :: viscous_lodi_S
           real(rkind), dimension(bc_size,ny-2*bc_size,ne), intent(in)    :: viscous_lodi_E
           real(rkind), dimension(bc_size,ny-2*bc_size,ne), intent(in)    :: viscous_lodi_W
           real(rkind), dimension(nx+1,ny,ne)             , intent(in)    :: flux_x
           real(rkind), dimension(nx,ny+1,ne)             , intent(in)    :: flux_y
           real(rkind), dimension(nx,ny,ne)               , intent(inout) :: timedev
           integer                                        , intent(in)    :: flow_user_N
           integer                                        , intent(in)    :: flow_user_S
           integer                                        , intent(in)    :: flow_user_E
           integer                                        , intent(in)    :: flow_user_W


           integer :: j
           integer :: j_offset
           logical :: side_y
           integer :: flow_y_user
           
           integer :: i
           integer :: i_offset
           logical :: side_x
           integer :: flow_x_user

           type(lodi_edge_inflow)  :: edge_inflow_bc
           type(lodi_edge_outflow) :: edge_outflow_bc


           !S layer (j=1)
           j           = 1
           j_offset    = 0
           side_y      = left
           flow_y_user = flow_user_S
           !gradient_y = gradient_y_y_oneside_L0

           call compute_timedev_y_layer_2ndorder(
     $          p_model,
     $          t,nodes,x_map,y_map,
     $          transverse_lodi_S,
     $          viscous_lodi_S,
     $          flow_user_W,
     $          flow_user_E,
     $          flux_x,
     $          gradient_y_y_oneside_L0,
     $          j,
     $          j_offset,
     $          side_y,
     $          flow_y_user,
     $          timedev)


           !S layer (j=2)
           j           = 2
           j_offset    = 0
           side_y      = left
           flow_y_user = flow_user_S
           !gradient_y = gradient_y_y_oneside_L1

           call compute_timedev_y_layer_2ndorder(
     $          p_model,
     $          t,nodes,x_map,y_map,
     $          transverse_lodi_S, 
     $          viscous_lodi_S,
     $          flow_user_W,
     $          flow_user_E,
     $          flux_x,
     $          gradient_y_y_oneside_L1,
     $          j,
     $          j_offset,
     $          side_y,
     $          flow_y_user,
     $          timedev)


           !W and E layers (j \in [3,ny-2])
           do j=bc_size+1, ny-bc_size

              !W layers
              i_offset    = 0
              side_x      = left
              flow_x_user = flow_user_W

              !W layer (i=1)
              i=1
              call compute_timedev_x_edge(
     $             p_model,
     $             t,nodes,x_map,y_map,i,j,
     $             flux_y,
     $             transverse_lodi_W(i-i_offset,j-bc_size,:),
     $             viscous_lodi_W(i,j-bc_size,:),
     $             side_x,
     $             gradient_x_x_oneside_L0,
     $             edge_inflow_bc,
     $             edge_outflow_bc,
     $             flow_x_user,
     $             timedev)
              
              !W layer (i=2)
              i=2
              call compute_timedev_x_edge(
     $             p_model,
     $             t,nodes,x_map,y_map,i,j,
     $             flux_y,
     $             transverse_lodi_W(i-i_offset,j-bc_size,:),
     $             viscous_lodi_W(i,j-bc_size,:),
     $             side_x,
     $             gradient_x_x_oneside_L1,
     $             edge_inflow_bc,
     $             edge_outflow_bc,
     $             flow_x_user,
     $             timedev)


              !E layers
              i_offset    = nx-bc_size
              side_x      = right
              flow_x_user = flow_user_E

              !E layer (i=nx-1)
              i=nx-1
              call compute_timedev_x_edge(
     $             p_model,
     $             t,nodes,x_map,y_map,i,j,
     $             flux_y,
     $             transverse_lodi_E(i-i_offset,j-bc_size,:),
     $             viscous_lodi_E(i-i_offset,j-bc_size,:),
     $             side_x,
     $             gradient_x_x_oneside_R1,
     $             edge_inflow_bc,
     $             edge_outflow_bc,
     $             flow_x_user,
     $             timedev)

              !E layer (i=nx)
              i=nx
              call compute_timedev_x_edge(
     $             p_model,
     $             t,nodes,x_map,y_map,i,j,
     $             flux_y,
     $             transverse_lodi_E(i-i_offset,j-bc_size,:),
     $             viscous_lodi_E(i-i_offset,j-bc_size,:),
     $             side_x,
     $             gradient_x_x_oneside_R0,
     $             edge_inflow_bc,
     $             edge_outflow_bc,
     $             flow_x_user,
     $             timedev)
              
           end do


           !N layer (j=ny-1)
           j           = ny-1
           j_offset    = ny-bc_size
           side_y      = right
           flow_y_user = flow_user_N
           !gradient_y = gradient_y_y_oneside_R1

           call compute_timedev_y_layer_2ndorder(
     $          p_model,
     $          t,nodes,x_map,y_map,
     $          transverse_lodi_N, 
     $          viscous_lodi_N,
     $          flow_user_W,
     $          flow_user_E,
     $          flux_x,
     $          gradient_y_y_oneside_R1,
     $          j,
     $          j_offset,
     $          side_y,
     $          flow_y_user,
     $          timedev)


           !N layer (j=ny)
           j           = ny
           j_offset    = ny-bc_size
           side_y      = right
           flow_y_user = flow_user_N
           !gradient_y = gradient_y_y_oneside_R0

           call compute_timedev_y_layer_2ndorder(
     $          p_model,
     $          t,nodes,x_map,y_map,
     $          transverse_lodi_N, 
     $          viscous_lodi_N,
     $          flow_user_W,
     $          flow_user_E,
     $          flux_x,
     $          gradient_y_y_oneside_R0,
     $          j,
     $          j_offset,
     $          side_y,
     $          flow_y_user,
     $          timedev)

        end subroutine compute_timedev_at_the_edges_2ndorder


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the boundary time derivatives
        !> from the LODI transverse and viscous vectors as well
        !> as the fluxes at the edges of the computational domain
        !
        !> @date
        !> 11_09_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> physical model
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> array of grid point data
        !
        !>@param x_map
        !> map of the x-coordinates
        !
        !>@param y_map
        !> map of the y-coordinates
        !
        !>@param transverse_lodi
        !> LODI transverse vector for the y-layer
        !
        !>@param viscous_lodi
        !> LODI viscous vector for the y-layer
        !
        !>@param flow_user_W
        !> user flow configuration choice for the W layer
        !
        !>@param flow_user_E
        !> user flow configuration choice for the E layer
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param gradient_y
        !> gradient procedure along the y-direction
        !
        !>@param j
        !> index identifying the y-layer computed
        !
        !>@param j_offset
        !> index needed to match the data from the transverse_lodi
        !> and viscous_lodi tables with the nodes tables
        !
        !>@param side_y
        !> boolean identifying the spatial location of the b.c. (left or right)
        !
        !>@param flow_y_user
        !> user flow configuration choice for the y layer
        !
        !>@param timedev
        !> time derivatives
        !-------------------------------------------------------------
        subroutine compute_timedev_y_layer_2ndorder(
     $     p_model,
     $     t,nodes,x_map,y_map,
     $     transverse_lodi, 
     $     viscous_lodi,
     $     flow_user_W,
     $     flow_user_E,
     $     flux_x,
     $     gradient_y,
     $     j,
     $     j_offset,
     $     side_y,
     $     flow_y_user,
     $     timedev)

          implicit none

          type(pmodel_eq)                                , intent(in)    :: p_model
          real(rkind)                                    , intent(in)    :: t
          real(rkind), dimension(nx,ny,ne)               , intent(in)    :: nodes
          real(rkind), dimension(nx)                     , intent(in)    :: x_map
          real(rkind), dimension(ny)                     , intent(in)    :: y_map
          real(rkind), dimension(nx-2*bc_size,bc_size,ne), intent(in)    :: transverse_lodi
          real(rkind), dimension(nx-2*bc_size,bc_size,ne), intent(in)    :: viscous_lodi
          integer                                        , intent(in)    :: flow_user_W
          integer                                        , intent(in)    :: flow_user_E
          real(rkind), dimension(nx+1,ny,ne)             , intent(in)    :: flux_x
          procedure(gradient_y_proc)                                     :: gradient_y
          integer(ikind)                                 , intent(in)    :: j
          integer(ikind)                                 , intent(in)    :: j_offset
          logical                                        , intent(in)    :: side_y
          integer                                        , intent(in)    :: flow_y_user
          real(rkind), dimension(nx,ny,ne)               , intent(inout) :: timedev

          logical :: side_x
          integer :: flow_x_user
          integer :: i

          type(lodi_edge_inflow)            :: edge_inflow_bc
          type(lodi_edge_outflow)           :: edge_outflow_bc
          type(lodi_corner_inflow_inflow)   :: corner_inflow_inflow_bc
          type(lodi_corner_inflow_outflow)  :: corner_inflow_outflow_bc
          type(lodi_corner_outflow_outflow) :: corner_outflow_outflow_bc


          !(S or N) W corner (i={1,2},j)
          !------------------------------
          side_x = left
          flow_x_user = flow_user_W

          i=1
          call compute_timedev_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x, side_y,
     $         gradient_x_x_oneside_L0,
     $         gradient_y,
     $         corner_inflow_inflow_bc,
     $         corner_inflow_outflow_bc,
     $         corner_outflow_outflow_bc,
     $         flow_x_user, flow_y_user,
     $         timedev)

          i=2
          call compute_timedev_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x, side_y,
     $         gradient_x_x_oneside_L1,
     $         gradient_y,
     $         corner_inflow_inflow_bc,
     $         corner_inflow_outflow_bc,
     $         corner_outflow_outflow_bc,
     $         flow_x_user, flow_y_user,
     $         timedev)


          !(N or S) layer (i\in[3,nx-3],j)
          !------------------------------
          do i=bc_size+1, nx-bc_size
             call compute_timedev_y_edge(
     $            p_model,
     $            t,nodes,x_map,y_map,i,j,
     $            flux_x,
     $            transverse_lodi(i-bc_size,j-j_offset,:),
     $            viscous_lodi(i-bc_size,j-j_offset,:),
     $            side_y,
     $            gradient_y,
     $            edge_inflow_bc,
     $            edge_outflow_bc,
     $            flow_y_user,
     $            timedev)
          end do


          !(N or S) E corner (i={nx-1,nx},j)
          !------------------------------
          side_x = right
          flow_x_user = flow_user_E

          i=nx-1
          call compute_timedev_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x, side_y,
     $         gradient_x_x_oneside_R1,
     $         gradient_y,
     $         corner_inflow_inflow_bc,
     $         corner_inflow_outflow_bc,
     $         corner_outflow_outflow_bc,
     $         flow_x_user, flow_y_user,
     $         timedev)

          i=nx
          call compute_timedev_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x, side_y,
     $         gradient_x_x_oneside_R0,
     $         gradient_y,
     $         corner_inflow_inflow_bc,
     $         corner_inflow_outflow_bc,
     $         corner_outflow_outflow_bc,
     $         flow_x_user, flow_y_user,
     $         timedev)

        end subroutine compute_timedev_y_layer_2ndorder

      end module bc_operators_class
