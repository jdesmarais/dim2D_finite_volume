      !> @file
      !> class encapsulating subroutines to apply open boundary
      !> conditions at the edge of the computational domain using
      !> hedstrom boundary conditions with special treatment for the
      !> corners
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to apply open boundary
      !> conditions at the edge of the computational domain using
      !> hedstrom boundary conditions with special treatment for the
      !> corners
      !
      !> @date
      !> 06_08_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_class

        use bc_operators_default_class, only :
     $     bc_operators_default

        use hedstrom_xy_module, only :
     $       compute_timedev_xlayer,
     $       compute_timedev_ylayer

        use hedstrom_ncoords_module, only :
     $       compute_timedev_corner_ncoords

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

        use openbc_operators_module, only :
     $       compute_fluxes_at_the_edges_2ndorder,
     $       incoming_left,
     $       incoming_right

        use pmodel_eq_class, only :
     $       pmodel_eq

        use parameters_constant, only :
     $       bc_timedev_choice,
     $       x_direction,
     $       y_direction

        use parameters_input, only :
     $       nx,ny,ne,bc_size

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

        use sd_operators_fd_ncoords_module, only :
     $       gradient_n1_oneside_L0,
     $       gradient_n1_oneside_L1,
     $       gradient_n1_oneside_R1,
     $       gradient_n1_oneside_R0,
     $       gradient_n2_oneside_L0,
     $       gradient_n2_oneside_L1,
     $       gradient_n2_oneside_R1,
     $       gradient_n2_oneside_R0

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
        public :: bc_operators



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

          contains

          procedure,   pass :: ini
          procedure, nopass :: apply_bc_on_timedev => apply_bc_on_timedev_2ndorder_corners

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
        !> 04_08_2014 - initial version - J.L. Desmarais
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

        end subroutine ini


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
        !> 06_08_2014 - initial version - J.L. Desmarais
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
        subroutine apply_bc_on_timedev_2ndorder_corners(
     $    nodes,dx,dy,
     $    p_model,
     $    flux_x,flux_y,
     $    timedev)
        
          implicit none
        
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev


          type(sd_operators_x_oneside_L0) :: s_x_L0
          type(sd_operators_x_oneside_L1) :: s_x_L1
          type(sd_operators_x_oneside_R1) :: s_x_R1
          type(sd_operators_x_oneside_R0) :: s_x_R0
          type(sd_operators_y_oneside_L0) :: s_y_L0
          type(sd_operators_y_oneside_L1) :: s_y_L1
          type(sd_operators_y_oneside_R1) :: s_y_R1
          type(sd_operators_y_oneside_R0) :: s_y_R0


          integer(ikind) :: i,j


          !compute the fluxes at the edge of the computational
          !domain
          call compute_fluxes_at_the_edges_2ndorder(
     $         nodes, dx, dy,
     $         s_x_L0, s_x_L1, s_x_R1, s_x_R0,
     $         s_y_L0, s_y_L1, s_y_R1, s_y_R0,
     $         p_model,
     $         flux_x, flux_y)


          !apply the boundary conditions on the south layer
          j=1
          call compute_timedev_corner_ncoords(
     $         nodes, [1,2], j, dx, dy, p_model,
     $         gradient_n2_oneside_L0,
     $         gradient_n2_oneside_L0,
     $         incoming_left,
     $         y_direction,
     $         timedev)

          call compute_timedev_ylayer(
     $         nodes, j, dx, dy, p_model, flux_x,
     $         gradient_y_y_oneside_L0, incoming_left,
     $         timedev)

          call compute_timedev_corner_ncoords(
     $         nodes, [nx-1,nx], j, dx, dy, p_model,
     $         gradient_n1_oneside_R0,
     $         gradient_n1_oneside_R0,
     $         incoming_right,
     $         x_direction,
     $         timedev)

          j=2
          call compute_timedev_corner_ncoords(
     $         nodes, [1,2], j, dx, dy, p_model,
     $         gradient_n2_oneside_L0,
     $         gradient_n2_oneside_L1,
     $         incoming_left,
     $         y_direction,
     $         timedev)

          call compute_timedev_ylayer(
     $         nodes, j, dx, dy, p_model, flux_x,
     $         gradient_y_y_oneside_L1, incoming_left,
     $         timedev)

          call compute_timedev_corner_ncoords(
     $         nodes, [nx-1,nx], j, dx, dy, p_model,
     $         gradient_n1_oneside_R1,
     $         gradient_n1_oneside_R0,
     $         incoming_right,
     $         x_direction,
     $         timedev)


          !apply the boundary conditions on the west and east
          !layers
          do j=bc_size+1, ny-bc_size

             i=1
             call compute_timedev_xlayer(
     $            nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_L0, incoming_left,
     $            timedev)

             i=bc_size
             call compute_timedev_xlayer(
     $            nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_L1, incoming_left,
     $            timedev)

             i=nx-1
             call compute_timedev_xlayer(
     $            nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_R1, incoming_right,
     $            timedev)

             i=nx
             call compute_timedev_xlayer(
     $            nodes, i,j, dx,dy, p_model,  flux_y,
     $            gradient_x_x_oneside_R0, incoming_right,
     $            timedev)

          end do


          !apply the boundary conditions on the north layer
          j=ny-1
          call compute_timedev_corner_ncoords(
     $         nodes, [1,2], j, dx, dy, p_model,
     $         gradient_n1_oneside_L0,
     $         gradient_n1_oneside_L1,
     $         incoming_left,
     $         x_direction,
     $         timedev)

          call compute_timedev_ylayer(
     $         nodes, j, dx, dy, p_model, flux_x,
     $         gradient_y_y_oneside_R1, incoming_right,
     $         timedev)

          call compute_timedev_corner_ncoords(
     $         nodes, [nx-1,nx], j, dx, dy, p_model,
     $         gradient_n2_oneside_R1,
     $         gradient_n2_oneside_R0,
     $         incoming_right,
     $         y_direction,
     $         timedev)

          j=ny
          call compute_timedev_corner_ncoords(
     $         nodes, [1,2], j, dx, dy, p_model,
     $         gradient_n1_oneside_L0,
     $         gradient_n1_oneside_L0,
     $         incoming_left,
     $         x_direction,
     $         timedev)

          call compute_timedev_ylayer(
     $         nodes, j, dx, dy, p_model, flux_x,
     $         gradient_y_y_oneside_R0, incoming_right,
     $         timedev)

          call compute_timedev_corner_ncoords(
     $         nodes, [nx-1,nx], j, dx, dy, p_model,
     $         gradient_n2_oneside_R0,
     $         gradient_n2_oneside_R0,
     $         incoming_right,
     $         y_direction,
     $         timedev)
        
        end subroutine apply_bc_on_timedev_2ndorder_corners

      end module bc_operators_class
