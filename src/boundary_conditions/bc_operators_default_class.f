      module bc_operators_default_class

        use bc_operators_abstract_class, only :
     $       bc_operators_abstract

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        implicit none

        private
        public :: bc_operators_default


        !> @class bc_operators_default
        !> class encapsulating default subroutines for the
        !> application of boundary conditions
        !
        !> @param ini
        !> default initialization: nothing initialized
        !
        !> @param apply_bc_on_nodes
        !> default application of the boundary conditions on
        !> the nodes : program stops
        !
        !> @param apply_bc_on_fluxes
        !> default application of the boundary conditions on
        !> the fluxes : program stops
        !
        !> @param apply_bc_on_timedev
        !> default application of the boundary conditions on
        !> the time derivatives : program stops
        !
        !> @param apply_bc_on_timedev_x_edge
        !> default application of the boundary conditions on
        !> the time derivatives for E and W edges: program stops
        !
        !> @param apply_bc_on_timedev_y_edge
        !> default application of the boundary conditions on
        !> the time derivatives for N and S edges: program stops
        !
        !> @param apply_bc_on_timedev_xy_corner
        !> default application of the boundary conditions on
        !> the time derivatives for corners: program stops
        !---------------------------------------------------------------
        type, abstract, extends(bc_operators_abstract) :: bc_operators_default

          contains

          procedure,   pass :: apply_bc_on_nodes
          procedure, nopass :: apply_bc_on_fluxes
          procedure,   pass :: apply_bc_on_timedev

          procedure,   pass :: apply_bc_on_timedev_x_edge
          procedure,   pass :: apply_bc_on_timedev_y_edge
          procedure,   pass :: apply_bc_on_timedev_xy_corner

        end type bc_operators_default


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default subroutine applying the boundary conditions
        !> on the grid point governing variables
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary conditions
        !
        !>@param nodes
        !> object encapsulating the main variables
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes(this,nodes)

          implicit none

          class(bc_operators_default)     , intent(in)    :: this
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes

          real(rkind) :: node_s          
          integer :: bcx_type_s

          stop 'bc_operator%apply_bc_on_nodes() not implemented'

          node_s     = nodes(1,1,1)
          bcx_type_s = this%bcx_type
          

        end subroutine apply_bc_on_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default subroutine applying the boundary conditions
        !> on the grid point fluxes
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !
        !>@param dy
        !> grid size along the y-axis
        !
        !>@param s
        !> space discretization operators
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !--------------------------------------------------------------
        subroutine apply_bc_on_fluxes(nodes,dx,dy,s,flux_x,flux_y)

          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(sd_operators)                , intent(in)    :: s
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y

          real(rkind) :: node,flux,dx_s,dy_s
          integer :: bc_s

          stop 'bc_operator%apply_bc_on_fluxes() not implemented'

          node=nodes(1,1,1)
          dx_s = dx
          dy_s = dy
          bc_s = s%get_bc_size()

          flux=flux_x(1,1,1)
          flux=flux_y(1,1,1)

        end subroutine apply_bc_on_fluxes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> default subroutine applying the boundary conditions
        !> on the grid point time derivatives
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !
        !>@param dy
        !> grid size along the y-axis
        !
        !>@param s
        !> space discretization operators
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param timedev
        !> time derivatives of the grid points
        !--------------------------------------------------------------
        subroutine apply_bc_on_timedev(
     $     this,
     $     p_model,
     $     t,nodes,x_map,y_map,
     $     flux_x,flux_y,
     $     timedev)

          implicit none
           
          class(bc_operators_default)       , intent(in)    :: this
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev

          real(rkind) :: node,flux,dx_s,dy_s,timedev_s,t_s
          integer     :: neq,bc_s

          stop 'bc_operator%apply_bc_on_time_dev() not implemented'

          !to prevent unused param warnings
          node=nodes(1,1,1)
          dx_s = x_map(2)-x_map(1)
          dy_s = y_map(2)-y_map(1)
          t_s  = t
          neq  = p_model%get_eq_nb()
          flux=flux_x(1,1,1)
          flux=flux_y(1,1,1)
          timedev_s = timedev(1,1,1)
          bc_s = this%bcx_type

        end subroutine apply_bc_on_timedev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> and x edge: W_edge or E_edge
        !
        !> @date
        !> 21_10_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param t
        !> simulation time for boundary conditions depending
        !> on time
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !
        !>@param dy
        !> grid size along the y-axis
        !
        !>@param i
        !> grid point index along the x-axis
        !
        !>@param j
        !> grid point index along the y-axis
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param side_x
        !> edge side to determine the boundary normal vector
        !
        !>@param gradient_x
        !> procedure to compute the gradient along the x-direction
        !> at (i,j)
        !
        !>@param timedev
        !> time derivatives of the grid points
        !--------------------------------------------------------------
        function apply_bc_on_timedev_x_edge(
     $     this, p_model, t, nodes, dx, dy, i, j, flux_y, side_x, gradient_x)
     $     result(timedev)

          implicit none

          class(bc_operators_default)  , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind), dimension(:,:,:), intent(in) :: flux_y
          logical                      , intent(in) :: side_x
          procedure(gradient_x_proc)                :: gradient_x
          real(rkind), dimension(ne)                :: timedev

          real(rkind)       :: flux,dy_s,t_s
          integer           :: bc_s
          logical           :: side_s

          print '(''bc_operators_default_class'')'
          print '(''apply_bc_on_timedev_x_edge'')'
          stop 'function not implemented'

          !to prevent unused param warnings
          t_s  = t
          dy_s = dy
          flux=flux_y(1,1,1)
          side_s = side_x
          timedev = p_model%compute_x_gradient(nodes,i,j,gradient_x,dx)
          bc_s = this%bcx_type
          

        end function apply_bc_on_timedev_x_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> an y edge: N_edge or S_edge
        !
        !> @date
        !> 21_10_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param t
        !> simulation time for boundary conditions depending
        !> on time
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !
        !>@param dy
        !> grid size along the y-axis
        !
        !>@param i
        !> grid point index along the x-axis
        !
        !>@param j
        !> grid point index along the y-axis
        !
        !>@param flux_x
        !> fluxes along the y-direction
        !
        !>@param side_y
        !> edge side to determine the boundary normal vector
        !
        !>@param gradient_y
        !> procedure to compute the gradient along the y-direction
        !> at (i,j)
        !
        !>@param timedev
        !> time derivatives of the grid points
        !--------------------------------------------------------------
        function apply_bc_on_timedev_y_edge(
     $     this, p_model, t, nodes, dx, dy, i, j, flux_x, side_y, gradient_y)
     $     result(timedev)

          implicit none

          class(bc_operators_default)  , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind), dimension(:,:,:), intent(in) :: flux_x
          logical                      , intent(in) :: side_y
          procedure(gradient_y_proc)                :: gradient_y
          real(rkind), dimension(ne)                :: timedev

          real(rkind)       :: flux,dx_s,t_s
          integer           :: bc_s
          logical           :: side_s

          print '(''bc_operators_default_class'')'
          print '(''apply_bc_on_timedev_y_edge'')'
          stop 'function not implemented'

          !to prevent unused param warnings
          t_s  = t
          dx_s = dx
          flux=flux_x(1,1,1)
          side_s = side_y
          timedev = p_model%compute_y_gradient(nodes,i,j,gradient_y,dy)
          bc_s = this%bcx_type

        end function apply_bc_on_timedev_y_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> a corner: SE_corner, SW_corner, NE_corner, NW_corner
        !
        !> @date
        !> 21_10_2014 - initial version - J.L. Desmarais
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param t
        !> simulation time for boundary conditions depending
        !> on time
        !
        !>@param nodes
        !> object encapsulating the main variables
        !
        !>@param dx
        !> grid size along the x-axis
        !
        !>@param dy
        !> grid size along the y-axis
        !
        !>@param i
        !> grid point index along the x-axis
        !
        !>@param j
        !> grid point index along the y-axis
        !
        !>@param flux_x
        !> fluxes along the y-direction
        !
        !>@param side_y
        !> edge side to determine the boundary normal vector
        !
        !>@param gradient_y
        !> procedure to compute the gradient along the y-direction
        !> at (i,j)
        !
        !>@param timedev
        !> time derivatives of the grid points
        !--------------------------------------------------------------
        function apply_bc_on_timedev_xy_corner(
     $     this, p_model, t, nodes, dx, dy, i, j, side_x, side_y, gradient_x, gradient_y)
     $     result(timedev)

          implicit none

          class(bc_operators_default)  , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          logical                      , intent(in) :: side_x
          logical                      , intent(in) :: side_y
          procedure(gradient_x_proc)                :: gradient_x
          procedure(gradient_y_proc)                :: gradient_y
          real(rkind), dimension(ne)                :: timedev

          real(rkind)       :: t_s
          integer           :: bc_s
          logical           :: side_x_s, side_y_s

          print '(''bc_operators_default_class'')'
          print '(''apply_bc_on_timedev_xy_corner'')'
          stop 'function not implemented'

          !to prevent unused param warnings
          t_s  = t
          side_x_s = side_x
          side_y_s = side_y          
          timedev = p_model%compute_x_gradient(nodes,i,j,gradient_x,dx)
          timedev = p_model%compute_y_gradient(nodes,i,j,gradient_y,dy)
          bc_s = this%bcx_type

        end function apply_bc_on_timedev_xy_corner

      end module bc_operators_default_class
