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

          procedure,   pass :: compute_fluxes_for_bc_x_edge
          procedure,   pass :: compute_fluxes_for_bc_y_edge
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
        !> subroutine computing the fluxes at the edge of the
        !> computational domain in the y-direction so that
        !> the time derivatives for an edge in the x-direction
        !> can be computed
        !
        !> @date
        !> 22_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param j_min
        !> index min along the y-direction corresponding
        !> to the beginning of the edge layer computed
        !
        !>@param j_max
        !> index max along the y-direction corresponding
        !> to the end of the edge layer computed
        !
        !>@param i
        !> index along the x-direction positioning the
        !> the edge boundary layer
        !
        !>@param edge_card_coord
        !> cardinal coordinate identifying the type of
        !> edge boundary layer
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !-------------------------------------------------------------
        subroutine compute_fluxes_for_bc_x_edge(
     $     this,
     $     p_model,
     $     nodes,
     $     s_x_L0, s_x_L1,
     $     s_x_R1, s_x_R0,
     $     dx, dy,
     $     j_min, j_max, i,
     $     edge_card_coord,
     $     flux_y)
        
          implicit none
        
          class(bc_operators_default)       , intent(in)    :: this
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(:,:,:)     , intent(in)    :: nodes
          type(sd_operators_x_oneside_L0)   , intent(in)    :: s_x_L0
          type(sd_operators_x_oneside_L1)   , intent(in)    :: s_x_L1
          type(sd_operators_x_oneside_R1)   , intent(in)    :: s_x_R1
          type(sd_operators_x_oneside_R0)   , intent(in)    :: s_x_R0
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          integer(ikind)                    , intent(in)    :: j_min
          integer(ikind)                    , intent(in)    :: j_max
          integer(ikind)                    , intent(in)    :: i
          integer                           , intent(in)    :: edge_card_coord
          real(rkind), dimension(:,:,:)     , intent(inout) :: flux_y
        
          integer           :: bc_s
          integer           :: neq
          real(rkind)       :: node
          real(rkind)       :: dx_s
          real(rkind)       :: dy_s
          integer(ikind)    :: n_s
          real(rkind)       :: flux

          print '(''bc_operators_default_class'')'
          print '(''compute_fluxes_for_bc_x_edge'')'
          stop 'function not implemented'

          !to prevent unused param warnings
          bc_s = this%bcx_type
          neq  = p_model%get_eq_nb()
          node = nodes(1,1,1)
          bc_s = s_x_L0%get_bc_size()
          bc_s = s_x_L1%get_bc_size()
          bc_s = s_x_R1%get_bc_size()
          bc_s = s_x_R0%get_bc_size()
          dx_s = dx
          dy_s = dy
          n_s  = j_min+j_max+i+edge_card_coord
          flux=flux_y(1,1,1)


        end subroutine compute_fluxes_for_bc_x_edge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the fluxes at the edge of the
        !> computational domain in the x-direction so that
        !> the time derivatives for an edge in the y-direction
        !> can be computed
        !
        !> @date
        !> 22_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@param p_model
        !> object encapsulating the physical model
        !
        !>@param nodes
        !> array containing the grid point data
        !
        !>@param dx
        !> space step along the x-direction
        !
        !>@param dy
        !> space step along the y-direction
        !
        !>@param i_min
        !> index min along the x-direction corresponding
        !> to the beginning of the edge layer computed
        !
        !>@param i_max
        !> index max along the x-direction corresponding
        !> to the end of the edge layer computed
        !
        !>@param j
        !> index along the y-direction positioning the
        !> the edge boundary layer
        !
        !>@param edge_card_coord
        !> cardinal coordinate identifying the type of
        !> edge boundary layer
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !-------------------------------------------------------------
        subroutine compute_fluxes_for_bc_y_edge(
     $     this,
     $     p_model,
     $     nodes,
     $     s_y_L0, s_y_L1,
     $     s_y_R1, s_y_R0,
     $     dx, dy,
     $     i_min, i_max, j,
     $     edge_card_coord,
     $     flux_x)
        
          implicit none
        
          class(bc_operators_default)       , intent(in)    :: this
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(:,:,:)     , intent(in)    :: nodes
          type(sd_operators_y_oneside_L0)   , intent(in)    :: s_y_L0
          type(sd_operators_y_oneside_L1)   , intent(in)    :: s_y_L1
          type(sd_operators_y_oneside_R1)   , intent(in)    :: s_y_R1
          type(sd_operators_y_oneside_R0)   , intent(in)    :: s_y_R0
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          integer(ikind)                    , intent(in)    :: i_min
          integer(ikind)                    , intent(in)    :: i_max
          integer(ikind)                    , intent(in)    :: j
          integer                           , intent(in)    :: edge_card_coord
          real(rkind), dimension(:,:,:)     , intent(inout) :: flux_x


          integer           :: bc_s
          integer           :: neq
          real(rkind)       :: node
          real(rkind)       :: dx_s
          real(rkind)       :: dy_s
          integer(ikind)    :: n_s
          real(rkind)       :: flux

          print '(''bc_operators_default_class'')'
          print '(''compute_fluxes_for_bc_y_edge'')'
          stop 'function not implemented'

          !to prevent unused param warnings
          bc_s = this%bcx_type
          neq  = p_model%get_eq_nb()
          node = nodes(1,1,1)
          bc_s = s_y_L0%get_bc_size()
          bc_s = s_y_L1%get_bc_size()
          bc_s = s_y_R1%get_bc_size()
          bc_s = s_y_R0%get_bc_size()
          dx_s = dx
          dy_s = dy
          n_s  = i_min+i_max+j+edge_card_coord
          flux=flux_x(1,1,1)

        
        end subroutine compute_fluxes_for_bc_y_edge


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
     $     this,
     $     p_model, t,
     $     nodes, x_map, y_map, i,j,
     $     flux_y, side_x, gradient_x)
     $     result(timedev)

          implicit none

          class(bc_operators_default)  , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind), dimension(:,:,:), intent(in) :: flux_y
          logical                      , intent(in) :: side_x
          procedure(gradient_x_proc)                :: gradient_x
          real(rkind), dimension(ne)                :: timedev

          real(rkind)       :: dx
          real(rkind)       :: dy
          real(rkind)       :: flux,t_s
          integer           :: bc_s
          logical           :: side_s

          print '(''bc_operators_default_class'')'
          print '(''apply_bc_on_timedev_x_edge'')'
          stop 'function not implemented'

          !to prevent unused param warnings
          dx   = x_map(2)-x_map(1)
          dy   = y_map(2)-y_map(1)
          t_s  = t
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
     $     this,
     $     p_model, t,
     $     nodes, x_map, y_map, i,j,
     $     flux_x, side_y, gradient_y)
     $     result(timedev)

          implicit none

          class(bc_operators_default)  , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind), dimension(:,:,:), intent(in) :: flux_x
          logical                      , intent(in) :: side_y
          procedure(gradient_y_proc)                :: gradient_y
          real(rkind), dimension(ne)                :: timedev

          real(rkind)       :: dx
          real(rkind)       :: dy
          real(rkind)       :: flux,t_s
          integer           :: bc_s
          logical           :: side_s

          print '(''bc_operators_default_class'')'
          print '(''apply_bc_on_timedev_y_edge'')'
          stop 'function not implemented'

          !to prevent unused param warnings
          dx   = x_map(2)-x_map(1)
          dy   = y_map(2)-y_map(1)
          t_s  = t
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
     $     this,
     $     p_model, t,
     $     nodes, x_map, y_map, i,j,
     $     side_x, side_y,
     $     gradient_x, gradient_y)
     $     result(timedev)

          implicit none

          class(bc_operators_default)  , intent(in) :: this
          type(pmodel_eq)              , intent(in) :: p_model
          real(rkind)                  , intent(in) :: t
          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind), dimension(:)    , intent(in) :: x_map
          real(rkind), dimension(:)    , intent(in) :: y_map
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          logical                      , intent(in) :: side_x
          logical                      , intent(in) :: side_y
          procedure(gradient_x_proc)                :: gradient_x
          procedure(gradient_y_proc)                :: gradient_y
          real(rkind), dimension(ne)                :: timedev

          real(rkind)       :: dx
          real(rkind)       :: dy
          real(rkind)       :: t_s
          integer           :: bc_s
          logical           :: side_x_s, side_y_s

          print '(''bc_operators_default_class'')'
          print '(''apply_bc_on_timedev_xy_corner'')'
          stop 'function not implemented'

          !to prevent unused param warnings
          dx   = x_map(2)-x_map(1)
          dy   = y_map(2)-y_map(1)
          t_s  = t
          side_x_s = side_x
          side_y_s = side_y          
          timedev = p_model%compute_x_gradient(nodes,i,j,gradient_x,dx)
          timedev = p_model%compute_y_gradient(nodes,i,j,gradient_y,dy)
          bc_s = this%bcx_type

        end function apply_bc_on_timedev_xy_corner

      end module bc_operators_default_class
