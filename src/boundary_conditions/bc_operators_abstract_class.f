      !> @file
      !> abstract class encapsulating subroutine interfaces
      !> to apply boundary conditions at the edge of the
      !> computational domain
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutine interfaces to compute
      !> the gridpoints and/or the fluxes at the edge of the
      !> computational domain
      !
      !> @date
      !> 24_09_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_abstract_class

        use interface_primary, only :
     $     gradient_x_proc,
     $     gradient_y_proc

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
        public :: bc_operators_abstract


        !> @class bc_operators_abstract
        !> abstract class encapsulating subroutine interfaces
        !> to apply boundary conditions at the edge of the
        !> computational domain on the nodes and/or the fluxes
        !
        !>@param ini
        !> interface to apply the initial conditions for the
        !> boundary conditions
        !
        !>@param apply_bc_on_nodes
        !> interface to apply the boundary conditions along
        !> the x and y directions on the nodes at the edge of
        !> the computational domain
        !
        !>@param apply_bc_on_fluxes
        !> interface to apply the boundary conditions along
        !> the x and y directions on the fluxes at the edge of
        !> the computational domain
        !
        !>@param apply_bc_on_timedev
        !> interface to apply the boundary conditions along
        !> the x and y directions on the time derivatives at
        !> the edge of the computational domain
        !---------------------------------------------------------------
        type, abstract :: bc_operators_abstract

          integer :: bcx_type
          integer :: bcy_type

          contains

          procedure,   pass :: get_bcx_type
          procedure,   pass :: get_bcy_type

          procedure(ini_proc)      ,   pass, deferred :: ini
          procedure(nodes_proc)    ,   pass, deferred :: apply_bc_on_nodes
          procedure(fluxes_proc)   , nopass, deferred :: apply_bc_on_fluxes
          procedure(tdev_proc)     ,   pass, deferred :: apply_bc_on_timedev
          procedure(tdev_x_edge)   ,   pass, deferred :: apply_bc_on_timedev_x_edge
          procedure(tdev_y_edge)   ,   pass, deferred :: apply_bc_on_timedev_y_edge
          procedure(tdev_xy_corner),   pass, deferred :: apply_bc_on_timedev_xy_corner

        end type bc_operators_abstract


        abstract interface
           
           !> @author
           !> Julien L. Desmarais
           !
           !> @brief
           !> subroutine initializing the main attributes
           !> of the boundary conditions
           !
           !> @date
           !> 24_09_2013 - initial version - J.L. Desmarais
           !
           !>@param this
           !> abstract boundary conditions
           !
           !>@param s
           !> spatial discretisation operators
           !
           !>@param p_model
           !> physical model
           !-------------------------------------------------------------
           subroutine ini_proc(this,p_model)
        
             import bc_operators_abstract
             import pmodel_eq

             class(bc_operators_abstract), intent(inout) :: this
             type(pmodel_eq)             , intent(in)    :: p_model


           end subroutine ini_proc


           !> @author
           !> Julien L. Desmarais
           !
           !> @brief
           !> subroutine applying the boundary conditions
           !> along the x and y directions at the edge of the
           !> computational domain
           !
           !> @date
           !> 24_09_2013 - initial version - J.L. Desmarais
           !
           !>@param this
           !> abstract boundary conditions
           !
           !>@param f_used
           !> object encapsulating the main variables
           !
           !>@param s
           !> space discretization operators
           !-------------------------------------------------------------
           subroutine nodes_proc(this,nodes)
           
             import bc_operators_abstract
             import nx,ny,ne
             import rkind
           
             class(bc_operators_abstract)    , intent(in)    :: this
             real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes

           end subroutine nodes_proc

      
           !> @author
           !> Julien L. Desmarais
           !
           !> @brief
           !> subroutine applying the boundary conditions
           !> on the fluxes along the x directions at the
           !> edge of the computational domain
           !
           !> @date
           !> 24_09_2013 - initial version - J.L. Desmarais
           !
           !>@param this
           !> abstract boundary conditions
           !
           !>@param f_used
           !> object encapsulating the main variables
           !
           !>@param s
           !> space discretization operators
           !
           !>@param flux_x
           !> fluxes along the x-direction
           !
           !>@param flux_y
           !> fluxes along the y-direction
           !-------------------------------------------------------------
           subroutine fluxes_proc(nodes,dx,dy,s,flux_x,flux_y)
           
             import sd_operators
             import rkind
             import nx,ny,ne
           
             real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
             real(rkind)                       , intent(in)    :: dx
             real(rkind)                       , intent(in)    :: dy
             type(sd_operators)                , intent(in)    :: s
             real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
             real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
           
           end subroutine fluxes_proc


           !> @author
           !> Julien L. Desmarais
           !
           !> @brief
           !> subroutine applying the boundary conditions
           !> on the time derivatives along the x directions
           !> at the edge of the computational domain
           !
           !> @date
           !> 01_08_2014 - initial version - J.L. Desmarais
           !
           !>@param this
           !> abstract boundary conditions
           !
           !>@param f_used
           !> object encapsulating the main variables
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
           !> time derivatives
           !-------------------------------------------------------------
           subroutine tdev_proc(
     $       this,
     $       p_model,
     $       t,nodes,x_map,y_map,
     $       flux_x,flux_y,
     $       timedev)
           
             import bc_operators_abstract
             import nx,ny,ne
             import pmodel_eq
             import rkind
           
             class(bc_operators_abstract)      , intent(in)    :: this
             type(pmodel_eq)                   , intent(in)    :: p_model
             real(rkind)                       , intent(in)    :: t
             real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
             real(rkind), dimension(nx)        , intent(in)    :: x_map
             real(rkind), dimension(ny)        , intent(in)    :: y_map
             real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
             real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
             real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev
           
           end subroutine tdev_proc

      
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
           function tdev_x_edge(
     $        this,
     $        p_model, t,
     $        nodes, x_map, y_map, i,j,
     $        flux_y,
     $        side_x,
     $        gradient_x)
     $        result(timedev)
           
             import bc_operators_abstract
             import gradient_x_proc
             import ikind
             import ne
             import pmodel_eq
             import rkind
           
             class(bc_operators_abstract) , intent(in) :: this
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
           
           end function tdev_x_edge


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
           function tdev_y_edge(
     $        this,
     $        p_model, t,
     $        nodes, x_map, y_map, i,j,
     $        flux_x, side_y, gradient_y)
     $        result(timedev)
           
             import bc_operators_abstract
             import gradient_y_proc
             import ne
             import pmodel_eq
             import ikind
             import rkind
           
             class(bc_operators_abstract) , intent(in) :: this
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
           
           end function tdev_y_edge
           
           
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
           function tdev_xy_corner(
     $        this,
     $        p_model, t,
     $        nodes, x_map, y_map, i,j,
     $        side_x, side_y,
     $        gradient_x, gradient_y)
     $        result(timedev)
           
             import bc_operators_abstract
             import gradient_x_proc
             import gradient_y_proc
             import ikind
             import ne
             import pmodel_eq
             import rkind
           
             class(bc_operators_abstract) , intent(in) :: this
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
           
           end function tdev_xy_corner

        end interface


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the boundary condition type along the x-axis
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@return type
        !> type of boundary condition along the x-axis
        !-------------------------------------------------------------
        function get_bcx_type(this) result(type)

          implicit none

          class(bc_operators_abstract), intent(in) :: this
          integer :: type

          type = this%bcx_type

        end function get_bcx_type


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the boundary condition type along the y-axis
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@return type
        !> type of boundary condition along the y-axis
        !-------------------------------------------------------------
        function get_bcy_type(this) result(type)

          implicit none

          class(bc_operators_abstract), intent(in) :: this
          integer :: type

          type = this%bcy_type

        end function get_bcy_type        

      end module bc_operators_abstract_class
