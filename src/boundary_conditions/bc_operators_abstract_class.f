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

        use sd_operators_class, only : sd_operators
        use parameters_input  , only : nx,ny,ne
        use parameters_kind   , only : rkind
        use pmodel_eq_class   , only : pmodel_eq


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

          procedure(ini_proc)    ,   pass, deferred :: ini
          procedure(nodes_proc)  ,   pass, deferred :: apply_bc_on_nodes
          procedure(fluxes_proc) , nopass, deferred :: apply_bc_on_fluxes
          procedure(timedev_proc), nopass, deferred :: apply_bc_on_timedev

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
           subroutine timedev_proc(
     $       nodes,dx,dy,
     $       s,p_model,
     $       flux_x,flux_y,
     $       timedev)
           
             import nx,ny,ne
             import pmodel_eq
             import rkind
             import sd_operators
           
             real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
             real(rkind)                       , intent(in)    :: dx
             real(rkind)                       , intent(in)    :: dy
             type(sd_operators)                , intent(in)    :: s
             type(pmodel_eq)                   , intent(in)    :: p_model
             real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
             real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
             real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev
           
           end subroutine timedev_proc

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
