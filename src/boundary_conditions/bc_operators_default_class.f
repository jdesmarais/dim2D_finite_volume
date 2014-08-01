      module bc_operators_default_class

        use bc_operators_abstract_class, only : bc_operators_abstract
        use parameters_input           , only : nx,ny,ne
        use parameters_kind            , only : ikind, rkind
        use pmodel_eq_class            , only : pmodel_eq
        use sd_operators_class         , only : sd_operators

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
        !---------------------------------------------------------------
        type, abstract, extends(bc_operators_abstract) :: bc_operators_default

          contains

          procedure,   pass :: apply_bc_on_nodes
          procedure, nopass :: apply_bc_on_fluxes
          procedure, nopass :: apply_bc_on_timedev

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

          stop 'bc_operator%ini() not implemented'

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
     $     nodes,dx,dy,
     $     s,p_model,
     $     flux_x,flux_y,
     $     timedev)

          implicit none
           
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(sd_operators)                , intent(in)    :: s
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev

          real(rkind) :: node,flux,dx_s,dy_s,timedev_s
          integer     :: bc_s,neq

          stop 'bc_operator%apply_bc_on_time_dev() not implemented'

          node=nodes(1,1,1)
          dx_s = dx
          dy_s = dy
          bc_s = s%get_bc_size()
          neq  = p_model%get_eq_nb()
          flux=flux_x(1,1,1)
          flux=flux_y(1,1,1)
          timedev_s = timedev(1,1,1)

        end subroutine apply_bc_on_timedev

      end module bc_operators_default_class
