      module bc_operators_default_class

        use bc_operators_abstract_class, only :
     $       bc_operators_abstract

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
        !> @param apply_bc_on_nodes_nopt
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
        !> @param apply_bc_on_timedev_nopt
        !> default application of the boundary conditions on
        !> the time derivatives : program stops
        !---------------------------------------------------------------
        type, abstract, extends(bc_operators_abstract) :: bc_operators_default

          contains

          procedure, nopass :: apply_bc_on_nodes
          procedure,   pass :: apply_bc_on_nodes_nopt

          procedure, nopass :: apply_bc_on_fluxes

          procedure,   pass :: apply_bc_on_timedev
          procedure,   pass :: apply_bc_on_timedev_nopt

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
        subroutine apply_bc_on_nodes(
     $       bc_section,
     $       t,x_map,y_map,nodes_tmp,
     $       p_model,
     $       nodes)

          implicit none

          integer    , dimension(4)       , intent(in)    :: bc_section
          real(rkind)                     , intent(in)    :: t
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes_tmp
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes

          real(rkind) :: node_s

          stop 'bc_operator%apply_bc_on_nodes() not implemented'

          node_s  = nodes(1,1,1)+t+x_map(1)+y_map(1)+nodes_tmp(1,1,1)+
     $              bc_section(1)+p_model%get_eq_nb()
          

        end subroutine apply_bc_on_nodes


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
        subroutine apply_bc_on_nodes_nopt(this,nodes,bc_sections)

          implicit none

          class(bc_operators_default)                , intent(in)    :: this
          real(rkind)   , dimension(nx,ny,ne)        , intent(inout) :: nodes
          integer(ikind), dimension(:,:), allocatable, intent(in)    :: bc_sections

          real(rkind)           :: node_s          
          integer, dimension(4) :: bc_type

          stop 'bc_operator%apply_bc_on_nodes_opt() not implemented'

          node_s  = nodes(1,1,1)

          if(allocated(bc_sections)) then
             bc_type = this%bc_type
          end if

        end subroutine apply_bc_on_nodes_nopt


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
        subroutine apply_bc_on_fluxes(
     $     bc_section,
     $     t,x_map,y_map,nodes,s,
     $     flux_x,flux_y)

          implicit none

          integer    , dimension(4)         , intent(in)    :: bc_section
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(sd_operators)                , intent(in)    :: s
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y

          real(rkind) :: node,flux
          integer :: bc_s

          stop 'bc_operator%apply_bc_on_fluxes() not implemented'

          node=nodes(1,1,1)+x_map(1)+y_map(1)+t
          bc_s = s%get_bc_size() + bc_section(1)

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
     $     t,x_map,y_map,nodes,
     $     p_model,
     $     flux_x,flux_y,
     $     timedev)

          implicit none
           
          class(bc_operators_default)       , intent(in)    :: this
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev

          real(rkind)           :: node,flux,dx_s,dy_s,timedev_s,t_s
          integer               :: neq
          integer, dimension(4) :: bc_s

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
          bc_s = this%bc_type

        end subroutine apply_bc_on_timedev


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
        !
        !>@param bc_sections
        !> identification of the boundary layers
        !--------------------------------------------------------------
        subroutine apply_bc_on_timedev_nopt(
     $     this,
     $     t,
     $     bf_alignment,
     $     bf_grdpts_id,
     $     bf_x_map,
     $     bf_y_map,
     $     bf_nodes,
     $     interior_nodes,
     $     p_model,
     $     flux_x, flux_y,
     $     bc_sections,
     $     timedev)
        
          implicit none
          
          class(bc_operators_default)                    , intent(in)    :: this
          real(rkind)                                    , intent(in)    :: t
          integer(ikind), dimension(2,2)                 , intent(in)    :: bf_alignment
          integer       , dimension(:,:)                 , intent(in)    :: bf_grdpts_id
          real(rkind)   , dimension(:)                   , intent(in)    :: bf_x_map
          real(rkind)   , dimension(:)                   , intent(in)    :: bf_y_map
          real(rkind)   , dimension(:,:,:)               , intent(in)    :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne)            , intent(in)    :: interior_nodes
          type(pmodel_eq)                                , intent(in)    :: p_model
          real(rkind)   , dimension(:,:,:)               , intent(inout) :: flux_x
          real(rkind)   , dimension(:,:,:)               , intent(inout) :: flux_y
          integer(ikind), dimension(:,:)    , allocatable, intent(in)    :: bc_sections
          real(rkind)   , dimension(:,:,:)               , intent(inout) :: timedev

          real(rkind)           :: node,flux,dx_s,dy_s,timedev_s,t_s
          integer               :: neq
          integer, dimension(4) :: bc_s
          integer               :: bc_sections_s

          stop 'bc_operator%apply_bc_on_time_dev() not implemented'

          !to prevent unused param warnings
          node=bf_nodes(1,1,1)+interior_nodes(1,1,bf_grdpts_id(1,1))
          dx_s = bf_x_map(2)-bf_x_map(1)
          dy_s = bf_y_map(2)-bf_y_map(1)
          t_s  = t
          neq  = p_model%get_eq_nb() + bf_alignment(1,1)
          flux=flux_x(1,1,1)
          flux=flux_y(1,1,1)
          timedev_s = timedev(1,1,1)
          bc_s = this%bc_type
          if(allocated(bc_sections)) then
             bc_sections_s=bc_sections(1,1)
          end if

        end subroutine apply_bc_on_timedev_nopt

      end module bc_operators_default_class
