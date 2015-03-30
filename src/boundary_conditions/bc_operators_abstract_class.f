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

          integer, dimension(4) :: bc_type

          contains

          procedure                ,   pass           :: get_bc_type

          procedure(ini_proc)      ,   pass, deferred :: ini
          procedure(nodes_proc)    ,   pass, deferred :: apply_bc_on_nodes
          procedure(fluxes_proc)   , nopass, deferred :: apply_bc_on_fluxes
          procedure(tdev_proc)     ,   pass, deferred :: apply_bc_on_timedev
          procedure(tdev_nopt_proc),   pass, deferred :: apply_bc_on_timedev_nopt

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
     $       t,
     $       x_map,
     $       y_map,
     $       nodes,
     $       p_model,
     $       flux_x,
     $       flux_y,
     $       timedev)
           
             import bc_operators_abstract
             import nx,ny,ne
             import pmodel_eq
             import rkind
           
             class(bc_operators_abstract)      , intent(in)    :: this
             real(rkind)                       , intent(in)    :: t
             real(rkind), dimension(nx)        , intent(in)    :: x_map
             real(rkind), dimension(ny)        , intent(in)    :: y_map
             real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
             type(pmodel_eq)                   , intent(in)    :: p_model
             real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
             real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
             real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev
           
           end subroutine tdev_proc

      
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
           !>@param p_model
           !> object encapsulating the physical model
           !
           !>@param t
           !> simulation time
           !
           !>@param nodes
           !> gridpoint data
           !
           !>@param x_map
           !> coordinates along the x-direction
           !
           !>@param y_map
           !> coordinates along the y-direction
           !
           !>@param flux_x
           !> fluxes along the x-direction
           !
           !>@param flux_y
           !> fluxes along the y-direction
           !
           !>@param timedev
           !> time derivatives
           !
           !>@param bc_sections
           !> identification of the boundary layers
           !-------------------------------------------------------------
           subroutine tdev_nopt_proc(
     $       this,
     $       t,
     $       bf_alignment,
     $       bf_grdpts_id,
     $       bf_x_map,
     $       bf_y_map,
     $       bf_nodes,
     $       interior_nodes,
     $       p_model,
     $       flux_x, flux_y,
     $       bc_sections,
     $       timedev)
           
             import bc_operators_abstract
             import nx,ny,ne
             import pmodel_eq
             import ikind
             import rkind
           
             class(bc_operators_abstract)                   , intent(in)    :: this
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

           end subroutine tdev_nopt_proc

        end interface


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the types of boundary conditions
        !
        !> @date
        !> 01_08_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
        !
        !>@return type
        !> type of boundary conditions
        !-------------------------------------------------------------
        function get_bc_type(this) result(type)

          implicit none

          class(bc_operators_abstract), intent(in) :: this
          integer, dimension(4)                    :: type

          type = this%bc_type

        end function get_bc_type        

      end module bc_operators_abstract_class
