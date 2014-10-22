      !> @file
      !> abstract class encapsulating subroutines to compute
      !> the time derivatives of the main variables from the
      !> governing equations using the space discretisation
      !> operators and the physical model
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> the time derivatives of the main variables
      !> in the governing equations
      !
      !> @date
      !> 13_08_2013 - initial version               - J.L. Desmarais
      !> 14_07_2014 - interface for erymanthianboar - J.L. Desmarais
      !-----------------------------------------------------------------
      module td_operators_abstract_class

        use bc_operators_class, only : bc_operators
        use sd_operators_class, only : sd_operators
        use pmodel_eq_class   , only : pmodel_eq
        use parameters_input  , only : nx,ny,ne
        use parameters_kind   , only : ikind,rkind

        implicit none


        private
        public :: td_operators_abstract


        !>@class td_operators_abstract
        !> abstract class encapsulating operators to compute
        !> the time derivatives of the main variables from 
        !> the governing equations
        !>
        !>@param compute_time_dev
        !> compute the time derivatives
        !
        !>@param compute_time_dev_nopt
        !> compute the time derivatives without array dimension
        !> optimizations
        !---------------------------------------------------------------
        type, abstract :: td_operators_abstract

          contains
          procedure(time_proc)     , nopass, deferred :: compute_time_dev
          procedure(time_proc_nopt), nopass, deferred :: compute_time_dev_nopt

        end type td_operators_abstract


        abstract interface

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the time derivatives using the
          !> space discretisation operators and the physical model
          !
          !> @date
          !> 13_08_2013 - initial version - J.L. Desmarais
          !
          !>@param nodes
          !> array with the grid point data
          !
          !>@param dx
          !> grid step along the x-axis
          !
          !>@param dy
          !> grid step along the y-axis
          !
          !>@param s
          !> space discretization operators
          !
          !>@param p_model
          !> physical model
          !
          !>@param bc_used
          !> boundary conditions
          !
          !>@param time_dev
          !> time derivatives
          !--------------------------------------------------------------
          function time_proc(t,nodes,x_map,y_map,s,p_model,bc_used)
     $       result(time_dev)

            import bc_operators
            import sd_operators
            import pmodel_eq
            import rkind
            import nx,ny,ne

            real(rkind)                     , intent(in)   :: t
            real(rkind), dimension(nx,ny,ne), intent(in)   :: nodes
            real(rkind), dimension(nx)      , intent(in)   :: x_map
            real(rkind), dimension(ny)      , intent(in)   :: y_map
            type(sd_operators)              , intent(in)   :: s
            type(pmodel_eq)                 , intent(in)   :: p_model
            type(bc_operators)              , intent(in)   :: bc_used
            real(rkind), dimension(nx,ny,ne)               :: time_dev

          end function time_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the time derivatives using the
          !> space discretisation operators and the physical model
          !> without optimizations for the size of the arrays passed
          !> as arguments
          !
          !> @date
          !> 14_07_2014 - initial version - J.L. Desmarais
          !
          !>@param nodes
          !> array with the grid point data
          !
          !>@param dx
          !> grid step along the x-axis
          !
          !>@param dy
          !> grid step along the y-axis
          !
          !>@param s
          !> space discretization operators
          !
          !>@param p_model
          !> physical model
          !
          !>@param bc_used
          !> boundary conditions
          !
          !>@param time_dev
          !> time derivatives
          !
          !>@param grdpts_id
          !> array containing the role of the grid points (interior_pt,
          !> bc_interior_pt, bc_pt, no_pt)
          !
          !>@param bc_sections
          !> array identifying the boundary layers
          !
          !>@param x_borders
          !> array containing the limits of the computed grid points in
          !> the x-direction
          !
          !>@param y_borders
          !> array containing the limits of the computed grid points in
          !> the y-direction
          !--------------------------------------------------------------
          subroutine time_proc_nopt(
     $      t,nodes,x_map,y_map,
     $      s,p_model,bc_used,
     $      time_dev,
     $      grdpts_id,
     $      bc_sections,
     $      x_borders, y_borders)

            import bc_operators
            import ikind
            import pmodel_eq
            import rkind
            import sd_operators

            real(rkind)                                  , intent(in)    :: t
            real(rkind)   , dimension(:,:,:)             , intent(in)    :: nodes
            real(rkind)   , dimension(:)                 , intent(in)    :: x_map
            real(rkind)   , dimension(:)                 , intent(in)    :: y_map
            type(sd_operators)                           , intent(in)    :: s
            type(pmodel_eq)                              , intent(in)    :: p_model
            type(bc_operators)                           , intent(in)    :: bc_used
            real(rkind)   , dimension(:,:,:)             , intent(out)   :: time_dev
            integer       , dimension(:,:)               , intent(in)    :: grdpts_id
            integer       , dimension(:,:)  , allocatable, intent(inout) :: bc_sections
            integer(ikind), dimension(2)                 , intent(in)    :: x_borders
            integer(ikind), dimension(2)                 , intent(in)    :: y_borders

          end subroutine time_proc_nopt

        end interface

      end module td_operators_abstract_class
