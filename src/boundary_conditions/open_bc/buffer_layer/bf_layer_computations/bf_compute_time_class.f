      !> @file
      !> class encapsulating the main temporary tables for the time 
      !> integration of the buffer layer
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the main temporary tables for the time 
      !> integration of the buffer layer
      !
      !> @date
      !> 16_07_2014 - initial version         - J.L. Desmarais
      !> 15_10_2014 - interface modifications - J.L. Desmarais
      !> (to have a unique sd_operators, p_model... shared between the
      !>  field and the buffer layer objects)
      !-----------------------------------------------------------------
      module bf_compute_time_class

        use bc_operators_class, only :
     $       bc_operators

        use bf_compute_basic_class, only :
     $       bf_compute_basic

        use interface_integration_step, only :
     $       timeInt_step_nopt

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use td_operators_class, only :
     $       td_operators


        implicit none

        private
        public :: bf_compute_time


        !>@class bf_compute_time
        !> bf_compute_basic augmented with time integration functions
        !
        !>@param compute_time_dev
        !> compute the time_dev attribute
        !
        !>@param compute_integration_step
        !> compute the nodes and nodes_tmp using the integration
        !> procedure
        !---------------------------------------------------------------
        type, extends(bf_compute_basic) :: bf_compute_time

          contains

          procedure,   pass :: compute_time_dev
          procedure,   pass :: compute_integration_step

        end type bf_compute_time

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the grid points of the
        !> buffer layer
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute_time object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param td_operators_used
        !> object encapsulating the functions computing the time
        !> derivatives of the grid points
        !
        !>@param t
        !> time
        !
        !>@param nodes
        !> grid points of the buffer layer whose time derivatives
        !> are computed
        !
        !>@param x_map
        !> coordinates along the x-direction of the buffer layer
        !> grid points
        !
        !>@param y_map
        !> coordinates along the y-direction of the buffer layer
        !> grid points
        !
        !>@param s
        !> object encapsulating the space discretisation methods
        !
        !>@param p_model
        !> object encapsulating the governing equations of the
        !> physical model
        !
        !>@param bc_used
        !> object encapsulating the boundary conditions
        !
        !>@param grdpts_id
        !> identification of the role of the grid points
        !
        !>@param x_borders
        !> border indices along the x-direction for the time 
        !> integration
        !
        !>@param y_borders
        !> border indices along the y-direction for the time 
        !> integration
        !--------------------------------------------------------------
        subroutine compute_time_dev(
     $     this,
     $     td_operators_used,
     $     t, nodes, x_map, y_map,
     $     s,
     $     p_model,bc_used,
     $     bf_alignment,
     $     grdpts_id,
     $     interior_nodes,
     $     bc_sections,
     $     x_borders, y_borders)

          implicit none

          class(bf_compute_time)                     , intent(inout) :: this
          type(td_operators)                         , intent(in)    :: td_operators_used
          real(rkind)                                , intent(in)    :: t
          real(rkind), dimension(:,:,:)              , intent(in)    :: nodes
          real(rkind), dimension(:)                  , intent(in)    :: x_map
          real(rkind), dimension(:)                  , intent(in)    :: y_map
          type(sd_operators)                         , intent(in)    :: s
          type(pmodel_eq)                            , intent(in)    :: p_model
          type(bc_operators)                         , intent(in)    :: bc_used
          integer(ikind), dimension(2,2)             , intent(in)    :: bf_alignment
          integer       , dimension(:,:)             , intent(in)    :: grdpts_id
          real(rkind)   , dimension(nx,ny,ne)        , intent(in)    :: interior_nodes
          integer(ikind), dimension(:,:), allocatable, intent(in)    :: bc_sections
          integer(ikind), dimension(2)               , intent(in)    :: x_borders
          integer(ikind), dimension(2)               , intent(in)    :: y_borders
          
          call td_operators_used%compute_time_dev_nopt(
     $         t,nodes,x_map,y_map,
     $         s,p_model,bc_used,
     $         this%time_dev,
     $         bf_alignment,
     $         grdpts_id,
     $         interior_nodes,
     $         bc_sections,
     $         x_borders=x_borders,
     $         y_borders=y_borders)

        end subroutine compute_time_dev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the integration step
        !
        !> @date
        !> 27_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_compute object encapsulating the main
        !> tables extending the interior domain
        !
        !>@param grdpts_id
        !> identification of the role of the grid points
        !
        !>@param nodes
        !> grid points of the buffer layer whose time derivatives
        !> are computed
        !
        !>@param dt
        !> time
        !
        !>@param x_borders
        !> border indices along the x-direction for the time 
        !> integration
        !
        !>@param y_borders
        !> border indices along the y-direction for the time 
        !> integration
        !
        !>@param integration_step_nopt
        !> function determining the procedure for the time integration
        !--------------------------------------------------------------
        subroutine compute_integration_step(
     $     this,
     $     grdpts_id, nodes, dt,
     $     x_borders, y_borders,
     $     integration_step_nopt)

          implicit none

          class(bf_compute_time)       , intent(inout) :: this
          integer    , dimension(:,:)  , intent(in)    :: grdpts_id
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind)                  , intent(in)    :: dt
          integer(ikind), dimension(2) , intent(in)    :: x_borders
          integer(ikind), dimension(2) , intent(in)    :: y_borders
          procedure(timeInt_step_nopt)                 :: integration_step_nopt

          call integration_step_nopt(
     $         nodes,
     $         dt,
     $         this%nodes_tmp,
     $         this%time_dev,
     $         grdpts_id,
     $         x_borders=x_borders,
     $         y_borders=y_borders)

        end subroutine compute_integration_step


      end module bf_compute_time_class
