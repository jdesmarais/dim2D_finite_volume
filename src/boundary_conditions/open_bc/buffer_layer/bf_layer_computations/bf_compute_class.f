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
      module bf_compute_class

        use bc_operators_class        , only : bc_operators
        use interface_integration_step, only : timeInt_step_nopt
        use parameters_input          , only : nx,ny,ne
        use parameters_kind           , only : ikind, rkind
        use pmodel_eq_class           , only : pmodel_eq
        use sd_operators_class        , only : sd_operators
        use td_operators_class        , only : td_operators


        implicit none

        private
        public :: bf_compute


        !>@class bf_compute
        !> class encapsulating the main temporary tables for the time 
        !> integration of the buffer layer
        !
        !>@param sd_operators_used
        !> space discretization operators
        !
        !>@param pmodel_eq_used
        !> physical model
        !
        !>@param td_operators_used
        !> time discretization operators
        !
        !>@param nodes_tmp
        !> temporary array whose size is the same as the array containing
        !> the grid points integrated in time
        !
        !>@param time_dev
        !> temporary array containing the time derivatives and whose size
        !> is the same as the array containing the grid points integrated
        !> in time
        !
        !>@param dx
        !> space step along the x-axis
        !
        !>@param dy
        !> space step along the y-axis
        !
        !>@param ini
        !> initialize the attributes dx and dy
        !
        !>@param allocate_tables
        !> allocate the nodes_tmp and time_dev attributes
        !
        !>@param deallocate_tables
        !> deallocate the nodes_tmp and time_dev attributes
        !
        !>@param compute_time_dev
        !> compute the time_dev attribute
        !
        !>@param compute_integration_step
        !> compute the nodes and nodes_tmp using the integration
        !> procedure
        !---------------------------------------------------------------
        type :: bf_compute

          integer    , dimension(:,:)  , allocatable, private :: grdpts_id_tmp
          real(rkind), dimension(:,:,:), allocatable, private :: nodes_tmp
          real(rkind), dimension(:,:,:), allocatable, private :: time_dev

          contains

          procedure, pass :: allocate_tables
          procedure, pass :: deallocate_tables

          procedure, pass :: compute_time_dev
          procedure, pass :: compute_integration_step

          procedure, pass :: get_time_dev !only for tests

        end type bf_compute

        contains


        !allocate the nodes_tmp and time_dev tables
        subroutine allocate_tables(this, size_x, size_y)

          implicit none

          class(bf_compute), intent(inout) :: this
          integer(ikind)   , intent(in)    :: size_x
          integer(ikind)   , intent(in)    :: size_y

          allocate(this%grdpts_id_tmp(size_x,size_y))
          allocate(this%nodes_tmp(size_x,size_y,ne))
          allocate(this%time_dev(size_x,size_y,ne))

        end subroutine allocate_tables


        !allocate the nodes_tmp and time_dev tables
        subroutine deallocate_tables(this)

          implicit none

          class(bf_compute), intent(inout) :: this

          deallocate(this%grdpts_id_tmp)
          deallocate(this%nodes_tmp)
          deallocate(this%time_dev)

        end subroutine deallocate_tables


        !compute the time derivatives
        subroutine compute_time_dev(
     $     this,
     $     nodes,
     $     dx,
     $     dy,
     $     sd_operators_used,
     $     pmodel_eq_used,
     $     bc_operators_used,
     $     td_operators_used,
     $     grdpts_id)

          implicit none

          class(bf_compute)            , intent(inout) :: this
          real(rkind), dimension(:,:,:), intent(in)    :: nodes
          real(rkind)                  , intent(in)    :: dx
          real(rkind)                  , intent(in)    :: dy
          type(sd_operators)           , intent(in)    :: sd_operators_used
          type(pmodel_eq)              , intent(in)    :: pmodel_eq_used
          type(bc_operators)           , intent(in)    :: bc_operators_used
          type(td_operators)           , intent(in)    :: td_operators_used
          integer    , dimension(:,:)  , intent(in)    :: grdpts_id
          
          call td_operators_used%compute_time_dev_nopt(
     $         nodes,
     $         dx,
     $         dy,
     $         sd_operators_used,
     $         pmodel_eq_used,
     $         bc_operators_used,
     $         this%time_dev,
     $         grdpts_id)

        end subroutine compute_time_dev


        !compute the integration step
        subroutine compute_integration_step(
     $     this, grdpts_id, nodes, dt, integration_step_nopt)

          implicit none

          class(bf_compute)            , intent(inout) :: this
          integer    , dimension(:,:)  , intent(in)    :: grdpts_id
          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind)                  , intent(in)    :: dt
          procedure(timeInt_step_nopt) :: integration_step_nopt

          call integration_step_nopt(
     $         nodes, dt, this%nodes_tmp, this%time_dev, grdpts_id)

        end subroutine compute_integration_step


        !get the time_dev attribute
        subroutine get_time_dev(this, time_dev)

          implicit none

          class(bf_compute)                         , intent(in)  :: this
          real(rkind), dimension(:,:,:), allocatable, intent(out) :: time_dev


          if(allocated(this%time_dev)) then
             allocate(time_dev(
     $            size(this%time_dev,1),
     $            size(this%time_dev,2),
     $            size(this%time_dev,3)))

             time_dev = this%time_dev

          else
             print '(''bf_compute_class'')'
             print '(''get_time_dev'')'
             stop 'time dev not allocated'

          end if

        end subroutine get_time_dev        

      end module bf_compute_class
