      !> @file
      !> bf_mainlayer_sync augmented with time integration
      !> procedures
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_mainlayer_sync augmented with time integration
      !> procedures
      !
      !> @date
      ! 11_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_mainlayer_time_class

        use bc_operators_class, only :
     $       bc_operators

        use bf_layer_errors_module, only :
     $       error_mainlayer_id

        use bf_mainlayer_sync_class, only :
     $       bf_mainlayer_sync

        use bf_sublayer_class, only :
     $       bf_sublayer

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
        public :: bf_mainlayer_time
        
        
        !> @class bf_mainlayer
        !> bf_mainlayer_sync augmented with time integration
        !> procedures
        !
        !> @param initialize_before_timeInt
        !> allocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @param finalize_after_timeInt
        !> deallocate memory space for the intermediate
        !> variables needed to perform the time integration
        !
        !> @param compute_time_dev
        !> compute the time derivatives
        !
        !> @param compute_integration_step
        !> compute the integration step
        !---------------------------------------------------------------
        type, extends(bf_mainlayer_sync) :: bf_mainlayer_time

          contains

          procedure, pass :: initialize_before_timeInt
          procedure, pass :: finalize_after_timeInt
          procedure, pass :: compute_time_dev
          procedure, pass :: compute_integration_step

        end type bf_mainlayer_time


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the necessary attributes of the buffer layers
        !> before the time integration
        !>     - update the time integration borders
        !>     - update the localization of the bc_sections
        !>     - allocate the intermediate arrays for the integration
        !
        !> @date
        !> 07_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !--------------------------------------------------------------
        subroutine initialize_before_timeInt(this)

          implicit none

          class(bf_mainlayer_time), intent(inout) :: this


          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: i

          
          !initialize to the pointer to the first sublayer
          !written
          current_sublayer => this%head_sublayer


          !go through the chained list and allocate the 
          !temporary variables for the time integration
          do i=1, this%nb_sublayers

             !update the time integration borders
             call current_sublayer%update_integration_borders()

             !update the position of the bc_sections
             call current_sublayer%update_bc_sections()

             !allocate the temporary variables for the 
             !time integration
             call current_sublayer%allocate_before_timeInt()

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine initialize_before_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate the memory space for the intermediate
        !> variables needed to perform the time integration
        !> for each sublayer contained in this main layer
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !--------------------------------------------------------------
        subroutine finalize_after_timeInt(this)

          implicit none

          class(bf_mainlayer_time), intent(inout) :: this


          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: i

          
          !initialize to the pointer to the first sublayer
          !written
          current_sublayer => this%head_sublayer


          !go through the chained list and allocate the 
          !temporary variables for the time integration
          do i=1, this%nb_sublayers


             !allocate the temporary variables for the 
             !time integration
             call current_sublayer%deallocate_after_timeInt()

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine finalize_after_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the sublayers
        !> contained in this main layer
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !--------------------------------------------------------------
        subroutine compute_time_dev(
     $     this,
     $     td_operators_used,
     $     t,s,p_model,bc_used,
     $     interior_nodes)

          implicit none

          class(bf_mainlayer_time)        , intent(inout) :: this
          type(td_operators)              , intent(in)    :: td_operators_used
          real(rkind)                     , intent(in)    :: t
          type(sd_operators)              , intent(in)    :: s
          type(pmodel_eq)                 , intent(in)    :: p_model
          type(bc_operators)              , intent(in)    :: bc_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes

          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: i

          
          !initialize to the pointer to the first sublayer
          !written
          current_sublayer => this%head_sublayer


          !go through the chained list and allocate the 
          !temporary variables for the time integration
          do i=1, this%nb_sublayers


             !allocate the temporary variables for the 
             !time integration
             call current_sublayer%compute_time_dev(
     $            td_operators_used,
     $            t,s,p_model,bc_used,
     $            interior_nodes)

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine compute_time_dev



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives of the sublayers
        !> contained in this main layer
        !
        !> @date
        !> 17_07_2014 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating the double chained list of sublayers,
        !> pointers to the head and tail elements of the list and the
        !> total number of elements in the list
        !
        !> @param dt
        !> integration time step
        !
        !> @param integration_step_nopt
        !> procedure performing the time integration
        !
        !>@param full
        !> logical to enforce whether the full domain is computed
        !> discarding the x_borders and y_borders supplied
        !> (important for the first integration step when the
        !> previous integration step is saved temporary in another
        !> array)
        !--------------------------------------------------------------
        subroutine compute_integration_step(
     $     this, dt, integration_step_nopt)

          implicit none

          class(bf_mainlayer_time)    , intent(inout) :: this
          real(rkind)                 , intent(in)    :: dt
          procedure(timeInt_step_nopt)                :: integration_step_nopt

          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: i

          
          !initialize to the pointer to the first sublayer
          !written
          current_sublayer => this%head_sublayer


          !go through the chained list and allocate the 
          !temporary variables for the time integration
          do i=1, this%nb_sublayers


             !allocate the temporary variables for the 
             !time integration
             call current_sublayer%compute_integration_step(
     $            dt, integration_step_nopt)

             !get the next sublayer in the mainlayer
             current_sublayer => current_sublayer%get_next()

          end do

        end subroutine compute_integration_step        

      end module bf_mainlayer_time_class
