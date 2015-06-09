      !> @file
      !> bf_interface_dyn augmented with time integration procedures
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_interface_dyn augmented with time integration procedures
      !
      !> @date
      ! 11_03_2015 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_interface_time_class

        use bc_operators_gen_class, only :
     $       bc_operators_gen

        use bf_interface_dyn_class, only :
     $       bf_interface_dyn

        use bf_mainlayer_bc_sections_module, only :
     $       update_interior_bc_sections_from_mainlayers

        use bf_sublayer_class, only :
     $       bf_sublayer

        use interface_integration_step, only :
     $       timeInt_step_nopt

        use parameters_input, only :
     $       nx,ny,ne,
     $       obc_edge_overlap_ac

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
        public :: bf_interface_time

        
        !>@class bf_interface_time
        !> bf_interface_dyn augmented with time integration
        !> procedures
        !
        !> @param apply_initial_conditions
        !> apply the initial conditions of the physical model
        !> in the domain extension
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
        !------------------------------------------------------------
        type, extends(bf_interface_dyn) :: bf_interface_time

          contains
          
          procedure, pass :: apply_initial_conditions

          procedure, pass :: initialize_before_timeInt
          procedure, pass :: finalize_after_timeInt
          procedure, pass :: compute_time_dev
          procedure, pass :: compute_integration_step

        end type bf_interface_time


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> apply the initial conditions of the physical model
        !> in the domain extension
        !
        !> @date
        !> 26_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_interface_dyn with time integration procedures
        !
        !> @param p_model
        !> physical model
        !--------------------------------------------------------------
        subroutine apply_initial_conditions(this,p_model)

          implicit none

          class(bf_interface_time), intent(inout) :: this
          type(pmodel_eq)         , intent(in)    :: p_model

          integer :: k
          integer :: nb_sublayers
          integer :: m

          type(bf_sublayer), pointer :: bf_sublayer_ptr


          !> loop over the main layers and apply the initial
          !> conditions on each buffer layer of the domain
          !> extension
          do k=1,4

             nb_sublayers = this%mainlayer_pointers(k)%get_nb_sublayers()

             bf_sublayer_ptr => this%mainlayer_pointers(k)%get_head_sublayer()

             do m=1, nb_sublayers

                call bf_sublayer_ptr%apply_initial_conditions(p_model)

                bf_sublayer_ptr => bf_sublayer_ptr%get_next()

             end do
          
          end do

        end subroutine apply_initial_conditions


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the bc_sections for the interior and initialize
        !> the buffer layers for the time integration
        !
        !> @date
        !> 11_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_interface_dyn with time integration procedures
        !
        !> @param interior_bc_sections
        !> identification of the bc_sections in the interior domain
        !--------------------------------------------------------------
        subroutine initialize_before_timeInt(this,interior_bc_sections)
        
          implicit none

          class(bf_interface_time)                   , intent(inout) :: this
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: interior_bc_sections


          integer :: k
          integer :: m
          integer :: nb_sublayers

          type(bf_sublayer), pointer :: bf_sublayer_ptr
          
          
          !> update the bc_sections for the interior domain
          call update_interior_bc_sections_from_mainlayers(
     $         this%mainlayer_pointers,
     $         interior_bc_sections)


          !> ask the buffer layers to initialize the intermediate
          !> variables for the time integration
          do k=1, size(this%mainlayer_pointers,1)

             if(this%mainlayer_pointers(k)%associated_ptr()) then

                ! initialize by taking into account only data
                ! from the buffer layer
                call this%mainlayer_pointers(k)%initialize_before_timeInt()


                ! modify the bc_sections initialized by the buffer layer
                ! by overlapping edge and anti-corners using information
                ! from the neighboring buffer layers
                if(obc_edge_overlap_ac) then

                   bf_sublayer_ptr => this%mainlayer_pointers(k)%get_head_sublayer()
                   nb_sublayers = this%mainlayer_pointers(k)%get_nb_sublayers()

                   do m=1, nb_sublayers

                      call this%mainlayer_interfaces%update_bc_sections(bf_sublayer_ptr)
                      bf_sublayer_ptr => bf_sublayer_ptr%get_next()

                   end do

                end if

             end if
             
          end do

        end subroutine initialize_before_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> finalize the buffer layers after the time integration
        !
        !> @date
        !> 11_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_interface_dyn with time integration procedures
        !--------------------------------------------------------------
        subroutine finalize_after_timeInt(this)
        
          implicit none

          class(bf_interface_time), intent(inout) :: this

          
          integer :: k
          

          !> ask the buffer layers to finalize the intermediate
          !> variables after the time integration
          do k=1, size(this%mainlayer_pointers,1)

             call this%mainlayer_pointers(k)%finalize_after_timeInt()
             
          end do

        end subroutine finalize_after_timeInt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives in the buffer layers
        !
        !> @date
        !> 11_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_interface_dyn with time integration procedures
        !--------------------------------------------------------------
        subroutine compute_time_dev(
     $     this,
     $     td_operators_used,
     $     t,s,p_model,bc_used,
     $     interior_nodes)

          implicit none

          class(bf_interface_time)        , intent(inout) :: this
          type(td_operators)              , intent(in)    :: td_operators_used
          real(rkind)                     , intent(in)    :: t
          type(sd_operators)              , intent(in)    :: s
          type(pmodel_eq)                 , intent(in)    :: p_model
          type(bc_operators_gen)          , intent(in)    :: bc_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes

          
          integer :: k


          !> ask the buffer layers to compute their time
          !> derivatives
          do k=1, size(this%mainlayer_pointers,1)
          
             if(this%mainlayer_pointers(k)%associated_ptr()) then
                
                call this%mainlayer_pointers(k)%compute_time_dev(
     $               td_operators_used,
     $               t,s,p_model,bc_used,
     $               interior_nodes)

             end if
          end do

        end subroutine compute_time_dev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives in the buffer layers
        !
        !> @date
        !> 11_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> bf_interface_dyn with time integration procedures
        !--------------------------------------------------------------
        subroutine compute_integration_step(
     $     this, dt, integration_step_nopt,
     $     interior_nodes)

          implicit none

          class(bf_interface_time)        , intent(inout) :: this
          real(rkind)                     , intent(in)    :: dt
          procedure(timeInt_step_nopt)                    :: integration_step_nopt
          real(rkind), dimension(nx,ny,ne), intent(inout) :: interior_nodes

          
          integer :: k


          !> ask the buffer layers to compute the
          !> time integration step
          do k=1, size(this%mainlayer_pointers,1)
          
             if(this%mainlayer_pointers(k)%associated_ptr()) then
                
                call this%mainlayer_pointers(k)%compute_integration_step(
     $               dt, integration_step_nopt)

             end if
          end do


          !> synchronize the nodes between the interior
          !> domain and the buffer layers
          call this%sync_nodes(interior_nodes)          

        end subroutine compute_integration_step          

      end module bf_interface_time_class
