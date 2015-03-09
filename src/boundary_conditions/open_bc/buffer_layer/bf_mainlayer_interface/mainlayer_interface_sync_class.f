      !> @file
      !> mainlayer_interface_basic augmented with synchronization
      !> methods
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> mainlayer_interface_basic augmented with synchronization
      !> methods
      !
      !> @date
      ! 09_03_2015 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module mainlayer_interface_sync_class

        use mainlayer_interface_basic_class, only :
     $     mainlayer_interface_basic

        use bf_layer_errors_module, only :
     $       error_mainlayer_interface_type

        use parameters_bf_layer, only :
     $       NE_interface_type,
     $       NW_interface_type,
     $       SE_interface_type,
     $       SW_interface_type

        implicit none

        private
        public :: mainlayer_interface_sync


        !>@class mainlayer_interface_sync
        !> mainlayer_interface_basic augmented with synchronization
        !> methods
        !
        !>@param sync_nodes_at_mainlayer_interface
        !> synchronize the nodes at the interface if there are
        !> two buffer layers sharing grid points
        !--------------------------------------------------------------
        type, extends(mainlayer_interface_basic) :: mainlayer_interface_sync

          contains

          procedure, pass :: sync_nodes_at_mainlayer_interfaces

        end type mainlayer_interface_sync


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> synchronize the nodes at the mainlayer interface
        !
        !> @date
        !> 09_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating links to buffer layers at the edge
        !> between different main layers
        !------------------------------------------------------------
        subroutine sync_nodes_at_mainlayer_interfaces(this)

          implicit none

          class(mainlayer_interface_sync), intent(inout) :: this

          
          !NW interface
          !------------------------------------------------------------
          if(associated(this%NW_interface_N_ptr).and.
     $       associated(this%NW_interface_W_ptr)) then

             call this%NW_interface_N_ptr%sync_nodes_with_neighbor1(
     $            this%NW_interface_W_ptr)

          end if


          !NE interface
          !------------------------------------------------------------
          if(associated(this%NE_interface_N_ptr).and.
     $       associated(this%NE_interface_E_ptr)) then

             call this%NE_interface_N_ptr%sync_nodes_with_neighbor2(
     $            this%NE_interface_E_ptr)

          end if


          !SW interface
          !------------------------------------------------------------
          if(associated(this%SW_interface_S_ptr).and.
     $       associated(this%SW_interface_W_ptr)) then

             call this%SW_interface_S_ptr%sync_nodes_with_neighbor1(
     $            this%SW_interface_W_ptr)

          end if


          !SE interface
          !------------------------------------------------------------
          if(associated(this%SE_interface_S_ptr).and.
     $       associated(this%SE_interface_E_ptr)) then

             call this%SE_interface_S_ptr%sync_nodes_with_neighbor2(
     $            this%SE_interface_E_ptr)

          end if


        end subroutine sync_nodes_at_mainlayer_interfaces

      end module mainlayer_interface_sync_class
