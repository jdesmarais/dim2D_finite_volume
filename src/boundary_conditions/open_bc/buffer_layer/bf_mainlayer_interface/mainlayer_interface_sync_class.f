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

        use bf_layer_errors_module, only :
     $       error_mainlayer_id,
     $       error_mainlayer_interface_type

        use bf_sublayer_class, only :
     $       bf_sublayer

        use mainlayer_interface_basic_class, only :
     $       mainlayer_interface_basic

        use parameters_bf_layer, only :
     $       NE_interface_type,
     $       NW_interface_type,
     $       SE_interface_type,
     $       SW_interface_type

        use parameters_constant, only :
     $       N,S,E,W

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
        !
        !>@param copy_from_mainlayer_neighbors
        !> copy the nodes and the grdpts_id from the buffer layers
        !> in the mainlayers sharing grid points with the current
        !> buffer layer to the buffer layer
        !
        !>@param copy_to_mainlayer_neighbors
        !> copythe nodes and the grdpts_id of the current buffer
        !> layer in the buffer layers in the mainlayers sharing
        !> grid points with the current buffer layer
        !--------------------------------------------------------------
        type, extends(mainlayer_interface_basic) :: mainlayer_interface_sync

          contains

          procedure, pass :: sync_nodes_at_mainlayer_interfaces

          procedure, pass :: copy_from_mainlayer_neighbors
          procedure, pass :: copy_to_mainlayer_neighbors

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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the nodes and the grdpts_id from the buffer layers
        !> in the mainlayers sharing grid points with the current
        !> buffer layer to the buffer layer
        !
        !> @date
        !> 10_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating links to buffer layers at the edge
        !> between different main layers
        !
        !> @param added_sublayer
        !> current buffer layer updated by its neighbors
        !------------------------------------------------------------
        subroutine copy_from_mainlayer_neighbors(this,added_sublayer)

          implicit none

          class(mainlayer_interface_sync), intent(in)    :: this
          type(bf_sublayer), pointer     , intent(inout) :: added_sublayer

          
          integer :: bf_mainlayer_id


          bf_mainlayer_id = added_sublayer%get_localization()


          select case(bf_mainlayer_id)
            case(N)
               if(added_sublayer%can_exchange_with_neighbor1()) then
                  call added_sublayer%copy_from_neighbor1(this%NW_interface_W_ptr)
               end if

               if(added_sublayer%can_exchange_with_neighbor2()) then
                  call added_sublayer%copy_from_neighbor2(this%NE_interface_E_ptr)
               end if

            case(S)
               if(added_sublayer%can_exchange_with_neighbor1()) then
                  call added_sublayer%copy_from_neighbor1(this%SW_interface_W_ptr)
               end if

               if(added_sublayer%can_exchange_with_neighbor2()) then
                  call added_sublayer%copy_from_neighbor2(this%SE_interface_E_ptr)
               end if

            case(E)
               if(added_sublayer%can_exchange_with_neighbor1()) then
                  call added_sublayer%copy_from_neighbor1(this%SE_interface_S_ptr)
               end if

               if(added_sublayer%can_exchange_with_neighbor2()) then
                  call added_sublayer%copy_from_neighbor2(this%NE_interface_N_ptr)
               end if

            case(W)
               if(added_sublayer%can_exchange_with_neighbor1()) then
                  call added_sublayer%copy_from_neighbor1(this%SW_interface_S_ptr)
               end if

               if(added_sublayer%can_exchange_with_neighbor2()) then
                  call added_sublayer%copy_from_neighbor2(this%NW_interface_N_ptr)
               end if

            case default
               call error_mainlayer_id(
     $              'mainlayer_interface_sync_class',
     $              'copy_from_mainlayer_neighbors',
     $              bf_mainlayer_id)

          end select

        end subroutine copy_from_mainlayer_neighbors


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the nodes and the grdpts_id from the buffer layers
        !> in the mainlayers sharing grid points with the current
        !> buffer layer to the buffer layer
        !
        !> @date
        !> 10_03_2015 - initial version - J.L. Desmarais
        !
        !> @param this
        !> object encapsulating links to buffer layers at the edge
        !> between different main layers
        !
        !> @param added_sublayer
        !> current buffer layer updated by its neighbors
        !------------------------------------------------------------
        subroutine copy_to_mainlayer_neighbors(this,added_sublayer)

          implicit none

          class(mainlayer_interface_sync), intent(inout) :: this
          type(bf_sublayer), pointer     , intent(in)    :: added_sublayer

          
          integer :: bf_mainlayer_id


          bf_mainlayer_id = added_sublayer%get_localization()


          select case(bf_mainlayer_id)
            case(N)
               if(added_sublayer%can_exchange_with_neighbor1()) then
                  call added_sublayer%copy_to_neighbor1(this%NW_interface_W_ptr)
               end if

               if(added_sublayer%can_exchange_with_neighbor2()) then
                  call added_sublayer%copy_to_neighbor2(this%NE_interface_E_ptr)
               end if

            case(S)
               if(added_sublayer%can_exchange_with_neighbor1()) then
                  call added_sublayer%copy_to_neighbor1(this%SW_interface_W_ptr)
               end if

               if(added_sublayer%can_exchange_with_neighbor2()) then
                  call added_sublayer%copy_to_neighbor2(this%SE_interface_E_ptr)
               end if

            case(E)
               if(added_sublayer%can_exchange_with_neighbor1()) then
                  call added_sublayer%copy_to_neighbor1(this%SE_interface_S_ptr)
               end if

               if(added_sublayer%can_exchange_with_neighbor2()) then
                  call added_sublayer%copy_to_neighbor2(this%NE_interface_N_ptr)
               end if

            case(W)
               if(added_sublayer%can_exchange_with_neighbor1()) then
                  call added_sublayer%copy_to_neighbor1(this%SW_interface_S_ptr)
               end if

               if(added_sublayer%can_exchange_with_neighbor2()) then
                  call added_sublayer%copy_to_neighbor2(this%NW_interface_N_ptr)
               end if

            case default
               call error_mainlayer_id(
     $              'mainlayer_interface_sync_class',
     $              'copy_to_mainlayer_neighbors',
     $              bf_mainlayer_id)

          end select

        end subroutine copy_to_mainlayer_neighbors

      end module mainlayer_interface_sync_class
