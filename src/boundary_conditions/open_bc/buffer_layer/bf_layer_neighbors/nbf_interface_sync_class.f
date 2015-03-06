      !> @file
      !> nbf_interface_basic augmented with procedures for
      !> the synchronization of the buffer layers with the
      !> interior
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> nbf_interface_basic augmented with procedures for
      !> the synchronization of the buffer layers with the
      !> interior
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module nbf_interface_sync_class
      
        use nbf_interface_basic_class, only :
     $       nbf_interface_basic

        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,bc_size,debug
        
        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none


        private
        public :: nbf_interface_sync

        
        !>@class nbf_interface_sync
        !> nbf_interface_basic augmented with procedures for
        !> the synchronization of the buffer layers with the
        !> interior
        !
        !>@param update_grdpts_from_neighbors
        !> ask all the buffer layers that have grid points in common
        !> with the current main layer to send data to the current
        !> buffer layer
        !
        !>@param update_neighbor_grdpts
        !> ask all the main layers that have grid points in common
        !> with the current main layer to receive data from the current
        !> buffer layer
        !
        !>@param sync_interface_nodes
        !> synchronize the nodes located at the interface between
        !> buffer main layers
        !--------------------------------------------------------------
        type, extends(nbf_interface_basic) :: nbf_interface_sync

          contains

          procedure, pass :: update_grdpts_from_neighbors
          procedure, pass :: update_neighbor_grdpts
          procedure, pass :: sync_interface_nodes

        end type nbf_interface_sync

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> ask all the buffer layers that have grid points in common
        !> with the current main layer to send data to the current
        !> buffer layer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_sublayer
        !> bf_sublayer object updated from its references
        !--------------------------------------------------------------
        subroutine update_grdpts_from_neighbors(
     $     this, nbf_sublayer)

          implicit none

          class(nbf_interface_sync), intent(in)    :: this
          type(bf_sublayer)        , intent(inout) :: nbf_sublayer

          integer :: mainlayer_id
          
          mainlayer_id = nbf_sublayer%get_localization()


          !if the current buffer layer shares gridpoints with the
          !neighboring layers of type 1, the neighboring buffer
          !layers send data to the current buffer
          !layer
          if(nbf_sublayer%can_exchange_with_neighbor1()) then

             call this%nbf_links(mainlayer_id,1)%copy_from_neighbors_to_bf_layer(
     $            1, nbf_sublayer)
          end if
          

          !if the current buffer layer shares gridpoints with the
          !neighboring layers of type 2, the neighboring buffer
          !layers send data to the current buffer
          !layer
          if(nbf_sublayer%can_exchange_with_neighbor2()) then
             call this%nbf_links(mainlayer_id,2)%copy_from_neighbors_to_bf_layer(
     $            2, nbf_sublayer)
          end if

        end subroutine update_grdpts_from_neighbors


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> ask all the main layers that have grid points in common
        !> with the current main layer to receive data from the current
        !> buffer layer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_sublayer
        !> bf_sublayer object updated from its references
        !--------------------------------------------------------------
        subroutine update_neighbor_grdpts(this, nbf_sublayer)

          implicit none

          class(nbf_interface_sync), intent(inout) :: this
          type(bf_sublayer)        , intent(in)    :: nbf_sublayer

          integer :: mainlayer_id          


          mainlayer_id = nbf_sublayer%get_localization()


          !if the current buffer layer shares gridpoints with the
          !neighboring layers of type 1, the neighboring buffer
          !layers are updated with data from the current buffer
          !layer
          if(nbf_sublayer%can_exchange_with_neighbor1()) then
             call this%nbf_links(mainlayer_id,1)%copy_to_neighbors_from_bf_layer(
     $            1, nbf_sublayer)
          end if


          !if the current buffer layer shares gridpoints with the
          !neighboring layers of type 2, the neighboring buffer
          !layers are updated with data from the current buffer
          !layer
          if(nbf_sublayer%can_exchange_with_neighbor2()) then
             call this%nbf_links(mainlayer_id,2)%copy_to_neighbors_from_bf_layer(
     $            2, nbf_sublayer)
          end if

        end subroutine update_neighbor_grdpts


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> synchronize the nodes located at the interface between
        !> buffer main layer
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !--------------------------------------------------------------
        subroutine sync_interface_nodes(this)

          implicit none

          class(nbf_interface_sync), intent(inout) :: this


          !synchronize the nodes at the SW interface
          !-----------------------------------------
          !nbf_links(W,1) are references to buffer
          !layers of the S layer that may have grid
          !points in common with the W main layer
          !-----------------------------------------
          !The W layer is a neighbor of type 1 for
          !the S layer
          !-----------------------------------------
          call this%nbf_links(W,1)%sync_nodes_with_neighbor1(
     $         this%nbf_links(S,1))


          !synchronize the nodes at the SE interface
          !-----------------------------------------
          !nbf_links(E,1) are references to buffer
          !layers of the S layer that may have grid
          !points in common with the E main layer
          !-----------------------------------------
          !The E layer is a neighbor of type 2 for
          !the S layer
          !-----------------------------------------
          call this%nbf_links(E,1)%sync_nodes_with_neighbor2(
     $         this%nbf_links(S,2))


          !synchronize the nodes at the NW interface
          !-----------------------------------------
          !nbf_links(W,2) are references to buffer
          !layers of the N layer that may have grid
          !points in common with the W main layer
          !-----------------------------------------
          !The W layer is a neighbor of type 1 for
          !the N layer
          !-----------------------------------------
          call this%nbf_links(W,2)%sync_nodes_with_neighbor1(
     $         this%nbf_links(N,1))


          !synchronize the nodes at the NE interface
          !-----------------------------------------
          !nbf_links(E,2) are references to buffer
          !layers of the N layer that may have grid
          !points in common with the E main layer
          !-----------------------------------------
          !The E layer is a neighbor of type 2 for
          !the N layer
          !-----------------------------------------
          call this%nbf_links(E,2)%sync_nodes_with_neighbor2(
     $         this%nbf_links(N,2))

        end subroutine sync_interface_nodes

      end module nbf_interface_sync_class
