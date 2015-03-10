      !> @file
      !> bf_interface_print augmented with synchronization procedures
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> bf_interface_print augmented with synchronization procedures
      !
      !> @date
      ! 10_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_interface_sync_class

        use bf_interface_print_class, only :
     $     bf_interface_print

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       rkind

        implicit none

        private
        public :: bf_interface_sync


        !> @class bf_interface_sync
        !> bf_interface_print augmented with synchronization procedures
        !
        !> @param sync_nodes
        !> synchronize nodes between the interior and the buffer layers
        !> as well as at the interface between main layers
        !------------------------------------------------------------
        type, extends(bf_interface_print) :: bf_interface_sync

          contains

          procedure, pass :: sync_nodes

        end type bf_interface_sync


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> synchronize the nodes between the interior and the
        !> buffer layers and at the interface between the main
        !> buffer layers
        !
        !> @date
        !> 10_03_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        !> the data between them
        !
        !>@param interior_nodes
        !> grid points from the interior domain
        !--------------------------------------------------------------!
        subroutine sync_nodes(this,interior_nodes)

          implicit none

          class(bf_interface_sync)        , intent(inout) :: this
          real(rkind), dimension(nx,ny,ne), intent(inout) :: interior_nodes

          integer, dimension(4) :: exch_order
          integer               :: k


          !synchronize the nodes between the interior domain and
          !the main buffer layers
          exch_order = [E,W,N,S]

          do k=1, size(exch_order,1)
             
             if(this%mainlayer_pointers(exch_order(k))%associated_ptr()) then
               
               call this%mainlayer_pointers(exch_order(k))%sync_nodes_with_interior(
     $              interior_nodes)

            end if

          end do


          !synchronize the nodes at the mainlayer interfaces
          call this%mainlayer_interfaces%sync_nodes_at_mainlayer_interfaces()

        end subroutine sync_nodes

      end module bf_interface_sync_class
