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
     $      bf_interface_print

        use bf_sublayer_class, only :
     $       bf_sublayer

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
          procedure, pass :: get_neighbor_sublayer

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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the reference to the neighboring buffer layer of the
        !> buffer layer passed as argument
        !
        !> @date
        !> 25_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        !> the data between them
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying the main layer to which the
        !> bf_sublayer is belonging
        !
        !>@param neighbor_type
        !> type of neighbor (1 or 2) for the neihgboring buffer layer
        !
        !>@return bf_sublayer_ptr
        !> reference to the neighboring buffer layer
        !--------------------------------------------------------------!
        function get_neighbor_sublayer(
     $     this,
     $     mainlayer_id,
     $     neighbor_type)
     $     result(bf_sublayer_ptr)

          implicit none

          class(bf_interface_sync)  , intent(in) :: this
          integer                   , intent(in) :: mainlayer_id
          integer                   , intent(in) :: neighbor_type
          type(bf_sublayer), pointer             :: bf_sublayer_ptr


          bf_sublayer_ptr => this%mainlayer_interfaces%get_neighbor_sublayer_ptr(
     $         mainlayer_id, neighbor_type)


        end function get_neighbor_sublayer
          

      end module bf_interface_sync_class
