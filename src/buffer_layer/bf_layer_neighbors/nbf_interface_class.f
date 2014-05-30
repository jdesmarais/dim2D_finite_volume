      module nbf_interface_class
      
        use bf_sublayer_class  , only : bf_sublayer
        use nbf_list_class     , only : nbf_list
        use parameters_constant, only : N,S,E,W
        use parameters_input   , only : debug
        use parameters_bf_layer, only : align_N, align_S,
     $                                  align_E, align_W

        implicit none


        private
        public :: nbf_interface

        
        !< object saving the links between buffer layers
        !> that are in different main layers but share a
        !> border along the x-direction
        !> the links are saved in the array nbf_links
        !> ex: nbf_links(N,1) : links to the neighbor1
        !>                      of the north main layer
        !>     nbf_links(N,2) : links to the neighbor2
        !>                      of the north main layer
        type nbf_interface

          type(nbf_list), dimension(4,2) :: nbf_links

          contains

          procedure, pass :: ini

          procedure, pass :: link_neighbor1_to_bf_sublayer
          procedure, pass :: link_neighbor2_to_bf_sublayer
          procedure, pass :: update_link_from_neighbor1_to_bf_sublayer
          procedure, pass :: update_link_from_neighbor2_to_bf_sublayer
          procedure, pass :: remove_link_from_neighbor1_to_bf_sublayer
          procedure, pass :: remove_link_from_neighbor2_to_bf_sublayer

          procedure, pass :: update_grdpts_from_neighbors
          procedure, pass :: update_neighbor_grdpts

          procedure, pass :: print_on_screen

        end type nbf_interface


        contains

        subroutine ini(this)
        
          implicit none

          class(nbf_interface), intent(inout) :: this

          integer :: i,j

          do j=1, size(this%nbf_links,2)
             do i=1, size(this%nbf_links,1)
                call this%nbf_links(i,j)%ini()
             end do
          end do

        end subroutine ini


        !< a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers, i.e. this
        !> buffer layer share some gridpoints with other layers
        !> each main layer sharing grid points with this buffer layer
        !> is informed that this buffer layer now exists and should be
        !> considered when updating information
        subroutine link_neighbor1_to_bf_sublayer(this, nbf_sublayer)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer

          integer :: neighbor1_id
          integer :: neighbor_index
          
          call nbf_sublayer%get_neighbor1_id(neighbor1_id, neighbor_index)
          call this%nbf_links(neighbor1_id,neighbor_index)%add_link_in_list(
     $         nbf_sublayer)

        end subroutine link_neighbor1_to_bf_sublayer


        !< a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers, i.e. this
        !> buffer layer share some gridpoints with other layers
        !> each main layer sharing grid points with this buffer layer
        !> is informed that this buffer layer now exists and should be
        !> considered when updating information
        subroutine link_neighbor2_to_bf_sublayer(this, nbf_sublayer)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer

          integer :: neighbor2_id
          integer :: neighbor_index
          
          call nbf_sublayer%get_neighbor2_id(neighbor2_id, neighbor_index)
          call this%nbf_links(neighbor2_id,neighbor_index)%add_link_in_list(
     $         nbf_sublayer)

        end subroutine link_neighbor2_to_bf_sublayer


        !< update the link from : neighbor1 -> nbf_sublayer1
        !>                 to   : neighbor1 -> nbf_sublayer2
        subroutine update_link_from_neighbor1_to_bf_sublayer(
     $     this, nbf_sublayer1, nbf_sublayer2)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer1
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer2
          

          integer :: neighbor1_id
          integer :: neighbor_index
          
          call nbf_sublayer1%get_neighbor1_id(neighbor1_id, neighbor_index)
          call this%nbf_links(neighbor1_id,neighbor_index)%update_link_in_list(
     $         nbf_sublayer1, nbf_sublayer2)

        end subroutine update_link_from_neighbor1_to_bf_sublayer


        !< update the link from : neighbor2 -> nbf_sublayer1
        !>                 to   : neighbor2 -> nbf_sublayer2
        subroutine update_link_from_neighbor2_to_bf_sublayer(
     $     this, nbf_sublayer1, nbf_sublayer2)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer1
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer2
          

          integer :: neighbor2_id
          integer :: neighbor_index
          
          call nbf_sublayer1%get_neighbor2_id(neighbor2_id, neighbor_index)
          call this%nbf_links(neighbor2_id,neighbor_index)%update_link_in_list(
     $         nbf_sublayer1, nbf_sublayer2)

        end subroutine update_link_from_neighbor2_to_bf_sublayer
      

        !< remove the link existing from neighbor1 -> nbf_sublayer
        subroutine remove_link_from_neighbor1_to_bf_sublayer(
     $     this, nbf_sublayer)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer
          

          integer :: neighbor1_id
          integer :: neighbor_index
          
          call nbf_sublayer%get_neighbor1_id(neighbor1_id, neighbor_index)
          call this%nbf_links(neighbor1_id,neighbor_index)%remove_link_from_list(
     $         nbf_sublayer)

        end subroutine remove_link_from_neighbor1_to_bf_sublayer


        !< remove the link existing from neighbor2 -> nbf_sublayer
        subroutine remove_link_from_neighbor2_to_bf_sublayer(
     $     this, nbf_sublayer)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer
          

          integer :: neighbor2_id
          integer :: neighbor_index
          
          call nbf_sublayer%get_neighbor2_id(neighbor2_id, neighbor_index)
          call this%nbf_links(neighbor2_id,neighbor_index)%remove_link_from_list(
     $         nbf_sublayer)

        end subroutine remove_link_from_neighbor2_to_bf_sublayer


        !> we ask all the main layers that have grid points in common
        !> with the current main layer to send data to the current
        !> buffer layer
        subroutine update_grdpts_from_neighbors(this, nbf_sublayer)

          implicit none

          class(nbf_interface), intent(in)    :: this
          type(bf_sublayer)   , intent(inout) :: nbf_sublayer

          integer :: mainlayer_id
          
          mainlayer_id = nbf_sublayer%get_localization()


          !if the current buffer layer shares gridpoints with the
          !neighboring layers of type 1, the neighboring buffer
          !layers send data to the current buffer
          !layer
          if(nbf_sublayer%shares_grdpts_with_neighbor1()) then
             call this%nbf_links(mainlayer_id,1)%copy_from_neighbors_to_bf_layer(
     $            1, nbf_sublayer)
          end if


          !if the current buffer layer shares gridpoints with the
          !neighboring layers of type 2, the neighboring buffer
          !layers send data to the current buffer
          !layer
          if(nbf_sublayer%shares_grdpts_with_neighbor2()) then
             call this%nbf_links(mainlayer_id,2)%copy_from_neighbors_to_bf_layer(
     $            2, nbf_sublayer)
          end if

        end subroutine update_grdpts_from_neighbors


        !< we ask all the main layers that have grid points in common
        !> with the current main layer to receive data from the current
        !> buffer layer
        subroutine update_neighbor_grdpts(this, nbf_sublayer)

          implicit none

          class(nbf_interface), intent(inout) :: this
          type(bf_sublayer)   , intent(in)    :: nbf_sublayer

          integer :: mainlayer_id          


          mainlayer_id = nbf_sublayer%get_localization()


          !if the current buffer layer shares gridpoints with the
          !neighboring layers of type 1, the neighboring buffer
          !layers are updated with data from the current buffer
          !layer
          if(nbf_sublayer%shares_grdpts_with_neighbor1()) then
             call this%nbf_links(mainlayer_id,1)%copy_to_neighbors_from_bf_layer(
     $            1, nbf_sublayer)
          end if


          !if the current buffer layer shares gridpoints with the
          !neighboring layers of type 2, the neighboring buffer
          !layers are updated with data from the current buffer
          !layer
          if(nbf_sublayer%shares_grdpts_with_neighbor2()) then
             call this%nbf_links(mainlayer_id,2)%copy_to_neighbors_from_bf_layer(
     $            2, nbf_sublayer)
          end if

        end subroutine update_neighbor_grdpts


        subroutine print_on_screen(this)

          implicit none

          class(nbf_interface), intent(in) :: this

          integer     , dimension(4,2) :: neighbors
          character(1), dimension(4)   :: bf_layer_char
          integer                      :: i,j

          neighbors(N,1) = W
          neighbors(N,2) = E
          neighbors(S,1) = W
          neighbors(S,2) = E
          neighbors(E,1) = S
          neighbors(E,2) = N
          neighbors(W,1) = S
          neighbors(W,2) = N

          bf_layer_char = ['N','S','E','W']          

          do j=1, size(this%nbf_links,2)
             do i=1, size(this%nbf_links,1)
                print '(A1,'' --> '',A1)',
     $               bf_layer_char(neighbors(i,j)),
     $               bf_layer_char(i)
                call this%nbf_links(i,j)%print_on_screen()
             end do
          end do

        end subroutine print_on_screen


      end module nbf_interface_class
