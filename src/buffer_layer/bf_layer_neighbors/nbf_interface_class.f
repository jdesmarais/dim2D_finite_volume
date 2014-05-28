      module nbf_interface_class
      
        use nbf_list_class     , only : nbf_list
        use parameters_input   , only : debug
        use parameters_bf_layer, only : align_N, align_S,
     $                                  align_E, align_W

        implicit none


        private
        public :: nbf_interface

        
        !< object saving the links between buffer layers
        !> that are in different main layers but share a
        !> border along the x-direction
        type nbf_interface

          type(nbf_list), dimension(8) :: nbf_links

          contains

          procedure, pass, private :: get_index
          procedure, pass          :: add_links_to_neighbor1_bf_layer
          procedure, pass          :: add_links_to_neighbor2_bf_layer
          procedure, pass          :: update_links_to_neighbor1_bf_layer
          procedure, pass          :: update_links_to_neighbor2_bf_layer
          procedure, pass          :: remove_links_to_neighbor1_bf_layer
          procedure, pass          :: remove_links_to_neighbor2_bf_layer

        end type nbf_interface


        contains


        !< from the coordinates identifying the main layer
        !> to which the buffer layer belongs and the index
        !> identifying its neighbor, the index where the links
        !> between the main layer and its neighbors is determined
        !> this index should be used for the nbf_links table
        function get_index(mainlayer_id, neighbor_index)

          implicit none

          integer, intent(in) :: mainlayer_id
          integer, intent(in) :: neighbor_index
          integer             :: get_index

          get_index = (mainlayer_id-1)*2 + neighbor_index

        end function get_index


        !< a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers, i.e. this
        !> buffer layer share some gridpoints with other layers
        !> each main layer sharing grid points with this buffer layer
        !> is informed that this buffer layer has been updated or should
        !> send newer grid points to other main layers
        subroutine add_links_to_neighbor1_bf_layer(this, nbf_sublayer)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer
          

          integer :: mainlayer_id
          integer :: neighbor1_id
          integer :: neighbor2_id
          integer :: neighbor_index
          integer :: nbf_index

          
          !get the cardinal coordinate of the buffer layer
          mainlayer_id = nbf_sublayer%get_localization()

          !find the coordinates identifying the neighboring buffer layers
          !with which this buffer layer may be sharing grid points
          call get_neighbor_coords(
     $       mainlayer_id,
     $       neighbor1_id,
     $       neighbor2_id,
     $       neighbor_index)

          !get the index corresponding to current_bf_layer -> neighbor1
          nbf_index = get_index(mainlayer_id,1)
          !add links for current_bf_layer -> neighbor1
          call this%nbf_links(nbf_index)%add_link_in_list(nbf_sublayer)

          !get the index corresponding to neighbor1 -> current_bf_layer
          nbf_index = get_index(neighbor1_id, neighbor_index)
          !add links for neighbor1 -> current_bf_layer
          call this%nbf_linfs(nbf_index)%add_link_in_list(nbf_sublayer)


        end subroutine add_links_to_neighbor1_bf_layer


        subroutine add_links_to_neighbor2_bf_layer(this, nbf_sublayer)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer
          

          integer :: mainlayer_id
          integer :: neighbor1_id
          integer :: neighbor2_id
          integer :: neighbor_index
          integer :: nbf_index


          !get the cardinal coordinate of the buffer layer
          mainlayer_id = nbf_sublayer%get_localization()

          !find the coordinates identifying the neighboring buffer layers
          !with which this buffer layer may be sharing grid points
          call get_neighbor_coords(
     $       mainlayer_id,
     $       neighbor1_id,
     $       neighbor2_id,
     $       neighbor_index)

          !get the index corresponding to current_bf_layer -> neighbor2
          nbf_index = get_index(mainlayer_id,2)
          !add links for current_bf_layer -> neighbor2
          call this%nbf_links(nbf_index)%add_link_in_list(nbf_sublayer)

          !get the index corresponding to neighbor2 -> current_bf_layer
          nbf_index = get_index(neighbor2_id, neighbor_index)
          !add links for neighbor2 -> current_bf_layer
          call this%nbf_linfs(nbf_index)%add_link_in_list(nbf_sublayer)

        end subroutine add_links_to_neighbor2_bf_layer


        !< a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers, i.e. this
        !> buffer layer share some gridpoints with other layers
        !> this bf layer has been updated (merge) and pointers to
        !> this buffer layer should be updated as they will be no
        !> longer exist. They are replaced by another link
        subroutine update_links_to_neighbor1_bf_layer(this,
     $     nbf_sublayer1, nbf_sublayer2)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer1
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer2
          

          integer :: mainlayer_id, mainlayer_id2
          integer :: neighbor1_id
          integer :: neighbor2_id
          integer :: neighbor_index
          integer :: nbf_index
          
          !get the cardinal coordinate of the buffer layer
          mainlayer_id = nbf_sublayer1%get_localization()
          if(debug) then
             mainlayer_id2 = nbf_sublayer2%get_localization()
             if(mainlayer_id.ne.mainlayer_id2) then
                call error_diff_mainlayer_id(
     $               'nbf_interafce_class.f',
     $               'update_links_to_neighbor2_bf_layer',
     $               mainlayer_id,
     $               mainlayer_id2)
             end if
          end if


          !find the coordinates identifying the neighboring buffer layers
          !with which this buffer layer may be sharing grid points
          call get_neighbor_coords(
     $       mainlayer_id,
     $       neighbor1_id,
     $       neighbor2_id,
     $       neighbor_index)


          !get the index corresponding to nbf_sublayer1 -> neighbor1
          nbf_index = get_index(mainlayer_id,1)
          !update links from: nbf_sublayer1 -> neighbor1
          !             to  : nbf_sublayer2 -> neighbor1
          call this%nbf_links(nbf_index)%update_link_in_list(
     $         nbf_sublayer1, nbf_sublayer2)

          !get the index corresponding to neighbor1 -> nbf_sublayer1
          nbf_index = get_index(neighbor1_id, neighbor_index)
          !update links from: neighbor1 -> nbf_sublayer1
          !             to  : neighbor1 -> nbf_sublayer2
          call this%nbf_linfs(nbf_index)%update_link_in_list(
     $         nbf_sublayer1, nbf_sublayer2)

        end subroutine update_links_to_neighbor1_bf_layer


        !< a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers, i.e. this
        !> buffer layer share some gridpoints with other layers
        !> this bf layer has been updated (merge) and pointers to
        !> this buffer layer should be updated as they will be no
        !> longer exist. They are replaced by another link
        subroutine update_links_to_neighbor2_bf_layer(this,
     $     nbf_sublayer1, nbf_sublayer2)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer1
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer2
          

          integer :: mainlayer_id, mainlayer_id2
          integer :: neighbor1_id
          integer :: neighbor2_id
          integer :: neighbor_index
          integer :: nbf_index
          
          !get the cardinal coordinate of the buffer layer
          mainlayer_id = nbf_sublayer1%get_localization()
          if(debug) then
             mainlayer_id2 = nbf_sublayer2%get_localization()
             if(mainlayer_id.ne.mainlayer_id2) then
                call error_diff_mainlayer_id(
     $               'nbf_interafce_class.f',
     $               'update_links_to_neighbor2_bf_layer',
     $               mainlayer_id,
     $               mainlayer_id2)
             end if
          end if


          !find the coordinates identifying the neighboring buffer layers
          !with which this buffer layer may be sharing grid points
          call get_neighbor_coords(
     $       mainlayer_id,
     $       neighbor1_id,
     $       neighbor2_id,
     $       neighbor_index)


          !get the index corresponding to nbf_sublayer1 -> neighbor1
          nbf_index = get_index(mainlayer_id,2)
          !update links from: nbf_sublayer1 -> neighbor1
          !             to  : nbf_sublayer2 -> neighbor1
          call this%nbf_links(nbf_index)%update_link_in_list(
     $         nbf_sublayer1, nbf_sublayer2)

          !get the index corresponding to neighbor1 -> nbf_sublayer1
          nbf_index = get_index(neighbor2_id, neighbor_index)
          !update links from: neighbor1 -> nbf_sublayer1
          !             to  : neighbor1 -> nbf_sublayer2
          call this%nbf_linfs(nbf_index)%update_link_in_list(
     $         nbf_sublayer1, nbf_sublayer2)

        end subroutine update_links_to_neighbor2_bf_layer


        subroutine remove_links_to_neighbor1_bf_layer(this, nbf_sublayer)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer
          

          integer :: mainlayer_id
          integer :: neighbor1_id
          integer :: neighbor2_id
          integer :: neighbor_index
          integer :: nbf_index
          
          !get the cardinal coordinate of the buffer layer
          mainlayer_id = nbf_sublayer%get_localization()


          !find the coordinates identifying the neighboring buffer layers
          !with which this buffer layer may be sharing grid points
          call get_neighbor_coords(
     $       mainlayer_id,
     $       neighbor1_id,
     $       neighbor2_id,
     $       neighbor_index)


          !get the index corresponding to nbf_sublayer1 -> neighbor1
          nbf_index = get_index(mainlayer_id,1)
          !update links from: nbf_sublayer1 -> neighbor1
          !             to  : nbf_sublayer2 -> neighbor1
          call this%nbf_links(nbf_index)%remove_link_from_list(nbf_sublayer)

          !get the index corresponding to neighbor1 -> nbf_sublayer1
          nbf_index = get_index(neighbor1_id, neighbor_index)
          !update links from: neighbor1 -> nbf_sublayer1
          !             to  : neighbor1 -> nbf_sublayer2
          call this%nbf_linfs(nbf_index)%remove_link_from_list(nbf_sublayer)

        end subroutine remove_links_to_neighbor1_bf_layer


        subroutine remove_links_to_neighbor2_bf_layer(this,nbf_sublayer)

          implicit none

          class(nbf_interface)         , intent(inout) :: this
          type(bf_sublayer)   , pointer, intent(in)    :: nbf_sublayer
          

          integer :: mainlayer_id
          integer :: neighbor1_id
          integer :: neighbor2_id
          integer :: neighbor_index
          integer :: nbf_index
          
          !get the cardinal coordinate of the buffer layer
          mainlayer_id = nbf_sublayer%get_localization()


          !find the coordinates identifying the neighboring buffer layers
          !with which this buffer layer may be sharing grid points
          call get_neighbor_coords(
     $       mainlayer_id,
     $       neighbor1_id,
     $       neighbor2_id,
     $       neighbor_index)


          !get the index corresponding to nbf_sublayer -> neighbor1
          nbf_index = get_index(mainlayer_id,2)
          !update links from: nbf_sublayer -> neighbor1
          !             to  : nbf_sublayer -> neighbor1
          call this%nbf_links(nbf_index)%remove_link_from_list(nbf_sublayer)

          !get the index corresponding to neighbor1 -> nbf_sublayer
          nbf_index = get_index(neighbor2_id, neighbor_index)
          !update links from: neighbor1 -> nbf_sublayer
          !             to  : neighbor1 -> nbf_sublayer
          call this%nbf_linfs(nbf_index)%remove_link_from_list(nbf_sublayer)

        end subroutine remove_links_to_neighbor2_bf_layer

      end module nbf_interface_class
