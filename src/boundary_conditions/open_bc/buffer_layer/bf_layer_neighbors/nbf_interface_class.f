      !> @file
      !> module implementing the object encapsulating links
      !> to buffer layers at the edge between different main layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the object encapsulating links
      !> to buffer layers at the edge between different main layers
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module nbf_interface_class
      
        use bf_sublayer_class  , only : bf_sublayer
        use nbf_list_class     , only : nbf_list
        use parameters_constant, only : N,S,E,W
        use parameters_input   , only : debug
        use parameters_bf_layer, only : align_N, align_S,
     $                                  align_E, align_W
        use sbf_list_class     , only : sbf_list

        implicit none


        private
        public :: nbf_interface

        
        !>@class nbf_interface
        !> object encapsulting links to buffer layers at the edge between
        !> different main layers
        !
        !>@param nbf_links
        !> array of pointers to the buffer layers at the edge of different
        !> main layers
        !> ex: nbf_links(N,1) : links to the buffer layers considered
        !>                      neighbors of type 1 by the north main
        !>                      layer
        !>     nbf_links(N,2) : links to the buffer layers considered
        !>                      neighbors of type 2 by the north main
        !>                      layer
        !
        !>@param ini
        !> initialize the array of links by initializing the
        !> nbf_list containing the links
        !
        !>@param link_neighbor1_to_bf_sublayer
        !> a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers. Each main
        !> layer sharing grid points with this buffer layer as
        !> neighbor of type 1 is informed that this buffer layer now
        !> exists and should be considered when information are updated
        !> between buffer layers
        !
        !>@param link_neighbor2_to_bf_sublayer
        !> a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers. Each main
        !> layer sharing grid points with this buffer layer as
        !> neighbor of type 2 is informed that this buffer layer now
        !> exists and should be considered when information are updated
        !> between buffer layers
        !
        !>@param update_link_from_neighbor1_to_bf_sublayer
        !> the links that were refering to nbf_sublayer1 are changed into
        !> links to nbf_sublayer2. Only buffer layers considered neighbors
        !> of type 1 by nbf_sublayer1 have their links updated
        !
        !>@param update_link_from_neighbor2_to_bf_sublayer
        !> the links that were refering to nbf_sublayer1 are changed into
        !> links to nbf_sublayer2. Only buffer layers considered neighbors
        !> of type 2 by nbf_sublayer1 have their links updated
        !
        !>@param remove_link_from_neighbor1_to_bf_sublayer
        !> remove the links existing from neighbor1 to nbf_sublayer
        !
        !>@param remove_link_from_neighbor2_to_bf_sublayer
        !> remove the links existing from neighbor2 to nbf_sublayer
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
        !
        !>@param get_nbf_layers_sharing_grdpts_with
        !> add to the list of sublayer pointers the neighboring 
        !> buffer layers that shares grid points in the
        !> x-direction with the current buffer layer
        !
        !>@param bf_layer_depends_on_neighbors
        !> test whether the bf_sublayer is sharing grid points with
        !> its neighboring buffer layers
        !
        !>@param does_a_neighbor_remains
        !> test whether one of the bf_sublayer neighbors is remaining
        !
        !>@param print_on_screen
        !> print the links between bf_sublayers on screen
        !--------------------------------------------------------------
        type nbf_interface

          type(nbf_list), dimension(4,2), private :: nbf_links

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

          procedure, pass :: sync_interface_nodes

          procedure, pass :: get_nbf_layers_sharing_grdpts_with
          procedure, pass :: bf_layer_depends_on_neighbors
          procedure, pass :: does_a_neighbor_remains

          procedure, pass :: print_on_screen

        end type nbf_interface


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the array of links by initializing the
        !> nbf_list containing the links
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers. Each main
        !> layer sharing grid points with this buffer layer as
        !> neighbor of type 1 is informed that this buffer layer now
        !> exists and should be considered when information are updated
        !> between buffer layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !--------------------------------------------------------------        
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> a buffer layer has been identified as a buffer layer
        !> located at the interface between main layers. Each main
        !> layer sharing grid points with this buffer layer as
        !> neighbor of type 2 is informed that this buffer layer now
        !> exists and should be considered when information are updated
        !> between buffer layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !--------------------------------------------------------------        
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> the links that were refering to nbf_sublayer1 are changed into
        !> links to nbf_sublayer2. Only buffer layers considered neighbors
        !> of type 1 by nbf_sublayer1 have their links updated
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_sublayer1
        !> bf_sublayer reference before the update
        !
        !>@param nbf_sublayer2
        !> bf_sublayer reference after the update
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> the links that were refering to nbf_sublayer1 are changed into
        !> links to nbf_sublayer2. Only buffer layers considered neighbors
        !> of type 2 by nbf_sublayer1 have their links updated
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_sublayer1
        !> bf_sublayer reference before the update
        !
        !>@param nbf_sublayer2
        !> bf_sublayer reference after the update
        !--------------------------------------------------------------
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
      

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the links existing from neighbor1 to nbf_sublayer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_sublayer
        !> bf_sublayer reference removed
        !--------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the links existing from neighbor2 to nbf_sublayer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_sublayer
        !> bf_sublayer reference removed
        !--------------------------------------------------------------
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

          class(nbf_interface), intent(in)    :: this
          type(bf_sublayer)   , intent(inout) :: nbf_sublayer

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

          class(nbf_interface), intent(inout) :: this
          type(bf_sublayer)   , intent(in)    :: nbf_sublayer

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

          class(nbf_interface), intent(inout) :: this


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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add to the list of sublayer pointers the neighboring 
        !> buffer layers that shares grid points in the
        !> x-direction with the current buffer layer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_type
        !> type of the neighboring bf_sublayer investigated
        !
        !>@param bf_sublayer_i
        !> reference to the bf_sublayer whose neighbors are investigated
        !
        !>@param bf_sublayer_list
        !> list of the bf_sublayer objects that share grid points in the
        !> x-direction with the current buffer layer
        !
        !>@param bf_mainlayer_id
        !> cardinal coordinate of the buffer layer investigated
        !--------------------------------------------------------------
        subroutine get_nbf_layers_sharing_grdpts_with(
     $     this,
     $     nbf_type,
     $     bf_sublayer_i,
     $     bf_sublayer_list,
     $     bf_mainlayer_id)

          implicit none

          class(nbf_interface)      , intent(in)    :: this
          integer                   , intent(in)    :: nbf_type
          type(bf_sublayer), pointer, intent(in)    :: bf_sublayer_i
          type(sbf_list)            , intent(inout) :: bf_sublayer_list
          integer         , optional, intent(in)    :: bf_mainlayer_id


          integer :: mainlayer_id


          if(present(bf_mainlayer_id)) then
             mainlayer_id = bf_mainlayer_id
          else
             mainlayer_id = bf_sublayer_i%get_localization()
          end if

          call this%nbf_links(mainlayer_id,nbf_type)%get_nbf_layers_sharing_grdpts_with(
     $         bf_sublayer_i, bf_sublayer_list)

        end subroutine get_nbf_layers_sharing_grdpts_with


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> test whether the bf_sublayer is sharing grid points with
        !> its neighboring buffer layers
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_type
        !> type of the neighboring bf_sublayer investigated
        !
        !>@param bf_mainlayer_id
        !> cardinal coordinate of the buffer layer investigated
        !
        !>@param dependent
        !> logical stating whether the buffer layer is sharing
        !> grid points with the neighboring buffer layers
        !--------------------------------------------------------------
        function bf_layer_depends_on_neighbors(
     $     this, nbf_type, bf_sublayer_i, bf_mainlayer_id)
     $     result(dependent)

          implicit none

          class(nbf_interface)      , intent(in) :: this
          integer                   , intent(in) :: nbf_type
          type(bf_sublayer), pointer, intent(in) :: bf_sublayer_i
          integer         , optional, intent(in) :: bf_mainlayer_id
          logical                                :: dependent

          integer :: mainlayer_id


          if(present(bf_mainlayer_id)) then
             mainlayer_id = bf_mainlayer_id
          else
             mainlayer_id = bf_sublayer_i%get_localization()
          end if

          dependent = this%nbf_links(mainlayer_id,nbf_type)%bf_layer_depends_on_neighbors(
     $         bf_sublayer_i)

        end function bf_layer_depends_on_neighbors


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> test whether one of the bf_sublayer neighbors is remaining
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !
        !>@param nbf_type
        !> type of the neighboring bf_sublayer investigated
        !
        !>@param bf_sublayer_id
        !> bf_sublayer 
        !
        !>@param bf_mainlayer_id
        !> cardinal coordinate of the buffer layer investigated
        !
        !>@param a_neighbor_remains
        !> logical stating whether the buffer layer cannot be removed
        !> because a neighboring buffer layer should remain
        !--------------------------------------------------------------
        function does_a_neighbor_remains(
     $     this, nbf_type, bf_sublayer_i, bf_mainlayer_id)
     $     result(a_neighbor_remains)

          implicit none

          class(nbf_interface)      , intent(in)    :: this
          integer                   , intent(in)    :: nbf_type
          type(bf_sublayer), pointer, intent(in)    :: bf_sublayer_i
          integer         , optional, intent(in)    :: bf_mainlayer_id
          logical                                   :: a_neighbor_remains


          integer :: mainlayer_id


          if(present(bf_mainlayer_id)) then
             mainlayer_id = bf_mainlayer_id
          else
             mainlayer_id = bf_sublayer_i%get_localization()
          end if

          a_neighbor_remains = this%nbf_links(mainlayer_id,nbf_type)%does_a_neighbor_remains(
     $         bf_sublayer_i)

        end function does_a_neighbor_remains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the links between bf_sublayers on screen
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_interface object encapsulting links to buffer
        !> layers at the edge between different main layers
        !--------------------------------------------------------------
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
