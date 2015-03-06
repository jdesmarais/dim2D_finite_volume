      !> @file
      !> object encapsulating links to buffer layers at the
      !> edge between different main layers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> object encapsulating links to buffer layers at the
      !> edge between different main layers
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module nbf_interface_basic_class

        use bf_sublayer_class, only :
     $       bf_sublayer

        use nbf_element_class, only :
     $       nbf_element

        use nbf_list_class, only :
     $       nbf_list

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none


        private
        public :: nbf_interface_basic

        
        !>@class nbf_interface_basic
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
        !>@param print_on_screen
        !> print the links between bf_sublayers on screen
        !--------------------------------------------------------------
        type nbf_interface_basic

          type(nbf_list), dimension(4,2) :: nbf_links

          contains

          procedure, pass :: ini

          procedure, pass :: link_neighbor1_to_bf_sublayer
          procedure, pass :: link_neighbor2_to_bf_sublayer
          procedure, pass :: update_link_from_neighbor1_to_bf_sublayer
          procedure, pass :: update_link_from_neighbor2_to_bf_sublayer
          procedure, pass :: remove_link_from_neighbor1_to_bf_sublayer
          procedure, pass :: remove_link_from_neighbor2_to_bf_sublayer

          procedure, pass :: print_on_screen

        end type nbf_interface_basic


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

          class(nbf_interface_basic), intent(inout) :: this

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

          class(nbf_interface_basic), intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: nbf_sublayer

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

          class(nbf_interface_basic), intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: nbf_sublayer

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

          class(nbf_interface_basic), intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: nbf_sublayer1
          type(bf_sublayer), pointer, intent(in)    :: nbf_sublayer2
          

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

          class(nbf_interface_basic), intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: nbf_sublayer1
          type(bf_sublayer), pointer, intent(in)    :: nbf_sublayer2
          

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

          class(nbf_interface_basic), intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: nbf_sublayer
          

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

          class(nbf_interface_basic), intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: nbf_sublayer
          

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

          class(nbf_interface_basic), intent(in) :: this

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

      end module nbf_interface_basic_class
