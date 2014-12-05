      !> @file
      !> module implementing a doubled chained list
      !> of bf_sublayer pointers
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing a doubled chained list
      !> of bf_sublayer pointers
      !
      !> @date
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module nbf_list_class
      
        use bf_interior_bc_sections_module, only :
     $       determine_interior_bc_sections
      
        use bf_layer_errors_module, only :
     $       error_neighbor_index

        use bf_sublayer_class, only :
     $       bf_sublayer

        use nbf_element_class, only :
     $       nbf_element

        use parameters_input, only :
     $       ne,bc_size

        use parameters_kind, only :
     $       ikind,rkind

        use sbf_list_class, only :
     $       sbf_list

        implicit none

        private
        public :: nbf_list


        !>@class nbf_list
        !> double chained list saving ordered references
        !> to bf_sublayer objects. The references are ordered
        !> with increasing bf_alignment(1,1)
        !
        !>@param head
        !> pointer to the first bf_sublayer reference of the list
        !
        !>@param tail
        !> pointer to the last bf_sublayer reference of the list
        !
        !>@param nb_elements
        !> number of references to bf_sublayer stored in the list
        !        
        !>@param ini
        !> initialize the list by initializing the
        !> number of element to 0 and nullifying the head
        !> and tail attributes
        !
        !>@param get_head
        !> get the head attribute
        !
        !>@param get_tail
        !> get the tail attribute
        !
        !>@param get_nb_elements
        !> get the nb_elements attribute
        !
        !>@param add_link_in_list
        !> add an element corresponding to the specific link
        !> in the chained list
        !
        !>@param update_link_in_list
        !> update the element linked to bf_sublayer1
        !> by modifying its reference to bf_sublayer2
        !
        !>@param remove_link_from_list
        !> remove the element corresponding to a specific link
        !> from the chained list
        !
        !>@param copy_from_neighbors_to_bf_layer
        !> copy the common layers from the buffer layers
        !> saved in the chained list as neighbors to the
        !> bf_layer object passed as argument
        !
        !>@param copy_to_neighbors_from_bf_layer
        !> copy the common layers to the buffer layers
        !> saved in the chained list as neighbors
        !> from the bf_layer object passed as argument
        !
        !>@param copy_from_neighbors1_to_bf_layer
        !> copy the common layers from the buffer layers
        !> saved in the chained list as neighbors of type 1
        !> to the bf_layer object passed as argument
        !
        !>@param copy_from_neighbors2_to_bf_layer
        !> copy the common layers from the buffer layers
        !> saved in the chained list as neighbors of type 2
        !> to the bf_layer object passed as argument
        !
        !>@param copy_to_neighbors1_from_bf_layer
        !> copy the common layers to the buffer layers
        !> saved in the chained list as neighbors of type 1
        !> from the bf_layer object passed as argument
        !
        !>@param copy_to_neighbors2_from_bf_layer
        !> copy the common layers to the buffer layers
        !> saved in the chained list as neighbors of type 2
        !> from the bf_layer object passed as argument
        !
        !>@param get_nbf_layers_sharing_grdpts_with
        !> add to the list of sublayer pointers the buffer layers
        !> in the current list that share grid points with the bf layer
        !> given
        !
        !>@param bf_layer_depends_on_neighbors
        !> check whether the buffer layer passed as argument has grdpts
        !> in common with the buffer layers contained in the list
        !
        !>@param does_a_neighbor_remains
        !> check whether a bf_sublayer can be removed considering the
        !> neighboring buffer layers and if their removal has been
        !> confirmed
        !
        !>@param sync_nodes_with_neighbor1
        !> synchronize the nodes at the interface between buffer main
        !> layer of type neighbor1
        !
        !>@param sync_nodes_with_neighbor2
        !> synchronize the nodes at the interface between buffer main
        !> layer of type neighbor2
        !
        !>@param update_bc_sections
        !> update the bc_sections of the buffer layer by checking the
        !> grid points shared with the neighboring buffer layers
        !
        !>@param get_data_for_newgrdpt
        !> get the data needed for the computation of the new grid point
        !
        !>@param get_grdpts_id_part
        !> get the grid points needed to determine whether a suspicious
        !> bc_interior_pt should be updated to the interior_pt status
        !> at t
        !
        !>@param add_element
        !> add an element in the chained list ensuring that
        !> the doubled chained list element are ordered
        !
        !>@param remove_element
        !> remove an element from the chained list
        !
        !>@param print_on_screen
        !> print the alignment of the buffer layer
        !> referenced in the chained list
        !---------------------------------------------------------------
        type :: nbf_list

          type(nbf_element), pointer, private :: head
          type(nbf_element), pointer, private :: tail
          integer                   , private :: nb_elements

          contains

          procedure, pass :: ini

          procedure, pass :: get_head
          procedure, pass :: get_tail
          procedure, pass :: get_nb_elements

          procedure, pass :: add_link_in_list
          procedure, pass :: update_link_in_list
          procedure, pass :: remove_link_from_list

          procedure, pass :: copy_from_neighbors_to_bf_layer
          procedure, pass :: copy_to_neighbors_from_bf_layer
          procedure, pass :: copy_from_neighbors1_to_bf_layer
          procedure, pass :: copy_from_neighbors2_to_bf_layer
          procedure, pass :: copy_to_neighbors1_from_bf_layer
          procedure, pass :: copy_to_neighbors2_from_bf_layer

          procedure, pass :: get_nbf_layers_sharing_grdpts_with
          procedure, pass :: bf_layer_depends_on_neighbors
          procedure, pass :: does_a_neighbor_remains

          procedure, pass :: sync_nodes_with_neighbor1
          procedure, pass :: sync_nodes_with_neighbor2
          procedure, pass :: update_bc_sections

          procedure, pass :: get_data_for_newgrdpt
          procedure, pass :: get_grdpts_id_part

          procedure, pass, private :: add_element
          procedure, pass, private :: remove_element

          procedure, pass :: print_on_screen

        end type nbf_list

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the list by initializing the
        !> number of element to 0 and nullifying the head
        !> and tail attributes
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(nbf_list), intent(inout) :: this
          
          nullify(this%head)
          nullify(this%tail)
          this%nb_elements = 0
          
        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the head attribute
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param get_head
        !> reference to the first element of the doubled chained
        !> list
        !--------------------------------------------------------------
        function get_head(this)

          implicit none

          class(nbf_list), intent(in) :: this
          type(nbf_element), pointer  :: get_head
          
          get_head => this%head

        end function get_head


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the tail attribute
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param get_head
        !> reference to the last element of the doubled chained
        !> list
        !--------------------------------------------------------------
        function get_tail(this)

          implicit none

          class(nbf_list), intent(in) :: this
          type(nbf_element), pointer  :: get_tail
          
          get_tail => this%tail

        end function get_tail


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the nb_elements attribute
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param get_nb_elements
        !> nb_element attribute
        !--------------------------------------------------------------
        function get_nb_elements(this)

          implicit none

          class(nbf_list), intent(in) :: this
          integer                     :: get_nb_elements
          
          get_nb_elements = this%nb_elements

        end function get_nb_elements


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add an element corresponding to the specific link
        !> in the chained list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param added_bf_sublayer_ptr
        !> reference to the bf_sublayer added in the doubled chained
        !> list
        !--------------------------------------------------------------
        subroutine add_link_in_list(this, added_bf_sublayer_ptr)

          implicit none

          class(nbf_list)           , intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: added_bf_sublayer_ptr

          type(nbf_element), pointer :: new_element

          allocate(new_element)
          call new_element%ini(added_bf_sublayer_ptr)

          call this%add_element(new_element)

        end subroutine add_link_in_list


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the element linked to bf_sublayer1
        !> by modifying its reference to bf_sublayer2
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param bf_sublayer1
        !> bf_sublayer reference to be updated
        !
        !>@param bf_sublayer2
        !> bf_sublayer reference after update
        !--------------------------------------------------------------
        subroutine update_link_in_list(this, bf_sublayer1, bf_sublayer2)

          implicit none
          
          class(nbf_list)           , intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: bf_sublayer1
          type(bf_sublayer), pointer, intent(in)    :: bf_sublayer2

          type(nbf_element), pointer :: current_element
          integer :: i

          if(associated(this%head)) then
             current_element => this%head

             do i=1, this%nb_elements
                if(current_element%refers_to(bf_sublayer1)) then
                   call current_element%set_ptr(bf_sublayer2)
                   exit
                end if
                current_element => current_element%get_next()
             end do

          end if

        end subroutine update_link_in_list


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the element corresponding to a specific link
        !> from the chained list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param rm_bf_sublayer_ptr
        !> bf_sublayer reference identifying the element removed
        !--------------------------------------------------------------
        subroutine remove_link_from_list(this, rm_bf_sublayer_ptr)

          implicit none
          
          class(nbf_list)           , intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: rm_bf_sublayer_ptr

          type(nbf_element), pointer :: current_element
          integer :: i

          if(associated(this%head)) then
             current_element => this%head

             do i=1, this%nb_elements
                if(current_element%refers_to(rm_bf_sublayer_ptr)) then
                   call this%remove_element(current_element)
                   exit
                end if
                current_element => current_element%get_next()
             end do

          end if

        end subroutine remove_link_from_list


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layers from the buffer layers
        !> saved in the chained list as neighbors to the
        !> bf_layer object passed as argument
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param neighbor_index
        !> index identifying the type of neighbors saved in the
        !> current doubled chained list
        !
        !>@param bf_exchanged
        !> bf_layer object exchanging data with the elements
        !> of the doubled chained list
        !--------------------------------------------------------------
        subroutine copy_from_neighbors_to_bf_layer(
     $     this, neighbor_index, bf_exchanged)

          implicit none

          class(nbf_list)  , intent(in)    :: this
          integer          , intent(in)    :: neighbor_index
          type(bf_sublayer), intent(inout) :: bf_exchanged

          
          select case(neighbor_index)
            case(1)
               call copy_from_neighbors1_to_bf_layer(
     $              this, bf_exchanged)
            case(2)
               call copy_from_neighbors2_to_bf_layer(
     $              this, bf_exchanged)
            case default
               call error_neighbor_index(
     $              'nbf_list_class.f',
     $              'copy_from_neighbors_to_bf_layer',
     $              neighbor_index)
          end select

        end subroutine copy_from_neighbors_to_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layers from the buffer layers
        !> saved in the chained list as neighbors of type 1
        !> to the bf_layer object passed as argument
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param bf_exchanged
        !> bf_layer object exchanging data with the elements
        !> of the doubled chained list
        !--------------------------------------------------------------
        subroutine copy_from_neighbors1_to_bf_layer(
     $     this, bf_exchanged)

          implicit none

          class(nbf_list)  , intent(in)    :: this
          type(bf_sublayer), intent(inout) :: bf_exchanged

          type(nbf_element), pointer :: current_element
          integer :: i

          if(bf_exchanged%can_exchange_with_neighbor1()) then
             if(associated(this%head)) then
                current_element => this%head
                
                do i=1, this%nb_elements
                   call current_element%copy_from_neighbor1_to(bf_exchanged)
                   current_element => current_element%get_next()
                end do
             end if
          end if
             
        end subroutine copy_from_neighbors1_to_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layers from the buffer layers
        !> saved in the chained list as neighbors of type 2
        !> to the bf_layer object passed as argument
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param bf_exchanged
        !> bf_layer object exchanging data with the elements
        !> of the doubled chained list
        !--------------------------------------------------------------
        subroutine copy_from_neighbors2_to_bf_layer(
     $     this, bf_exchanged)

          implicit none

          class(nbf_list)  , intent(in)    :: this
          type(bf_sublayer), intent(inout) :: bf_exchanged

          type(nbf_element), pointer :: current_element
          integer :: i


          if(bf_exchanged%can_exchange_with_neighbor2()) then
             if(associated(this%head)) then
                current_element => this%head

                do i=1, this%nb_elements
                   call current_element%copy_from_neighbor2_to(bf_exchanged)
                   current_element => current_element%get_next()
                end do
             end if
          end if

        end subroutine copy_from_neighbors2_to_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layers to the buffer layers
        !> saved in the chained list as neighbors
        !> from the bf_layer object passed as argument
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param neighbor_index
        !> index identifying the type of neighbors saved in the
        !> current doubled chained list
        !
        !>@param bf_exchanged
        !> bf_layer object exchanging data with the elements
        !> of the doubled chained list
        !--------------------------------------------------------------
        subroutine copy_to_neighbors_from_bf_layer(
     $     this, neighbor_index, bf_exchanged)

          implicit none

          class(nbf_list)  , intent(inout) :: this
          integer          , intent(in)    :: neighbor_index
          type(bf_sublayer), intent(in)    :: bf_exchanged


          select case(neighbor_index)
            case(1)
               call copy_to_neighbors1_from_bf_layer(
     $              this, bf_exchanged)
            case(2)
               call copy_to_neighbors2_from_bf_layer(
     $              this, bf_exchanged)
            case default
               call error_neighbor_index(
     $              'nbf_list_class.f',
     $              'copy_to_neighbors_from_bf_layer',
     $              neighbor_index)
          end select

        end subroutine copy_to_neighbors_from_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layers to the buffer layers
        !> saved in the chained list as neighbors of type 1
        !> from the bf_layer object passed as argument
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param bf_exchanged
        !> bf_layer object exchanging data with the elements
        !> of the doubled chained list
        !--------------------------------------------------------------
        subroutine copy_to_neighbors1_from_bf_layer(
     $     this, bf_exchanged)

          implicit none

          class(nbf_list)  , intent(inout) :: this
          type(bf_sublayer), intent(in)    :: bf_exchanged

          type(nbf_element), pointer :: current_element
          integer :: i

          if(bf_exchanged%can_exchange_with_neighbor1()) then

             if(associated(this%head)) then
                
                current_element => this%head
                
                do i=1, this%nb_elements
                   call current_element%copy_to_neighbor1_from(bf_exchanged)
                   current_element => current_element%get_next()
                end do

             end if
          end if

        end subroutine copy_to_neighbors1_from_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> copy the common layers to the buffer layers
        !> saved in the chained list as neighbors of type 2
        !> from the bf_layer object passed as argument
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param bf_exchanged
        !> bf_layer object exchanging data with the elements
        !> of the doubled chained list
        !--------------------------------------------------------------
        subroutine copy_to_neighbors2_from_bf_layer(
     $     this, bf_exchanged)

          implicit none

          class(nbf_list)  , intent(inout) :: this
          type(bf_sublayer), intent(in)    :: bf_exchanged

          type(nbf_element), pointer :: current_element
          integer :: i

          if(bf_exchanged%can_exchange_with_neighbor2()) then

             if(associated(this%head)) then
                
                current_element => this%head
                
                do i=1, this%nb_elements
                   call current_element%copy_to_neighbor2_from(bf_exchanged)
                   current_element => current_element%get_next()
                end do
                
             end if
          end if

        end subroutine copy_to_neighbors2_from_bf_layer


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add to the list of sublayer pointers the buffer layers
        !> in the current list that share grid points with the bf layer
        !> given
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param bf_sublayer_i
        !> bf_layer object that may be exchanging data with the elements
        !> of the double chained list
        !
        !>@param bf_sublayer_list
        !> list of bf_layer object from the current doubled chained list
        !> that are exchanging data with bf_sublayer_i
        !--------------------------------------------------------------
        subroutine get_nbf_layers_sharing_grdpts_with(
     $     this, bf_sublayer_i, bf_sublayer_list)

          implicit none

          class(nbf_list)           , intent(in)    :: this
          type(bf_sublayer), pointer, intent(in)    :: bf_sublayer_i
          type(sbf_list)            , intent(inout) :: bf_sublayer_list

          type(nbf_element), pointer :: current_element
          integer                    :: k

          current_element => this%head

          do k=1, this%nb_elements
             if(current_element%shares_grdpts_along_x_dir_with(bf_sublayer_i)) then
                call bf_sublayer_list%add_ele(current_element%get_ptr())
             end if
             current_element => current_element%get_next()
          end do

        end subroutine get_nbf_layers_sharing_grdpts_with


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the buffer layer passed as argument has grdpts
        !> in common with the buffer layers contained in the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param bf_sublayer_i
        !> bf_layer object that may be exchanging data with the elements
        !> of the double chained list
        !
        !>@return dependent
        !> logical stating whether the bf_sublayer object passed as
        !> argument has grid points in common with the elements of the
        !> doubled chained list
        !--------------------------------------------------------------
        function bf_layer_depends_on_neighbors(
     $     this, bf_sublayer_i) result(dependent)

          implicit none

          class(nbf_list)           , intent(in)    :: this
          type(bf_sublayer), pointer, intent(in)    :: bf_sublayer_i
          logical                                   :: dependent

          type(nbf_element), pointer :: current_element
          integer                    :: k

          dependent = .false.

          current_element => this%head

          do k=1, this%nb_elements
             if(current_element%shares_grdpts_along_x_dir_with(bf_sublayer_i)) then
                dependent = .true.
                exit
             end if
             current_element => current_element%get_next()
          end do

        end function bf_layer_depends_on_neighbors


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether a bf_sublayer can be removed considering the
        !> neighboring buffer layers and if their removal has been
        !> confirmed
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param bf_sublayer_i
        !> bf_layer object that may be exchanging data with the elements
        !> of the double chained list
        !
        !>@return a_neighbor_remains
        !> logical stating whether a neighboring bufer layer should remain
        !--------------------------------------------------------------
        function does_a_neighbor_remains(this, bf_sublayer_i) result(a_neighbor_remains)

          implicit none

          class(nbf_list)           , intent(in)    :: this
          type(bf_sublayer), pointer, intent(in)    :: bf_sublayer_i
          logical                                   :: a_neighbor_remains

          type(nbf_element), pointer :: current_element
          integer                    :: k

          a_neighbor_remains = .false.

          current_element => this%head

          do k=1, this%nb_elements
             if(current_element%shares_grdpts_along_x_dir_with(bf_sublayer_i)) then
                if(current_element%get_remain_status()) then
                   a_neighbor_remains = .true.
                   exit
                end if
             end if
             current_element => current_element%get_next()
          end do

        end function does_a_neighbor_remains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> synchronize nodes at the interface between two buffer main
        !> layers, the second main layer is a neighbor of type 1 for
        !> the first buffer layer
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param nbf_list_interface
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !--------------------------------------------------------------
        subroutine sync_nodes_with_neighbor1(this, this2)

          implicit none

          class(nbf_list), intent(inout) :: this
          class(nbf_list), intent(inout) :: this2

          type(nbf_element), pointer :: nbf_ele1
          type(nbf_element), pointer :: nbf_ele2
          integer                    :: k1, k2

          nbf_ele1 => this%head
          do k1=1, this%nb_elements

             nbf_ele2 => this2%head
             do k2=1, this2%nb_elements

                call nbf_ele1%sync_nodes_with_neighbor1(nbf_ele2)

                nbf_ele2 => nbf_ele2%get_next()
             end do

             nbf_ele1 => nbf_ele1%get_next()
          end do

        end subroutine sync_nodes_with_neighbor1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> synchronize nodes at the interface between two buffer main
        !> layers, the second main layer is a neighbor of type 2 for
        !> the first buffer layer
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param nbf_list_interface
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !--------------------------------------------------------------
        subroutine sync_nodes_with_neighbor2(this, this2)

          implicit none

          class(nbf_list), intent(inout) :: this
          class(nbf_list), intent(inout) :: this2

          type(nbf_element), pointer :: nbf_ele1
          type(nbf_element), pointer :: nbf_ele2
          integer                    :: k1, k2

          nbf_ele1 => this%head
          do k1=1, this%nb_elements

             nbf_ele2 => this2%head
             do k2=1, this2%nb_elements

                call nbf_ele1%sync_nodes_with_neighbor2(nbf_ele2)

                nbf_ele2 => nbf_ele2%get_next()
             end do

             nbf_ele1 => nbf_ele1%get_next()
          end do

        end subroutine sync_nodes_with_neighbor2

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the bc_sections of the buffer layer by
        !> checking the grid points shared with the
        !> neighboring buffer layers
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param nbf_list_interface
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !--------------------------------------------------------------
        subroutine update_bc_sections(
     $     this,
     $     interior_inf,
     $     interior_sup,
     $     nb_bc_sections,
     $     bc_sections,
     $     min_initialized,
     $     max_initialized,
     $     no_bf_common_with_bf_layer)
        
          implicit none

          class(nbf_list)                            , intent(in)    :: this
          integer(ikind)                             , intent(in)    :: interior_inf
          integer(ikind)                             , intent(in)    :: interior_sup
          integer(ikind)                             , intent(inout) :: nb_bc_sections
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: bc_sections
          logical                                    , intent(inout) :: min_initialized
          logical                                    , intent(inout) :: max_initialized
          logical                                    , intent(inout) :: no_bf_common_with_bf_layer

          type(nbf_element), pointer   :: current_element
          type(bf_sublayer), pointer   :: bf_sublayer_ptr
          integer(ikind), dimension(2) :: bf_alignment
          integer                      :: i

          current_element => this%get_head()

          do i=1, this%get_nb_elements()

             bf_sublayer_ptr => current_element%get_ptr()
             bf_alignment(1) = bf_sublayer_ptr%get_alignment(1,1)
             bf_alignment(2) = bf_sublayer_ptr%get_alignment(1,2)

             !update the bc_sections
             call determine_interior_bc_sections(
     $            bf_alignment,
     $            interior_inf,
     $            interior_sup,
     $            nb_bc_sections,
     $            bc_sections,
     $            min_initialized,
     $            max_initialized,
     $            no_bf_common_with_bf_layer)

             current_element => current_element%get_next()
          end do
          

        end subroutine update_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the bc_sections of the buffer layer by
        !> checking the grid points shared with the
        !> neighboring buffer layers
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param nbf_list_interface
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !--------------------------------------------------------------
        subroutine get_data_for_newgrdpt(
     $     this,
     $     tmp_grdpts_id0,
     $     tmp_nodes0,
     $     tmp_nodes1,
     $     gen_borders)
        
          implicit none

          class(nbf_list)                                               , intent(in)    :: this
          integer        , dimension(2*(bc_size+1)+1,2*(bc_size+1)+1)   , intent(inout) :: tmp_grdpts_id0
          real(rkind)    , dimension(2*(bc_size+1)+1,2*(bc_size+1)+1,ne), intent(inout) :: tmp_nodes0
          real(rkind)    , dimension(2*(bc_size+1)+1,2*(bc_size+1)+1,ne), intent(inout) :: tmp_nodes1
          integer(ikind) , dimension(2,2)                               , intent(in)    :: gen_borders


          type(nbf_element), pointer   :: current_element
          type(bf_sublayer), pointer   :: bf_sublayer_ptr
          integer                      :: i

          current_element => this%get_head()

          do i=1, this%get_nb_elements()

             bf_sublayer_ptr => current_element%get_ptr()

             !get the data needed for the new grdpt
             call bf_sublayer_ptr%get_data_for_newgrdpt(
     $            tmp_grdpts_id0,
     $            tmp_nodes0,
     $            tmp_nodes1,
     $            gen_borders)

             current_element => current_element%get_next()
          end do


        end subroutine get_data_for_newgrdpt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> update the bc_sections of the buffer layer by
        !> checking the grid points shared with the
        !> neighboring buffer layers
        !
        !> @date
        !> 30_10_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param tmp_grdpts_id1
        !> temporary array where the identity of the grid points
        !> needed to determine whether the suspicious bc_interior_pt
        !> should be turned into an interior_pt are stored at t
        !
        !>@param gen_borders
        !> general coordinates identifying the SW and NE borders of the
        !> grid-points asked
        !--------------------------------------------------------------
        subroutine get_grdpts_id_part(
     $     this,
     $     tmp_grdpts_id1,
     $     gen_borders)
        
          implicit none

          class(nbf_list)                                            , intent(in)    :: this
          integer        , dimension(2*(bc_size+1)+1,2*(bc_size+1)+1), intent(inout) :: tmp_grdpts_id1
          integer(ikind) , dimension(2,2)                            , intent(in)    :: gen_borders


          type(nbf_element), pointer   :: current_element
          type(bf_sublayer), pointer   :: bf_sublayer_ptr
          integer                      :: i

          current_element => this%get_head()

          do i=1, this%get_nb_elements()

             bf_sublayer_ptr => current_element%get_ptr()

             !get the data needed for the new grdpt
             call bf_sublayer_ptr%get_grdpts_id_part(
     $            tmp_grdpts_id1,
     $            gen_borders)

             current_element => current_element%get_next()

          end do

        end subroutine get_grdpts_id_part


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add an element in the chained list ensuring that
        !> the doubled chained list element are ordered
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param new_element
        !> nbf_element added in the list
        !--------------------------------------------------------------
        subroutine add_element(this, new_element)

          implicit none

          class(nbf_list)           , intent(inout) :: this
          type(nbf_element), pointer, intent(inout) :: new_element

          type(nbf_element), pointer :: current_element
          type(nbf_element), pointer :: current_element_prev
          integer :: i          

          select case(this%nb_elements)
            case(0)
               this%head => new_element
               this%tail => this%head

            case(1)               
               if(this%head%is_before(new_element)) then
                  call this%head%set_next(new_element)
                  call new_element%set_prev(this%head)
                  this%tail => new_element

               else                
                  call new_element%set_next(this%head)
                  call this%head%set_prev(new_element)
                  call this%head%nullify_next()
                  this%head => new_element

               end if
                
            case default
             
               current_element => this%head

               do i=1, this%nb_elements-1
                  
                  if(current_element%is_before(new_element)) then
                     current_element => current_element%get_next()
                  else
                     exit
                  end if

               end do


               if(i.ne.this%nb_elements) then
                  if(associated(current_element%get_prev())) then

                     current_element_prev => current_element%get_prev()
                     call insert_element(
     $                    new_element,
     $                    current_element_prev,
     $                    current_element)

                  else

                     this%head=> new_element
                     call new_element%set_next(current_element)
                     call current_element%set_prev(new_element)

                  end if

               else

                  if(current_element%is_before(new_element)) then

                     call current_element%set_next(new_element)
                     call new_element%set_prev(current_element)
                     this%tail => new_element

                  else

                     current_element_prev => current_element%get_prev()
                     call insert_element(
     $                    new_element,
     $                    current_element_prev,
     $                    current_element)

                  end if
               end if

          end select

          !update the number of elements in the chained list
          this%nb_elements = this%nb_elements+1

        end subroutine add_element


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove an element from the chained list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param rm_element
        !> nbf_element removed from the list
        !--------------------------------------------------------------
        subroutine remove_element(this, rm_element)

          implicit none

          class(nbf_list)           , intent(inout) :: this
          type(nbf_element), pointer, intent(inout) :: rm_element

          type(nbf_element), pointer :: prev
          type(nbf_element), pointer :: next
          

          !reconnect the list elements
          if(associated(rm_element%get_next())) then
             next => rm_element%get_next()

             if(associated(rm_element%get_prev())) then
                prev => rm_element%get_prev()
                call prev%set_next(next)
                call next%set_prev(prev)

             else
                call next%nullify_prev()
                this%head => next

             end if
          else

             if(associated(rm_element%get_prev())) then
                prev => rm_element%get_prev()
                call prev%nullify_next()
                this%tail => prev

             else
                nullify(this%tail)
                nullify(this%head)

             end if
          end if


          !destroy element
          call rm_element%nullify_prev()
          call rm_element%nullify_next()
          call rm_element%nullify_ptr()
          deallocate(rm_element)


          !there is one less element in the list
          this%nb_elements = this%nb_elements - 1

        end subroutine remove_element


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> insert element between two elements of the chained list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !
        !>@param prev_element
        !> element before the element inserted
        !
        !>@param next_element
        !> element after the element inserted
        !--------------------------------------------------------------
        subroutine insert_element(
     $     inserted_element,
     $     prev_element,
     $     next_element)

          implicit none

          type(nbf_element), pointer, intent(inout) :: inserted_element
          type(nbf_element), pointer, intent(inout) :: prev_element
          type(nbf_element), pointer, intent(inout) :: next_element

          call prev_element%set_next(inserted_element)
          call inserted_element%set_prev(prev_element)

          call inserted_element%set_next(next_element)
          call next_element%set_prev(inserted_element)

        end subroutine insert_element


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the alignment of the buffer layer
        !> referenced in the chained list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> nbf_list object implementing a doubled chained
        !> list of bf_sublayer references
        !--------------------------------------------------------------
        subroutine print_on_screen(this)

          implicit none

          class(nbf_list), intent(in) :: this

          type(nbf_element), pointer :: current_element
          type(bf_sublayer), pointer :: bf_sublayer_ptr
          integer(ikind), dimension(2,2) :: alignment
          character(1)  , dimension(4)   :: bf_layer_char

          integer :: i

          bf_layer_char = ['N','S','E','W']

          current_element => this%get_head()

          do i=1, this%get_nb_elements()
             bf_sublayer_ptr => current_element%get_ptr()
             alignment       = bf_sublayer_ptr%get_alignment_tab()
             
             print '(A1,'': ('',I3,I3,'')  ('',I3,I3,'')'')',
     $            bf_layer_char(bf_sublayer_ptr%get_localization()),
     $            alignment(1,1), alignment(1,2),
     $            alignment(2,1), alignment(2,2)

             current_element => current_element%get_next()
          end do
          print '()'

        end subroutine print_on_screen

      end module nbf_list_class
