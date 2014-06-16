      module nbf_list_class

        use bf_layer_errors_module, only : error_neighbor_index
        use bf_sublayer_class     , only : bf_sublayer
        use nbf_element_class     , only : nbf_element
        use parameters_kind       , only : ikind
        use sbf_list_class        , only : sbf_list


        implicit none

        private
        public :: nbf_list


        !< double chained list saving ordered references
        !> to buffer sublayers. The references are ordered
        !> with increasing bf_alignment(1,1)
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
          
          procedure, pass, private :: add_element
          procedure, pass, private :: remove_element

          procedure, pass :: print_on_screen

        end type nbf_list


        contains

        !< initialize the chained list
        subroutine ini(this)

          implicit none

          class(nbf_list), intent(inout) :: this
          
          nullify(this%head)
          nullify(this%tail)
          this%nb_elements = 0
          
        end subroutine ini


        function get_head(this)

          implicit none

          class(nbf_list), intent(in) :: this
          type(nbf_element), pointer  :: get_head
          
          get_head => this%head

        end function get_head


        function get_tail(this)

          implicit none

          class(nbf_list), intent(in) :: this
          type(nbf_element), pointer  :: get_tail
          
          get_tail => this%tail

        end function get_tail


        function get_nb_elements(this)

          implicit none

          class(nbf_list), intent(in) :: this
          integer                     :: get_nb_elements
          
          get_nb_elements = this%nb_elements

        end function get_nb_elements


        !< add an element corresponding to the specific link
        !> in the chained list
        subroutine add_link_in_list(this, added_bf_sublayer_ptr)

          implicit none

          class(nbf_list)           , intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: added_bf_sublayer_ptr

          type(nbf_element), pointer :: new_element

          allocate(new_element)
          call new_element%ini(added_bf_sublayer_ptr)

          call this%add_element(new_element)

        end subroutine add_link_in_list


        !< update the element corresponding to a specific link
        !> in the chained list
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


        !< remove the element corresponding to a specific link
        !> in the chained list
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


        !< copy the layers common between the buffer layers
        !> saved in the chained list as neighbors1 to the
        !> buffer layer passed as argument
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


        !< copy the layers common between the buffer layers
        !> saved in the chained list as neighbors to the
        !> buffer layer passed as argument
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


        !< add to the list of sublayer pointers the buffer layers
        !> in the current list that share grid points with the bf layer
        !> given
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


        !< check whether the buffer layer passed as argument has grdpts
        !> in common with the buffe rlayers contained in the list
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


        !< add element in the chained list and its position
        !> is such that the chained list is ordered
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


        !< remove element from the chained list
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


        !< insert element between two elements of the chained list
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


        !< print the alignment of the buffer layer
        !> referenced by the list
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
