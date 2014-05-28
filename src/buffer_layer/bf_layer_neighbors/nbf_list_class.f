      module nbf_list_class

        use bf_sublayer_class, only : bf_sublayer
        use nbf_element_class, only : nbf_element

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

          procedure, pass :: add_link_in_list
          procedure, pass :: udpate_link_in_list
          procedure, pass :: remove_link_from_list

          procedure, pass :: copy_from_neighbors_to_bf_layer
          procedure, pass :: copy_to_neighbors_from_bf_layer

          procedure, pass, private :: add_element
          procedure, pass, private :: remove_element

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

          type(nbf_element), pointer :: rm_element


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

          integer :: i

          if(associated(this%head)) then

             current_element => this%head

             if(neighbor_index.eq.1) then

                do i=1, this%nb_elements
                   call current_element%copy_from_neighbor1_to(bf_exchanged)
                   current_element => current_element%get_next()
                end do

             else

                if(neighbor_index.eq.2) then
                   do i=1, this%nb_elements
                      call current_element%copy_from_neighbor2_to(bf_exchanged)
                      current_element => current_element%get_next()
                   end do
                else
                   print '(''nbf_list_class'')'
                   print '(''copy_from_neighbors_to_bf_layer'')'
                   print '(''neighbor_index not recognized'')'
                   print '(''neighbor_index: '',I2)', neighbor_index
                   stop 'verify neighbor_index (should be 1 or 2)'
                end if

             end if                
          end if

        end subroutine copy_from_neighbors_to_bf_layer


        !< copy the layers common between the buffer layers
        !> saved in the chained list as neighbors to the
        !> buffer layer passed as argument
        subroutine copy_to_neighbors_from_bf_layer(
     $     this, neighbor_index, bf_exchanged)

          implicit none

          class(nbf_list)  , intent(inout) :: this
          integer          , intent(in)    :: neighbor_index
          type(bf_sublayer), intent(in)    :: bf_exchanged

          integer :: i

          if(associated(this%head)) then

             current_element => this%head

             if(neighbor_index.eq.1) then

                do i=1, this%nb_elements
                   call current_element%copy_to_neighbor1_to(bf_exchanged)
                   current_element => current_element%get_next()
                end do

             else

                if(neighbor_index.eq.2) then
                   do i=1, this%nb_elements
                      call current_element%copy_to_neighbor2_to(bf_exchanged)
                      current_element => current_element%get_next()
                   end do
                else
                   print '(''nbf_list_class'')'
                   print '(''copy_from_neighbors_to_bf_layer'')'
                   print '(''neighbor_index not recognized'')'
                   print '(''neighbor_index: '',I2)', neighbor_index
                   stop 'verify neighbor_index (should be 1 or 2)'
                end if

             end if                
          end if

        end subroutine copy_to_neighbors_from_bf_layer


        !< add element in the chained list and its position
        !> is such that the chained list is ordered
        subroutine add_element(this, new_element)

          implicit none

          class(nbf_list)           , intent(inout) :: this
          type(nbf_element), pointer, intent(in)    :: new_element

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
                  this%head => added_bf_sublayer_ptr

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
                     this%tail_sublayer => new_element

                  else

                     current_sublayer_prev => current_sublayer%get_prev()
                     call insert_element(
     $                    new_element,
     $                    current_sublayer_prev,
     $                    current_sublayer)

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

      end module nbf_list_class
