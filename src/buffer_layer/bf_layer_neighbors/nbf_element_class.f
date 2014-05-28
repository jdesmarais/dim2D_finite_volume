      module nbf_element_class

        use bf_layer_class     , only : bf_layer
        use bf_sublayer_class  , only : bf_sublayer
        use parameters_constant, only : x_direction, y_direction
        use parameters_kind    , only : ikind

        implicit none

        private
        public : nbf_element


        !< element of a double chained list to save a pointer
        !> to a buffer sublayer
        type :: nbf_element

          type(nbf_element), pointer, private :: prev
          type(nbf_element), pointer, private :: next
          type(bf_sublayer), pointer, private :: ptr

          contains

          procedure, pass :: ini
          procedure, pass :: set_prev
          procedure, pass :: set_next
          procedure, pass :: get_prev
          procedure, pass :: get_next
          procedure, pass :: nullify_prev
          procedure, pass :: nullify_next
          procedure, pass :: set_element
          procedure, pass :: is_before
          procedure, pass :: refers_to
          procedure, pass :: can_exchange_with
          procedure, pass :: copy_from_neighbor1_to
          procedure, pass :: copy_from_neighbor2_to
          procedure, pass :: copy_to_neighbor1_from
          procedure, pass :: copy_to_neighbor2_from

        end type nbf_element

        contains

        subroutine ini(this, added_bf_sublayer_ptr)

          implicit none

          class(nbf_element)         , intent(inout) :: this
          type(bf_sublayer) , pointer, intent(in)    :: added_bf_sublayer_ptr

          if(associated(added_bf_sublayer_ptr)) then
             this%ptr => added_bf_sublayer_ptr
          else
             nullify(this%ptr)
          end if

          nullify(this%prev)
          nullify(this%next)

        end subroutine ini


        subroutine set_prev(this, prev_ptr)

          implicit none

          class(nbf_element)        , intent(inout) :: this
          type(nbf_element), pointer, intent(in)    :: prev_ptr

          this%prev => prev_ptr

        end subroutine set_prev


        subroutine set_next(this, next_ptr)

          implicit none

          class(nbf_element)        , intent(inout) :: this
          type(nbf_element), pointer, intent(in)    :: next_ptr

          this%next => next_ptr

        end subroutine set_next


        function get_prev(this)

          implicit none

          class(nbf_element), intent(in) :: this
          type(nbf_element) , pointer    :: get_prev

          get_prev => this%prev

        end function get_prev


        subroutine get_next(this)

          implicit none

          class(nbf_element), intent(inout) :: this
          type(nbf_element) , pointer       :: get_next

          get_next => this%next

        end subroutine get_next


        subroutine nullify_prev(this)

          implicit none

          class(nbf_element), intent(inout) :: this

          nullify(this%prev)

        end subroutine nullify_prev


        subroutine nullify_next(this)

          implicit none

          class(nbf_element), intent(inout) :: this

          nullify(this%next)

        end subroutine nullify_next
      

        subroutine set_element(this, added_bf_sublayer_ptr)

          implicit none

          class(nbf_element)        , intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: added_bf_sublayer_ptr

          if(associated(added_bf_sublayer_ptr)) then
             this%ptr => added_bf_sublayer_ptr
          end if

        end subroutine set_element

      
        function is_before(this, element2)

          implicit none

          class(nbf_element), intent(in) :: this
          class(nbf_element), intent(in) :: element2
          logical                        :: is_before

          if(associated(this%ptr)) then
             if(associated(element2%ptr)) then
                
                is_before = 
     $               (element2%ptr%get_alignment(x_direction,1)).lt.
     $               (this%ptr%get_alignment(x_direction,1))

             else
                print '(''nbf_element_class'')'
                print '(''is_before'')'
                print '(''element2%ptr not associated'')'
                stop 'cannot compare elements'
             end if
          else
             print '(''nbf_element_class'')'
             print '(''is_before'')'
             print '(''this%ptr not associated'')'
             stop 'is element initialized ?'
          end if

        end function is_before


        function refers_to(this, added_bf_sublayer_ptr)

          implicit none

          class(nbf_element)        , intent(in) :: this
          type(bf_sublayer), pointer, intent(in) :: added_bf_sublayer_ptr
          logical                                :: refers_to

          if(associated(this%ptr)) then
             if(associated(added_bf_sublayer_ptr)) then
                
                refers_to = associated(this%ptr, added_bf_sublayer_ptr)

             else
                print '(''nbf_element_class'')'
                print '(''is_same_element'')'
                print '(''added_bf_sublayer_ptr not associated'')'
                stop 'cannot compare elements'
             end if
          else
             print '(''nbf_element_class'')'
             print '(''is_same_element'')'
             print '(''this%ptr not associated'')'
             stop 'is element initialized ?'
          end if

        end function refers_to


        !< check if the exchange between the element and the
        !> buffer sublayer are possible
        function can_exchange_with(this, bf_exchanged)

          implicit none

          class(nbf_element), intent(in) :: this
          class(bf_layer)   , intent(in) :: bf_exchanged
          logical                        :: can_exchange_with

          integer(ikind) :: min_border
          integer(ikind) :: max_border


          min_border = max(this%ptr%get_alignment(x_direction,1),
     $                     bf_exchanged%get_alignment(x_direction,1))
          max_border = min(this%ptr%get_alignment(x_direction,2),
     $                     bf_exchanged%get_alignment(x_direction,2))

          can_exchange_with = (max_border-min_border).gt.0

        end function can_exchange_with


        !< copy from neighbor1 saved in the element to the
        !> buffer layer passed as argument
        subroutine copy_from_neighbor1_to(this, bf_exchanged)
        
          implicit none

          class(nbf_element), intent(in)    :: this
          class(bf_layer)   , intent(inout) :: bf_exchanged

          if(this%can_exchange_with(bf_exchanged)) then
             call bf_exchanged%copy_from_neighbor1(this%ptr)
          end if

        end subroutine copy_from_neighbor1_to


        !< copy from neighbor2 saved in the element to the
        !> buffer layer passed as argument
        subroutine copy_from_neighbor2_to(this, bf_exchanged)
        
          implicit none

          class(nbf_element), intent(in)    :: this
          class(bf_layer)   , intent(inout) :: bf_exchanged

          if(this%can_exchange_with(bf_exchanged)) then
             call bf_exchanged%copy_from_neighbor2(this%ptr)
          end if

        end subroutine copy_from_neighbor2_to


        !< copy from the buffer layer passed as argument to
        !> the neighbor1 saved in the element
        subroutine copy_to_neighbor1_from(this, bf_exchanged)
        
          implicit none

          class(nbf_element), intent(inout) :: this
          class(bf_layer)   , intent(in)    :: bf_exchanged

          if(this%can_exchange_with(bf_exchanged)) then
             call bf_exchanged%copy_to_neighbor1(this%ptr)
          end if

        end subroutine copy_to_neighbor1_from


        !< copy from the buffer layer passed as argument to
        !> to the neighbor2 saved in the element
        subroutine copy_to_neighbor2_from(this, bf_exchanged)
        
          implicit none

          class(nbf_element), intent(inout) :: this
          class(bf_layer)   , intent(in)    :: bf_exchanged

          if(this%can_exchange_with(bf_exchanged)) then
             call bf_exchanged%copy_to_neighbor2(this%ptr)
          end if

        end subroutine copy_to_neighbor2_from

      end module nbf_element_class
