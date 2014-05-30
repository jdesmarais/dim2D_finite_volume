      module nbf_element_class

        use bf_layer_class        , only : bf_layer
        use bf_layer_errors_module, only : error_incompatible_neighbor
        use bf_sublayer_class     , only : bf_sublayer
        use parameters_bf_layer   , only : align_N, align_S
        use parameters_constant   , only : N,S,E,W,
     $                                     x_direction, y_direction
        use parameters_input      , only : bc_size
        use parameters_kind       , only : ikind

        implicit none

        private
        public :: nbf_element


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
          procedure, pass :: set_ptr
          procedure, pass :: get_prev
          procedure, pass :: get_next
          procedure, pass :: get_ptr
          procedure, pass :: nullify_prev
          procedure, pass :: nullify_next
          procedure, pass :: nullify_ptr
          procedure, pass :: is_before
          procedure, pass :: refers_to
          !procedure, pass :: can_exchange_with
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


        subroutine set_ptr(this, added_bf_sublayer_ptr)

          implicit none

          class(nbf_element)        , intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: added_bf_sublayer_ptr

          if(associated(added_bf_sublayer_ptr)) then
             this%ptr => added_bf_sublayer_ptr
          end if

        end subroutine set_ptr
      
        
        function get_prev(this)

          implicit none

          class(nbf_element), intent(in) :: this
          type(nbf_element) , pointer    :: get_prev

          get_prev => this%prev

        end function get_prev


        function get_next(this)

          implicit none

          class(nbf_element), intent(inout) :: this
          type(nbf_element) , pointer       :: get_next

          get_next => this%next

        end function get_next


        function get_ptr(this)

          implicit none

          class(nbf_element), intent(inout) :: this
          type(bf_sublayer) , pointer       :: get_ptr

          get_ptr => this%ptr

        end function get_ptr


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
      

        subroutine nullify_ptr(this)

          implicit none

          class(nbf_element), intent(inout) :: this

          nullify(this%ptr)

        end subroutine nullify_ptr

      
        function is_before(this, element2)

          implicit none

          class(nbf_element), intent(in) :: this
          class(nbf_element), intent(in) :: element2
          logical                        :: is_before

          if(associated(this%ptr)) then
             if(associated(element2%ptr)) then
                
                is_before = 
     $               (this%ptr%get_alignment(x_direction,1)).le.
     $               (element2%ptr%get_alignment(x_direction,1))


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
c$$$        !> buffer sublayer are possible
c$$$        function can_exchange_with(this, bf_exchanged)
c$$$
c$$$          implicit none
c$$$
c$$$          class(nbf_element), intent(in) :: this
c$$$          class(bf_layer)   , intent(in) :: bf_exchanged
c$$$          logical                        :: can_exchange_with
c$$$
c$$$          integer(ikind) :: min_border
c$$$          integer(ikind) :: max_border
c$$$          integer        :: mainlayer_id
c$$$          integer        :: n_mainlayer_id
c$$$
c$$$
c$$$          !restrictions along the x-direction
c$$$          min_border = max(this%ptr%get_alignment(x_direction,1),
c$$$     $                     bf_exchanged%get_alignment(x_direction,1))
c$$$          max_border = min(this%ptr%get_alignment(x_direction,2),
c$$$     $                     bf_exchanged%get_alignment(x_direction,2))
c$$$
c$$$          can_exchange_with = (max_border-min_border+2*bc_size+1).gt.0
c$$$
c$$$
c$$$          !restrictions for E and W buffer layers
c$$$          !along the y-direction
c$$$          if(can_exchange_with) then
c$$$
c$$$             n_mainlayer_id = this%ptr%get_localization()
c$$$
c$$$             select case(n_mainlayer_id)
c$$$
c$$$               case(N)
c$$$                  mainlayer_id = bf_exchanged%get_localization()
c$$$                  select case(mainlayer_id)
c$$$                    case(E,W)
c$$$                       can_exchange_with =
c$$$     $                      bf_exchanged%get_alignment(y_direction,2).eq.(align_N-1)
c$$$                    
c$$$                    case default
c$$$                       call error_incompatible_neighbor(
c$$$     $                      'nbf_element_class.f',
c$$$     $                      'can_exchange_with',
c$$$     $                      mainlayer_id,
c$$$     $                      n_mainlayer_id)
c$$$                  end select
c$$$
c$$$               case(S)
c$$$                  mainlayer_id = bf_exchanged%get_localization()
c$$$                  select case(mainlayer_id)
c$$$                    case(E,W)
c$$$                       can_exchange_with =
c$$$     $                      bf_exchanged%get_alignment(y_direction,1).eq.(align_S+1)
c$$$                    case default
c$$$                       call error_incompatible_neighbor(
c$$$     $                      'nbf_element_class.f',
c$$$     $                      'can_exchange_with',
c$$$     $                      mainlayer_id,
c$$$     $                      n_mainlayer_id)
c$$$                  end select
c$$$
c$$$               case(E,W)
c$$$                  mainlayer_id = bf_exchanged%get_localization()
c$$$                  select case(mainlayer_id)
c$$$                    case(N)
c$$$                       can_exchange_with =
c$$$     $                      this%ptr%get_alignment(y_direction,2).eq.(align_N-1)
c$$$                    case(S)
c$$$                       can_exchange_with =
c$$$     $                      this%ptr%get_alignment(y_direction,1).eq.(align_S+1)
c$$$                    case default
c$$$                       call error_incompatible_neighbor(
c$$$     $                      'nbf_element_class.f',
c$$$     $                      'can_exchange_with',
c$$$     $                      mainlayer_id,
c$$$     $                      n_mainlayer_id)
c$$$                  end select
c$$$
c$$$             end select
c$$$
c$$$          end if
c$$$
c$$$        end function can_exchange_with


        !< copy from neighbor1 saved in the element to the
        !> buffer layer passed as argument
        subroutine copy_from_neighbor1_to(this, bf_exchanged)
        
          implicit none

          class(nbf_element), intent(in)    :: this
          class(bf_layer)   , intent(inout) :: bf_exchanged

          !if(this%can_exchange_with(bf_exchanged)) then
             call bf_exchanged%copy_from_neighbor1(this%ptr)
          !end if

        end subroutine copy_from_neighbor1_to


        !< copy from neighbor2 saved in the element to the
        !> buffer layer passed as argument
        subroutine copy_from_neighbor2_to(this, bf_exchanged)
        
          implicit none

          class(nbf_element), intent(in)    :: this
          class(bf_layer)   , intent(inout) :: bf_exchanged

          !if(this%can_exchange_with(bf_exchanged)) then
             call bf_exchanged%copy_from_neighbor2(this%ptr)
          !end if

        end subroutine copy_from_neighbor2_to


        !< copy from the buffer layer passed as argument to
        !> the neighbor1 saved in the element
        subroutine copy_to_neighbor1_from(this, bf_exchanged)
        
          implicit none

          class(nbf_element), intent(inout) :: this
          class(bf_layer)   , intent(in)    :: bf_exchanged

          !if(this%can_exchange_with(bf_exchanged)) then
             call bf_exchanged%copy_to_neighbor1(this%ptr)
          !end if

        end subroutine copy_to_neighbor1_from


        !< copy from the buffer layer passed as argument to
        !> to the neighbor2 saved in the element
        subroutine copy_to_neighbor2_from(this, bf_exchanged)
        
          implicit none

          class(nbf_element), intent(inout) :: this
          class(bf_layer)   , intent(in)    :: bf_exchanged

          !if(this%can_exchange_with(bf_exchanged)) then
             call bf_exchanged%copy_to_neighbor2(this%ptr)
          !end if

        end subroutine copy_to_neighbor2_from

      end module nbf_element_class
