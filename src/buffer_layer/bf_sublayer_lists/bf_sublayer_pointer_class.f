      module bf_sublayer_pointer_class

        use bf_sublayer_class, only : bf_sublayer

        implicit none

        private
        public :: bf_sublayer_pointer

        
        !< object containing a pointer reference to a
        !> bf_sublayer object
        type :: bf_sublayer_pointer

          type(bf_sublayer), pointer :: ptr

          contains

          procedure, pass :: get_ptr
          procedure, pass :: set_ptr

        end type bf_sublayer_pointer


        contains


        !< get the pointer to the buffer sublayer
        function get_ptr(this) result(ptr)
        
          implicit none

          class(bf_sublayer_pointer), intent(in) :: this
          type(bf_sublayer), pointer             :: ptr

          if(associated(this%ptr)) then
              ptr => this%ptr
          else
             nullify(ptr)
          end if

        end function get_ptr


        !< set the pointer to the buffer sublayer
        subroutine set_ptr(this, ptr)
        
          implicit none

          class(bf_sublayer_pointer), intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: ptr

          if(associated(ptr)) then
             this%ptr => ptr
          else
             nullify(this%ptr)
          end if

        end subroutine set_ptr

      end module bf_sublayer_pointer_class
