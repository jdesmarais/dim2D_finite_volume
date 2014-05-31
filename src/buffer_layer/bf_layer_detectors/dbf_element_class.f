      module dbf_element_class

        use parameters_kind, only : ikind


        implicit none

        private
        public :: dbf_element


        !< element of a double chained list to save the general
        !> coordinates of a detector
        type :: dbf_element

          integer(ikind), dimension(2), private :: coords
          type(dbf_element), pointer  , private :: prev
          type(dbf_element), pointer  , private :: next

          contains

          procedure, pass :: ini
          procedure, pass :: set_prev
          procedure, pass :: set_next
          procedure, pass :: get_prev
          procedure, pass :: get_next
          procedure, pass :: get_coords
          procedure, pass :: nullify_prev
          procedure, pass :: nullify_next

          procedure, pass :: print_on_screen

        end type dbf_element


        contains


        subroutine ini(this, coords)

          implicit none

          class(dbf_element)          , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: coords

          this%coords = coords
          nullify(this%prev)
          nullify(this%next)

        end subroutine ini


        subroutine set_prev(this, prev_ptr)

          implicit none

          class(dbf_element)        , intent(inout) :: this
          type(dbf_element), pointer, intent(in)    :: prev_ptr

          this%prev => prev_ptr

        end subroutine set_prev


        subroutine set_next(this, next_ptr)

          implicit none

          class(dbf_element)        , intent(inout) :: this
          type(dbf_element), pointer, intent(in)    :: next_ptr

          this%next => next_ptr

        end subroutine set_next
      
        
        function get_prev(this)

          implicit none

          class(dbf_element), intent(in) :: this
          type(dbf_element) , pointer    :: get_prev

          get_prev => this%prev

        end function get_prev


        function get_next(this)

          implicit none

          class(dbf_element), intent(in) :: this
          type(dbf_element) , pointer    :: get_next

          get_next => this%next

        end function get_next


        function get_coords(this)

          implicit none
          
          class(dbf_element), intent(in) :: this
          integer(ikind), dimension(2)   :: get_coords

          get_coords = this%coords

        end function get_coords


        subroutine nullify_prev(this)

          implicit none

          class(dbf_element), intent(inout) :: this

          nullify(this%prev)

        end subroutine nullify_prev


        subroutine nullify_next(this)

          implicit none

          class(dbf_element), intent(inout) :: this

          nullify(this%next)

        end subroutine nullify_next


        subroutine print_on_screen(this)

          implicit none

          class(dbf_element), intent(in) :: this


          print '(''+ coords: '',2I3)', this%coords
          if(associated(this%prev)) then
             print '(''+ prev: '')'
             print '(''  - coords: '',2I3)', this%prev%coords
             print '(''  - prev associated: '', L1)', associated(this%prev%prev)
             print '(''  - next associated: '', L1)', associated(this%prev%next)
          else
             print '(''+ prev: X'')'
          end if
          if(associated(this%next)) then
             print '(''+ next: '')'
             print '(''  - coords: '',2I3)', this%next%coords
             print '(''  - prev associated: '', L1)', associated(this%next%prev)
             print '(''  - next associated: '', L1)', associated(this%next%next)
          else
             print '(''+ next: X'')'
          end if

        end subroutine print_on_screen


      end module dbf_element_class
