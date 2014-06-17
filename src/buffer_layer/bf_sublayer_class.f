      !> @file
      !> module implementing the object encapsulating 
      !> a sublayer element, an element of a double
      !> chained list
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the sublyare element, an
      !> element of a double chained list saving a buffer
      !> layer
      !
      !> @date
      ! 11_04_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_sublayer_class

        use bf_layer_class     , only : bf_layer

        use parameters_constant, only : N,S,E,W,
     $                                  x_direction, y_direction,
     $                                  min_border
        use parameters_input   , only : nx,ny,ne
        use parameters_kind    , only : ikind, rkind

        implicit none

        private
        public :: bf_sublayer


        !> @class bf_sublayer
        !> class encapsulating the bf_layer/interior interface
        !> object but only its main attributes are implemented
        !>
        !> @param element
        !> buffer layer saved in the double chained list element
        !>
        !> @param next_bf_sublayer
        !> pointer of the next element of the double chained list
        !>
        !> @param prev_bf_sublayer
        !> pointer of the previous element of the double chained list
        !>
        !> @param ini
        !> nullify the pointers of the object
        !>
        !> @param link_next
        !> links the current sublayer element with the next sublayer
        !> element given
        !>
        !> @param link_prev
        !> links the current sublayer element with the previous sublayer
        !> element given
        !---------------------------------------------------------------
        type, extends(bf_layer) :: bf_sublayer

          type(bf_sublayer), pointer, private :: next
          type(bf_sublayer), pointer, private :: prev

          contains

          procedure, pass   :: ini

          procedure, pass   :: get_prev
          procedure, pass   :: get_next
          procedure, pass   :: set_prev
          procedure, pass   :: set_next
          procedure, pass   :: nullify_prev
          procedure, pass   :: nullify_next

          procedure, pass   :: remove

        end type bf_sublayer

        contains

        !< initialize the sublayer
        subroutine ini(this, localization)

          implicit none

          class(bf_sublayer), intent(inout) :: this
          integer(ikind)    , intent(in)    :: localization

          call this%bf_layer%ini(localization)

          nullify(this%prev)
          nullify(this%next)

        end subroutine ini


        !< access the previous list element
        function get_prev(this)

          implicit none

          class(bf_sublayer), intent(in) :: this
          type(bf_sublayer) , pointer    :: get_prev

          if(associated(this%prev)) then
             get_prev => this%prev
          else
             nullify(get_prev)
          end if

        end function get_prev


        !< access the next list element
        function get_next(this)

          implicit none

          class(bf_sublayer), intent(in) :: this
          type(bf_sublayer) , pointer    :: get_next

          if(associated(this%next)) then
             get_next => this%next
          else
             nullify(get_next)
          end if

        end function get_next


        !< link with the previous element in the chained list
        subroutine set_prev(this, bf_sublayer_prev)

          implicit none

          class(bf_sublayer)         , intent(inout) :: this
          type(bf_sublayer) , pointer, intent(in)    :: bf_sublayer_prev

          this%prev => bf_sublayer_prev

        end subroutine set_prev


        !< link with the next element in the chained list
        subroutine set_next(this, bf_sublayer_next)

          implicit none

          class(bf_sublayer)         , intent(inout) :: this
          type(bf_sublayer) , pointer, intent(in)    :: bf_sublayer_next

          this%next => bf_sublayer_next

        end subroutine set_next


        !< nullify the previous list element
        subroutine nullify_prev(this)

          implicit none

          class(bf_sublayer), intent(inout) :: this

          nullify(this%prev)

        end subroutine nullify_prev


        !< nullify the next list element
        subroutine nullify_next(this)

          implicit none

          class(bf_sublayer), intent(inout) :: this

          nullify(this%next)

        end subroutine nullify_next


        !< remove buffer layer
        subroutine remove(this)

          implicit none

          class(bf_sublayer), intent(inout) :: this

          
          if(associated(this%prev)) then
             if(associated(this%next)) then
                this%prev%next => this%next
                this%next%prev => this%prev
             else
                nullify(this%prev%next)
             end if
          else
             if(associated(this%next)) then
                nullify(this%next%prev)
             end if
          end if

        end subroutine remove        

      end module bf_sublayer_class
