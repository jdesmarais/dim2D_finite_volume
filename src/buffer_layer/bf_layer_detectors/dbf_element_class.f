      !> @file
      !> module implementing the temporary object saving the
      !> general coordinates of a detector
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the temporary object saving the
      !> general coordinates of a detector
      !
      !> @date
      ! 16_04_2014 - initial version      - J.L. Desmarais
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !----------------------------------------------------------------
      module dbf_element_class

        use parameters_kind, only : ikind

        implicit none

        private
        public :: dbf_element


        !>@class dbf_element
        !< element of a double chained list of  save the general
        !> coordinates of a detector
        !
        !>@param prev
        !> pointer to the previous element in the list
        !
        !>@param next
        !> pointer to the next element in the list
        !
        !>@param coords
        !> general coordinates of the detector
        !
        !>@param ini
        !> initialize the dbf_element object its links to
        !> previous and next elements in the list and setting
        !> coordinates of the detector
        !
        !>@param set_prev
        !> set the link to the previous element of the list
        !
        !>@param set_next
        !> set the link to the next element of the list
        !
        !>@param get_prev
        !> access the previous element of the list
        !
        !>@param get_next
        !> access the next element of the list
        !
        !>@param get_coords
        !> get the coords attribute
        !
        !>@param nullify_prev
        !> remove the link to the previous element of the list
        !
        !>@param nullify_next
        !> remove the link to the next element of the list
        !
        !>@param print_on_screen
        !> print the content of the list on screen (only for test)
        !---------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the dbf_element object its links to
        !> previous and next elements in the list and setting
        !> coordinates of the detector
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_element object encapsulating the detector coordinates
        !> and the previous and next elements in the doubled chained
        !> list
        !
        !>@param added_bf_sublayer_ptr
        !> pointer to the bf_sublayer stored in the element
        !--------------------------------------------------------------
        subroutine ini(this, coords)

          implicit none

          class(dbf_element)          , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: coords

          this%coords = coords
          nullify(this%prev)
          nullify(this%next)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the link to the previous element of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_element object encapsulating the detector coordinates
        !> and the previous and next elements in the doubled chained
        !> list
        !
        !>@param prev_ptr
        !> pointer to the previous element of the list
        !--------------------------------------------------------------
        subroutine set_prev(this, prev_ptr)

          implicit none

          class(dbf_element)        , intent(inout) :: this
          type(dbf_element), pointer, intent(in)    :: prev_ptr

          this%prev => prev_ptr

        end subroutine set_prev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the link to the next element of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_element object encapsulating the detector coordinates
        !> and the previous and next elements in the doubled chained
        !> list
        !
        !>@param next_ptr
        !> pointer to the next element of the list
        !--------------------------------------------------------------
        subroutine set_next(this, next_ptr)

          implicit none

          class(dbf_element)        , intent(inout) :: this
          type(dbf_element), pointer, intent(in)    :: next_ptr

          this%next => next_ptr

        end subroutine set_next
      
        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> access the previous element of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_element object encapsulating the detector coordinates
        !> and the previous and next elements in the doubled chained
        !> list
        !
        !>@return get_prev
        !> pointer to the previous element of the list
        !--------------------------------------------------------------
        function get_prev(this)

          implicit none

          class(dbf_element), intent(in) :: this
          type(dbf_element) , pointer    :: get_prev

          get_prev => this%prev

        end function get_prev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> access the next element of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_element object encapsulating the detector coordinates
        !> and the previous and next elements in the doubled chained
        !> list
        !
        !>@return get_next
        !> pointer to the next element of the list
        !--------------------------------------------------------------
        function get_next(this)

          implicit none

          class(dbf_element), intent(in) :: this
          type(dbf_element) , pointer    :: get_next

          get_next => this%next

        end function get_next


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the ptr attribute
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_element object encapsulating the detector coordinates
        !> and the previous and next elements in the doubled chained
        !> list
        !
        !>@return get_coords
        !> coords attribute
        !--------------------------------------------------------------
        function get_coords(this)

          implicit none
          
          class(dbf_element), intent(in) :: this
          integer(ikind), dimension(2)   :: get_coords

          get_coords = this%coords

        end function get_coords


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the link to the previous element of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_element object encapsulating the detector coordinates
        !> and the previous and next elements in the doubled chained
        !> list
        !--------------------------------------------------------------
        subroutine nullify_prev(this)

          implicit none

          class(dbf_element), intent(inout) :: this

          nullify(this%prev)

        end subroutine nullify_prev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the link to the next element of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_element object encapsulating the detector coordinates
        !> and the previous and next elements in the doubled chained
        !> list
        !--------------------------------------------------------------
        subroutine nullify_next(this)

          implicit none

          class(dbf_element), intent(inout) :: this

          nullify(this%next)

        end subroutine nullify_next


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the content of the element on screen
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_element object encapsulating the detector coordinates
        !> and the previous and next elements in the doubled chained
        !> list
        !--------------------------------------------------------------
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
