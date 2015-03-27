      !> @file
      !> icr_path augmented with references to the first and 
      !> next elements of a double chain list
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> icr_path augmented with references to the first and 
      !> next elements of a double chain list
      !
      !> @date
      ! 21_03_2015 - initial version      - J.L. Desmarais
      !-----------------------------------------------------------------
      module icr_path_chain_class

        use icr_path_class, only :
     $     icr_path

        implicit none

        private
        public :: icr_path_chain


        !> @class icr_path_chain
        !> icr_path augmented with references to the first and 
        !> next elements of a double chain list
        !
        !> @param next
        !> pointer of the next element of the double chained list
        !
        !> @param prev
        !> pointer of the previous element of the double chained list
        !
        !> @param ini
        !> initialize the bf_layer parent and nullify the pointers
        !> to the prev and next elements of the doubled chained list
        !
        !>@param get_prev
        !> access the previous element of the list
        !
        !>@param get_next
        !> access the next element of the list
        !
        !>@param set_prev
        !> set the link to the previous element of the list
        !
        !>@param set_next
        !> set the link to the next element of the list
        !
        !>@param nullify_prev
        !> remove the link to the previous element of the list
        !
        !>@param nullify_next
        !> remove the link to the next element of the list
        !
        !>@param remove
        !> remove the element of the list by removing the bf_layer
        !> parent object and the links to the previous and next
        !> elements in the list
        !---------------------------------------------------------------
        type, extends(icr_path) :: icr_path_chain

          type(icr_path_chain), pointer :: next
          type(icr_path_chain), pointer :: prev

          contains

          procedure, pass :: ini

          procedure, pass :: get_prev
          procedure, pass :: get_next
          procedure, pass :: set_prev
          procedure, pass :: set_next
          procedure, pass :: nullify_prev
          procedure, pass :: nullify_next

          procedure, pass :: remove

        end type icr_path_chain
        

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the icr_path parent and nullify the pointers
        !> to the prev and next elements of the doubled chained list
        !
        !> @date
        !> 21_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> icr_path augmented with references to the first and 
        !> next elements of a double chain list
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(icr_path_chain), intent(inout) :: this

          call this%icr_path%ini()

          nullify(this%prev)
          nullify(this%next)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> access the previous element of the list
        !
        !> @date
        !> 21_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> icr_path augmented with references to the first and 
        !> next elements of a double chain list
        !
        !>@return get_prev
        !> pointer to the previous element of the list
        !--------------------------------------------------------------
        function get_prev(this)

          implicit none

          class(icr_path_chain), intent(in) :: this
          type(icr_path_chain) , pointer    :: get_prev

          if(associated(this%prev)) then
             get_prev => this%prev
          else
             nullify(get_prev)
          end if

        end function get_prev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> access the next element of the list
        !
        !> @date
        !> 21_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> icr_path augmented with references to the first and 
        !> next elements of a double chain list
        !
        !>@return get_next
        !> pointer to the next element of the list
        !--------------------------------------------------------------
        function get_next(this)

          implicit none

          class(icr_path_chain), intent(in) :: this
          type(icr_path_chain) , pointer    :: get_next

          if(associated(this%next)) then
             get_next => this%next
          else
             nullify(get_next)
          end if

        end function get_next


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the link to the previous element of the list
        !
        !> @date
        !> 21_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> icr_path augmented with references to the first and 
        !> next elements of a double chain list
        !
        !>@param icr_path_prev
        !> pointer to the previous element of the list
        !--------------------------------------------------------------
        subroutine set_prev(this, icr_path_prev)

          implicit none

          class(icr_path_chain)         , intent(inout) :: this
          type(icr_path_chain) , pointer, intent(in)    :: icr_path_prev

          this%prev => icr_path_prev

        end subroutine set_prev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the link to the next element of the list
        !
        !> @date
        !> 21_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> icr_path augmented with references to the first and 
        !> next elements of a double chain list
        !
        !>@param icr_path_prev
        !> pointer to the next element of the list
        !--------------------------------------------------------------
        subroutine set_next(this, icr_path_next)

          implicit none

          class(icr_path_chain)         , intent(inout) :: this
          type(icr_path_chain) , pointer, intent(in)    :: icr_path_next

          this%next => icr_path_next

        end subroutine set_next


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the link to the previous element of the list
        !
        !> @date
        !> 21_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> icr_path augmented with references to the first and 
        !> next elements of a double chain list
        !--------------------------------------------------------------
        subroutine nullify_prev(this)

          implicit none

          class(icr_path_chain), intent(inout) :: this

          nullify(this%prev)

        end subroutine nullify_prev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the link to the next element of the list
        !
        !> @date
        !> 21_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> icr_path augmented with references to the first and 
        !> next elements of a double chain list
        !--------------------------------------------------------------
        subroutine nullify_next(this)

          implicit none

          class(icr_path_chain), intent(inout) :: this

          nullify(this%next)

        end subroutine nullify_next


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> remove the element of the list by removing the icr_path
        !> parent object and the links to the previous and next
        !> elements in the list
        !
        !> @date
        !> 21_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> icr_path augmented with references to the first and 
        !> next elements of a double chain list
        !--------------------------------------------------------------
        subroutine remove(this)

          implicit none

          class(icr_path_chain), intent(inout) :: this


          !remove icr_path
          call this%icr_path%remove()

          !nullify the links
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

      end module icr_path_chain_class
