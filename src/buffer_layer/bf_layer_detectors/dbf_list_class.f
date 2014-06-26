      !> @file
      !> module implementing the temporary object saving the
      !> general coordinates of a list of detectors
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the temporary object saving the
      !> general coordinates of a list of detectors
      !
      !> @date
      ! 16_04_2014 - initial version      - J.L. Desmarais
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !----------------------------------------------------------------
      module dbf_list_class

        use dbf_element_class, only : dbf_element
        use parameters_kind  , only : ikind

        implicit none


        !>@class dbf_list
        !> double chained list saving the general coordinates of
        !> a list of detectors
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
        !>@param add_to_list
        !> add the detector general coordinates to the list
        !
        !>@param add_element
        !> add a new element to the list
        !
        !>@param print_on_screen
        !> print the content of the list on screen
        !
        !>@param destroy
        !> deallocate the content of the list
        !---------------------------------------------------------------
        type :: dbf_list

          type(dbf_element), pointer :: head
          type(dbf_element), pointer :: tail
          integer                    :: nb_elements
         
          contains

          procedure, pass :: ini
          
          procedure, pass :: get_head
          procedure, pass :: get_tail
          procedure, pass :: get_nb_elements

          procedure, pass          :: add_to_list
          procedure, pass, private :: add_element

          procedure, pass :: print_on_screen

          procedure, pass :: destroy

        end type dbf_list


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
        !> dbf_list object implementing a doubled chained
        !> list storing the general coordinates of a list
        !> of detectors
        !--------------------------------------------------------------
        subroutine ini(this)

          implicit none

          class(dbf_list), intent(inout) :: this
          
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
        !> dbf_list object implementing a doubled chained
        !> list storing the general coordinates of a list
        !> of detectors
        !
        !>@param get_head
        !> reference to the first element of the doubled chained
        !> list
        !--------------------------------------------------------------
        function get_head(this)

          implicit none

          class(dbf_list), intent(in) :: this
          type(dbf_element), pointer  :: get_head
          
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
        !> dbf_list object implementing a doubled chained
        !> list storing the general coordinates of a list
        !> of detectors
        !
        !>@param get_head
        !> reference to the last element of the doubled chained
        !> list
        !--------------------------------------------------------------
        function get_tail(this)

          implicit none

          class(dbf_list), intent(in) :: this
          type(dbf_element), pointer  :: get_tail
          
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
        !> dbf_list object implementing a doubled chained
        !> list storing the general coordinates of a list
        !> of detectors
        !
        !>@param get_nb_elements
        !> nb_element attribute
        !--------------------------------------------------------------
        function get_nb_elements(this)

          implicit none

          class(dbf_list), intent(in) :: this
          integer                     :: get_nb_elements
          
          get_nb_elements = this%nb_elements

        end function get_nb_elements
        

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add the detector general coordinates to the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_list object implementing a doubled chained
        !> list storing the general coordinates of a list
        !> of detectors
        !
        !>@param coords
        !> general coordinates of the detector
        !--------------------------------------------------------------
        subroutine add_to_list(this, coords)

          implicit none

          class(dbf_list)             , intent(inout) :: this
          integer(ikind), dimension(2), intent(in)    :: coords

          type(dbf_element), pointer :: new_element

          !create the new element
          allocate(new_element)
          call new_element%ini(coords)

          !add the new element to the list
          call add_element(this, new_element)

        end subroutine add_to_list

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add a new element to the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_list object implementing a doubled chained
        !> list storing the general coordinates of a list
        !> of detectors
        !
        !>@param new_element
        !> dbf_element object added to the list
        !--------------------------------------------------------------
        subroutine add_element(this, new_element)

          implicit none

          class(dbf_list)           , intent(inout) :: this
          type(dbf_element), pointer, intent(inout) :: new_element

          !add the new element to the list
          if(this%nb_elements.eq.0) then
             this%head => new_element
             this%tail => new_element
          else
             call new_element%set_prev(this%tail)
             call this%tail%set_next(new_element)
             this%tail => new_element
          end if

          this%nb_elements = this%nb_elements+1

        end subroutine add_element


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> print the content of the list on screen
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_list object implementing a doubled chained
        !> list storing the general coordinates of a list
        !> of detectors
        !--------------------------------------------------------------
        subroutine print_on_screen(this)
        
          implicit none

          class(dbf_list), intent(in) :: this

          integer :: i
          type(dbf_element), pointer :: current_element

          if(this%nb_elements.gt.0) then
             current_element => this%head

             do i=1, this%nb_elements
                print '(''element '',I3)', i
                print '(''----------------'')'
                call current_element%print_on_screen()
                current_element => current_element%get_next()
                print '()'
             end do

          else
             print '(''empty list'')'
          end if

        end subroutine print_on_screen


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate the content of the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dbf_list object implementing a doubled chained
        !> list storing the general coordinates of a list
        !> of detectors
        !--------------------------------------------------------------
        subroutine destroy(this)

          implicit none

          class(dbf_list), intent(inout) :: this

          type(dbf_element), pointer :: current_element
          type(dbf_element), pointer :: next_element
          integer                    :: i


          if(this%nb_elements.gt.0) then

             current_element => this%head

             do i=1, this%nb_elements
                next_element => current_element%get_next()
                call current_element%nullify_prev()
                call current_element%nullify_next()
                deallocate(current_element)
                current_element => next_element
             end do
             
             nullify(this%head)
             nullify(this%tail)

          end if

        end subroutine destroy

      end module dbf_list_class
