      !> @file
      !> module implementing the dcr_list object. It encapsulates an array
      !> of pointers to bf_sublayer objects
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module implementing the dcr_list object. It encapsulates an array
      !> of pointers to bf_sublayer objects
      !
      !> @date
      ! 07_04_2014 - initial version      - J.L. Desmarais
      ! 27_06_2014 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
       module dcr_list_class

        use bf_sublayer_class, only :
     $      bf_sublayer

        use dcr_chain_class, only :
     $       dcr_chain

        implicit none

        private
        public :: dcr_list


        !>@class dcr_listx
        !> class encapsulating an array of pointers to bf_sublayer
        !> objects
        !
        !>@param list
        !> array of pointers to bf_sublayer objects
        !
        !>@param nb_ele
        !> number of references in the list attribute
        !
        !>@param ini
        !> initialize the number of elements to 0
        !> and allocate the array containing the pointers
        !> to bf_sublayer objects
        !
        !>@param get_nb_ele
        !> get the nb_ele attribute
        !
        !>@param get_ele
        !> get the bf_sublayer reference at index i
        !
        !>@param add_ele
        !> add element in the list only if the element
        !> is not present already
        !
        !>@param does_not_contain
        !> check that a reference is not present in the list
        !
        !>@param remove
        !> deallocate the content of the object
        !--------------------------------------------------------------
        type :: dcr_list

          type(dcr_chain), dimension(:), allocatable :: list
          integer :: nb_ele
                    
          contains

          procedure, pass :: ini
          procedure, pass :: get_nb_ele
          procedure, pass :: get_ele
          procedure, pass :: add_ele
          procedure, pass :: does_not_contain
          procedure, pass :: remove

        end type dcr_list


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the number of elements to 0
        !> and allocate the array containing the pointers
        !> to bf_sublayer
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dcr_list object encapsulating references to bf_sublayer
        !> objects
        !
        !>@param nb_ele_max
        !> maximum number of references stored in the list
        !--------------------------------------------------------------
        subroutine ini(this, nb_ele_max)

          implicit none

          class(dcr_list)          , intent(inout) :: this
          integer        , optional, intent(in)    :: nb_ele_max

          this%nb_ele = 0
          if(present(nb_ele_max)) then
             allocate(this%list(nb_ele_max))
          end if

        end subroutine ini

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the nb_ele attribute
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dcr_list object encapsulating references to bf_sublayer
        !> objects
        !
        !>@return get_nb_ele
        !> nb_ele attribute
        !--------------------------------------------------------------
        function get_nb_ele(this)

          implicit none

          class(dcr_list), intent(in) :: this
          integer                     :: get_nb_ele

          get_nb_ele = this%nb_ele

        end function get_nb_ele


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the bf_sublayer reference at index i
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dcr_list object encapsulating references to bf_sublayer
        !> objects
        !
        !>@param i
        !> index identifying the reference asked
        !
        !>@return ptr
        !> bf_sublayer reference
        !--------------------------------------------------------------
        function get_ele(this, i) result(ptr)

          implicit none

          class(dcr_list), intent(inout) :: this
          integer        , intent(in)    :: i
          type(bf_sublayer), pointer     :: ptr

          ptr => this%list(i)%get_ptr()

        end function get_ele


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> add element in the list only if the element
        !> is not present already
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dcr_list object encapsulating references to bf_sublayer
        !> objects
        !
        !>@param ptr
        !> bf_sublayer reference
        !--------------------------------------------------------------
        subroutine add_ele(this, ptr)

          implicit none

          class(dcr_list)           , intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: ptr

          if(does_not_contain(this,ptr)) then
             this%nb_ele = this%nb_ele+1
             call this%list(this%nb_ele)%set_ptr(ptr)
          end if
          
        end subroutine add_ele

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check that a reference is not present in the list
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dcr_list object encapsulating references to bf_sublayer
        !> objects
        !
        !>@param ptr
        !> bf_sublayer reference asked
        !
        !>@return does_not_contain
        !> logical stating whether the reference was found in the list
        !--------------------------------------------------------------
        function does_not_contain(this, ptr)

          implicit none

          class(dcr_list)           , intent(in) :: this
          type(bf_sublayer), pointer, intent(in) :: ptr
          logical                                :: does_not_contain

          integer :: k

          does_not_contain = .true.

          do k=1, this%nb_ele
             
             if(this%list(k)%associated_ptr(ptr)) then
                does_not_contain = .false.
                exit
             end if
             
          end do

        end function does_not_contain
        

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> deallocate the content of the object
        !
        !> @date
        !> 27_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> dcr_list object encapsulating references to bf_sublayer
        !> objects
        !--------------------------------------------------------------
        subroutine remove(this)
        
          implicit none

          class(dcr_list), intent(inout) :: this

          if(allocated(this%list)) then
             deallocate(this%list)
          end if

        end subroutine remove

      end module dcr_list_class
