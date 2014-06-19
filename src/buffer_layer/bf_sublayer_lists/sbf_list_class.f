      module sbf_list_class

        use bf_sublayer_class        , only : bf_sublayer
        use bf_sublayer_pointer_class, only : bf_sublayer_pointer

        implicit none

        private
        public :: sbf_list


        !< object containing an array of pointers to bf_sublayer objects
        type :: sbf_list

          type(bf_sublayer_pointer), dimension(:), allocatable :: list
          integer :: nb_ele
                    
          contains

          procedure, pass :: ini
          procedure, pass :: get_nb_ele
          procedure, pass :: get_ele
          procedure, pass :: add_ele
          procedure, pass :: does_not_contain
          procedure, pass :: destroy

        end type sbf_list


        contains

        
        !< initialize the number of elements to 0
        !> and allocate the array containing the pointers
        !> to bf_sublayer
        subroutine ini(this, nb_ele_max)

          implicit none

          class(sbf_list)          , intent(inout) :: this
          integer        , optional, intent(in)    :: nb_ele_max

          this%nb_ele = 0
          if(present(nb_ele_max)) then
             allocate(this%list(nb_ele_max))
          end if

        end subroutine ini

      
        function get_nb_ele(this)

          implicit none

          class(sbf_list), intent(in) :: this
          integer                     :: get_nb_ele

          get_nb_ele = this%nb_ele

        end function get_nb_ele


        function get_ele(this, i) result(ptr)

          implicit none

          class(sbf_list), intent(inout) :: this
          integer        , intent(in)    :: i
          type(bf_sublayer), pointer     :: ptr

          ptr => this%list(i)%get_ptr()

        end function get_ele


        !add element in the list only if the element
        !is not present already
        subroutine add_ele(this, ptr)

          implicit none

          class(sbf_list)           , intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: ptr

          if(does_not_contain(this,ptr)) then
             this%nb_ele = this%nb_ele+1
             call this%list(this%nb_ele)%set_ptr(ptr)
          end if
          
        end subroutine add_ele

        
        function does_not_contain(this, ptr)

          implicit none

          class(sbf_list)           , intent(in) :: this
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
        

        subroutine destroy(this)
        
          implicit none

          class(sbf_list), intent(inout) :: this

          if(allocated(this%list)) then
             deallocate(this%list)
          end if

        end subroutine destroy

      end module sbf_list_class
