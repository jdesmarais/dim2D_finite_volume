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

        end type sbf_list


        contains

        
        !< initialize the number of elements to 0
        !> and allocate the array containing the pointers
        !> to bf_sublayer
        subroutine ini(this, nb_ele_max)

          implicit none

          class(sbf_list), intent(inout) :: this
          integer        , intent(in)    :: nb_ele_max

          this%nb_ele = 0
          allocate(this%list(nb_ele_max))

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


        subroutine add_ele(this, ptr)

          implicit none

          class(sbf_list)           , intent(inout) :: this
          type(bf_sublayer), pointer, intent(in)    :: ptr

          this%nb_ele = this%nb_ele+1

          call this%list(this%nb_ele)%set_ptr(ptr)
          
        end subroutine add_ele
        

      end module sbf_list_class
