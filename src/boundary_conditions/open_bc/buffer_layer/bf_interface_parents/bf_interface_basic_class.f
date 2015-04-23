      !> @file
      !> object encapsulating the buffer layers around the interior domain
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> object encapsulating the buffer layers around the interior domain
      !
      !> @date
      ! 09_03_2015 - documentation update - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_interface_basic_class

        use netcdf

        use bf_mainlayer_class, only :
     $       bf_mainlayer

        use bf_mainlayer_pointer_class, only :
     $       bf_mainlayer_pointer

        use mainlayer_interface_icr_class, only :
     $       mainlayer_interface_icr

        use parameters_input, only :
     $       nx,ny

        use parameters_kind, only :
     $       rkind

        implicit none

        private
        public :: bf_interface_basic
        

        !>@class bf_interface_basic
        !> object encapsulating the buffer layers around
        !> the interior domain
        !
        !>@param mainlayer_pointers
        !> table with reference to the bf_mainlayers objects gathering
        !> the bf_sublayer corresponding to a specific cardinal point
        !
        !>@param mainlayer_interfaces
        !> object encapsulating references to the sublayers at the edge
        !> between the main layers and subroutines to synchronize the data
        !> between them
        !
        !>@param ini
        !> initialize the buffer/interior domain interface
        !---------------------------------------------------------------
        type :: bf_interface_basic

          type(bf_mainlayer_pointer), dimension(4) :: mainlayer_pointers
          type(mainlayer_interface_icr)            :: mainlayer_interfaces

          contains

          procedure, pass :: ini
          procedure, pass :: get_mainlayer_ptr

        end type bf_interface_basic

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the buffer/interior domain interface
        !
        !> @date
        !> 26_06_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        ! the data between them
        !--------------------------------------------------------------
        subroutine ini(this,interior_x_map,interior_y_map)

          implicit none

          class(bf_interface_basic) , intent(inout) :: this
          real(rkind), dimension(nx), intent(in)    :: interior_x_map
          real(rkind), dimension(ny), intent(in)    :: interior_y_map

          integer     :: i
          real(rkind) :: x_s
          real(rkind) :: y_s
          

          !initialize the references to the mainlayer pointers
          do i=1, size(this%mainlayer_pointers,1)
             call this%mainlayer_pointers(i)%ini()             
          end do

          !initialize the mainlayer interfaces
          call this%mainlayer_interfaces%ini()

          x_s = interior_x_map(1)
          y_s = interior_y_map(2)

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get a reference to the mainlayer corresponding to
        !> the cardinal coordinate passed
        !
        !> @date
        !> 25_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> bf_interface object encapsulating the buffer layers
        !> around the interior domain and subroutines to synchronize
        ! the data between them
        !
        !>@param mainlayer_id
        !> cardinal coordinate identifying the main layer asked
        !--------------------------------------------------------------
        function get_mainlayer_ptr(this,mainlayer_id)
     $     result(bf_mainlayer_ptr)

          implicit none

          class(bf_interface_basic), intent(in) :: this
          integer                  , intent(in) :: mainlayer_id
          type(bf_mainlayer), pointer           :: bf_mainlayer_ptr


          bf_mainlayer_ptr => this%mainlayer_pointers(mainlayer_id)%get_ptr()

        end function get_mainlayer_ptr
      

      end module bf_interface_basic_class
