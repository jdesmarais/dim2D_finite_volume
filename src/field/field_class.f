      !-------------------------------------------------------------------
      !
      ! MODULE: field_class
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the main tables for the variables and the
      !> coordinates
      !
      ! REVISION HISTORY:
      ! 07_08_2013 - initial version - J.L. Desmarais
      !-------------------------------------------------------------------
      module field_class

        use parameters_kind, only : ikind, rkind
        use surrogate_class, only : surrogate


        implicit none


        private
        public :: field


        !-----------------------------------------------------------------
        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> class encapsulating the main tables for the variables and the
        !> coordinates
        !
        ! REVISION HISTORY:
        ! 07_08_2013 - initial version - J.L. Desmarais
        !
        !>@attribute : nodes  : real, dimension(:,:,:), allocatable
        !>@attribute : x_map  : real, dimension(:)    , allocatable
        !>@attribute : y_map  : real, dimension(:)    , allocatable
        !>@attribute : dx     : real
        !>@attribute : dy     : real
        !-----------------------------------------------------------------
        type, extends(surrogate) :: field

          real(rkind), dimension(:,:,:), allocatable :: nodes
          real(rkind) ,dimension(:)    , allocatable :: x_map
          real(rkind) ,dimension(:)    , allocatable :: y_map
          real(rkind) :: dx
          real(rkind) :: dy

          contains

          procedure, pass, non_overridable :: allocate_tables

        end type field


        contains


        subroutine allocate_tables(this,nx,ny,ne)

          implicit none

          class(field)  , intent(inout) :: this
          integer(ikind), intent(in)    :: nx,ny
          integer       , intent(in)    :: ne

          if(allocated(this%nodes)) then
             deallocate(this%nodes)
          end if

          if(allocated(this%x_map)) then
             deallocate(this%x_map)
          end if

          if(allocated(this%y_map)) then
             deallocate(this%y_map)
          end if

          allocate(this%nodes(nx,ny,ne))
          allocate(this%x_map(nx))
          allocate(this%y_map(ny))

        end subroutine allocate_tables

      end module field_class
