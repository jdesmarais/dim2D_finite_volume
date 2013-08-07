      !> @file
      !> class encapsulating the main tables for the variables and the
      !> coordinates
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating the main tables for the variables and the
      !> coordinates
      !
      !> @date
      ! 07_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module field_class

        use parameters_kind, only : ikind, rkind
        use surrogate_class, only : surrogate


        implicit none


        private
        public :: field


        !> @class field
        !> class encapsulating the variables of the governing equations
        !> and the discretisation maps
        !>
        !> @param nodes
        !> variables computed during the simulation
        !> (ex: mass=nodes(:,:,1),
        !> momentum_x=nodes(:,:,2),
        !> momentum_y=nodes(:,:,3),
        !> energy=nodes(:,:,4))
        !>
        !> @param x_map
        !> discretisation map along the x-axis
        !>
        !> @param y_map
        !> discretisation map along the y-axis
        !>
        !> @param dx
        !> space step along the x-axis
        !>
        !> @param dy
        !> space step along the y-axis
        !---------------------------------------------------------------
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


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> give the boundary layer size
        !
        !> @date
        !> 07_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> field object allocated
        !
        !>@param nx
        !> extension along the x-axis
        !
        !>@param ny
        !> extension along the y-axis
        !
        !>@param ne
        !> number of conservative variables in the governing equations
        !---------------------------------------------------------------
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
