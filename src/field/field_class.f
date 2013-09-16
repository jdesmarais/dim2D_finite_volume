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

        use parameters_input, only : nx,ny,ne,x_min,x_max,y_min,y_max
        use parameters_kind , only : ikind, rkind
        use surrogate_class , only : surrogate

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
          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(nx)       :: x_map
          real(rkind), dimension(ny)       :: y_map
          real(rkind) :: dx
          real(rkind) :: dy

          contains

          procedure, pass :: ini_coordinates

        end type field


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to integrate the governing equations using
        !> the numerical scheme developed by C.W.Shu and S.Osher
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the main variables
        !
        !>@param x_min
        !> coordinate along the x-axis of the SW border
        !
        !>@param x_max
        !> coordinate along the x-axis of the NE border
        !
        !>@param y_min
        !> coordinate along the y-axis of the SW border
        !
        !>@param y_max
        !> coordinate along the y-axis of the NE border
        !
        !>@param bc_size
        !> size of the boundary layer
        !--------------------------------------------------------------
        subroutine ini_coordinates(this,bc_size)

          implicit none

          class(field), intent(inout) :: this
          integer     , intent(in)    :: bc_size


          integer(ikind) :: i,j


          !< initialize the space steps along the 
          !> x and y directions
          this%dx = (x_max-x_min)/(nx-1-2*bc_size)
          this%dy = (y_max-y_min)/(ny-1-2*bc_size)


          !< initialize the coordinates along the
          !> x-direction
          do i=1, nx
             this%x_map(i)=x_min + (i-1-bc_size)*this%dx
          end do


          !< initialize the coordinates along the
          !> y-direction
          do j=1, ny
             this%y_map(j)=y_min + (j-1-bc_size)*this%dy
          end do

        end subroutine ini_coordinates

      end module field_class
