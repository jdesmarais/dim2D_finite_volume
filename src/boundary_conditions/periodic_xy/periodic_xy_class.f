      !> @file
      !> class encapsulating subroutines to apply periodic
      !> boundary conditions in the x and y directions
      !> at the edge of the computational domain
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to apply periodic
      !> boundary conditions in the x and y directions
      !> at the edge of the computational domain
      !
      !> @date
      !> 13_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module periodic_xy_class

        use cg_operators_class, only : cg_operators
        use field_class       , only : field
        use parameters_kind   , only : ikind

        implicit none
        
        
        private
        public :: periodic_xy


        !> @class periodic_xy
        !> class encapsulating subroutines to apply periodic
        !> boundary conditions in the x and y directions
        !> at the edge of the computational domain
        !>
        !> @param apply_bc_on_nodes
        !> apply the periodic boundary conditions along the x and y
        !> directions at the edge of the computational domain
        !---------------------------------------------------------------
        type :: periodic_xy

          contains

          procedure, nopass, non_overridable :: apply_bc_on_nodes

        end type periodic_xy


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the periodic boundary conditions
        !> along the x and y directions at the edge of the
        !> computational domain
        !
        !> @date
        !> 13_08_2013 - initial version - J.L. Desmarais
        !
        !>@param f
        !> object encapsulating the main variables
        !
        !>@param s
        !> space discretization operators
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes(f,s)
        
          implicit none

          class(field)      , intent(inout) :: f
          type(cg_operators), intent(in)    :: s


          integer(ikind) :: i
          integer(ikind) :: j
          integer        :: k
          integer        :: bc_size
          integer        :: period_x
          integer        :: period_y


          bc_size  = s%get_bc_size()
          period_x = size(f%nodes,1)-2*bc_size
          period_y = size(f%nodes,2)-2*bc_size


          !<compute the east and west boundary layers
          !>without the north and south corners
          do k=1, size(f%nodes,3)
             do j=1+bc_size, size(f%nodes,2)-bc_size
                do i=1, bc_size

                   f%nodes(i,j,k)=f%nodes(i+period_x,j,k)
                   f%nodes(i+period_x+bc_size,j,k)=f%nodes(i+bc_size,j,k)
                   
                end do
             end do
          end do


          !<compute the south and north layers
          !>with the east and west corners
          do k=1, size(f%nodes,3)
             do j=1, bc_size
                do i=1, size(f%nodes,1)

                   f%nodes(i,j,k)=f%nodes(i,j+period_y,k)
                   f%nodes(i,j+period_y+bc_size,k)=f%nodes(i,j+bc_size,k)

                end do
             end do
          end do

        end subroutine apply_bc_on_nodes

      end module periodic_xy_class
