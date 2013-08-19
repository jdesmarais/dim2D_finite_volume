      !> @file
      !> class encapsulating subroutines to apply boundary conditions
      !> at the edge of the computational domain
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the 
      !> gridpoints at the edge of the computational domain
      !
      !> @date
      !> 13_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_class

        use cg_operators_class , only : cg_operators
        use parameters_constant, only : periodic_xy_choice
        use parameters_input   , only : nx,ny,ne,bc_choice
        use parameters_kind    , only : rkind,ikind
        use periodic_xy_class  , only : periodic_xy
        
        implicit none


        private
        public :: bc_operators


        !> @class bc_operators
        !> class encapsulating subroutines to apply
        !> boundary conditions in the x and y directions
        !> at the edge of the computational domain
        !>
        !> @param apply_bc_on_nodes
        !> apply the boundary conditions along the x and y
        !> directions at the edge of the computational domain
        !> for the field (ex:periodic_bc)
        !---------------------------------------------------------------
        type :: bc_operators

          contains

          procedure, nopass :: apply_bc_on_nodes

        end type bc_operators


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
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
        subroutine apply_bc_on_nodes(nodes,s)

          implicit none

          real(rkind), dimension(nx,ny,ne) :: nodes
          type(cg_operators), intent(in)    :: s

          type(periodic_xy) :: periodic_bc


          !<select the type of boundary conditions
          select case(bc_choice)

            case(periodic_xy_choice)

               !TAG INLINE
               !call periodic_bc%apply_bc_on_nodes(nodes,s)
               nodes(1,1,1)=1

            case default
               
               print '(''bc_class: apply_bc_on_nodes'')'
               stop 'boundary conditions not recognized'

          end select

        end subroutine apply_bc_on_nodes

      end module bc_operators_class
