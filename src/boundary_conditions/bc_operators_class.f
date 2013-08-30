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

        use cg_operators_class  , only : cg_operators
        use dim2d_eq_class      , only : dim2d_eq
        use field_class         , only : field
        use parameters_constant , only : periodic_xy_choice,
     $                                   reflection_xy_choice
        use parameters_input    , only : nx,ny,ne,bc_choice
        use parameters_kind     , only : rkind,ikind
        use periodic_xy_module  , only : apply_periodic_xy_on_nodes
        use reflection_xy_module, only : apply_reflection_xy_on_nodes
        
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
        subroutine apply_bc_on_nodes(field_used,p_model,s)

          implicit none

          class(field)      , intent(inout) :: field_used
          type(dim2d_eq)    , intent(in)    :: p_model
          type(cg_operators), intent(in)    :: s

          !<select the type of boundary conditions
          select case(bc_choice)

            case(periodic_xy_choice)

               !DEC$ FORCEINLINE RECURSIVE
               call apply_periodic_xy_on_nodes(
     $              field_used%nodes,s)


            case(reflection_xy_choice)

               !DEC$ FORCEINLINE RECURSIVE
               call apply_reflection_xy_on_nodes(
     $              field_used%nodes,p_model,s)
               

            case default
               
               print '(''bc_class: apply_bc_on_nodes'')'
               stop 'boundary conditions not recognized'

          end select

        end subroutine apply_bc_on_nodes

      end module bc_operators_class
