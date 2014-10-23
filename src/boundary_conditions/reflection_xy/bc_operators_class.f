      !> @file
      !> class encapsulating subroutines to apply reflection boundary
      !> conditions at the edge of the computational domain
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the 
      !> gridpoints at the edge of the computational domain
      !
      !> @date
      !> 24_09_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_class

        use bc_operators_default_class, only : bc_operators_default
        use pmodel_eq_class           , only : pmodel_eq
        use parameters_constant       , only : bc_nodes_choice
        use parameters_input          , only : nx,ny,ne,bc_size
        use parameters_kind           , only : rkind,ikind
        use reflection_xy_module      , only : reflection_x_prefactor,
     $                                         reflection_y_prefactor
        use sd_operators_class        , only : sd_operators
        
        implicit none


        private
        public :: bc_operators


        !> @class bc_operators
        !> class encapsulating subroutines to apply
        !> reflection boundary conditions in the x and
        !> y directions at the edge of the computational
        !> domain
        !>
        !> @param prefactor_x
        !> prefactor for the computation of the reflection
        !> boundary conditions along the x-direction
        !>
        !> @param prefactor_y
        !> prefactor for the computation of the reflection
        !> boundary conditions along the y-direction
        !>
        !> @param initialize
        !> initialize the prefactor_x and prefactor_y
        !> attributes of the boundary conditions
        !>
        !> @param apply_bc_on_nodes
        !> apply the reflection boundary conditions along
        !> the x and y directions at the edge of the
        !> computational domain for the field
        !>
        !> @param apply_bc_on_fluxes
        !> apply the reflection boundary conditions for
        !> the fluxes
        !---------------------------------------------------------------
        type, extends(bc_operators_default) :: bc_operators

          integer, dimension(ne) :: prefactor_x
          integer, dimension(ne) :: prefactor_y

          contains

          procedure,   pass :: ini
          procedure,   pass :: apply_bc_on_nodes

        end type bc_operators


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing the main attributes
        !> of the boundary conditions
        !
        !> @date
        !> 24_09_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary conditions initialized
        !
        !>@param p_model
        !> physical model to know the type of the main variables
        !--------------------------------------------------------------
        subroutine ini(this, p_model)
        
          implicit none

          class(bc_operators), intent(inout) :: this
          type(pmodel_eq)    , intent(in)    :: p_model

          this%prefactor_x = reflection_x_prefactor(p_model)
          this%prefactor_y = reflection_y_prefactor(p_model)

          this%bc_type = [
     $         bc_nodes_choice,
     $         bc_nodes_choice,
     $         bc_nodes_choice,
     $         bc_nodes_choice]

        end subroutine ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> along the x and y directions at the edge of the
        !> computational domain
        !
        !> @date
        !> 23_09_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary conditions
        !
        !>@param f
        !> object encapsulating the main variables
        !
        !>@param s
        !> space discretization operators
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes(this,nodes)

          implicit none

          class(bc_operators)             , intent(in)    :: this
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes


          integer(ikind)         :: i,j
          integer                :: k
          

          !< compute the reflection b.c. in E and W boundary layers
          do k=1,ne
             do j=1+bc_size, ny-bc_size
                !DEC$ IVDEP
                do i=1,bc_size
                   
                   nodes(i,j,k) = 
     $                  this%prefactor_x(k)*nodes(2*bc_size+1-i,j,k)
                   nodes(nx-bc_size+i,j,k) = 
     $                  this%prefactor_x(k)*nodes(nx-bc_size-i+1,j,k)
                   
                end do
             end do
          end do


          !< compute the reflection b.c. in N and S boundary layers
          do k=1, ne
             do j=1, bc_size
                !DEC$ IVDEP
                do i=1, nx
                   
                   nodes(i,j,k) = 
     $                  this%prefactor_y(k)*nodes(i,2*bc_size+1-j,k)
                   nodes(i,ny-bc_size+j,k) = 
     $                  this%prefactor_y(k)*nodes(i,ny-bc_size-j+1,k)
                   
                end do
             end do
          end do

        end subroutine apply_bc_on_nodes

      end module bc_operators_class
