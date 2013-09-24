      !> @file
      !> class encapsulating subroutines to apply periodic boundary
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
      !> 23_09_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_class

        use bc_abstract_class , only : bc_abstract
        use cg_operators_class, only : cg_operators
        use dim2d_eq_class    , only : dim2d_eq
        use field_class       , only : field
        use parameters_input  , only : nx,ny,ne
        use parameters_kind   , only : rkind,ikind
        
        implicit none


        private
        public :: bc_operators


        !> @class bc_operators
        !> class encapsulating subroutines to apply
        !> periodic boundary conditions in the x and
        !> y directions at the edge of the computational
        !> domain
        !>
        !> @param period_x
        !> period along the x-direction
        !>
        !> @param period_y
        !> period along the y-direction
        !> 
        !> @param initialize
        !> initialize the period_x and period_y
        !> attributes of the boundary conditions
        !>
        !> @param apply_bc_on_nodes
        !> apply the periodic boundary conditions along the x and y
        !> directions at the edge of the computational domain
        !> for the field
        !>
        !> @param apply_bc_on_fluxes
        !> apply the periodic boundary conditions for the fluxes
        !---------------------------------------------------------------
        type, extends(bc_abstract) :: bc_operators

          integer(ikind) :: period_x
          integer(ikind) :: period_y

          contains

          procedure, pass :: initialize
          procedure, pass :: apply_bc_on_nodes
          procedure, pass :: apply_bc_on_fluxes

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
        subroutine initialize(this, p_model,s)
        
          implicit none

          class(bc_operators), intent(inout) :: this
          type(dim2d_eq)     , intent(in)    :: p_model
          type(cg_operators) , intent(in)    :: s
          
          integer :: neq,bc_size

          neq     = p_model%get_eq_nb()
          bc_size = s%get_bc_size()

          this%period_x = nx-2*bc_size
          this%period_y = ny-2*bc_size

        end subroutine initialize


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> along the x and y directions at the edge of the
        !> computational domain
        !
        !> @date
        !> 24_09_2013 - initial version - J.L. Desmarais
        !
        !>@param f
        !> object encapsulating the main variables
        !
        !>@param p_model
        !> physical model
        !
        !>@param s
        !> space discretization operators
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes(this,f_used,p_model,s)

          implicit none

          class(bc_operators), intent(in)    :: this
          class(field)       , intent(inout) :: f_used
          type(dim2d_eq)     , intent(in)    :: p_model
          type(cg_operators) , intent(in)    :: s

          integer(ikind) :: i,j
          integer        :: neq,bc_size,k

          neq      = p_model%get_eq_nb()
          bc_size  = s%get_bc_size()

          !<compute the east and west boundary layers
          !>without the north and south corners
          do k=1, ne
             !DEC$ IVDEP
             do j=1+bc_size, ny-bc_size
                !DEC$ IVDEP
                do i=1, bc_size

                   f_used%nodes(i,j,k)=
     $                  f_used%nodes(i+this%period_x,j,k)
                   f_used%nodes(i+this%period_x+bc_size,j,k)=
     $                  f_used%nodes(i+bc_size,j,k)
                   
                end do
             end do
          end do


          !<compute the south and north layers
          !>with the east and west corners
          do k=1, ne
             !DEC$ IVDEP
             do j=1, bc_size
                !DEC$ IVDEP
                do i=1, nx

                   f_used%nodes(i,j,k)=
     $                  f_used%nodes(i,j+this%period_y,k)
                   f_used%nodes(i,j+this%period_y+bc_size,k)=
     $                  f_used%nodes(i,j+bc_size,k)

                end do
             end do
          end do

        end subroutine apply_bc_on_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes along the x directions at the
        !> edge of the computational domain
        !
        !> @date
        !> 24_09_2013 - initial version - J.L. Desmarais
        !
        !>@param this
        !> boundary conditions
        !
        !>@param f_used
        !> object encapsulating the main variables
        !
        !>@param s
        !> space discretization operators
        !
        !>@param flux_x
        !> flux along the x-direction
        !
        !>@param flux_y
        !> flux along the y-direction
        !--------------------------------------------------------------
        subroutine apply_bc_on_fluxes(this,f_used,s,flux_x,flux_y)

          implicit none

          class(bc_operators)               , intent(in)    :: this
          class(field)                      , intent(in)    :: f_used
          type(cg_operators)                , intent(in)    :: s
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y

          integer     :: period,bc_size
          real(rkind) :: node,flux

          stop 'periodic_xy: apply_bc_on_fluxes not implemented'

          period=this%period_x
          node=f_used%nodes(1,1,1)
          bc_size=s%get_bc_size()
          flux=flux_x(1,1,1)
          flux=flux_y(1,1,1)

        end subroutine apply_bc_on_fluxes

      end module bc_operators_class
