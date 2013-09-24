      !> @file
      !> class encapsulating subroutines to apply wall boundary
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

        use bc_abstract_class , only : bc_abstract
        use cg_operators_class, only : cg_operators
        use dim2d_eq_class    , only : dim2d_eq
        use field_class       , only : field
        use parameters_input  , only : nx,ny,ne
        use parameters_kind   , only : rkind,ikind
        use wall_xy_module    , only : wall_prefactor,
     $                                 wall_fx_momentum_x,
     $                                 wall_fx_momentum_y,
     $                                 wall_fy_momentum_x,
     $                                 wall_fy_momentum_y,
     $                                 wall_heat_flux
        
        implicit none


        private
        public :: bc_operators


        !> @class bc_operators
        !> class encapsulating subroutines to apply
        !> wall boundary conditions in the x and
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
        !> apply the wall boundary conditions along the x and y
        !> directions at the edge of the computational domain
        !> for the field
        !>
        !> @param apply_bc_on_fluxes
        !> apply the wall boundary conditions for the fluxes
        !---------------------------------------------------------------
        type, extends(bc_abstract) :: bc_operators

          integer, dimension(ne) :: prefactor

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
        !>@param s
        !> spatial discretisation operators
        !
        !>@param p_model
        !> physical model to know the type of the main variables
        !--------------------------------------------------------------
        subroutine initialize(this, s, p_model)
        
          implicit none

          class(bc_operators), intent(inout) :: this
          type(cg_operators) , intent(in)    :: s
          type(dim2d_eq)     , intent(in)    :: p_model

          
          integer :: bc_size

          bc_size = s%get_bc_size()

          this%prefactor = wall_prefactor(p_model)

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
        !>@param s
        !> space discretization operators
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes(this,f_used,s)

          implicit none

          class(bc_operators), intent(in)    :: this
          class(field)       , intent(inout) :: f_used
          type(cg_operators) , intent(in)    :: s


          integer(ikind) :: i,j
          integer        :: bc_size,k


          !< compute the size of the boundary layer
          bc_size = s%get_bc_size()


          !< compute the reflection b.c. in E and W boundary layers
          do k=1,ne
             do j=1+bc_size, ny-bc_size
                !DEC$ IVDEP
                do i=1,bc_size
                   
                   f_used%nodes(i,j,k) = 
     $                  this%prefactor(k)*f_used%nodes(2*bc_size+1-i,j,k)
                   f_used%nodes(nx-bc_size+i,j,k) = 
     $                  this%prefactor(k)*f_used%nodes(nx-bc_size-i+1,j,k)
                   
                end do
             end do
          end do


          !< compute the reflection b.c. in N and S boundary layers
          do k=1, ne
             do j=1, bc_size
                !DEC$ IVDEP
                do i=1, nx
                   
                   f_used%nodes(i,j,k) = 
     $                  this%prefactor(k)*f_used%nodes(i,2*bc_size+1-j,k)
                   f_used%nodes(i,ny-bc_size+j,k) = 
     $                  this%prefactor(k)*f_used%nodes(i,ny-bc_size-j+1,k)
                   
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

          integer(ikind) :: i,j
          integer        :: k, bc_size
          integer(ikind), dimension(2) :: id


          !< get the size of the boundary layer
          bc_size = s%get_bc_size()


          !< provide the x-indices modified
          id(1)=bc_size+1
          id(2)=nx+1-bc_size


          !< modify the fluxes along the x-direction
          !> at the E and W borders
          !> W border: i= bc_size+1
          !> E border: i= nx-bc_size+1
          do j=bc_size+1, ny-bc_size             
             do k=1,2

                !< no mass entering the system
                flux_x(id(k),j,1) = 0.0d0
                
                !< b.c. for the momentum along the x-direction
                flux_x(id(k),j,2) = wall_fx_momentum_x(f_used,s,id(k),j)

                !< b.c. for the momentum along the y-direction
                flux_x(id(k),j,3) = wall_fx_momentum_y(f_used,s,id(k),j)

                !< constant heat flux entering the system
                flux_x(id(k),j,4) = wall_heat_flux(f_used,id(k),j)

             end do
          end do


          !< provide the y-indices modified
          id(1)=bc_size+1
          id(2)=ny-bc_size+1

          
          !< modify the fluxes along the y-direction
          !> at the N and S borders
          !> S border: j= bc_size+1
          !> N border: j= ny-bc_size+1
          do i=bc_size+1, nx-bc_size
             do k=1,2

                !< no mass entering the system
                flux_y(i,id(k),1)= 0.0d0
                
                !< b.c. for the momentum along the x-direction
                flux_y(i,id(k),2)= wall_fy_momentum_x(f_used,s,i,id(k))

                !< b.c. for the momentum along the y-direction
                flux_y(i,id(k),3)= wall_fy_momentum_y(f_used,s,i,id(k))

                !< constant heat flux entering the system
                flux_y(i,id(k),4)= wall_heat_flux(f_used,i,id(k))

             end do
          end do

        end subroutine apply_bc_on_fluxes

      end module bc_operators_class
