      !> @file
      !> class encapsulating subroutines to apply wall boundary
      !> conditions and reflection at the edge of the
      !> computational domain
      !
      !                      wall b.c.
      !                          |
      !>                 //////////////////
      !>                 |               //
      !>                 |               //
      !>    reflection ->| computational //  <-  wall b.c.
      !>     b.c.        |    domain     //
      !>                 |               //
      !>                 |               //
      !                  //////////////////
      !                          |
      !                      wall b.c.
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the 
      !> gridpoints at the edge of the computational domain
      !
      !> @date
      !> 27_09_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_class

        use bc_operators_abstract_class, only : bc_operators_abstract
        use sd_operators_class         , only : sd_operators
        use pmodel_eq_class            , only : pmodel_eq
        use parameters_input           , only : nx,ny,ne,bc_size
        use parameters_kind            , only : rkind,ikind
        use reflection_xy_module       , only : reflection_x_prefactor
        use wall_xy_module             , only : wall_prefactor,
     $                                          compute_wall_flux_x,
     $                                          compute_wall_flux_y
        
        implicit none


        private
        public :: bc_operators


        !> @class bc_operators
        !> class encapsulating subroutines to apply
        !> wall and reflection boundary conditions in the
        !> x and y directions at the edge of the computational
        !> domain
        !>
        !> @param prefactor
        !> wall prefactor for the application of the nodes b.c.
        !>
        !> @param prefactor_x
        !> reflection prefactor for the application of the nodes
        !> b.c.
        !> 
        !> @param initialize
        !> initialize the prefactor and x_prefactor
        !> attributes of the boundary conditions
        !>
        !> @param apply_bc_on_nodes
        !> apply the wall and reflection boundary conditions
        !> along the x and y directions at the edge of the
        !> computational domain for the field
        !>
        !> @param apply_bc_on_fluxes
        !> apply the wall and reflection boundary conditions
        !> for the fluxes
        !---------------------------------------------------------------
        type, extends(bc_operators_abstract) :: bc_operators

          integer, dimension(ne) :: prefactor
          integer, dimension(ne) :: prefactor_x

          contains

          procedure,   pass :: ini
          procedure,   pass :: apply_bc_on_nodes
          procedure, nopass :: apply_bc_on_fluxes

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
        subroutine ini(this, p_model)
        
          implicit none

          class(bc_operators), intent(inout) :: this
          type(pmodel_eq)    , intent(in)    :: p_model

          
          this%prefactor   = wall_prefactor(p_model)
          this%prefactor_x = reflection_x_prefactor(p_model)

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
        !> 24_09_2013 - initial version - J.L. Desmarais
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


          integer(ikind) :: i,j
          integer        :: k


          !< compute the reflection b.c. in W boundary layers
          !> and the wall b.c. in the E boundary layers
          do k=1,ne
             do j=1+bc_size, ny-bc_size
                !DEC$ IVDEP
                do i=1,bc_size
                   
                   nodes(i,j,k) = 
     $                this%prefactor_x(k)*nodes(2*bc_size+1-i,j,k)
                   nodes(nx-bc_size+i,j,k) = 
     $                this%prefactor(k)*nodes(nx-bc_size-i+1,j,k)
                   
                end do
             end do
          end do


          !< compute the wall b.c. in N and S boundary layers
          do k=1, ne
             do j=1, bc_size
                !DEC$ IVDEP
                do i=1, nx
                   
                   nodes(i,j,k) = 
     $                  this%prefactor(k)*nodes(i,2*bc_size+1-j,k)
                   nodes(i,ny-bc_size+j,k) = 
     $                  this%prefactor(k)*nodes(i,ny-bc_size-j+1,k)
                   
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
        subroutine apply_bc_on_fluxes(nodes,dx,dy,s,flux_x,flux_y)

          implicit none

          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          real(rkind)                       , intent(in)    :: dx
          real(rkind)                       , intent(in)    :: dy
          type(sd_operators)                , intent(in)    :: s
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y

          integer        :: k
          integer(ikind), dimension(2) :: id


          !< provide the x-indices modified
          id(1)=nx+1-bc_size

          !< modify the fluxes along the x-direction
          !> only the E border: i=nx-bc_size+1
          !DEC$ FORCEINLINE RECURSIVE
          call compute_wall_flux_x(nodes,dx,dy,s,id(1),flux_x)

          !< provide the y-indices modified
          id(1)=bc_size+1
          id(2)=ny-bc_size+1

          !< modify the fluxes along the y-direction
          !> at the N and S borders
          !> S border: j= bc_size+1
          !> N border: j= ny-bc_size+1
          do k=1,2
             !DEC FORCEINLINE RECURSIVE
             call compute_wall_flux_y(nodes,dx,dy,s,id(k),flux_y)
          end do

        end subroutine apply_bc_on_fluxes

      end module bc_operators_class
