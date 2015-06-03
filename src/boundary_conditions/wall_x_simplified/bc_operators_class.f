      !> @file
      !> wall b.c. at bottom, fixed nodes/no fluxes
      !> for the other boundaries
      !
      !                    fixed nodes
      !                   + no fluxes
      !                 ________|________
      !                 |               |
      !                 |               |
      !  fixed nodes -> | computational | <-  fixed nodes
      ! + no fluxes     |   domain      |     + no fluxes
      !                 |               |
      !                 |               |
      !                 /////////////////
      !                         |
      !                     wall b.c.
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to apply
      !> wall boundary conditions at the bottom
      !> and simple fixed nodes+no-fluxes at the
      !> other boundaries
      !
      !> @date
      !> 03_06_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_class

        use bc_operators_default_class, only :
     $       bc_operators_default

        use dim2d_parameters, only :
     $       cv_r

        use dim2d_prim_module, only :
     $       temperature_eff

        use dim2d_state_eq_module, only :
     $       get_mass_density_vapor,
     $       get_mass_density_liquid

        use parameters_constant, only :
     $       N,S,E,W,
     $       bc_flux_and_node_choice,
     $       left

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size,
     $       debug_initialize_bc_nodes,
     $       debug_real

        use parameters_kind, only :
     $       rkind,
     $       ikind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use sd_operators_fd_module, only :
     $       gradient_x_interior,
     $       gradient_y_interior

        use wall_xy_equilibrium_module, only :
     $       compute_wall_S_ghost_cell,
     $       compute_wall_flux_y

        implicit none


        private
        public :: bc_operators


        !> @class bc_operators
        !> class encapsulating subroutines to apply
        !> wall and reflection boundary conditions in the
        !> x and y directions at the edge of the computational
        !> domain
        !> 
        !> @param ini
        !> initialize the type of boundary conditions
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
        type, extends(bc_operators_default) :: bc_operators

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

          integer :: i_s

          this%bc_type(N) = bc_flux_and_node_choice
          this%bc_type(S) = bc_flux_and_node_choice
          this%bc_type(E) = bc_flux_and_node_choice
          this%bc_type(W) = bc_flux_and_node_choice

          i_s = p_model%get_eq_nb()

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
        subroutine apply_bc_on_nodes(this,t,x_map,y_map,nodes_tmp,nodes)

          implicit none

          class(bc_operators)             , intent(in)    :: this
          real(rkind)                     , intent(in)    :: t
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes

          integer(ikind) :: i,j
          integer        :: k

          integer, dimension(4) :: prefactor_x
          integer, dimension(4) :: prefactor_y

          real(rkind) :: dx
          real(rkind) :: dy
          real(rkind) :: T_guess
          real(rkind) :: md_guess

          i = this%bc_type(N)


          prefactor_x = [ 1,-1, 1, 1]
          prefactor_y = [ 1, 1,-1, 1]

          
          ! to check whether the nodes corresponding to
          ! where the wall b.c are correctly initialized,
          ! they are first initialized with debug_real
          !--------------------------------------------------
          if(debug_initialize_bc_nodes) then
             do k=1,ne
                do j=1,bc_size
                   do i=bc_size+1,nx-bc_size
                      nodes(i,j,k) = debug_real
                   end do
                end do
             end do
          end if


          ! E and W layer: reflection
          !--------------------------------------------------
          do k=1,ne
             !DEC$ IVDEP
             do j=1+bc_size, ny-bc_size
                !DEC$ IVDEP
                do i=1,bc_size
                   
                   nodes(i,j,k)            = prefactor_x(k)*nodes( 2*bc_size-i+1,j,k)
                   nodes(nx-bc_size+i,j,k) = prefactor_x(k)*nodes(nx-bc_size-i+1,j,k)
                   
                end do
             end do
          end do


          ! SE and SW corners: reflection
          !--------------------------------------------------
          do k=1,ne
             !DEC$ IVDEP
             do j=1,bc_size
                !DEC$ IVDEP
                do i=1,bc_size
                   nodes(i,j,k) = prefactor_y(k)*nodes(i,2*bc_size-j+1,k)
                end do

                !DEC$ IVDEP
                do i=nx-bc_size+1,nx
                   nodes(i,j,k) = prefactor_y(k)*nodes(i,2*bc_size-j+1,k)
                end do
             end do
          end do


          ! S layer: wall b.c.
          !--------------------------------------------------
          j=bc_size
          !DEC$ IVDEP
          do i=1+bc_size, nx-bc_size
                   
             dx = x_map(2)-x_map(1)
             dy = y_map(2)-y_map(1)
             
             !temperature: temperature at the previous step
             T_guess = 3.0d0/(8.0d0*cv_r)*temperature_eff(
     $            nodes_tmp,i,j+1,
     $            dx,dy,
     $            gradient_x_interior,
     $            gradient_y_interior)

             !mass density: mass density at the previous step
             md_guess = nodes_tmp(i,j,1)

             call compute_wall_S_ghost_cell(
     $            i,j,
     $            t,
     $            x_map,
     $            y_map,
     $            nodes,
     $            T_guess,
     $            md_guess)

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
        !>@param t
        !> time
        !
        !>@param x_map
        !> x-coordinates
        !
        !>@param y_map
        !> y-coordinates
        !
        !>@param nodes
        !> governing variables
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
        subroutine apply_bc_on_fluxes(t,x_map,y_map,nodes,s,flux_x,flux_y)

          implicit none

          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(sd_operators)                , intent(in)    :: s
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y

          integer(ikind) :: i
          integer(ikind) :: j


          !> the fluxes for the S layer are computed
          !> using the wall boundary conditions
          j=bc_size+1
          do i=bc_size+1,nx+1-bc_size
             flux_y(i,j,:) = compute_wall_flux_y(
     $            nodes,s,i,j,
     $            t,x_map,y_map,
     $            left)
          end do

        end subroutine apply_bc_on_fluxes

      end module bc_operators_class
