      !> @file
      !> wall boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to apply
      !> wall boundary conditions
      !
      !> @date
      !> 08_06_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_wall_xy_class

        use bc_operators_default_class, only :
     $       bc_operators_default

        use dim2d_parameters, only :
     $       cv_r

        use errors_module, only :
     $       error_bc_section_type

        use dim2d_prim_module, only :
     $       temperature_eff

        use dim2d_state_eq_module, only :
     $       get_mass_density_vapor,
     $       get_mass_density_liquid

        use interface_primary, only :
     $       gradient_proc

        use parameters_constant, only :
     $       N,S,E,W,
     $       bc_flux_and_node_choice,
     $       left,right,
     $       sd_interior_type,
     $       sd_L0_type,
     $       sd_L1_type,
     $       sd_R1_type,
     $       sd_R0_type,
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type

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
     $       gradient_y_interior,
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_L1,
     $       gradient_x_x_oneside_R1,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_L1,
     $       gradient_y_y_oneside_R1,
     $       gradient_y_y_oneside_R0

        use sd_operators_x_oneside_L0_class, only :
     $       sd_operators_x_oneside_L0

        use sd_operators_x_oneside_L1_class, only :
     $       sd_operators_x_oneside_L1

        use sd_operators_x_oneside_R1_class, only :
     $       sd_operators_x_oneside_R1

        use sd_operators_x_oneside_R0_class, only :
     $       sd_operators_x_oneside_R0

        use wall_xy_equilibrium_module, only :
     $       compute_wall_N_ghost_cell,
     $       compute_wall_S_ghost_cell,
     $       compute_wall_E_ghost_cell,
     $       compute_wall_W_ghost_cell,
     $       compute_wall_flux_x,
     $       compute_wall_flux_y

        implicit none


        private
        public :: bc_operators_wall_xy


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
        type, extends(bc_operators_default) :: bc_operators_wall_xy

          contains

          procedure, nopass :: apply_bc_on_nodes
          procedure, nopass :: apply_bc_on_nodes_nopt

          procedure, nopass :: apply_bc_on_fluxes
          procedure, nopass :: apply_bc_on_fluxes_nopt

        end type bc_operators_wall_xy


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the nodes along the x and y directions
        !> at the edge of the computational domain
        !
        !> @date
        !> 08_06_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_section
        !> boundary section computed
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> coordinate map along the x-direction
        !
        !>@param y_map
        !> coordinate map along the y-direction
        !
        !>@param nodes_tmp
        !> governing variables at t-dt
        !
        !>@param p_model
        !> physical model
        !
        !>@param nodes
        !> governing variables at t
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes(
     $       bc_section,
     $       t,x_map,y_map,nodes_tmp,
     $       p_model,
     $       nodes)

          implicit none

          integer    , dimension(4)       , intent(in)    :: bc_section
          real(rkind)                     , intent(in)    :: t
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes_tmp
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes


          call apply_bc_on_nodes_nopt(
     $         [bc_section,0],
     $         t,x_map,y_map,nodes_tmp,
     $         p_model,
     $         nodes)

        end subroutine apply_bc_on_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the nodes for a specific boundary section
        !> on the sub-domain
        !
        !> @date
        !> 08_06_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_section
        !> boundary section computed on sub-domain
        !
        !>@param t
        !> time
        !
        !>@param x_map
        !> coordinate map along the x-direction on sub-domain
        !
        !>@param y_map
        !> coordinate map along the y-direction on sub-domain
        !
        !>@param nodes_tmp
        !> governing variables at t-dt on sub-domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param nodes
        !> governing variables at t on sub-domain
        !--------------------------------------------------------------
        subroutine apply_bc_on_nodes_nopt(
     $       bc_section,
     $       t,x_map,y_map,nodes_tmp,
     $       p_model,
     $       nodes)

          implicit none

          integer    , dimension(5)    , intent(in)    :: bc_section
          real(rkind)                  , intent(in)    :: t
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map
          real(rkind), dimension(:,:,:), intent(in)    :: nodes_tmp
          type(pmodel_eq)              , intent(in)    :: p_model
          real(rkind), dimension(:,:,:), intent(inout) :: nodes

          integer(ikind) :: i,j

          real(rkind) :: dx
          real(rkind) :: dy
          real(rkind) :: T_guess
          real(rkind) :: md_guess

          i = p_model%get_eq_nb()

          dx = x_map(2)-x_map(1)
          dy = y_map(2)-y_map(1)

          
          select case(bc_section(1))

            case(N_edge_type)
               call ini_debug_nodes(
     $              nodes,
     $              [bc_section(2),bc_section(4)],
     $              [bc_section(3),bc_section(3)+bc_size-1])

               j=bc_section(3)
               !DEC$ IVDEP
               do i=bc_section(2), bc_section(4)
                  
                  call compute_wall_xy_guess(
     $                 [i,j-1],
     $                 [i,j],
     $                 dx,dy,
     $                 nodes_tmp,
     $                 gradient_x_interior,
     $                 gradient_y_interior,
     $                 T_guess,
     $                 md_guess)

                  call compute_wall_N_ghost_cell(
     $                 i,j,
     $                 t,
     $                 x_map,
     $                 y_map,
     $                 nodes,
     $                 T_guess,
     $                 md_guess,
     $                 sd_interior_type)
                  
               end do


            case(S_edge_type)
               call ini_debug_nodes(
     $              nodes,
     $              [bc_section(2),bc_section(4)],
     $              [bc_section(3),bc_section(3)+bc_size-1])

               j=bc_section(3)+bc_size-1
               !DEC$ IVDEP
               do i=bc_section(2), bc_section(4)

                  call compute_wall_xy_guess(
     $                 [i,j+1],
     $                 [i,j],
     $                 dx,dy,
     $                 nodes_tmp,
     $                 gradient_x_interior,
     $                 gradient_y_interior,
     $                 T_guess,
     $                 md_guess)

                  call compute_wall_S_ghost_cell(
     $                 i,j,
     $                 t,
     $                 x_map,
     $                 y_map,
     $                 nodes,
     $                 T_guess,
     $                 md_guess,
     $                 sd_interior_type)
                  
               end do

               
            case(E_edge_type)
               call ini_debug_nodes(
     $              nodes,
     $              [bc_section(2),bc_section(2)+bc_size-1],
     $              [bc_section(3),bc_section(4)])

               i=bc_section(2)
               !DEC$ IVDEP
               do j=bc_section(3), bc_section(4)
                   
                  call compute_wall_xy_guess(
     $                 [i-1,j],
     $                 [i,j],
     $                 dx,dy,
     $                 nodes_tmp,
     $                 gradient_x_interior,
     $                 gradient_y_interior,
     $                 T_guess,
     $                 md_guess)

                  call compute_wall_E_ghost_cell(
     $                 i,j,
     $                 t,
     $                 x_map,
     $                 y_map,
     $                 nodes,
     $                 T_guess,
     $                 md_guess,
     $                 sd_interior_type)
                  
               end do


            case(W_edge_type)
               call ini_debug_nodes(
     $              nodes,
     $              [bc_section(2),bc_section(2)+bc_size-1],
     $              [bc_section(3),bc_section(4)])

               i=bc_section(2)+bc_size-1
               !DEC$ IVDEP
               do j=bc_section(3), bc_section(4)
                   
                  call compute_wall_xy_guess(
     $                 [i+1,j],
     $                 [i,j],
     $                 dx,dy,
     $                 nodes_tmp,
     $                 gradient_x_interior,
     $                 gradient_y_interior,
     $                 T_guess,
     $                 md_guess)

                  call compute_wall_W_ghost_cell(
     $                 i,j,
     $                 t,
     $                 x_map,
     $                 y_map,
     $                 nodes,
     $                 T_guess,
     $                 md_guess,
     $                 sd_interior_type)

               end do


            case(SW_corner_type)
               call ini_debug_nodes(
     $              nodes,
     $              [bc_section(2),bc_section(2)+bc_size-1],
     $              [bc_section(3),bc_section(3)+bc_size-1])

               j=bc_section(3)+bc_size-1

               i=bc_section(2)
               call compute_wall_xy_guess(
     $              [i,j+1],
     $              [i,j],
     $              dx,dy,
     $              nodes_tmp,
     $              gradient_x_x_oneside_L0,
     $              gradient_y_interior,
     $              T_guess,
     $              md_guess)

               call compute_wall_S_ghost_cell(
     $              i,j,
     $              t,
     $              x_map,
     $              y_map,
     $              nodes,
     $              T_guess,
     $              md_guess,
     $              sd_L0_type)

               i=bc_section(2)+1
               call compute_wall_xy_guess(
     $              [i,j+1],
     $              [i,j],
     $              dx,dy,
     $              nodes_tmp,
     $              gradient_x_x_oneside_L1,
     $              gradient_y_interior,
     $              T_guess,
     $              md_guess)

               call compute_wall_S_ghost_cell(
     $              i,j,
     $              t,
     $              x_map,
     $              y_map,
     $              nodes,
     $              T_guess,
     $              md_guess,
     $              sd_L1_type)
               
               
            case(SE_corner_type)
               call ini_debug_nodes(
     $              nodes,
     $              [bc_section(2),bc_section(2)+bc_size-1],
     $              [bc_section(3),bc_section(3)+bc_size-1])

               j=bc_section(3)+bc_size-1

               i=bc_section(2)
               call compute_wall_xy_guess(
     $              [i,j+1],
     $              [i,j],
     $              dx,dy,
     $              nodes_tmp,
     $              gradient_x_x_oneside_R1,
     $              gradient_y_interior,
     $              T_guess,
     $              md_guess)

               call compute_wall_S_ghost_cell(
     $              i,j,
     $              t,
     $              x_map,
     $              y_map,
     $              nodes,
     $              T_guess,
     $              md_guess,
     $              sd_R1_type)

               i=bc_section(2)+1
               call compute_wall_xy_guess(
     $              [i,j+1],
     $              [i,j],
     $              dx,dy,
     $              nodes_tmp,
     $              gradient_x_x_oneside_R0,
     $              gradient_y_interior,
     $              T_guess,
     $              md_guess)

               call compute_wall_S_ghost_cell(
     $              i,j,
     $              t,
     $              x_map,
     $              y_map,
     $              nodes,
     $              T_guess,
     $              md_guess,
     $              sd_R0_type)


            case(NW_corner_type)
               call ini_debug_nodes(
     $              nodes,
     $              [bc_section(2),bc_section(2)+bc_size-1],
     $              [bc_section(3),bc_section(3)+bc_size-1])

               j=bc_section(3)

               i=bc_section(2)
               call compute_wall_xy_guess(
     $              [i,j-1],
     $              [i,j],
     $              dx,dy,
     $              nodes_tmp,
     $              gradient_x_x_oneside_L0,
     $              gradient_y_interior,
     $              T_guess,
     $              md_guess)

               call compute_wall_N_ghost_cell(
     $              i,j,
     $              t,
     $              x_map,
     $              y_map,
     $              nodes,
     $              T_guess,
     $              md_guess,
     $              sd_L0_type)

               i=bc_section(2)+1
               call compute_wall_xy_guess(
     $              [i,j-1],
     $              [i,j],
     $              dx,dy,
     $              nodes_tmp,
     $              gradient_x_x_oneside_R1,
     $              gradient_y_interior,
     $              T_guess,
     $              md_guess)

               call compute_wall_N_ghost_cell(
     $              i,j,
     $              t,
     $              x_map,
     $              y_map,
     $              nodes,
     $              T_guess,
     $              md_guess,
     $              sd_L1_type)

            case(NE_corner_type)
               call ini_debug_nodes(
     $              nodes,
     $              [bc_section(2),bc_section(2)+bc_size-1],
     $              [bc_section(3),bc_section(3)+bc_size-1])

               j=bc_section(3)

               i=bc_section(2)
               call compute_wall_xy_guess(
     $              [i,j-1],
     $              [i,j],
     $              dx,dy,
     $              nodes_tmp,
     $              gradient_x_x_oneside_R1,
     $              gradient_y_interior,
     $              T_guess,
     $              md_guess)

               call compute_wall_N_ghost_cell(
     $              i,j,
     $              t,
     $              x_map,
     $              y_map,
     $              nodes,
     $              T_guess,
     $              md_guess,
     $              sd_R1_type)

               i=bc_section(2)+1
               call compute_wall_xy_guess(
     $              [i,j-1],
     $              [i,j],
     $              dx,dy,
     $              nodes_tmp,
     $              gradient_x_x_oneside_R0,
     $              gradient_y_interior,
     $              T_guess,
     $              md_guess)

               call compute_wall_N_ghost_cell(
     $              i,j,
     $              t,
     $              x_map,
     $              y_map,
     $              nodes,
     $              T_guess,
     $              md_guess,
     $              sd_R0_type)


            case default
               call error_bc_section_type(
     $              'bc_operators_wall_xy_class',
     $              'apply_bc_on_nodes',
     $              bc_section(1))
          end select

        end subroutine apply_bc_on_nodes_nopt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes along the x directions at the
        !> edge of the computational domain
        !
        !> @date
        !> 08_06_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_section
        !> boundary section computed
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
        subroutine apply_bc_on_fluxes(
     $     bc_section,
     $     t,x_map,y_map,nodes,s,
     $     flux_x,flux_y)

          implicit none

          integer    , dimension(4)         , intent(in)    :: bc_section
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(sd_operators)                , intent(in)    :: s
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y

          call apply_bc_on_fluxes_nopt(
     $         [bc_section,0],
     $         t,x_map,y_map,nodes,s,
     $         flux_x,flux_y)

        end subroutine apply_bc_on_fluxes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes along the x directions at the
        !> edge of the computational domain
        !
        !> @date
        !> 08_06_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_section
        !> boundary section computed
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
        subroutine apply_bc_on_fluxes_nopt(
     $     bc_section,
     $     t,x_map,y_map,nodes,s,
     $     flux_x,flux_y)

          implicit none

          integer    , dimension(5)    , intent(in)    :: bc_section
          real(rkind)                  , intent(in)    :: t
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map
          real(rkind), dimension(:,:,:), intent(in)    :: nodes
          type(sd_operators)           , intent(in)    :: s
          real(rkind), dimension(:,:,:), intent(inout) :: flux_x
          real(rkind), dimension(:,:,:), intent(inout) :: flux_y

          integer(ikind) :: i
          integer(ikind) :: j

          type(sd_operators_x_oneside_L0) :: s_x_oneside_L0
          type(sd_operators_x_oneside_L1) :: s_x_oneside_L1
          type(sd_operators_x_oneside_R1) :: s_x_oneside_R1
          type(sd_operators_x_oneside_R0) :: s_x_oneside_R0


          select case(bc_section(1))

            case(N_edge_type)

               j=bc_section(3)
               do i=bc_section(2),bc_section(4)
                  flux_y(i,j,:) = compute_wall_flux_y(
     $                 nodes,s,i,j,
     $                 t,x_map,y_map,
     $                 right,
     $                 gradient_x_interior)
               end do

            case(S_edge_type)

               j=bc_section(3)+bc_size
               do i=bc_section(2),bc_section(4)
                  flux_y(i,j,:) = compute_wall_flux_y(
     $                 nodes,s,i,j,
     $                 t,x_map,y_map,
     $                 left,
     $                 gradient_x_interior)
               end do

            case(E_edge_type)

               i=bc_section(2)
               do j=bc_section(3),bc_section(4)
                  flux_x(i,j,:) = compute_wall_flux_x(
     $                 nodes,s,i,j,
     $                 t,x_map,y_map,
     $                 right,
     $                 gradient_y_interior)
               end do

            case(W_edge_type)

               i=bc_section(2)+bc_size
               do j=bc_section(3),bc_section(4)
                  flux_x(i,j,:) = compute_wall_flux_x(
     $                 nodes,s,i,j,
     $                 t,x_map,y_map,
     $                 left,
     $                 gradient_y_interior)
               end do

            case(SW_corner_type)

               j=bc_section(3)+bc_size

               i=bc_section(2)
               flux_y(i,j,:) = compute_wall_flux_y(
     $              nodes,s_x_oneside_L0,i,j,
     $              t,x_map,y_map,
     $              left,
     $              gradient_x_x_oneside_L0)

               i=bc_section(2)+1
               flux_y(i,j,:) = compute_wall_flux_y(
     $              nodes,s_x_oneside_L1,i,j,
     $              t,x_map,y_map,
     $              left,
     $              gradient_x_interior)

            case(SE_corner_type)

               j=bc_section(3)+bc_size

               i=bc_section(2)
               flux_y(i,j,:) = compute_wall_flux_y(
     $              nodes,s_x_oneside_R1,i,j,
     $              t,x_map,y_map,
     $              left,
     $              gradient_x_interior)

               i=bc_section(2)+1
               flux_y(i,j,:) = compute_wall_flux_y(
     $              nodes,s_x_oneside_R0,i,j,
     $              t,x_map,y_map,
     $              left,
     $              gradient_x_x_oneside_R0)

            case(NW_corner_type)

               j=bc_section(3)

               i=bc_section(2)
               flux_y(i,j,:) = compute_wall_flux_y(
     $              nodes,s_x_oneside_L0,i,j,
     $              t,x_map,y_map,
     $              right,
     $              gradient_x_x_oneside_L0)

               i=bc_section(2)+1
               flux_y(i,j,:) = compute_wall_flux_y(
     $              nodes,s_x_oneside_L1,i,j,
     $              t,x_map,y_map,
     $              right,
     $              gradient_x_interior)

            case(NE_corner_type)

               j=bc_section(3)

               i=bc_section(2)
               flux_y(i,j,:) = compute_wall_flux_y(
     $              nodes,s_x_oneside_R1,i,j,
     $              t,x_map,y_map,
     $              right,
     $              gradient_x_interior)

               i=bc_section(2)+1
               flux_y(i,j,:) = compute_wall_flux_y(
     $              nodes,s_x_oneside_R0,i,j,
     $              t,x_map,y_map,
     $              right,
     $              gradient_x_x_oneside_R0)

            case default
               print '(''bc_operators_wall_xy_class'')'
               print '(''apply_bc_on_fluxes_nopt'')'
               print '(''bc_section not recognized'')'
               print '(''bc_section: '',I2)', bc_section(1)
               stop ''
               
          end select

        end subroutine apply_bc_on_fluxes_nopt


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the quantities that need to be guessed at
        !> the wall (T_guess and md_guess) when solving for
        !> the mass densities ensuring the wall equilibrium
        !> boundary conditions
        !
        !> @date
        !> 08_06_2015 - initial version - J.L. Desmarais
        !
        !>@param T_guess_index
        !> indices where the temperature is estimated
        !>
        !>@param md_guess_index
        !> indices where the mass-density is estimated
        !>
        !>@param nodes
        !> array with the grid point data
        !>
        !>@param gradient_x_proc
        !> gradient along the x-direction for computing the
        !> temperature
        !>
        !>@param gradient_y_proc
        !> gradient along the y-direction for computing the
        !> temperature
        !>
        !>@param dx
        !> space step along the x-direction
        !>
        !>@param dy
        !> space step along the y-direction
        !>
        !>@param T_guess
        !> temperature estimated at the wall
        !>
        !>@param md_guess
        !> mass density estimated at the wall
        !---------------------------------------------------------------
        subroutine compute_wall_xy_guess(
     $     T_guess_index,
     $     md_guess_index,
     $     dx,dy,
     $     nodes,
     $     gradient_x_proc,
     $     gradient_y_proc,
     $     T_guess,
     $     md_guess)
        
          implicit none

          integer(ikind), dimension(2)    , intent(in)  :: T_guess_index
          integer(ikind), dimension(2)    , intent(in)  :: md_guess_index
          real(rkind)                     , intent(in)  :: dx
          real(rkind)                     , intent(in)  :: dy
          real(rkind)   , dimension(:,:,:), intent(in)  :: nodes
          procedure(gradient_proc)                      :: gradient_x_proc
          procedure(gradient_proc)                      :: gradient_y_proc
          real(rkind)                     , intent(out) :: T_guess
          real(rkind)                     , intent(out) :: md_guess

          !temperature: temperature at the previous step
          T_guess = 3.0d0/(8.0d0*cv_r)*temperature_eff(
     $         nodes,T_guess_index(1),T_guess_index(2),
     $         dx,dy,
     $         gradient_x_proc,
     $         gradient_y_proc)
          
          !mass density: mass density at the previous step
          md_guess = nodes(md_guess_index(1),md_guess_index(2),1)

        end subroutine compute_wall_xy_guess


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the nodes for the debug
        !
        !> @date
        !> 08_06_2015 - initial version - J.L. Desmarais
        !
        !>@param x_borders
        !> [i_min,i_max]
        !
        !>@param y_borders
        !> [j_min,j_max]
        !
        !>@param nodes
        !> nodes initialized
        !--------------------------------------------------------------
        subroutine ini_debug_nodes(nodes,x_borders,y_borders)
         
          implicit none

          real(rkind)   , dimension(:,:,:), intent(inout) :: nodes
          integer(ikind), dimension(2)    , intent(in)    :: x_borders
          integer(ikind), dimension(2)    , intent(in)    :: y_borders


          integer(ikind) :: i,j
          integer        :: k

          if(debug_initialize_bc_nodes) then

             do k=1,ne
                do j=y_borders(1),y_borders(2)
                   do i=x_borders(1),x_borders(2)
                      nodes(i,j,k) = debug_real
                   end do
                end do
             end do

          end if

        end subroutine ini_debug_nodes

      end module bc_operators_wall_xy_class
