      !> @file
      !> generalized boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> object that encapulates the different types of boundary
      !> conditions implemented as relatives of bc_operators_abstract.
      !> The boundary conditions are then applied by choosing according
      !> to the user's choice
      !>
      !> @date
      !> 07_06_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module bc_operators_gen_class

        use bc_operators_hedstrom_xy_class, only :
     $       bc_operators_hedstrom_xy

        use bc_operators_reflection_x_class, only :
     $       bc_operators_reflection_x

        use bc_operators_reflection_y_class, only :
     $       bc_operators_reflection_y

        use bc_operators_wall_xy_class, only :
     $       bc_operators_wall_xy

        use errors_module, only :
     $       error_bc_choice

        use parameters_constant, only :
     $       bc_nodes_choice,
     $       bc_fluxes_choice,
     $       bc_timedev_choice,
     $       bc_flux_and_node_choice,
     $       
     $       reflection_x_choice,
     $       reflection_y_choice,
     $       wall_choice,
     $       hedstrom_choice,
     $       
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       bc_N_choice,
     $       bc_S_choice,
     $       bc_E_choice,
     $       bc_W_choice,
     $       bc_NW_choice,
     $       bc_NE_choice,
     $       bc_SW_choice,
     $       bc_SE_choice,
     $       bc_N_type_choice,
     $       bc_S_type_choice,
     $       bc_E_type_choice,
     $       bc_W_type_choice,
     $       bc_NW_type_choice,
     $       bc_NE_type_choice,
     $       bc_SW_type_choice,
     $       bc_SE_type_choice,
     $       x_min,x_max,
     $       y_min,y_max,
     $       bc_order1,
     $       bc_order2,
     $       bc_order3,
     $       bc_order4,
     $       bc_order5,
     $       bc_order6,
     $       bc_order7,
     $       bc_order8

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        implicit none


        integer, dimension(4), parameter :: SW_corner_section = [SW_corner_type, 1           , 1           , 0         ]
        integer, dimension(4), parameter :: S_edge_section    = [S_edge_type   , bc_size+1   , 1           , nx-bc_size]
        integer, dimension(4), parameter :: SE_corner_section = [SE_corner_type, nx-bc_size+1, 1           , 0         ]
        integer, dimension(4), parameter :: W_edge_section    = [W_edge_type   , 1           , bc_size+1   , ny-bc_size]
        integer, dimension(4), parameter :: E_edge_section    = [E_edge_type   , nx-bc_size+1, bc_size+1   , ny-bc_size]
        integer, dimension(4), parameter :: NW_corner_section = [NW_corner_type, 1           , ny-bc_size+1, 0         ]
        integer, dimension(4), parameter :: N_edge_section    = [N_edge_type   , bc_size+1   , ny-bc_size+1, nx-bc_size]
        integer, dimension(4), parameter :: NE_corner_section = [NE_corner_type, nx-bc_size+1, ny-bc_size+1, 0         ]

        logical :: apply_bc_on_nodes_N  = (bc_N_type_choice.eq.bc_nodes_choice).or.
     $                                    (bc_N_type_choice.eq.bc_flux_and_node_choice)
        logical :: apply_bc_on_nodes_S  = (bc_S_type_choice.eq.bc_nodes_choice).or.
     $                                    (bc_S_type_choice.eq.bc_flux_and_node_choice)
        logical :: apply_bc_on_nodes_E  = (bc_E_type_choice.eq.bc_nodes_choice).or.
     $                                    (bc_E_type_choice.eq.bc_flux_and_node_choice)
        logical :: apply_bc_on_nodes_W  = (bc_W_type_choice.eq.bc_nodes_choice).or.
     $                                    (bc_W_type_choice.eq.bc_flux_and_node_choice)

        logical :: apply_bc_on_nodes_NW = (bc_NW_type_choice.eq.bc_nodes_choice).or.
     $                                    (bc_NW_type_choice.eq.bc_flux_and_node_choice)
        logical :: apply_bc_on_nodes_NE = (bc_NE_type_choice.eq.bc_nodes_choice).or.
     $                                    (bc_NE_type_choice.eq.bc_flux_and_node_choice)
        logical :: apply_bc_on_nodes_SE = (bc_SE_type_choice.eq.bc_nodes_choice).or.
     $                                    (bc_SE_type_choice.eq.bc_flux_and_node_choice)
        logical :: apply_bc_on_nodes_SW = (bc_SW_type_choice.eq.bc_nodes_choice).or.
     $                                    (bc_SW_type_choice.eq.bc_flux_and_node_choice)

        logical :: apply_bc_on_fluxes_N = (bc_N_type_choice.eq.bc_flux_and_node_choice).or.
     $                                    (bc_N_type_choice.eq.bc_fluxes_choice)
        logical :: apply_bc_on_fluxes_S = (bc_S_type_choice.eq.bc_flux_and_node_choice).or.
     $                                    (bc_S_type_choice.eq.bc_fluxes_choice)
        logical :: apply_bc_on_fluxes_E = (bc_E_type_choice.eq.bc_flux_and_node_choice).or.
     $                                    (bc_E_type_choice.eq.bc_fluxes_choice)
        logical :: apply_bc_on_fluxes_W = (bc_W_type_choice.eq.bc_flux_and_node_choice).or.
     $                                    (bc_W_type_choice.eq.bc_fluxes_choice)

        logical :: apply_bc_on_fluxes_NW = (bc_NW_type_choice.eq.bc_flux_and_node_choice).or.
     $                                     (bc_NW_type_choice.eq.bc_fluxes_choice)
        logical :: apply_bc_on_fluxes_NE = (bc_NE_type_choice.eq.bc_flux_and_node_choice).or.
     $                                     (bc_NE_type_choice.eq.bc_fluxes_choice)
        logical :: apply_bc_on_fluxes_SE = (bc_SE_type_choice.eq.bc_flux_and_node_choice).or.
     $                                     (bc_SE_type_choice.eq.bc_fluxes_choice)
        logical :: apply_bc_on_fluxes_SW = (bc_SW_type_choice.eq.bc_flux_and_node_choice).or.
     $                                     (bc_SW_type_choice.eq.bc_fluxes_choice)

        logical :: apply_bc_on_timedev_N  = bc_N_type_choice.eq.bc_timedev_choice
        logical :: apply_bc_on_timedev_S  = bc_S_type_choice.eq.bc_timedev_choice
        logical :: apply_bc_on_timedev_E  = bc_E_type_choice.eq.bc_timedev_choice
        logical :: apply_bc_on_timedev_W  = bc_W_type_choice.eq.bc_timedev_choice

        logical :: apply_bc_on_timedev_NW = bc_NW_type_choice.eq.bc_timedev_choice
        logical :: apply_bc_on_timedev_NE = bc_NE_type_choice.eq.bc_timedev_choice
        logical :: apply_bc_on_timedev_SE = bc_SE_type_choice.eq.bc_timedev_choice
        logical :: apply_bc_on_timedev_SW = bc_SW_type_choice.eq.bc_timedev_choice


        private
        public :: bc_operators_gen


        !> @class bc_operators_gen
        !> object that encapulates the different types of boundary
        !> conditions implemented as relatives of bc_operators_abstract.
        !> The boundary conditions are then applied by choosing according
        !> to the user's choice
        !
        !>@param apply_bc_on_nodes
        !> interface to apply the boundary conditions along
        !> the x and y directions on the nodes at the edge of
        !> the computational domain
        !
        !>@param apply_bc_on_fluxes
        !> interface to apply the boundary conditions along
        !> the x and y directions on the fluxes at the edge of
        !> the computational domain
        !
        !>@param apply_bc_on_timedev
        !> interface to apply the boundary conditions along
        !> the x and y directions on the time derivatives at
        !> the edge of the computational domain
        !
        !>@param apply_bc_on_nodes_nopt
        !> interface to apply the boundary conditions along
        !> the x and y directions on the nodes at the edge of
        !> the computational domain for an adaptive computational
        !> domain
        !
        !>@param apply_bc_on_fluxes_nopt
        !> interface to apply the boundary conditions along
        !> the x and y directions on the fluxes at the edge of
        !> the computational domain for an adaptive computational
        !> domain
        !
        !>@param apply_bc_on_timedev_nopt
        !> interface to apply the boundary conditions along
        !> the x and y directions on the time derivatives at
        !> the edge of the computational domain for an adaptive
        !> computational domain
        !---------------------------------------------------------------
        type :: bc_operators_gen

          contains

          procedure, nopass :: apply_bc_on_nodes
          procedure, nopass :: apply_bc_on_fluxes
          procedure, nopass :: apply_bc_on_timedev

          procedure, nopass :: apply_bc_on_fluxes_nopt
          procedure, nopass :: apply_bc_on_timedev_nopt


        end type bc_operators_gen


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
        !> 07_06_2015 - initial version - J.L. Desmarais
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
        !>@param nodes
        !> governing variables at t
        !-------------------------------------------------------------
        subroutine apply_bc_on_nodes(
     $       t,x_map,y_map,
     $       nodes_tmp,p_model,
     $       nodes)
        
          implicit none
        
          real(rkind)                     , intent(in)    :: t
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes_tmp
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes

          integer, dimension(8) :: bc_nodes_order
          integer               :: k


          ! choose order in which the boundary conditions are applied
          call choose_bc_nodes_order(bc_nodes_order)


          ! apply the boundary conditions for each layer
          do k=1, size(bc_nodes_order,1)

             select case(bc_nodes_order(k))
             
               ! N layer
               case(N_edge_type)

                  if(apply_bc_on_nodes_N) then

                     call apply_bc_on_nodes_with_bc_operator(
     $                    N_edge_section,
     $                    bc_N_choice,
     $                    t,x_map,y_map,
     $                    nodes_tmp,p_model,
     $                    nodes)
                     
                  end if

               ! S layer
               case(S_edge_type)

                  if(apply_bc_on_nodes_S) then

                     call apply_bc_on_nodes_with_bc_operator(
     $                    S_edge_section,
     $                    bc_S_choice,
     $                    t,x_map,y_map,
     $                    nodes_tmp,p_model,
     $                    nodes)

                  end if

               ! E layer
               case(E_edge_type)

                  if(apply_bc_on_nodes_E) then

                     call apply_bc_on_nodes_with_bc_operator(
     $                    E_edge_section,
     $                    bc_E_choice,
     $                    t,x_map,y_map,
     $                    nodes_tmp,p_model,
     $                    nodes)

                  end if

               ! W layer
               case(W_edge_type)

                  if(apply_bc_on_nodes_W) then

                     call apply_bc_on_nodes_with_bc_operator(
     $                    W_edge_section,
     $                    bc_W_choice,
     $                    t,x_map,y_map,
     $                    nodes_tmp,p_model,
     $                    nodes)

                  end if

               ! NW corner
               case(NW_corner_type)

                  if(apply_bc_on_nodes_NW) then

                     call apply_bc_on_nodes_with_bc_operator(
     $                    NW_corner_section,
     $                    bc_NW_choice,
     $                    t,x_map,y_map,
     $                    nodes_tmp,p_model,
     $                    nodes)

                  end if
                     
               ! NE corner
               case(NE_corner_type)

                  if(apply_bc_on_nodes_NE) then

                     call apply_bc_on_nodes_with_bc_operator(
     $                    NE_corner_section,
     $                    bc_NE_choice,
     $                    t,x_map,y_map,
     $                    nodes_tmp,p_model,
     $                    nodes)
                     
                  end if

               ! SW corner
               case(SW_corner_type)

                  if(apply_bc_on_nodes_SW) then

                     call apply_bc_on_nodes_with_bc_operator(
     $                    SW_corner_section,
     $                    bc_SW_choice,
     $                    t,x_map,y_map,
     $                    nodes_tmp,p_model,
     $                    nodes)

                  end if

               ! SE corner
               case(SE_corner_type)

                  if(apply_bc_on_nodes_SE) then

                     call apply_bc_on_nodes_with_bc_operator(
     $                    SE_corner_section,
     $                    bc_SE_choice,
     $                    t,x_map,y_map,
     $                    nodes_tmp,p_model,
     $                    nodes)
                     
                  end if

             end select

          end do

        end subroutine apply_bc_on_nodes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the nodes using the bc_operators
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
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
        !>@param nodes
        !> governing variables at t
        !-------------------------------------------------------------
        subroutine apply_bc_on_nodes_with_bc_operator(
     $     bc_section, bc_choice,
     $     t,x_map,y_map,
     $     nodes_tmp,p_model,
     $     nodes)
        
          implicit none
        
          integer    , dimension(4)       , intent(in)    :: bc_section
          integer                         , intent(in)    :: bc_choice
          real(rkind)                     , intent(in)    :: t
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes_tmp
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes


          type(bc_operators_reflection_x) :: bc_reflection_x
          type(bc_operators_reflection_y) :: bc_reflection_y
          type(bc_operators_wall_xy)      :: bc_wall_xy
          

          select case(bc_choice)

            case(reflection_x_choice)

               call bc_reflection_x%apply_bc_on_nodes(
     $              bc_section,
     $              t,x_map,y_map,
     $              nodes_tmp,
     $              p_model,
     $              nodes)


            case(reflection_y_choice)

               call bc_reflection_y%apply_bc_on_nodes(
     $              bc_section,
     $              t,x_map,y_map,
     $              nodes_tmp,
     $              p_model,
     $              nodes)


            case(wall_choice)

               call bc_wall_xy%apply_bc_on_nodes(
     $              bc_section,
     $              t,x_map,y_map,
     $              nodes_tmp,
     $              p_model,
     $              nodes)

          
            case default

               call error_bc_choice(
     $              'bc_operators_gen_class',
     $              'apply_bc_on_nodes_with_bc_operator',
     $              bc_choice)

               
          end select

        end subroutine apply_bc_on_nodes_with_bc_operator


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes along the x directions at the
        !> edge of the computational domain
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
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
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !-------------------------------------------------------------
        subroutine apply_bc_on_fluxes(
     $     t,x_map,y_map,nodes,s,
     $     flux_x,flux_y)
        
          implicit none
        
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(sd_operators)                , intent(in)    :: s
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
        

          ! SW corner
          if(apply_bc_on_fluxes_SW) then             
             call apply_bc_on_fluxes_with_bc_operator(
     $            SW_corner_section,
     $            bc_SW_choice,
     $            t,x_map,y_map,nodes,
     $            s,flux_x,flux_y)
          end if


          ! S edge
          if(apply_bc_on_fluxes_S) then
             call apply_bc_on_fluxes_with_bc_operator(
     $            S_edge_section,
     $            bc_S_choice,
     $            t,x_map,y_map,nodes,
     $            s,flux_x,flux_y)
          end if


          ! SE corner
          if(apply_bc_on_fluxes_SE) then
             call apply_bc_on_fluxes_with_bc_operator(
     $            SE_corner_section,
     $            bc_SE_choice,
     $            t,x_map,y_map,nodes,
     $            s,flux_x,flux_y)
          end if


          ! W layer
          if(apply_bc_on_fluxes_W) then
             call apply_bc_on_fluxes_with_bc_operator(
     $            W_edge_section,
     $            bc_W_choice,
     $            t,x_map,y_map,nodes,
     $            s,flux_x,flux_y)
          end if


          ! E layer
          if(apply_bc_on_fluxes_E) then
             call apply_bc_on_fluxes_with_bc_operator(
     $            E_edge_section,
     $            bc_E_choice,
     $            t,x_map,y_map,nodes,
     $            s,flux_x,flux_y)
          end if


          ! NW corner
          if(apply_bc_on_fluxes_NW) then
             call apply_bc_on_fluxes_with_bc_operator(
     $            NW_corner_section,
     $            bc_NW_choice,
     $            t,x_map,y_map,nodes,
     $            s,flux_x,flux_y)
          end if
             
          ! N layer
          if(apply_bc_on_fluxes_N) then
             call apply_bc_on_fluxes_with_bc_operator(
     $            N_edge_section,
     $            bc_N_choice,
     $            t,x_map,y_map,nodes,
     $            s,flux_x,flux_y)
          end if

          ! NE corner
          if(apply_bc_on_fluxes_NE) then
             call apply_bc_on_fluxes_with_bc_operator(
     $            NE_corner_section,
     $            bc_NE_choice,
     $            t,x_map,y_map,nodes,
     $            s,flux_x,flux_y)             
          end if

        end subroutine apply_bc_on_fluxes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes using the bc_operators
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
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
        !>@param nodes
        !> governing variables at t
        !-------------------------------------------------------------
        subroutine apply_bc_on_fluxes_with_bc_operator(
     $     bc_section, bc_choice,
     $     t,x_map,y_map,nodes,
     $     s,flux_x,flux_y)
        
          implicit none
        
          integer    , dimension(4)         , intent(in)    :: bc_section
          integer                           , intent(in)    :: bc_choice
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(sd_operators)                , intent(in)    :: s
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y


          type(bc_operators_wall_xy) :: bc_wall_xy

          
          select case(bc_choice)

            case(wall_choice)

               ! there is here that the fluxes should be splitted
               ! into viscid and invisicd fluxes for the Yoo-Lodato
               ! b.c. using optional arguments in apply_bc_on_fluxes
               ! for the wall b.c.
               call bc_wall_xy%apply_bc_on_fluxes(
     $              bc_section,
     $              t,x_map,y_map,
     $              nodes,s,flux_x,flux_y)
          
            case default
               call error_bc_choice(
     $              'bc_operators_gen_class',
     $              'apply_bc_on_fluxes_with_bc_operator',
     $              bc_choice)
               
          end select

        end subroutine apply_bc_on_fluxes_with_bc_operator


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes along the x directions at the
        !> edge of the computational domain
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_sections
        !> boundary sections computed
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
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !-------------------------------------------------------------
        subroutine apply_bc_on_fluxes_nopt(
     $     bc_sections,
     $     t,x_map,y_map,nodes,s,
     $     flux_x,flux_y)
        
          implicit none
        
          integer(ikind), dimension(:,:), allocatable, intent(in)    :: bc_sections
          real(rkind)                                , intent(in)    :: t
          real(rkind)   , dimension(nx)              , intent(in)    :: x_map
          real(rkind)   , dimension(ny)              , intent(in)    :: y_map
          real(rkind)   , dimension(nx,ny,ne)        , intent(in)    :: nodes
          type(sd_operators)                         , intent(in)    :: s
          real(rkind)   , dimension(nx+1,ny,ne)      , intent(inout) :: flux_x
          real(rkind)   , dimension(nx,ny+1,ne)      , intent(inout) :: flux_y
        

          integer :: k
          integer :: bc_cardinal_choice
          integer :: bc_cardinal_type_choice

          
          if(allocated(bc_sections)) then
             do k=1, size(bc_sections,2)

                !1) determine which cardinal boundary operator
                !   should be used to compute the bc_section
                !   (N,S,E,W) ?
                call determine_cardinal_bc_operator(
     $               x_map(bc_sections(2,k)),
     $               y_map(bc_sections(3,k)),
     $               bc_cardinal_choice,
     $               bc_cardinal_type_choice)

                !2) apply the boundary operator corresponding to
                !   the cardinal coordinate to the bc_section
                if((bc_cardinal_type_choice.eq.bc_fluxes_choice).or.
     $             (bc_cardinal_type_choice.eq.bc_flux_and_node_choice)) then

                   call apply_bc_on_fluxes_nopt_with_bc_operator(
     $                  bc_sections(:,k),
     $                  bc_cardinal_choice,
     $                  t,x_map,y_map,nodes,
     $                  s,flux_x,flux_y)

                end if

             end do
          end if

        end subroutine apply_bc_on_fluxes_nopt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes for a specific boundary section
        !> using the bc_operators on a sub-domain
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_section
        !> boundary section
        !
        !>@param bc_choice
        !> bc_operator chosen
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
        !>@param nodes
        !> governing variables
        !
        !>@param s
        !> space discretization operator
        !
        !>@param flux_x
        !> fluxes in the x-direction
        !
        !>@param flux_y
        !> fluxes in the y-direction
        !-------------------------------------------------------------
        subroutine apply_bc_on_fluxes_nopt_with_bc_operator(
     $     bc_section, bc_choice,
     $     t,x_map,y_map,nodes,
     $     s,flux_x,flux_y)
        
          implicit none
        
          integer    , dimension(5)         , intent(in)    :: bc_section
          integer                           , intent(in)    :: bc_choice
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(sd_operators)                , intent(in)    :: s
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y


          type(bc_operators_wall_xy) :: bc_wall_xy

          
          select case(bc_choice)

            case(wall_choice)

               ! there is here that the fluxes should be splitted
               ! into viscid and invisicd fluxes for the Yoo-Lodato
               ! b.c. using optional arguments in apply_bc_on_fluxes
               ! for the wall b.c.
               call bc_wall_xy%apply_bc_on_fluxes_nopt(
     $              bc_section,
     $              t,x_map,y_map,nodes,s,
     $              flux_x,flux_y)
          
            case default
               call error_bc_choice(
     $              'bc_operators_gen_class',
     $              'apply_bc_on_fluxes_with_bc_operator',
     $              bc_choice)
               
          end select

        end subroutine apply_bc_on_fluxes_nopt_with_bc_operator


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes along the x directions at the
        !> edge of the computational domain
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> abstract boundary conditions
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
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev(
     $     t,x_map,y_map,nodes,
     $     p_model,
     $     flux_x,flux_y,timedev)
        
          implicit none
        
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev
        

          ! SW corner
          if(apply_bc_on_timedev_SW) then
             call apply_bc_on_timedev_with_bc_operator(
     $            SW_corner_section,
     $            bc_SW_choice,
     $            t,x_map,y_map,nodes,
     $            p_model,
     $            flux_x, flux_y, timedev)
          end if

          ! S layer
          if(apply_bc_on_timedev_S) then
             call apply_bc_on_timedev_with_bc_operator(
     $            S_edge_section,
     $            bc_S_choice,
     $            t,x_map,y_map,nodes,
     $            p_model,
     $            flux_x, flux_y, timedev)
          end if

          ! SE corner
          if(apply_bc_on_timedev_SE) then
             call apply_bc_on_timedev_with_bc_operator(
     $            SE_corner_section,
     $            bc_SE_choice,
     $            t,x_map,y_map,nodes,
     $            p_model,
     $            flux_x, flux_y, timedev)
          end if

          ! W layer
          if(apply_bc_on_timedev_W) then
             call apply_bc_on_timedev_with_bc_operator(
     $            W_edge_section,
     $            bc_W_choice,
     $            t,x_map,y_map,nodes,
     $            p_model,
     $            flux_x, flux_y, timedev)
          end if

          ! E layer
          if(apply_bc_on_timedev_E) then
             call apply_bc_on_timedev_with_bc_operator(
     $            E_edge_section,
     $            bc_E_choice,
     $            t,x_map,y_map,nodes,
     $            p_model,
     $            flux_x, flux_y, timedev)
          end if


          ! NW corner
          if(apply_bc_on_timedev_NW) then
             call apply_bc_on_timedev_with_bc_operator(
     $            NW_corner_section,
     $            bc_NW_choice,
     $            t,x_map,y_map,nodes,
     $            p_model,
     $            flux_x, flux_y, timedev)
          end if

          ! N layer
          if(apply_bc_on_timedev_N) then
             call apply_bc_on_timedev_with_bc_operator(
     $            N_edge_section,
     $            bc_N_choice,
     $            t,x_map,y_map,nodes,
     $            p_model,
     $            flux_x, flux_y, timedev)
          end if
             
          ! NE corner
          if(apply_bc_on_timedev_NE) then
             call apply_bc_on_timedev_with_bc_operator(
     $            NE_corner_section,
     $            bc_NE_choice,
     $            t,x_map,y_map,nodes,
     $            p_model,
     $            flux_x, flux_y, timedev)
          end if

        end subroutine apply_bc_on_timedev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the fluxes using the bc_operators
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
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
        !>@param nodes
        !> governing variables at t
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev_with_bc_operator(
     $     bc_section, bc_choice,
     $     t,x_map,y_map,nodes,
     $     p_model,
     $     flux_x, flux_y, timedev)
        
          implicit none
        
          integer    , dimension(4)         , intent(in)    :: bc_section
          integer                           , intent(in)    :: bc_choice
          real(rkind)                       , intent(in)    :: t
          real(rkind), dimension(nx)        , intent(in)    :: x_map
          real(rkind), dimension(ny)        , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne)  , intent(in)    :: nodes
          type(pmodel_eq)                   , intent(in)    :: p_model
          real(rkind), dimension(nx+1,ny,ne), intent(inout) :: flux_x
          real(rkind), dimension(nx,ny+1,ne), intent(inout) :: flux_y
          real(rkind), dimension(nx,ny,ne)  , intent(inout) :: timedev


          type(bc_operators_hedstrom_xy) :: bc_hedstrom_xy

          
          select case(bc_choice)

            case(hedstrom_choice)

               ! there is here that the fluxes should be splitted
               ! into viscid and invisicd fluxes for the Yoo-Lodato
               ! b.c. using optional arguments in apply_bc_on_fluxes
               ! for the wall b.c.
               call bc_hedstrom_xy%apply_bc_on_timedev(
     $              bc_section,
     $              t,x_map,y_map,nodes,
     $              p_model,
     $              flux_x,flux_y,timedev)
          
            case default
               call error_bc_choice(
     $              'bc_operators_gen_class',
     $              'apply_bc_on_timedev_with_bc_operator',
     $              bc_choice)
               
          end select

        end subroutine apply_bc_on_timedev_with_bc_operator


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the time derivatives for specific boundary
        !> sections on the sub-domain
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_sections
        !> sections computed for the boundary conditions
        !> of the sub-domain
        !
        !>@param bf_alignment
        !> relative alignment of the sub-domain compared
        !> to the interior domain
        !
        !>@param bf_grdpts_id
        !> identification of the grid-points in the
        !> sub-domain
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
        !>@param interior_nodes
        !> grid-points of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param flux_x
        !> fluxes along the x-direction
        !
        !>@param flux_y
        !> fluxes along the y-direction
        !
        !>@param timedev
        !> tim derivatives of the sub-domain
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev_nopt(
     $     bc_sections,
     $     bf_alignment, bf_grdpts_id,
     $     t,x_map,y_map,nodes,
     $     interior_nodes,
     $     p_model,flux_x,flux_y,
     $     timedev)
        
          implicit none
        
          integer(ikind), dimension(:,:), allocatable, intent(in)    :: bc_sections
          integer(ikind), dimension(2,2)             , intent(in)    :: bf_alignment
          integer       , dimension(:,:)             , intent(in)    :: bf_grdpts_id
          real(rkind)                                , intent(in)    :: t
          real(rkind)   , dimension(:)               , intent(in)    :: x_map
          real(rkind)   , dimension(:)               , intent(in)    :: y_map
          real(rkind)   , dimension(:,:,:)           , intent(in)    :: nodes
          real(rkind)   , dimension(nx,ny,ne)        , intent(in)    :: interior_nodes
          type(pmodel_eq)                            , intent(in)    :: p_model
          real(rkind)   , dimension(:,:,:)           , intent(inout) :: flux_x
          real(rkind)   , dimension(:,:,:)           , intent(inout) :: flux_y
          real(rkind)   , dimension(:,:,:)           , intent(inout) :: timedev
        
          
          integer :: k
          integer :: bc_cardinal_choice
          integer :: bc_cardinal_type_choice

          
          if(allocated(bc_sections)) then
             do k=1, size(bc_sections,2)

                !1) determine which cardinal boundary operator
                !   should be used to compute the bc_section
                !   (N,S,E,W) ?
                call determine_cardinal_bc_operator(
     $               x_map(bc_sections(2,k)),
     $               y_map(bc_sections(3,k)),
     $               bc_cardinal_choice,
     $               bc_cardinal_type_choice)

                !2) apply the boundary operator corresponding to
                !   the cardinal coordinate to the bc_section
                if(bc_cardinal_type_choice.eq.bc_timedev_choice) then

                   call apply_bc_on_timedev_nopt_with_bc_operator(
     $                  bc_sections(:,k), bc_cardinal_choice,
     $                  bf_alignment, bf_grdpts_id,
     $                  t,x_map,y_map,nodes,
     $                  interior_nodes,
     $                  p_model,
     $                  flux_x, flux_y, timedev)

                end if

             end do
          end if

        end subroutine apply_bc_on_timedev_nopt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine applying the boundary conditions
        !> on the time derivatives for a sub-domain
        !> using the bc_operators
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
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
        !>@param nodes
        !> governing variables at t
        !-------------------------------------------------------------
        subroutine apply_bc_on_timedev_nopt_with_bc_operator(
     $     bc_section, bc_choice,
     $     bf_alignment, bf_grdpts_id,
     $     t,x_map,y_map,nodes,
     $     interior_nodes,
     $     p_model,
     $     flux_x, flux_y, timedev)
        
          implicit none
        
          integer    , dimension(5)       , intent(in)    :: bc_section
          integer                         , intent(in)    :: bc_choice
          integer(ikind), dimension(2,2)  , intent(in)    :: bf_alignment
          integer       , dimension(:,:)  , intent(in)    :: bf_grdpts_id
          real(rkind)                     , intent(in)    :: t
          real(rkind), dimension(:)       , intent(in)    :: x_map
          real(rkind), dimension(:)       , intent(in)    :: y_map
          real(rkind), dimension(:,:,:)   , intent(in)    :: nodes
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind), dimension(:,:,:)   , intent(inout) :: flux_x
          real(rkind), dimension(:,:,:)   , intent(inout) :: flux_y
          real(rkind), dimension(:,:,:)   , intent(inout) :: timedev


          type(bc_operators_hedstrom_xy) :: bc_hedstrom_xy

          
          select case(bc_choice)

            case(hedstrom_choice)

               ! there is here that the fluxes should be splitted
               ! into viscid and invisicd fluxes for the Yoo-Lodato
               ! b.c. using optional arguments in apply_bc_on_fluxes
               ! for the wall b.c.
               call bc_hedstrom_xy%apply_bc_on_timedev_nopt(
     $              bc_section,
     $              bf_alignment, bf_grdpts_id,
     $              t,x_map,y_map,nodes,
     $              interior_nodes,
     $              p_model,
     $              flux_x,flux_y,
     $              timedev)
          
            case default
               call error_bc_choice(
     $              'bc_operators_gen_class',
     $              'apply_bc_on_timedev_nopt_with_bc_operator',
     $              bc_choice)
               
          end select

        end subroutine apply_bc_on_timedev_nopt_with_bc_operator


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the order in which the boundary conditions
        !> are applied
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_nodes_order
        !> order of the boundary conditions
        !--------------------------------------------------------------
        subroutine choose_bc_nodes_order(bc_nodes_order)

          implicit none

          integer, dimension(8), intent(out) :: bc_nodes_order

          bc_nodes_order = [
     $         bc_order1,
     $         bc_order2,
     $         bc_order3,
     $         bc_order4,
     $         bc_order5,
     $         bc_order6,
     $         bc_order7,
     $         bc_order8]

        end subroutine choose_bc_nodes_order


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the order in which the boundary conditions
        !> are applied
        !
        !> @date
        !> 07_06_2015 - initial version - J.L. Desmarais
        !
        !>@param bc_nodes_order
        !> order of the boundary conditions
        !--------------------------------------------------------------
        subroutine determine_cardinal_bc_operator(
     $     bc_section_x_min,
     $     bc_section_y_min,
     $     bc_cardinal_choice,
     $     bc_cardinal_type_choice)

          implicit none

          real(rkind), intent(in)  :: bc_section_x_min
          real(rkind), intent(in)  :: bc_section_y_min
          integer    , intent(out) :: bc_cardinal_choice
          integer    , intent(out) :: bc_cardinal_type_choice


          if(bc_section_y_min.lt.y_min) then
             if(bc_section_x_min.lt.x_min) then
                bc_cardinal_choice      = bc_SW_choice
                bc_cardinal_type_choice = bc_SW_type_choice

             else
                if(bc_section_x_min.gt.x_max) then
                   bc_cardinal_choice      = bc_SE_choice
                   bc_cardinal_type_choice = bc_SE_type_choice

                else
                   bc_cardinal_choice      = bc_S_choice
                   bc_cardinal_type_choice = bc_S_type_choice

                end if
             end if

          else
             if(bc_section_y_min.gt.y_max) then
                if(bc_section_x_min.lt.x_min) then
                   bc_cardinal_choice      = bc_NW_choice
                   bc_cardinal_type_choice = bc_NW_type_choice

                else
                   if(bc_section_x_min.gt.x_max) then
                      bc_cardinal_choice      = bc_NE_choice
                      bc_cardinal_type_choice = bc_NE_type_choice
                      
                   else
                      bc_cardinal_choice      = bc_N_choice
                      bc_cardinal_type_choice = bc_N_type_choice
                      
                   end if
                end if

             else
                if(bc_section_x_min.lt.x_min) then
                   bc_cardinal_choice      = bc_W_choice
                   bc_cardinal_type_choice = bc_W_type_choice

                else
                   if(bc_section_x_min.gt.x_max) then
                      bc_cardinal_choice      = bc_E_choice
                      bc_cardinal_type_choice = bc_E_type_choice
                      
                   else
                      print '(''bc_operators_gen_class'')'
                      print '(''determine_cardinal_bc_operator'')'
                      print '(''problem determine cardinal bc_operator'')'
                      print '(''x_min: '', F6.3)', x_min
                      print '(''x_max: '', F6.3)', x_max
                      print '(''y_min: '', F6.3)', y_min
                      print '(''y_max: '', F6.3)', y_max
                      print '(''bc_section_x_min: '', F6.3)', bc_section_x_min
                      print '(''bc_section_y_min: '', F6.3)', bc_section_y_min
                      stop ''                      
                   end if
                end if
             end if
          end if

        end subroutine determine_cardinal_bc_operator        

      end module bc_operators_gen_class
