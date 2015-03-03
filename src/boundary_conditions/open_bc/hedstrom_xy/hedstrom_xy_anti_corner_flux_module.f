      !> @file
      !> module implementing the subroutine computing the time
      !> derivatives for the anti-corner
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> computation of time derivatives for anti-corner using
      !> hedstrom_xy b.c. and (x,y) fluxes
      !
      !> @date
      !> 01_08_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module hedstrom_xy_anti_corner_flux_module

        use bc_operators_nopt_module, only :
     $       compute_edge_N,
     $       compute_edge_S,
     $       compute_edge_E,
     $       compute_edge_W,
     $       are_grdpts_needed_for_flux_x,
     $       are_grdpts_needed_for_flux_y,
     $       combine_grdpts_to_compute_fluxes

        use bf_layer_bc_sections_overlap_module, only :
     $       determine_corner_or_anti_corner_grdpts_computed

        use bf_layer_extract_module, only :
     $       get_bf_layer_match_table

        use hedstrom_xy_module, only :
     $       compute_timedev_x_edge_local,
     $       compute_timedev_y_edge_local,
     $       compute_timedev_corner_local

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right
        
        use parameters_bf_layer, only : 
     $       SE_edge_type,
     $       SW_edge_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       NE_corner_type,
     $       NW_corner_type

        use parameters_constant, only :
     $       bc_timedev_choice

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_R0

        use sd_operators_x_oneside_L1_class, only :
     $       sd_operators_x_oneside_L1

        use sd_operators_x_oneside_R1_class, only :
     $       sd_operators_x_oneside_R1

        use sd_operators_y_oneside_L1_class, only :
     $       sd_operators_y_oneside_L1

        use sd_operators_y_oneside_R1_class, only :
     $       sd_operators_y_oneside_R1


        implicit none


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> an xy edge: NE_edge, NW_edge, SE_edge, SW_edge
        !> using the fluxes
        !
        !> @date
        !> 26_01_2014 - initial version - J.L. Desmarais
        !
        !> @param p_model
        !> object encapsulating the physical model
        !
        !> @param t
        !> simulation time for boundary conditions depending
        !> on time
        !
        !> @param interior_nodes
        !> nodes from the interior computational domain
        !
        !> @param bf_alignment
        !> relative position of the buffer layer compared
        !> to the interiro domain
        !
        !> @param nodes
        !> object encapsulating the main variables
        !
        !> @param x_map
        !> coordinates along the x-direction
        !
        !> @param y_map
        !> coordinates along the y-direction
        !
        !> @param flux_x
        !> fluxes along the x-direction
        !
        !> @param flux_y
        !> fluxes along the y-direction
        !
        !> @param s_x_L1
        !> space discretization operator with one grid point
        !> on the left side (x-direction)
        !     
        !> @param s_x_R1
        !> space discretization operator with one grid point
        !> on the right side (x-direction)
        !
        !> @param s_y_L1
        !> space discretization operator with one grid point
        !> on the left side (y-direction)
        !
        !> @param s_y_R1
        !> space discretization operator with one grid point
        !> on the right side (y-direction)
        !
        !> @param dx
        !> space step along the x-direction
        !
        !> @param dy
        !> space step along the y-direction
        !
        !> @param bc_section
        !> type of edge (NE_edge_type, NW_edge_type, SE_edge_type,
        !> SW_edge_type) and localization of the edge
        !
        !> @param timedev
        !> time derivatives of the grid points
        !--------------------------------------------------------------
        subroutine compute_timedev_anti_corner_with_fluxes(
     $       p_model,
     $       t,nodes,x_map,y_map,
     $       flux_x, flux_y,
     $       timedev,
     $       s_x_L1, s_x_R1,
     $       s_y_L1, s_y_R1,
     $       dx, dy,
     $       bc_section,
     $       interior_nodes,
     $       bf_alignment)

          implicit none
        
          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)                        , intent(in)    :: t
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: nodes
          real(rkind)   , dimension(:)       , intent(in)    :: x_map
          real(rkind)   , dimension(:)       , intent(in)    :: y_map
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_x
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_y
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: timedev
          type(sd_operators_x_oneside_L1)    , intent(in)    :: s_x_L1
          type(sd_operators_x_oneside_R1)    , intent(in)    :: s_x_R1
          type(sd_operators_y_oneside_L1)    , intent(in)    :: s_y_L1
          type(sd_operators_y_oneside_R1)    , intent(in)    :: s_y_R1
          real(rkind)                        , intent(in)    :: dx
          real(rkind)                        , intent(in)    :: dy
          integer       , dimension(4)       , intent(in)    :: bc_section
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment

          integer(ikind) :: i_min,j_min
          integer(ikind) :: i,j
          logical        :: compute_edge

          logical :: compute_point1
          logical :: compute_point2
          logical :: compute_point3
          logical :: compute_point4


          i_min = bc_section(2)
          j_min = bc_section(3)

          
          call determine_edge_points_computed(
     $         bc_section(4),
     $         compute_point1,
     $         compute_point2,
     $         compute_point3,
     $         compute_point4)


          select case(bc_section(1))

            !  ___ ___
            ! |   |BBB|
            ! |___|BBB|
            ! |   |   |  NE_edge
            ! |___|___|
            !------------
            case(NE_edge_type)
               
               compute_edge =
     $              compute_edge_N(j_min,y_map,bc_timedev_choice).and.
     $              compute_edge_E(i_min,x_map,bc_timedev_choice)
               
               if(compute_edge) then

                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |CCC|   |  NE_edge(1,1): like NE_corner(1,1)
                  ! |CCC|___|
                  !------------
                  if(compute_point1) then

                     i=i_min
                     j=j_min
                     
                     timedev(i,j,:) = compute_timedev_corner_local(
     $                    p_model,
     $                    t, x_map, y_map, nodes,
     $                    dx,dy, i,j,
     $                    incoming_right, incoming_right,
     $                    gradient_x_x_oneside_R0, gradient_y_y_oneside_R0)

                  end if


                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |   |CCC|  NE_edge(2,1): like N_edge
                  ! |___|CCC|
                  !------------
                  if(compute_point2) then

                     i=i_min+1
                     j=j_min


                     ! compute fluxes N_edge
                     call compute_flux_x_anti_corner(
     $                    bf_alignment, nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_y_R1,
     $                    p_model,
     $                    i,j,
     $                    flux_x)
                     
                     ! deduce the time derivatives
                     timedev(i,j,:) = compute_timedev_y_edge_local(
     $                    p_model,
     $                    t,x_map,y_map,
     $                    nodes, dx,dy, i,j,
     $                    flux_x,
     $                    incoming_right, gradient_y_y_oneside_R0)
                     
                  end if


                  !  ___ ___
                  ! |CCC|   |
                  ! |CCC|___|
                  ! |   |   |  NE_edge(2,1): like E_edge
                  ! |___|___|
                  !------------
                  if(compute_point3) then

                     i=i_min
                     j=j_min+1
                     
                     ! compute fluxes E_edge
                     call compute_flux_y_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_x_R1,
     $                    p_model,
     $                    i,j,
     $                    flux_y)
                     
                     ! deduce the time derivatives
                     timedev(i,j,:) = compute_timedev_x_edge_local(
     $                    p_model,
     $                    t,x_map,y_map,
     $                    nodes, dx,dy, i,j,
     $                    flux_y,
     $                    incoming_right, gradient_x_x_oneside_R0)

                  end if


                  !  ___ ___
                  ! |   |CCC|
                  ! |___|CCC|
                  ! |   |   |  NE_edge(2,2): like NE_corner(2,2)
                  ! |___|___|
                  !------------
                  if(compute_point4) then

                     i=i_min+1
                     j=j_min+1

                     timedev(i,j,:) = compute_timedev_corner_local(
     $                    p_model,
     $                    t, x_map, y_map, nodes,
     $                    dx,dy, i,j,
     $                    incoming_right, incoming_right,
     $                    gradient_x_x_oneside_R0, gradient_y_y_oneside_R0)

                  end if
                     
               end if


            !  ___ ___
            ! |BBB|   |
            ! |BBB|___|
            ! |   |   |  NW_edge
            ! |___|___|
            !------------
            case(NW_edge_type)

               compute_edge =
     $              compute_edge_N(j_min,y_map,bc_timedev_choice).and.
     $              compute_edge_W(i_min,x_map,bc_timedev_choice)
               
               if(compute_edge) then

                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |CCC|   |  NW_edge(1,1): like N_edge
                  ! |CCC|___|
                  !------------
                  if(compute_point1) then

                     i=i_min
                     j=j_min

                     ! compute the x-fluxes
                     call compute_flux_x_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_y_R1,
     $                    p_model,
     $                    i,j,
     $                    flux_x)

                     ! deduce the time derivatives
                     timedev(i,j,:) = compute_timedev_y_edge_local(
     $                    p_model,
     $                    t,x_map,y_map,
     $                    nodes, dx,dy, i,j,
     $                    flux_x,
     $                    incoming_right, gradient_y_y_oneside_R0)

                  end if


                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |   |CCC|  NW_edge(2,1): like NW_corner(2,1)
                  ! |___|CCC|
                  !------------
                  if(compute_point2) then

                     i=i_min+1
                     j=j_min

                     timedev(i,j,:) = compute_timedev_corner_local(
     $                    p_model,
     $                    t, x_map, y_map, nodes,
     $                    dx,dy, i,j,
     $                    incoming_left, incoming_right,
     $                    gradient_x_x_oneside_L0, gradient_y_y_oneside_R0)

                  end if


                  !  ___ ___
                  ! |CCC|   |
                  ! |CCC|___|
                  ! |   |   |  NW_edge(1,2): like NW_corner(1,2)
                  ! |___|___|
                  !------------
                  if(compute_point3) then

                     i=i_min
                     j=j_min+1

                     timedev(i,j,:) = compute_timedev_corner_local(
     $                    p_model,
     $                    t, x_map, y_map, nodes,
     $                    dx,dy, i,j,
     $                    incoming_left, incoming_right,
     $                    gradient_x_x_oneside_L0, gradient_y_y_oneside_R0)

                  end if


                  !  ___ ___
                  ! |   |CCC|
                  ! |___|CCC|
                  ! |   |   |  NW_edge(2,2): like W_edge
                  ! |___|___|
                  !------------
                  if(compute_point4) then

                     i=i_min+1
                     j=j_min+1


                     ! compute the y-fluxes
                     call compute_flux_y_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_x_L1,
     $                    p_model,
     $                    i,j,
     $                    flux_y)

                     ! deduce the time derivatives
                     timedev(i,j,:) = compute_timedev_x_edge_local(
     $                    p_model,
     $                    t,x_map,y_map,
     $                    nodes, dx,dy, i,j,
     $                    flux_y,
     $                    incoming_left, gradient_x_x_oneside_L0)

                  end if

               end if


            !  ___ ___
            ! |   |   |
            ! |___|___|
            ! |CCC|   |  SW_edge
            ! |CCC|___|
            !------------
            case(SW_edge_type)

               compute_edge =
     $              compute_edge_S(j_min,y_map,bc_timedev_choice).and.
     $              compute_edge_W(i_min,x_map,bc_timedev_choice)
               
               if(compute_edge) then

                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |CCC|   |  SW_edge(1,1): like SW_corner(1,1)
                  ! |CCC|___|
                  !------------
                  if(compute_point1) then

                     i=i_min
                     j=j_min
                     
                     timedev(i,j,:) = compute_timedev_corner_local(
     $                    p_model,
     $                    t, x_map, y_map, nodes,
     $                    dx,dy, i,j,
     $                    incoming_left, incoming_left,
     $                    gradient_x_x_oneside_L0, gradient_y_y_oneside_L0)

                  end if
                     

                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |   |CCC|  SW_edge(2,1): like W_edge
                  ! |___|CCC|
                  !------------
                  if(compute_point2) then

                     i=i_min+1
                     j=j_min
                     
                     ! compute the y-fluxes
                     call compute_flux_y_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_x_L1,
     $                    p_model,
     $                    i,j,
     $                    flux_y)

                     ! deduce the time derivatives
                     timedev(i,j,:) = compute_timedev_x_edge_local(
     $                    p_model,
     $                    t,x_map,y_map,
     $                    nodes, dx,dy, i,j,
     $                    flux_y,
     $                    incoming_left, gradient_x_x_oneside_L0)

                  end if


                  !  ___ ___
                  ! |CCC|   |
                  ! |CCC|___|
                  ! |   |   |  SW_edge(2,1): like S_edge
                  ! |___|___|
                  !------------
                  if(compute_point3) then

                     i=i_min
                     j=j_min+1

                     ! compute the x-fluxes
                     call compute_flux_x_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_y_L1,
     $                    p_model,
     $                    i,j,
     $                    flux_x)
                     
                     ! deduce the time derivatives
                     timedev(i,j,:) = compute_timedev_y_edge_local(
     $                    p_model,
     $                    t,x_map,y_map,
     $                    nodes, dx,dy, i,j,
     $                    flux_x,
     $                    incoming_left, gradient_y_y_oneside_L0)

                  end if


                  !  ___ ___
                  ! |   |CCC|
                  ! |___|CCC|
                  ! |   |   |  SW_edge(2,2): like SW_corner(2,2)
                  ! |___|___|
                  !------------
                  if(compute_point4) then

                     i=i_min+1
                     j=j_min+1
                     
                     timedev(i,j,:) = compute_timedev_corner_local(
     $                    p_model,
     $                    t, x_map, y_map, nodes,
     $                    dx,dy, i,j,
     $                    incoming_left, incoming_left,
     $                    gradient_x_x_oneside_L0, gradient_y_y_oneside_L0)

                  end if

               end if


            !  ___ ___
            ! |   |   |
            ! |___|___|
            ! |   |BBB  SE_edge
            ! |___|BBB|
            !------------
            case(SE_edge_type)

               compute_edge =
     $              compute_edge_S(j_min,y_map,bc_timedev_choice).and.
     $              compute_edge_E(i_min,x_map,bc_timedev_choice)
               
               if(compute_edge) then

                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |CCC|   |  SE_edge(1,1): like E_edge
                  ! |CCC|___|
                  !------------
                  if(compute_point1) then

                     i=i_min
                     j=j_min
                     
                     ! compute the y-fluxes
                     call compute_flux_y_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_x_R1,
     $                    p_model,
     $                    i,j,
     $                    flux_y)

                     ! deduce the time derivatives
                     timedev(i,j,:) = compute_timedev_x_edge_local(
     $                    p_model,
     $                    t,x_map,y_map,
     $                    nodes, dx,dy, i,j,
     $                    flux_y,
     $                    incoming_right, gradient_x_x_oneside_R0)

                  end if
                  

                  !  ___ ___
                  ! |   |   |
                  ! |___|___|
                  ! |   |CCC|  SE_edge(2,1): like SE_corner(2,1)
                  ! |___|CCC|
                  !------------
                  if(compute_point2) then

                     i=i_min+1
                     j=j_min
                     
                     timedev(i,j,:) = compute_timedev_corner_local(
     $                    p_model,
     $                    t, x_map, y_map, nodes,
     $                    dx,dy, i,j,
     $                    incoming_right, incoming_left,
     $                    gradient_x_x_oneside_R0, gradient_y_y_oneside_L0)

                  end if
                     

                  !  ___ ___
                  ! |CCC|   |
                  ! |CCC|___|
                  ! |   |   |  SE_edge(2,1): like SE_corner(1,2)
                  ! |___|___|
                  !------------
                  if(compute_point3) then
                     
                     i=i_min
                     j=j_min+1
                     
                     timedev(i,j,:) = compute_timedev_corner_local(
     $                    p_model,
     $                    t, x_map, y_map, nodes,
     $                    dx,dy, i,j,
     $                    incoming_right, incoming_left,
     $                    gradient_x_x_oneside_R0, gradient_y_y_oneside_L0)

                  end if
                     
                     
                  !  ___ ___
                  ! |   |CCC|
                  ! |___|CCC|
                  ! |   |   |  SE_edge(2,2): like S_edge
                  ! |___|___|
                  !------------
                  if(compute_point4) then

                     i=i_min+1
                     j=j_min+1
                     
                     ! compute the x-fluxes
                     call compute_flux_x_anti_corner(
     $                    bf_alignment,
     $                    nodes,
     $                    interior_nodes,
     $                    dx,dy,
     $                    s_y_L1,
     $                    p_model,
     $                    i,j,
     $                    flux_x)
                     
                     ! deduce the time derivatives
                     timedev(i,j,:) = compute_timedev_y_edge_local(
     $                    p_model,
     $                    t,x_map,y_map,
     $                    nodes, dx,dy, i,j,
     $                    flux_x,
     $                    incoming_left, gradient_y_y_oneside_L0)

                  end if

               end if

          end select

        end subroutine compute_timedev_anti_corner_with_fluxes


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the fluxes in the x-direction needed to compute
        !> the time derivative of an anti-corner
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bf_alignment
        !> relative position of the buffer layer compared to
        !> the interior domain
        !
        !> @param bf_nodes
        !> nodes of the buffer layer
        !
        !> @param interior_nodes
        !> nodes of the interior domain
        !
        !> @param dx
        !> space step along the x-direction
        !
        !> @param dy
        !> space step along the y-direction
        !
        !> @param sd_used
        !> space discretization operator
        !
        !>@param p_model
        !> physical model
        !
        !>@param i
        !> x-index of the anti-corner in the buffer layer
        !
        !>@param j
        !> y-index of the anti-corner in the buffer layer
        !
        !>@param flux_x
        !> fluxes along the x-direction for the buffer layer
        !--------------------------------------------------------------
        subroutine compute_flux_x_anti_corner(
     $     bf_alignment,
     $     bf_nodes,
     $     interior_nodes,
     $     dx,dy,
     $     sd_used,
     $     p_model,
     $     i,j,
     $     flux_x)

          implicit none

          integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          real(rkind)                        , intent(in)    :: dx
          real(rkind)                        , intent(in)    :: dy
          class(sd_operators)                , intent(in)    :: sd_used
          type(pmodel_eq)                    , intent(in)    :: p_model
          integer(ikind)                     , intent(in)    :: i
          integer(ikind)                     , intent(in)    :: j
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_x


          logical                        :: grdpts_needed
          integer(ikind), dimension(2,2) :: border_coords
          integer(ikind), dimension(2)   :: cpt_coords

          real(rkind), dimension(:,:,:), allocatable :: tmp_nodes

          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2,2) :: gen_coords


          ! determine whether there are enough grid points
          grdpts_needed = are_grdpts_needed_for_flux_x(
     $         p_model,
     $         sd_used%get_operator_type(),
     $         i,j,
     $         size(bf_nodes,1),size(bf_nodes,2),
     $         border_coords,
     $         cpt_coords)


          ! if there are not enough grid points they should be
          ! extracted and the fluxes are computed from these
          ! temporary grid points
          if(grdpts_needed) then
             
             ! allocate space for the temporary gridpoints
             ! extracted
             allocate(tmp_nodes(
     $            border_coords(1,2)-border_coords(1,1)+1,
     $            border_coords(2,2)-border_coords(2,1)+1,
     $            ne))

             ! compute the general coordinates identifying the
             ! the borders of the gridpoints extracted
             match_table = get_bf_layer_match_table(
     $            bf_alignment)
             
             gen_coords(1,1) = border_coords(1,1) + match_table(1)
             gen_coords(1,2) = border_coords(1,2) + match_table(1)
             gen_coords(2,1) = border_coords(2,1) + match_table(2)
             gen_coords(2,2) = border_coords(2,2) + match_table(2)
             

             ! extract the grid points from the current nodes of
             ! the buffer layer and the interior domain
             call combine_grdpts_to_compute_fluxes(
     $            bf_alignment, bf_nodes,
     $            interior_nodes,
     $            gen_coords,
     $            tmp_nodes)


             !compute the x-fluxes
             flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $            tmp_nodes,dx,dy,
     $            cpt_coords(1),cpt_coords(2),
     $            sd_used)
             
             flux_x(i+1,j,:) = p_model%compute_flux_x_oneside(
     $            tmp_nodes,dx,dy,
     $            cpt_coords(1)+1,cpt_coords(2),
     $            sd_used)

             deallocate(tmp_nodes)


          ! otherwise the fluxes are directly computed from the
          ! existing nodes
          else

             flux_x(i,j,:) = p_model%compute_flux_x_oneside(
     $            bf_nodes,dx,dy,
     $            i,j,
     $            sd_used)
             
             flux_x(i+1,j,:) = p_model%compute_flux_x_oneside(
     $            bf_nodes,dx,dy,
     $            i+1,j,
     $            sd_used)

          end if

        end subroutine compute_flux_x_anti_corner



        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the fluxes in the x-direction needed to compute
        !> the time derivative of an anti-corner
        !
        !> @date
        !> 28_01_2015 - initial version - J.L. Desmarais
        !
        !> @param bf_alignment
        !> relative position of the buffer layer compared to
        !> the interior domain
        !
        !> @param bf_nodes
        !> nodes of the buffer layer
        !
        !> @param interior_nodes
        !> nodes of the interior domain
        !
        !> @param dx
        !> space step along the x-direction
        !
        !> @param dy
        !> space step along the y-direction
        !
        !> @param sd_used
        !> space discretization operator
        !
        !>@param p_model
        !> physical model
        !
        !>@param i
        !> x-index of the anti-corner in the buffer layer
        !
        !>@param j
        !> y-index of the anti-corner in the buffer layer
        !
        !>@param flux_y
        !> fluxes along the y-direction for the buffer layer
        !--------------------------------------------------------------
        subroutine compute_flux_y_anti_corner(
     $     bf_alignment,
     $     bf_nodes,
     $     interior_nodes,
     $     dx,dy,
     $     sd_used,
     $     p_model,
     $     i,j,
     $     flux_y)

          implicit none

          integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          real(rkind)                        , intent(in)    :: dx
          real(rkind)                        , intent(in)    :: dy
          class(sd_operators)                , intent(in)    :: sd_used
          type(pmodel_eq)                    , intent(in)    :: p_model
          integer(ikind)                     , intent(in)    :: i
          integer(ikind)                     , intent(in)    :: j
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: flux_y


          logical                        :: grdpts_needed
          integer(ikind), dimension(2,2) :: border_coords
          integer(ikind), dimension(2)   :: cpt_coords

          real(rkind), dimension(:,:,:), allocatable :: tmp_nodes

          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2,2) :: gen_coords


          ! determine whether there are enough grid points
          grdpts_needed = are_grdpts_needed_for_flux_y(
     $         p_model,
     $         sd_used%get_operator_type(),
     $         i,j,
     $         size(bf_nodes,1),size(bf_nodes,2),
     $         border_coords,
     $         cpt_coords)


          ! if there are not enough grid points they should be
          ! extracted and the fluxes are computed from these
          ! temporary grid points
          if(grdpts_needed) then
             
             ! allocate space for the temporary gridpoints
             ! extracted
             allocate(tmp_nodes(
     $            border_coords(1,2)-border_coords(1,1)+1,
     $            border_coords(2,2)-border_coords(2,1)+1,
     $            ne))

             ! compute the general coordinates identifying the
             ! the borders of the gridpoints extracted
             match_table = get_bf_layer_match_table(
     $            bf_alignment)
             
             gen_coords(1,1) = border_coords(1,1) + match_table(1)
             gen_coords(1,2) = border_coords(1,2) + match_table(1)
             gen_coords(2,1) = border_coords(2,1) + match_table(2)
             gen_coords(2,2) = border_coords(2,2) + match_table(2)
             

             ! extract the grid points from the current nodes of
             ! the buffer layer and the interior domain
             call combine_grdpts_to_compute_fluxes(
     $            bf_alignment, bf_nodes,
     $            interior_nodes,
     $            gen_coords,
     $            tmp_nodes)


             !compute the y-fluxes
             flux_y(i,j,:) = p_model%compute_flux_y_oneside(
     $            tmp_nodes,dx,dy,
     $            cpt_coords(1),cpt_coords(2),
     $            sd_used)
             
             flux_y(i,j+1,:) = p_model%compute_flux_y_oneside(
     $            tmp_nodes,dx,dy,
     $            cpt_coords(1),cpt_coords(2)+1,
     $            sd_used)

             deallocate(tmp_nodes)


          ! otherwise the fluxes are directly computed from the
          ! existing nodes
          else

             flux_y(i,j,:)   = p_model%compute_flux_y_oneside(
     $            bf_nodes,dx,dy,
     $            i,j,
     $            sd_used)
             
             flux_y(i,j+1,:) = p_model%compute_flux_y_oneside(
     $            bf_nodes,dx,dy,
     $            i,j+1,
     $            sd_used)

          end if

        end subroutine compute_flux_y_anti_corner        

      end module hedstrom_xy_anti_corner_flux_module
