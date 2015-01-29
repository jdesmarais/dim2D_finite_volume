      module hedstrom_xy_diag_fluxes_module

        use bc_operators_nopt_module, only :
     $       compute_edge_N,
     $       compute_edge_S,
     $       compute_edge_E,
     $       compute_edge_W,
     $       are_grdpts_needed_for_flux_x,
     $       are_grdpts_needed_for_flux_y,
     $       combine_grdpts_to_compute_fluxes

        use bf_layer_bc_procedure_module, only :
     $       SE_edge_type,
     $       SW_edge_type,
     $       NE_edge_type,
     $       NW_edge_type

        use bf_layer_bc_sections_class, only :
     $       determine_edge_points_computed

        use bf_layer_sync_module, only :
     $       get_bf_layer_match_table

        use hedstrom_xy_corners_module, only :
     $       compute_n_timedev_with_openbc_local

        use interface_primary, only :
     $       gradient_x_proc,
     $       gradient_y_proc

        use openbc_operators_module, only :
     $       incoming_proc,
     $       incoming_left,
     $       incoming_right

        use parameters_constant, only :
     $       bc_timedev_choice,
     $       left,right,
     $       n1_direction,
     $       n2_direction

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

        use sd_operators_n1_oneside_L0_class, only :
     $       sd_operators_n1_oneside_L0

        use sd_operators_n1_oneside_L1_class, only :
     $       sd_operators_n1_oneside_L1

        use sd_operators_n1_oneside_R1_class, only :
     $       sd_operators_n1_oneside_R1

        use sd_operators_n1_oneside_R0_class, only :
     $       sd_operators_n1_oneside_R0

        use sd_operators_n2_oneside_L0_class, only :
     $       sd_operators_n2_oneside_L0

        use sd_operators_n2_oneside_L1_class, only :
     $       sd_operators_n2_oneside_L1

        use sd_operators_n2_oneside_R1_class, only :
     $       sd_operators_n2_oneside_R1

        use sd_operators_n2_oneside_R0_class, only :
     $       sd_operators_n2_oneside_R0


        implicit none

        private
        public :: compute_timedev_anti_corner_with_diag_fluxes



        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the time derivatives at (i,j) resulting
        !> of the application of the boundary condition on
        !> an xy edge: NE_edge, NW_edge, SE_edge, SW_edge
        !> using diagonal fluxes
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
        subroutine compute_timedev_anti_corner_with_diag_fluxes(
     $     p_model,
     $     interior_nodes,
     $     bf_alignment,
     $     nodes, x_map, y_map,
     $     dx, dy,
     $     bc_section,
     $     timedev)
        
          implicit none
        
          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)   , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2)     , intent(in)    :: bf_alignment
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: nodes
          real(rkind)   , dimension(:)       , intent(in)    :: x_map
          real(rkind)   , dimension(:)       , intent(in)    :: y_map
          real(rkind)                        , intent(in)    :: dx
          real(rkind)                        , intent(in)    :: dy
          integer       , dimension(4)       , intent(in)    :: bc_section
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: timedev

          
          integer(ikind) :: i_min,j_min
          integer(ikind) :: i,j
          logical        :: compute_edge

          logical :: compute_point1
          logical :: compute_point2
          logical :: compute_point3
          logical :: compute_point4

          real(rkind) :: dn

          type(sd_operators_n1_oneside_L0) :: sd_n1_L0
          type(sd_operators_n1_oneside_L1) :: sd_n1_L1
          type(sd_operators_n1_oneside_R1) :: sd_n1_R1
          type(sd_operators_n1_oneside_R0) :: sd_n1_R0

          type(sd_operators_n2_oneside_L0) :: sd_n2_L0
          type(sd_operators_n2_oneside_L1) :: sd_n2_L1
          type(sd_operators_n2_oneside_R1) :: sd_n2_R1
          type(sd_operators_n2_oneside_R0) :: sd_n2_R0


          i_min = bc_section(2)
          j_min = bc_section(3)

          
          call determine_edge_points_computed(
     $         bc_section(4),
     $         compute_point1,
     $         compute_point2,
     $         compute_point3,
     $         compute_point4)

          dn = Sqrt(dx**2+dy**2)

          select case(bc_section(1))

            !  ___ ___
            ! |   |CCC|
            ! |___|CCC|
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
                     
                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n2_R1,
     $                    p_model,
     $                    n2_direction,
     $                    incoming_right,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_R0)

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


                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n2_R1,
     $                    p_model,
     $                    n2_direction,
     $                    incoming_right,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_R0)
                     
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
                     
                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n2_R1,
     $                    p_model,
     $                    n2_direction,
     $                    incoming_right,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_R0)

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

                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n2_R0,
     $                    p_model,
     $                    n2_direction,
     $                    incoming_right,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_R0)

                  end if
                     
               end if


            !  ___ ___
            ! |CCC|   |
            ! |CCC|___|
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

                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n1_L1,
     $                    p_model,
     $                    n1_direction,
     $                    incoming_left,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_R0)

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

                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n1_L1,
     $                    p_model,
     $                    n1_direction,
     $                    incoming_left,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_R0)

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

                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n1_L0,
     $                    p_model,
     $                    n1_direction,
     $                    incoming_left,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_R0)

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

                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment, 
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n1_L1,
     $                    p_model,
     $                    n1_direction,
     $                    incoming_left,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_R0)

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
                     
                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n2_L0,
     $                    p_model,
     $                    n2_direction,
     $                    incoming_left,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_L0)

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
                     
                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n2_L1,
     $                    p_model,
     $                    n2_direction,
     $                    incoming_left,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_L0)

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

                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n2_L1,
     $                    p_model,
     $                    n2_direction,
     $                    incoming_left,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_L0)

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
                     
                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n2_L1,
     $                    p_model,
     $                    n2_direction,
     $                    incoming_left,
     $                    gradient_x_x_oneside_L0,
     $                    gradient_y_y_oneside_L0)

                  end if

               end if


            !  ___ ___
            ! |   |   |
            ! |___|___|
            ! |   |CCC|  SE_edge
            ! |___|CCC|
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
                     
                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n1_R1,
     $                    p_model,
     $                    n1_direction,
     $                    incoming_right,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_L0)

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
                     
                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n1_R0,
     $                    p_model,
     $                    n1_direction,
     $                    incoming_right,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_L0)

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
                     
                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n1_R1,
     $                    p_model,
     $                    n1_direction,
     $                    incoming_right,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_L0)

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
                     
                     timedev(i,j,:) = compute_time_dev_openbc(
     $                    interior_nodes,
     $                    bf_alignment,
     $                    nodes,
     $                    i,j,
     $                    dx,dy,dn,
     $                    sd_n1_R1,
     $                    p_model,
     $                    n1_direction,
     $                    incoming_right,
     $                    gradient_x_x_oneside_R0,
     $                    gradient_y_y_oneside_L0)

                  end if

               end if

          end select

        end subroutine compute_timedev_anti_corner_with_diag_fluxes



        function compute_time_dev_openbc(
     $     interior_nodes,
     $     bf_alignment,
     $     bf_nodes,
     $     i,j,
     $     dx,dy,dn,
     $     sd_used,
     $     p_model,
     $     outward_dir,
     $     incoming_wave,
     $     gradient_x,
     $     gradient_y)
     $     result(timedev)

          implicit none

          real(rkind)   , dimension(nx,ny,ne), intent(in) :: interior_nodes
          integer(ikind), dimension(2,2)     , intent(in) :: bf_alignment
          real(rkind)   , dimension(:,:,:)   , intent(in) :: bf_nodes
          integer(ikind)                     , intent(in) :: i
          integer(ikind)                     , intent(in) :: j
          real(rkind)                        , intent(in) :: dx
          real(rkind)                        , intent(in) :: dy
          real(rkind)                        , intent(in) :: dn
          class(sd_operators)                , intent(in) :: sd_used
          type(pmodel_eq)                    , intent(in) :: p_model
          integer                            , intent(in) :: outward_dir
          procedure(incoming_proc)                        :: incoming_wave
          procedure(gradient_x_proc)                      :: gradient_x
          procedure(gradient_y_proc)                      :: gradient_y
          real(rkind)   , dimension(ne)                   :: timedev


          logical                        :: grdpts_needed
          integer(ikind), dimension(2,2) :: border_coords
          integer(ikind), dimension(2)   :: cpt_coords

          real(rkind), dimension(:,:,:), allocatable :: tmp_nodes

          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2,2) :: gen_coords

          
          real(rkind), dimension(ne) :: flux_diag1
          real(rkind), dimension(ne) :: flux_diag2

          integer :: k

          integer(ikind) :: i_c
          integer(ikind) :: j_c          


          !1) compute the contribution of the outward direction
          !   to the time derivatives
          timedev = compute_n_timedev_with_openbc_local(
     $         bf_nodes, i,j,
     $         p_model, dx,dy,
     $         gradient_x,
     $         gradient_y,
     $         incoming_wave,
     $         outward_dir)

          
          !2) compute the contribution of the direction normal
          !   to the outward direction to the time derivatives

          !2.1) determine whether there are enough grid points
          !     to compute the fluxes
          !---------------------------------------------------
          select case(outward_dir)

            case(n1_direction)

               ! check if there are enough grid points to compute
               ! the fluxes
               !-------------------------------------------------
               grdpts_needed = are_grdpts_needed_for_flux_y(
     $              p_model,
     $              sd_used%get_operator_type(),
     $              i,j,
     $              size(bf_nodes,1),size(bf_nodes,2),
     $              border_coords,
     $              cpt_coords)


            case(n2_direction)

               ! check if there are enough grid points to compute
               ! the fluxes
               !-------------------------------------------------
               grdpts_needed = are_grdpts_needed_for_flux_x(
     $              p_model,
     $              sd_used%get_operator_type(),
     $              i,j,
     $              size(bf_nodes,1),size(bf_nodes,2),
     $              border_coords,
     $              cpt_coords)

          end select


          !2.2) if grid points are needed, combine nodes from the
          !     interior and the buffer layer and compute the
          !     fluxes using the tmp_nodes
          !------------------------------------------------------
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

             ! compute the fluxes using the tmp_nodes
             i_c = cpt_coords(1)
             j_c = cpt_coords(2)

             select case(outward_dir)

               case(n1_direction)

                  flux_diag1 = p_model%compute_flux_y_oneside(
     $                 tmp_nodes,dn,dn,i_c  ,j_c  ,sd_used)
             
                  flux_diag2 = p_model%compute_flux_y_oneside(
     $                 tmp_nodes,dn,dn,i_c+1,j_c+1,sd_used)

               case(n2_direction)

                  flux_diag1 = p_model%compute_flux_x_oneside(
     $                 tmp_nodes,dn,dn,i_c  ,j_c  ,sd_used)
             
                  flux_diag2 = p_model%compute_flux_x_oneside(
     $                 tmp_nodes,dn,dn,i_c+1,j_c-1,sd_used)

             end select

          !2.3) if grid points are not needed, compute directly
          !     the fluxes using the bf_nodes
          !------------------------------------------------------
          else
             
             ! compute the fluxes using the bf_nodes
             select case(outward_dir)

               case(n1_direction)

                  flux_diag1 = p_model%compute_flux_y_oneside(
     $                 tmp_nodes,dn,dn,i  ,j  ,sd_used)
             
                  flux_diag2 = p_model%compute_flux_y_oneside(
     $                 tmp_nodes,dn,dn,i+1,j+1,sd_used)

               case(n2_direction)

                  flux_diag1 = p_model%compute_flux_x_oneside(
     $                 tmp_nodes,dn,dn,i  ,j  ,sd_used)
             
                  flux_diag2 = p_model%compute_flux_x_oneside(
     $                 tmp_nodes,dn,dn,i+1,j-1,sd_used)

             end select
             

          end if

          !2.4) add the contribution of the diagonal transverse
          !     fluxes to the time derivatives
          !------------------------------------------------------
          do k=1,ne
             timedev(k) = timedev(k) + 1.0d0/dn*(flux_diag1(k)-flux_diag2(k))
          end do

        end function compute_time_dev_openbc

      end module hedstrom_xy_diag_fluxes_module
