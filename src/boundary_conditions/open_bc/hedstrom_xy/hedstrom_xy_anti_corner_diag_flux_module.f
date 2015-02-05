      module hedstrom_xy_anti_corner_diag_flux_module

        use bc_operators_nopt_module, only :
     $       compute_edge_N,
     $       compute_edge_S,
     $       compute_edge_E,
     $       compute_edge_W,
     $       are_grdpts_needed_for_flux_x,
     $       are_grdpts_needed_for_flux_y,
     $       combine_grdpts_to_compute_fluxes

        use bf_layer_bc_sections_class, only :
     $       determine_edge_points_computed

        use bf_layer_sync_module, only :
     $       get_bf_layer_match_table

        use hedstrom_xy_corners_module, only :
     $       compute_timedev_with_openbc

        use interface_primary, only :
     $       gradient_proc

        use n_coords_module, only :
     $       get_dn

        use openbc_operators_module, only :
     $       incoming_proc,
     $       incoming_left,
     $       incoming_right

        use parameters_bf_layer, only :
     $       SE_edge_type,
     $       SW_edge_type,
     $       NE_edge_type,
     $       NW_edge_type

        use parameters_constant, only :
     $       bc_timedev_choice,
     $       left,right,
     $       n1_direction,
     $       n2_direction

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind        

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use sd_operators_fd_ncoords_module, only :
     $       gradient_n1_oneside_L0,
     $       gradient_n1_oneside_R0,
     $       gradient_n2_oneside_L0,
     $       gradient_n2_oneside_R0

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
     $       p_model,
     $       t,nodes,x_map,y_map,
     $       timedev,
     $       dx,dy,
     $       bc_section,
     $       interior_nodes,
     $       bf_alignment)
        
          implicit none
        
          type(pmodel_eq)                    , intent(in)    :: p_model
          real(rkind)                        , intent(in)    :: t
          real(rkind)   , dimension(:,:,:)   , intent(in)    :: nodes
          real(rkind)   , dimension(:)       , intent(in)    :: x_map
          real(rkind)   , dimension(:)       , intent(in)    :: y_map
          real(rkind)   , dimension(:,:,:)   , intent(inout) :: timedev
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

          real(rkind) :: dn
          
          real(rkind) :: x,y

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

          dn =  get_dn(dx,dy)

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

                     x = x_map(i)
                     y = y_map(j)
                     
                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n2_direction,
     $                    gradient_n2_oneside_R0,dn,
     $                    incoming_right,
     $                    sd_n2_R1,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)

                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n2_direction,
     $                    gradient_n2_oneside_R0,dn,
     $                    incoming_right,
     $                    sd_n2_R1,
     $                    interior_nodes,
     $                    bf_alignment)
                     
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
                     
                     x = x_map(i)
                     y = y_map(j)

                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n2_direction,
     $                    gradient_n2_oneside_R0,dn,
     $                    incoming_right,
     $                    sd_n2_R1,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)

                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n2_direction,
     $                    gradient_n2_oneside_R0,dn,
     $                    incoming_right,
     $                    sd_n2_R0,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)

                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n1_direction,
     $                    gradient_n1_oneside_L0,dn,
     $                    incoming_left,
     $                    sd_n1_L1,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)

                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n1_direction,
     $                    gradient_n1_oneside_L0,dn,
     $                    incoming_left,
     $                    sd_n1_L1,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)

                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n1_direction,
     $                    gradient_n1_oneside_L0,dn,
     $                    incoming_left,
     $                    sd_n1_L0,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)

                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n1_direction,
     $                    gradient_n1_oneside_L0,dn,
     $                    incoming_left,
     $                    sd_n1_L1,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)
                     
                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n1_direction,
     $                    gradient_n2_oneside_L0,dn,
     $                    incoming_left,
     $                    sd_n2_L0,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)
                     
                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n1_direction,
     $                    gradient_n2_oneside_L0,dn,
     $                    incoming_left,
     $                    sd_n2_L1,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)

                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n1_direction,
     $                    gradient_n2_oneside_L0,dn,
     $                    incoming_left,
     $                    sd_n2_L1,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)
                     
                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n1_direction,
     $                    gradient_n2_oneside_L0,dn,
     $                    incoming_left,
     $                    sd_n2_L1,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)
                     
                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n1_direction,
     $                    gradient_n1_oneside_R0,dn,
     $                    incoming_right,
     $                    sd_n1_R1,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)
                     
                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n1_direction,
     $                    gradient_n1_oneside_R0,dn,
     $                    incoming_right,
     $                    sd_n1_R0,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)
                     
                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n1_direction,
     $                    gradient_n1_oneside_R0,dn,
     $                    incoming_right,
     $                    sd_n1_R1,
     $                    interior_nodes,
     $                    bf_alignment)

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

                     x = x_map(i)
                     y = y_map(j)

                     timedev(i,j,:) = compute_timedev_local(
     $                    t,x,y,
     $                    nodes,i,j,
     $                    p_model,
     $                    n1_direction,
     $                    gradient_n1_oneside_R0,dn,
     $                    incoming_right,
     $                    sd_n1_R1,
     $                    interior_nodes,
     $                    bf_alignment)

                  end if

               end if

          end select

        end subroutine compute_timedev_anti_corner_with_diag_fluxes


        function compute_timedev_local(
     $     t,x,y,
     $     bf_nodes,i,j,
     $     p_model,
     $     outward_dir,
     $     gradient_n,dn,
     $     incoming_wave,
     $     sd_used,
     $     interior_nodes,
     $     bf_alignment)
     $     result(timedev)

          implicit none

          real(rkind)                        , intent(in) :: t
          real(rkind)                        , intent(in) :: x
          real(rkind)                        , intent(in) :: y
          real(rkind)   , dimension(:,:,:)   , intent(in) :: bf_nodes
          integer(ikind)                     , intent(in) :: i
          integer(ikind)                     , intent(in) :: j
          type(pmodel_eq)                    , intent(in) :: p_model
          integer                            , intent(in) :: outward_dir
          procedure(gradient_proc)                        :: gradient_n
          real(rkind)                        , intent(in) :: dn
          procedure(incoming_proc)                        :: incoming_wave
          class(sd_operators)                , intent(in) :: sd_used
          real(rkind)   , dimension(nx,ny,ne), intent(in) :: interior_nodes
          integer(ikind), dimension(2,2)     , intent(in) :: bf_alignment
          real(rkind)   , dimension(ne)                   :: timedev


          logical                        :: grdpts_needed
          integer(ikind), dimension(2,2) :: border_coords
          integer(ikind), dimension(2)   :: cpt_coords

          real(rkind), dimension(:,:,:), allocatable :: tmp_nodes
          real(rkind), dimension(:,:,:), allocatable :: tmp_nodes_n

          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2,2) :: gen_coords

          
          real(rkind), dimension(ne) :: flux_diag1
          real(rkind), dimension(ne) :: flux_diag2
          real(rkind), dimension(ne) :: timedev_n

          integer :: k

          integer(ikind) :: i_c
          integer(ikind) :: j_c
          integer(ikind) :: i_l
          integer(ikind) :: j_l          

          
          !1) determine whether there are enough grid points
          !   to compute the time derivatives
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


          ! allocate space for the temporary gridpoints
          ! extracted
          allocate(tmp_nodes(
     $         border_coords(1,2)-border_coords(1,1)+1,
     $         border_coords(2,2)-border_coords(2,1)+1,
     $         ne))


          !2) if grid points are needed, combine nodes from the
          !     interior and the buffer layer
          !------------------------------------------------------
          if(grdpts_needed) then             

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

          !3) otherwise simply use the nodes from bf_nodes
          !------------------------------------------------------
          else

             tmp_nodes = bf_nodes(
     $            border_coords(1,1):border_coords(1,2),
     $            border_coords(2,1):border_coords(2,2),
     $            :)

          end if

          
          !4) compute the contribution of the outgoing waves
          !   to the time derivatives
          !------------------------------------------------------
          ! determine the central pt computed
          i_c = cpt_coords(1)
          j_c = cpt_coords(2)

          ! compute the contribution of the outward direction
          ! to the time derivatives with the tmp_nodes
          timedev_n = compute_timedev_with_openbc(
     $         t,x,y,
     $         tmp_nodes, i_c,j_c,
     $         p_model,
     $         outward_dir,
     $         gradient_n, dn,
     $         incoming_wave)


          !5) compute the contribution of the transverse diagonal
          !   fluxes to the time derivatives
          !------------------------------------------------------
          ! convert the nodes into nodes_n
          ! (momentum_x,momentum_y) -> (momentum_n1,momentum_n2)
          allocate(tmp_nodes_n(size(tmp_nodes,1),size(tmp_nodes,2),ne))

          do j_l=1,size(tmp_nodes_n,2)
             do i_l=1, size(tmp_nodes_n,1)
                tmp_nodes_n(i_l,j_l,:) = p_model%compute_xy_to_n_var(tmp_nodes(i_l,j_l,:))
             end do
          end do


          ! compute the fluxes using the tmp_nodes
          select case(outward_dir)

            case(n1_direction)

               flux_diag1 = p_model%compute_flux_y_oneside(
     $              tmp_nodes_n,dn,dn, i_c  , j_c  , sd_used)
          
               flux_diag2 = p_model%compute_flux_y_oneside(
     $              tmp_nodes_n,dn,dn, i_c+1, j_c+1, sd_used)

            case(n2_direction)

               flux_diag1 = p_model%compute_flux_x_oneside(
     $              tmp_nodes_n,dn,dn, i_c  , j_c  , sd_used)
          
               flux_diag2 = p_model%compute_flux_x_oneside(
     $              tmp_nodes_n,dn,dn, i_c+1, j_c-1, sd_used)

          end select

          deallocate(tmp_nodes)
          deallocate(tmp_nodes_n)
          

          !6) add the contribution of the diagonal
          !   transverse fluxes to the time
          !   derivatives
          !------------------------------------------------------
          do k=1,ne

             timedev_n(k) = timedev_n(k) +
     $                      1.0d0/dn*(flux_diag1(k)-flux_diag2(k))

          end do


          !7) get the time derivatives in the (x,y)
          !   reference frame
          !------------------------------------------------------
          timedev = p_model%compute_n_to_xy_var(timedev_n)

        end function compute_timedev_local

      end module hedstrom_xy_anti_corner_diag_flux_module
