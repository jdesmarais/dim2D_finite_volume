      !> @file
      !> analyse of the bc_sections of the interior domain
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> analyse of the bc_sections of the interior domain
      !
      !> @date
      !> 23_03_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module icr_interior_bc_section_module

        use bf_interface_coords_class, only :
     $       bf_interface_coords

        use bf_layer_bc_sections_overlap_module, only :
     $       determine_edge_grdpts_computed,
     $       determine_corner_or_anti_corner_grdpts_computed

        use bf_layer_errors_module, only :
     $       error_bc_section_type

        use icr_interface_class, only :
     $       icr_interface

        use parameters_bf_layer, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $       
     $       dct_icr_distance,
     $       N_overlap,
     $       S_overlap,
     $       E_overlap,
     $       W_overlap,
     $       NS_overlap,
     $       EW_overlap,
     $       
     $       cptnot_type

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        private
        public ::
     $       analyze_interior_square_xy,
     $       analyze_interior_square_bounds,
     $       analyze_interior_xy,
     $       analyze_interior_edge_x,
     $       analyze_interior_edge_y


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of a corner-like or 
        !> anti-corner-like bc_section and stage the
        !> activated bc_interior_pt for update
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param bc_section
        !> boundary section
        !--------------------------------------------------------------
        subroutine analyze_interior_square_xy(
     $       bf_interface_used,
     $       icr_interface_used,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes,
     $       p_model,
     $       bc_section)

          implicit none

          class(bf_interface_coords)      , intent(in)    :: bf_interface_used
          type(icr_interface)             , intent(inout) :: icr_interface_used
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model
          integer(ikind), dimension(5)    , intent(in)    :: bc_section


          integer(ikind), dimension(2,2,2) :: analyzed_grdpts_bounds
          integer(ikind), dimension(2,3)   :: activated_grdpts
          integer                          :: nb_activated_grdpts
          logical                          :: no_activation


          !> determine the edges for the analyze of the grid-points
          !> and the bc_interior_pt activated
          call analyze_interior_square_bounds(
     $         bc_section,
     $         analyzed_grdpts_bounds,
     $         activated_grdpts,
     $         nb_activated_grdpts,
     $         no_activation)
          

          !> check whether the bc_interior_pt are indeed activated
          call analyze_interior_xy(
     $         bf_interface_used,
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         analyzed_grdpts_bounds,
     $         activated_grdpts,
     $         nb_activated_grdpts,
     $         no_activation)

        end subroutine analyze_interior_square_xy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of a corner-like or 
        !> anti-corner-like bc_section and stage the
        !> activated bc_interior_pt for update
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param analyzed_grdpts_bounds
        !> bounds for the grid-points to be analyzed
        !
        !>@param activated_grdpts
        !> general coordinates of the activated bc_interior_pt if
        !> the grid-points analyzed are activated
        !
        !>@param nb_activated_grdpts
        !> number of bc_interior_pt activated
        !--------------------------------------------------------------
        subroutine analyze_interior_square_bounds(
     $     bc_section,
     $     analyzed_grdpts_bounds,
     $     activated_grdpts,
     $     nb_activated_grdpts,
     $     no_activation)

          implicit none

          integer(ikind), dimension(5)    , intent(in)  :: bc_section
          integer(ikind), dimension(2,2,2), intent(out) :: analyzed_grdpts_bounds
          integer(ikind), dimension(2,3)  , intent(out) :: activated_grdpts
          integer                         , intent(out) :: nb_activated_grdpts
          logical                         , intent(out) :: no_activation


          integer(ikind)               :: i_min
          integer(ikind)               :: j_min
          integer       , dimension(4) :: compute_point
          integer(ikind), dimension(2) :: pt


          no_activation = .false.


          !> determine which grid-points are computed in the 
          !> square (corner or anti-corner)
          call determine_corner_or_anti_corner_grdpts_computed(
     $         bc_section(4),
     $         bc_section(5),
     $         compute_point)

          
          i_min = bc_section(2)
          j_min = bc_section(3)


          !> determine the bounds for the activation and the
          !> bc_interior_pt activated
          select case(bc_section(1))
            case(NE_corner_type)

               if(compute_point(1).eq.cptnot_type) then
                  no_activation = .true.

               else
                  nb_activated_grdpts = 1

                  pt = [i_min,j_min]
                  activated_grdpts(:,1) = pt

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)-dct_icr_distance-1,
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)-dct_icr_distance-1,
     $                 
     $                 pt(1)-dct_icr_distance-1,
     $                 pt(2)-dct_icr_distance,
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)-dct_icr_distance/),
     $                 (/2,2,2/))

               end if


            case(NW_corner_type)

               if(compute_point(2).eq.cptnot_type) then
                  no_activation = .true.

               else
                  nb_activated_grdpts = 1

                  pt = [i_min+1,j_min]
                  activated_grdpts(:,1) = pt

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)-dct_icr_distance-1,
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)-dct_icr_distance-1,
     $                 
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)-dct_icr_distance,
     $                 pt(1)+dct_icr_distance+1,
     $                 pt(2)-dct_icr_distance/),
     $                 (/2,2,2/))

               end if


            case(SE_corner_type)

               if(compute_point(3).eq.cptnot_type) then
                  no_activation = .true.

               else
                  nb_activated_grdpts = 1

                  pt = [i_min,j_min+1]
                  activated_grdpts(:,1) = pt

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)-dct_icr_distance-1,
     $                 pt(2)+dct_icr_distance,
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)+dct_icr_distance,
     $                 
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)+dct_icr_distance+1,
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)+dct_icr_distance+1/),
     $                 (/2,2,2/))

               end if

            case(SW_corner_type)

               if(compute_point(4).eq.cptnot_type) then
                  no_activation = .true.

               else
                  nb_activated_grdpts = 1

                  pt = [i_min+1,j_min+1]
                  activated_grdpts(:,1) = pt

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)+dct_icr_distance,
     $                 pt(1)+dct_icr_distance+1,
     $                 pt(2)+dct_icr_distance,
     $                 
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)+dct_icr_distance+1,
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)+dct_icr_distance+1/),
     $                 (/2,2,2/))

               end if

            case(NE_edge_type)
               
               !determine the grid-points potentially activated
               nb_activated_grdpts = 0

               if(compute_point(1).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min,j_min]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(2).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min+1,j_min]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(3).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min,j_min+1]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               !determine the bounds for the analyze of the activation
               if(nb_activated_grdpts.eq.0) then

                  no_activation = .true.

               else

                  pt = [i_min,j_min]

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)-dct_icr_distance,
     $                 pt(1)+2,
     $                 pt(2)-dct_icr_distance,
     $                 
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)-dct_icr_distance+1,
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)+2/),
     $                 (/2,2,2/))

               end if

            case(NW_edge_type)

               !determine the grid-points potentially activated
               nb_activated_grdpts = 0

               if(compute_point(1).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min,j_min]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(2).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min+1,j_min]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(4).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min+1,j_min+1]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               !determine the bounds for the analyze of the activation
               if(nb_activated_grdpts.eq.0) then

                  no_activation = .true.

               else

                  pt =[i_min+1,j_min]

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)-2,
     $                 pt(2)-dct_icr_distance,
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)-dct_icr_distance,
     $                 
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)-dct_icr_distance+1,
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)+2/),
     $                 (/2,2,2/))

               end if

            case(SE_edge_type)

               !determine the grid-points potentially activated
               nb_activated_grdpts = 0

               if(compute_point(1).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min,j_min]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(3).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min,j_min+1]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(4).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min+1,j_min+1]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               !determine the bounds for the analyze of the activation
               if(nb_activated_grdpts.eq.0) then

                  no_activation = .true.

               else

                  pt =[i_min,j_min+1]

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)-2,
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)+dct_icr_distance-1,
     $                 
     $                 pt(1)-dct_icr_distance,
     $                 pt(2)+dct_icr_distance,
     $                 pt(1)+2,
     $                 pt(2)+dct_icr_distance/),
     $                 (/2,2,2/))

               end if

            case(SW_edge_type)

               !determine the grid-points potentially activated
               nb_activated_grdpts = 0

               if(compute_point(2).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min+1,j_min]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(3).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min,j_min+1]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               if(compute_point(4).ne.cptnot_type) then
                  nb_activated_grdpts=nb_activated_grdpts+1
                  pt = [i_min+1,j_min+1]
                  activated_grdpts(:,nb_activated_grdpts) = pt
               end if

               !determine the bounds for the analyze of the activation
               if(nb_activated_grdpts.eq.0) then

                  no_activation = .true.

               else

                  pt =[i_min+1,j_min+1]

                  analyzed_grdpts_bounds = reshape((/
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)-2,
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)+dct_icr_distance-1,
     $                 
     $                 pt(1)-2,
     $                 pt(2)+dct_icr_distance,
     $                 pt(1)+dct_icr_distance,
     $                 pt(2)+dct_icr_distance/),
     $                 (/2,2,2/))

               end if
               
          end select

        end subroutine analyze_interior_square_bounds


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of a corner-like or 
        !> anti-corner-like bc_section and stage the
        !> activated bc_interior_pt for update
        !
        !> @date
        !> 24_03_2015 - initial version - J.L. Desmarais
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param analyzed_grdpts_bounds
        !> bounds for the grid-points to be analyzed
        !
        !>@param activated_grdpts
        !> general coordinates of the activated bc_interior_pt if
        !> the grid-points analyzed are activated
        !
        !>@param nb_activated_grdpts
        !> number of bc_interior_pt activated
        !--------------------------------------------------------------
        subroutine analyze_interior_xy(
     $     bf_interface_used,
     $     icr_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     p_model,
     $     analyzed_grdpts_bounds,
     $     activated_grdpts,
     $     nb_activated_grdpts,
     $     no_activation)

          implicit none

          class(bf_interface_coords)      , intent(in)    :: bf_interface_used
          type(icr_interface)             , intent(inout) :: icr_interface_used
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model
          integer(ikind), dimension(2,2,2), intent(in)    :: analyzed_grdpts_bounds
          integer(ikind), dimension(2,3)  , intent(in)    :: activated_grdpts
          integer                         , intent(in)    :: nb_activated_grdpts
          logical                         , intent(in)    :: no_activation

          logical        :: node_activated
          integer(ikind) :: i,j
          integer        :: k


          if(.not.no_activation) then

             node_activated = .false.
   
             !there are two loops of grid-points that should be analyzed
             !and trigger the activation of the bc_interior_pt
             do k=1,2
                do j=analyzed_grdpts_bounds(2,1,k),analyzed_grdpts_bounds(2,2,k)
                   do i=analyzed_grdpts_bounds(1,1,k),analyzed_grdpts_bounds(1,2,k)
   
                      node_activated = is_interior_node_activated(
     $                     [i,j],
     $                     interior_x_map,
     $                     interior_y_map,
     $                     interior_nodes,
     $                     p_model)
   
                      if(node_activated) then
                         exit
                      end if
   
                   end do
   
                   if(node_activated) then
                      exit
                   end if
   
                end do
   
                if(node_activated) then
                   exit
                end if
   
             end do
   
             !if only one node is activated in the grid-points that are
             !analyzed, the activated bc_interior_pt are staged
             if(node_activated) then

                do k=1,nb_activated_grdpts
   
                   call icr_interface_used%stage(
     $                  activated_grdpts(:,k),
     $                  bf_interface_used)
                   
                end do
                
             end if

          end if

        end subroutine analyze_interior_xy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of an E_edge or W_edge bc_section
        !> and stage the activated bc_interior_pt for update
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param bc_section
        !> boundary section
        !--------------------------------------------------------------
        subroutine analyze_interior_edge_x(
     $       bf_interface_used,
     $       icr_interface_used,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes,
     $       p_model,
     $       bc_section)

          implicit none

          class(bf_interface_coords)      , intent(in)    :: bf_interface_used
          type(icr_interface)             , intent(inout) :: icr_interface_used
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model
          integer(ikind), dimension(5)    , intent(in)    :: bc_section

          logical       , dimension(2) :: compute_edge
          integer                      :: i_analyzed
          integer                      :: dir_analyzed
          logical                      :: analyze
          integer(ikind)               :: loc_i_analyzed
          integer(ikind), dimension(2) :: loc_central_coords
          logical                      :: node_activated
          integer(ikind)               :: last_j_added
          integer(ikind)               :: j
          

          call determine_edge_grdpts_computed(
     $         bc_section(5),
     $         compute_edge)
          

          !determine :
          ! - i_analyzed: the position of the bc_interior_pt on the edge
          ! - dir_analyzed: left or right direction
          ! - analyze : whether the integration overlap the bc_interior_pt
          select case(bc_section(1))

            case(E_edge_type)

               i_analyzed   =  bc_section(2)
               dir_analyzed = -1

               analyze = compute_edge(1)

            case(W_edge_type)

               i_analyzed   =  bc_section(2)+1
               dir_analyzed = +1

               analyze = compute_edge(2)

            case default
               call error_bc_section_type(
     $              'icr_interface_bc_section_module',
     $              'analyze_interior_edge_x',
     $              bc_section(1))

          end select


          ! if the time integration borders do not overlap the
          ! bc_interior_pt in the bc_section, the grid-points
          ! are analyzed
          if(analyze) then

             loc_i_analyzed =
     $            i_analyzed +
     $            dir_analyzed*dct_icr_distance

             ! loop over the bc_interior_pt in the x-direction
             j=bc_section(3)
             last_j_added=j-3
             do j=bc_section(3), bc_section(4)

                loc_central_coords = [loc_i_analyzed,j]

                !check whether the node from which the grid-point
                !depends is activated
                node_activated = is_interior_node_activated(
     $               loc_central_coords,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               p_model)

                ! if the node is activated, the current grid-point
                ! and its nearest neighbors are staged
                if(node_activated) then
                   
                   if(j.ge.(last_j_added+3)) then
                   
                     ! grid-point [i,j-1] is staged
                      if(j.gt.bc_section(3)) then

                         call icr_interface_used%stage(
     $                        [i_analyzed,j-1],
     $                        bf_interface_used)
                         
                      end if
                   end if

                   if(j.ge.(last_j_added+2)) then

                     ! grid-point [i,j] is staged
                      call icr_interface_used%stage(
     $                     [i_analyzed,j],
     $                     bf_interface_used)

                   end if

                   last_j_added = j

                   ! grid-point [i,j+1] is staged
                   if(j.lt.bc_section(4)) then

                      call icr_interface_used%stage(
     $                     [i_analyzed,j+1],
     $                     bf_interface_used)

                   end if
                      
                end if

             end do

          end if

        end subroutine analyze_interior_edge_x


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> analyze the grid-points of an N_edge or S_edge bc_section
        !> and stage the activated bc_interior_pt for update
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param icr_interface_used
        !> object collecting the activated bc_interior_pt
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !>@param bc_section
        !> boundary section
        !--------------------------------------------------------------
        subroutine analyze_interior_edge_y(
     $       bf_interface_used,
     $       icr_interface_used,
     $       interior_x_map,
     $       interior_y_map,
     $       interior_nodes,
     $       p_model,
     $       bc_section)

          implicit none

          class(bf_interface_coords)      , intent(in)    :: bf_interface_used
          type(icr_interface)             , intent(inout) :: icr_interface_used
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes
          type(pmodel_eq)                 , intent(in)    :: p_model
          integer(ikind), dimension(5)    , intent(in)    :: bc_section

          logical       , dimension(2) :: compute_edge
          integer                      :: j_analyzed
          integer                      :: dir_analyzed
          logical                      :: analyze
          integer(ikind)               :: loc_j_analyzed
          integer(ikind), dimension(2) :: loc_central_coords
          logical                      :: node_activated
          integer(ikind)               :: last_i_added
          integer(ikind)               :: i
          
          
          call determine_edge_grdpts_computed(
     $         bc_section(5),
     $         compute_edge)


          !determine :
          ! - j_analyzed: the position of the bc_interior_pt on the edge
          ! - dir_analyzed: top or bottom direction
          ! - analyze : whether the integration overlap the bc_interior_pt
          select case(bc_section(1))

            case(N_edge_type)

               j_analyzed   =  bc_section(3)
               dir_analyzed = -1
               analyze      = compute_edge(1)

            case(S_edge_type)

               j_analyzed   =  bc_section(3)+1
               dir_analyzed = +1
               analyze      = compute_edge(2)

            case default
               call error_bc_section_type(
     $              'icr_interface_bc_section_module',
     $              'analyze_interior_edge_y',
     $              bc_section(1))

          end select


          ! if the time integration borders do not overlap the
          ! bc_interior_pt in the bc_section, the grid-points
          ! are analyzed
          if(analyze) then

             loc_j_analyzed =
     $            j_analyzed +
     $            dir_analyzed*dct_icr_distance

             ! loop over the bc_interior_pt in the x-direction
             i=bc_section(2)
             last_i_added=i-3
             do i=bc_section(2),bc_section(4)

                loc_central_coords = [i,loc_j_analyzed]

                !check whether the node from which the grid-point
                !depends is activated
                node_activated = is_interior_node_activated(
     $               loc_central_coords,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               p_model)

                ! if the node is activated, the current grid-point
                ! and its nearest neighbors are staged
                if(node_activated) then
                   
                   if(i.ge.(last_i_added+3)) then
                   
                     ! grid-point [i-1,j] is staged
                      if(i.gt.bc_section(2)) then

                         call icr_interface_used%stage(
     $                        [i-1,j_analyzed],
     $                        bf_interface_used)
                         
                      end if
                   end if

                   if(i.ge.(last_i_added+2)) then

                     ! grid-point [i,j] is staged
                      call icr_interface_used%stage(
     $                     [i,j_analyzed],
     $                     bf_interface_used)

                   end if

                   last_i_added = i

                   ! grid-point [i+1,j] is staged
                   if(i.lt.bc_section(4)) then

                      call icr_interface_used%stage(
     $                     [i+1,j_analyzed],
     $                     bf_interface_used)

                   end if
                      
                end if

             end do

          end if

        end subroutine analyze_interior_edge_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether the interior grid-point leads to undermined
        !> open boundary conditions
        !
        !> @date
        !> 23_03_2015 - initial version - J.L. Desmarais
        !
        !>@param loc_central_coords
        !> general coordinates of the grid-point checked
        !
        !>@param interior_x_map
        !> x-coordinates of the interior domain
        !
        !>@param interior_y_map
        !> y-coordinates of the interior domain
        !
        !>@param interior_nodes
        !> nodes of the interior domain
        !
        !>@param p_model
        !> physical model
        !
        !> @return
        !> whether the grid-point is activated
        !--------------------------------------------------------------
        function is_interior_node_activated(
     $     loc_central_coords,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     p_model)
     $     result(node_activated)

          implicit none

          integer(ikind) , dimension(2)       , intent(in) :: loc_central_coords
          real(rkind)    , dimension(nx)      , intent(in) :: interior_x_map
          real(rkind)    , dimension(ny)      , intent(in) :: interior_y_map
          real(rkind)    , dimension(nx,ny,ne), intent(in) :: interior_nodes
          type(pmodel_eq)                     , intent(in) :: p_model
          logical                                          :: node_activated

          integer(ikind) :: i_min
          integer(ikind) :: i_max
          integer(ikind) :: j_min
          integer(ikind) :: j_max

          i_min = loc_central_coords(1)-1
          i_max = loc_central_coords(1)+1
          j_min = loc_central_coords(2)-1
          j_max = loc_central_coords(2)+1

          node_activated = p_model%are_openbc_undermined(
     $         interior_x_map(i_min:i_max),
     $         interior_y_map(j_min:j_max),
     $         interior_nodes(i_min:i_max,j_min:j_max,:))

        end function is_interior_node_activated

      end module icr_interior_bc_section_module
