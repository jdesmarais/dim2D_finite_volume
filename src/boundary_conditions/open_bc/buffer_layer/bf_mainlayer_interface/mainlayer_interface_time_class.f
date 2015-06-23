      !> @file
      !> mainlayer_interface_dyn enhanced with procedures enabling
      !> the reorganization of the bc_sections initialized by
      !> the buffer layer
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> mainlayer_interface_dyn enhanced with procedures enabling
      !> the reorganization of the bc_sections initialized by
      !> the buffer layer
      !
      !> @date
      ! 10_04_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mainlayer_interface_time_class

        use bf_layer_bc_sections_class, only :
     $       bf_layer_bc_sections

        use bf_layer_bc_sections_merge_module, only :
     $       reallocate_bc_sections_for_merge,
     $       update_corner_for_merge,
     $       update_anticorner_for_merge,
     $       test_grdpts_id_config,
     $       get_edge_test_param,
     $       get_anticorner_test_param

        use bf_layer_class, only :
     $       bf_layer

        use bf_layer_extract_module, only :
     $       get_grdpts_id_from_interior

        use bf_sublayer_class, only :
     $       bf_sublayer

        use mainlayer_interface_dyn_class, only :
     $       mainlayer_interface_dyn

        use parameters_bf_layer, only :
     $       no_overlap

        use parameters_constant, only :
     $       left,
     $       right,
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        
        private
        public :: mainlayer_interface_time


        !> @class mainlayer_interface_dyn
        !> mainlayer_interface_dyn enhanced with procedures enabling
        !> the reorganization of the bc_sections initialized in the 
        !> buffer layer
        !
        !> @param update_bc_sections
        !> reorganize the bc_sections to overlap anti-corner by edges
        !> where it is possible...
        !------------------------------------------------------------
        type, extends(mainlayer_interface_dyn) :: mainlayer_interface_time

          contains

          procedure, pass :: update_bc_sections
          procedure, pass :: extract_grdpts_id_for_merge

        end type mainlayer_interface_time


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> reorganize the bc_sections of the buffer layer to match 
        !> differently the behavior at the anti-corners
        !
        !> @date
        !> 14_03_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> mainlayer_interface_dyn enhanced with procedures enabling
        !> the reorganization of the bc_sections initialized by
        !> the buffer layer
        !
        !>@param bf_layer_used
        !> buffer layer whose bc_sections are updated
        !--------------------------------------------------------------
        subroutine update_bc_sections(
     $       this,
     $       bf_layer_used)

          implicit none

          class(mainlayer_interface_time), intent(inout) :: this
          class(bf_layer)                , intent(inout) :: bf_layer_used


          integer(ikind), dimension(:,:), allocatable :: bc_sections
          logical       , dimension(2) :: side
          integer                      :: i
          integer                      :: nb_ele_removed

          type(bf_layer_bc_sections) :: bf_layer_bc_sections_used

          integer(ikind), dimension(2) :: x_borders
          integer(ikind), dimension(2) :: y_borders


          ! extract the bc_sections from the buffer layer
          call bf_layer_used%get_bc_sections(bc_sections)


          ! loop over the elements of the bc_sections and 
          ! update the elements
          if(allocated(bc_sections)) then


             ! reinitialize the overlaps in the bc_sections
             call reinitialize_overlaps(bc_sections)


             ! check whether edge can overlap anti-corners
             nb_ele_removed = 0

             side = [left,right]
             
             do i=1, size(bc_sections,2)

                select case(bc_sections(1,i))

                  ! if the bc_section is an edge, we check whether
                  ! this edge can be merged with neighboring
                  ! anti-corners and overlap corners
                  case(N_edge_type,S_edge_type,E_edge_type,W_edge_type)
                  
                     call overlap_anticorners_with_edges(
     $                    this,
     $                    bc_sections,i,
     $                    side,
     $                    bf_layer_used,
     $                    nb_ele_removed)                  

                end select
                  

             end do


             ! if the bc_sections were modified, some
             ! anti-corners should be removed and the
             ! bc_sections is reallocated
             !
             ! if the bc_sections were modified, the
             ! overlap of the bc_sections by the 
             ! integration borders should be updated
             if(nb_ele_removed.ne.0) then

                call reallocate_bc_sections_for_merge(
     $               bc_sections,
     $               nb_ele_removed)

             end if


             call bf_layer_bc_sections_used%bubble_sort(bc_sections)

             x_borders = bf_layer_used%get_x_borders()
             y_borders = bf_layer_used%get_y_borders()

             call bf_layer_bc_sections_used%check_overlaps(
     $            bc_sections,
     $            x_borders,
     $            y_borders)
             

             ! turn the anti-corners into corners
             if(allocated(bc_sections)) then
                do i=1, size(bc_sections,2)

                   select case(bc_sections(1,i))

                     ! if the bc_section is an anti-corner, we check
                     ! whether this anti-corner should behave as a
                     ! corner
                     case(NE_edge_type,NW_edge_type,SE_edge_type,SW_edge_type)

                        call identify_anticorners_as_corners(
     $                       this,
     $                       bc_sections(:,i),
     $                       side,
     $                       bf_layer_used)

                   end select

                end do
             end if


             ! set the bc_sections inside the buffer layer
             call bf_layer_used%set_bc_sections(bc_sections)


          end if

        end subroutine update_bc_sections


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> overlap the anti-corners located at the crossing between
        !> an edge and a corner by the edge
        !
        !> @date
        !> 13_04_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> mainlayer_interface_dyn enhanced with procedures enabling
        !> the reorganization of the bc_sections initialized by
        !> the buffer layer
        !
        !>@param bc_sections
        !> boundary sections identifying the edge of the computational
        !> domain
        !
        !>@param i
        !> index of the boundary section investigated
        !
        !>@param side
        !> convention for the identification of the left and right
        !
        !>@param bf_layer_used
        !> buffer layer whose bc_sections are investigated
        !
        !>@param nb_ele_removed
        !> number of elements removed when the 
        !--------------------------------------------------------------
        subroutine overlap_anticorners_with_edges(
     $     this,
     $     bc_sections,i,
     $     side,
     $     bf_layer_used,
     $     nb_ele_removed)

          implicit none

          class(mainlayer_interface_time), intent(in)    :: this
          integer(ikind), dimension(:,:) , intent(inout) :: bc_sections
          integer                        , intent(in)    :: i
          logical       , dimension(2)   , intent(in)    :: side
          type(bf_layer)                 , intent(in)    :: bf_layer_used
          integer                        , intent(inout) :: nb_ele_removed
          
          integer(ikind), dimension(2,2) :: grdpts_ex_borders
          integer(ikind), dimension(2,2) :: test_merge_loc_borders
          integer(ikind), dimension(4,4) :: test_merge_array
          integer(ikind), dimension(2,2) :: test_over_loc_borders
          integer(ikind), dimension(4,4) :: test_over_array
          integer(ikind)                 :: merge_anticorner_type
          integer(ikind), dimension(2)   :: merge_anticorner_position
          integer(ikind)                 :: over_corner_type
          integer(ikind), dimension(2)   :: over_corner_position
          integer                        :: over_corner_overlap
          integer(ikind), dimension(3)   :: edge_new_position

          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2,2) :: gen_coords
          integer       , dimension(4,4) :: grdptsid_side
          integer                        :: size_x
          integer                        :: size_y
          logical                        :: modify_edge

          integer :: j


          ! loop over the potential sides of the edge
          do j=1,2

             ! determine which grid points should be
             ! checked on the side of the bc_section
             call get_edge_test_param(
     $            bc_sections(:,i),
     $            side(j),
     $            
     $            grdpts_ex_borders,
     $            
     $            test_merge_loc_borders,
     $            test_merge_array,
     $            
     $            test_over_loc_borders,
     $            test_over_array,
     $            
     $            merge_anticorner_type,
     $            merge_anticorner_position,
     $            
     $            over_corner_type,
     $            over_corner_position,
     $            over_corner_overlap,
     $            
     $            edge_new_position)


             ! convert the borders asked expressed
             ! in the local coordinates of the buffer
             ! layer, into coordinates in the general
             ! reference frame
             match_table = bf_layer_used%get_general_to_local_coord_tab()
             gen_coords(1,1) = grdpts_ex_borders(1,1) + match_table(1)
             gen_coords(1,2) = grdpts_ex_borders(1,2) + match_table(1)
             gen_coords(2,1) = grdpts_ex_borders(2,1) + match_table(2)
             gen_coords(2,2) = grdpts_ex_borders(2,2) + match_table(2)

             size_x = gen_coords(1,2)-gen_coords(1,1)+1
             size_y = gen_coords(2,2)-gen_coords(2,1)+1


             ! extract the grid points on the side
             ! of the edge bc_section to determine
             ! whether it could overlap an anti-corner
             call this%extract_grdpts_id_for_merge(
     $            bf_layer_used,
     $            gen_coords,
     $            grdptsid_side(1:size_x,1:size_y))

             ! test whether the grdpts_id on the side
             ! leads to the merge with an anti-corner
             modify_edge = test_grdpts_id_config(
     $            grdptsid_side(
     $               test_merge_loc_borders(1,1):test_merge_loc_borders(1,2),
     $               test_merge_loc_borders(2,1):test_merge_loc_borders(2,2)),
     $            test_merge_array(
     $               test_merge_loc_borders(1,1):test_merge_loc_borders(1,2),
     $               test_merge_loc_borders(2,1):test_merge_loc_borders(2,2)))

             ! if the configuration of the grdpts_id on the side
             ! of the edge leads to the merge with an anti-corner
             ! the bc_sections are updated
             if(modify_edge) then

                call update_anticorner_for_merge(
     $               merge_anticorner_type,
     $               merge_anticorner_position,
     $               bc_sections,
     $               nb_ele_removed)

                bc_sections(2:4,i) = edge_new_position

             ! otherwise, we check whether the configuration of the
             ! grdpts_id on the side of the edge leads to the merge
             ! with an anti-corner and the overlap of a corner
             else

                ! if the configuration of the grdpts_id on the side
                ! of the edge leads to the merge with an anti-corner
                ! and the overlap of a corner, the bc_sections should
                ! be updated
                modify_edge = test_grdpts_id_config(
     $               grdptsid_side(
     $                  test_over_loc_borders(1,1):test_over_loc_borders(1,2),
     $                  test_over_loc_borders(2,1):test_over_loc_borders(2,2)),
     $               test_over_array(
     $                  test_over_loc_borders(1,1):test_over_loc_borders(1,2),
     $                  test_over_loc_borders(2,1):test_over_loc_borders(2,2)))

                if(modify_edge) then

                   call update_corner_for_merge(
     $                  merge_anticorner_type,
     $                  merge_anticorner_position,
     $                  over_corner_type,
     $                  over_corner_position,
     $                  over_corner_overlap,
     $                  bc_sections,
     $                  nb_ele_removed)

                   bc_sections(2:4,i) = edge_new_position

                end if

             end if

          end do

        end subroutine overlap_anticorners_with_edges


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether an anti-corner is surrounded by corners
        !> and so whether it should compute its grid points as if
        !> it was a corner
        !
        !> @date
        !> 13_04_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> mainlayer_interface_dyn enhanced with procedures enabling
        !> the reorganization of the bc_sections initialized by
        !> the buffer layer
        !
        !>@param bc_section
        !> boundary section investigated
        !
        !>@param side
        !> convention for the identification of the left and right
        !
        !>@param bf_layer_used
        !> buffer layer whose bc_sections are investigated
        !--------------------------------------------------------------
        subroutine identify_anticorners_as_corners(
     $     this,
     $     bc_section,
     $     side,
     $     bf_layer_used)

          implicit none

          class(mainlayer_interface_time), intent(inout) :: this
          integer(ikind), dimension(5)   , intent(inout) :: bc_section
          logical       , dimension(2)   , intent(in)    :: side
          class(bf_layer)                , intent(in)    :: bf_layer_used

          
          integer :: j

          integer(ikind), dimension(2,2) :: grdpts_ex_borders
          integer(ikind), dimension(2,2) :: test1_loc_borders
          integer       , dimension(2,2) :: test1_array
          integer       , dimension(2,2) :: test2_loc_borders
          integer       , dimension(2,2) :: test2_array
          integer                        :: new_anticorner_type

          logical :: modify_anticorner

          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2,2) :: gen_coords
          integer(ikind)                 :: size_x
          integer(ikind)                 :: size_y

          integer       , dimension(2,2) :: grdptsid_side


          do j=1,2

             ! determine the parameters for testing
             ! whether the anti-corner should be computed
             ! as a corner
             call get_anticorner_test_param(
     $            bc_section,
     $            side(j),
     $       
     $            grdpts_ex_borders,
     $            
     $            test1_loc_borders,
     $            test1_array,
     $            
     $            test2_loc_borders,
     $            test2_array,
     $            
     $            new_anticorner_type)

             ! convert the borders asked expressed
             ! in the local coordinates of the buffer
             ! layer, into coordinates in the general
             ! reference frame
             match_table = bf_layer_used%get_general_to_local_coord_tab()
             gen_coords(1,1) = grdpts_ex_borders(1,1) + match_table(1)
             gen_coords(1,2) = grdpts_ex_borders(1,2) + match_table(1)
             gen_coords(2,1) = grdpts_ex_borders(2,1) + match_table(2)
             gen_coords(2,2) = grdpts_ex_borders(2,2) + match_table(2)

             size_x = gen_coords(1,2)-gen_coords(1,1)+1
             size_y = gen_coords(2,2)-gen_coords(2,1)+1


             ! extract the grid points on the side
             ! of the edge bc_section to determine
             ! whether it could overlap an anti-corner
             call this%extract_grdpts_id_for_merge(
     $            bf_layer_used,
     $            gen_coords,
     $            grdptsid_side(1:size_x,1:size_y))

             ! test whether the grdpts_id on the side
             ! leads to considering the anti-corner as
             ! a corner
             modify_anticorner = test_grdpts_id_config(
     $            grdptsid_side(
     $               test1_loc_borders(1,1):test1_loc_borders(1,2),
     $               test1_loc_borders(2,1):test1_loc_borders(2,2)),
     $            test1_array(
     $               test1_loc_borders(1,1):test1_loc_borders(1,2),
     $               test1_loc_borders(2,1):test1_loc_borders(2,2)))

             if(.not.modify_anticorner) then

                modify_anticorner = test_grdpts_id_config(
     $            grdptsid_side(
     $               test2_loc_borders(1,1):test2_loc_borders(1,2),
     $               test2_loc_borders(2,1):test2_loc_borders(2,2)),
     $            test2_array(
     $               test2_loc_borders(1,1):test2_loc_borders(1,2),
     $               test2_loc_borders(2,1):test2_loc_borders(2,2)))

                if(.not.modify_anticorner) then
                   exit
                end if

             end if

          end do


          if(modify_anticorner) then
             bc_section(1) = new_anticorner_type
          end if


        end subroutine identify_anticorners_as_corners


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> extract the grdpts_id to check the configuration of the
        !> grid-points at the side of the edge-like bc_section
        !
        !> @date
        !> 10_04_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> mainlayer_interface_dyn enhanced with procedures enabling
        !> the reorganization of the bc_sections initialized by
        !> the buffer layer
        !
        !>@param bf_layer_used
        !> buffer layer whose bc_sections are investigated
        !
        !>@param gen_coords
        !> coordinates expressed in the general frame identifying the
        !> SW and NE corners of the grid-points extracted
        !
        !>@param grdpts_id
        !> configuration of the grid-points
        !--------------------------------------------------------------
        subroutine extract_grdpts_id_for_merge(
     $     this,
     $     bf_layer_used,
     $     gen_coords,
     $     grdpts_id)

           implicit none

           class(mainlayer_interface_time), intent(in)  :: this
           type(bf_layer)                 , intent(in)  :: bf_layer_used
           integer(ikind), dimension(2,2) , intent(in)  :: gen_coords
           integer       , dimension(:,:) , intent(out) :: grdpts_id


           type(bf_sublayer), pointer :: bf_neighbor_ptr


           ! extract the grdpts_id from the interior
           call get_grdpts_id_from_interior(
     $          grdpts_id,
     $          gen_coords)

           
           ! extract the grdpts_id from the neighboring buffer
           ! layer of type 1 if any
           if(bf_layer_used%can_exchange_with_neighbor1()) then

              bf_neighbor_ptr => this%get_neighbor_sublayer_ptr(
     $             bf_layer_used%get_localization(),1)

              call bf_neighbor_ptr%extract_grdpts_id(
     $             grdpts_id,
     $             gen_coords)

           end if


           ! extract the grdpts_id from the neighboring buffer
           ! layer of type 2 if any
           if(bf_layer_used%can_exchange_with_neighbor2()) then

              bf_neighbor_ptr => this%get_neighbor_sublayer_ptr(
     $             bf_layer_used%get_localization(),2)

              call bf_neighbor_ptr%extract_grdpts_id(
     $             grdpts_id,
     $             gen_coords)

           end if


           ! extract from the current buffer layer
           call bf_layer_used%extract_grdpts_id(
     $          grdpts_id,
     $          gen_coords)

        end subroutine extract_grdpts_id_for_merge


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> reinitialize the overlaps of the bc_sections
        !
        !> @date
        !> 11_04_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> mainlayer_interface_dyn enhanced with procedures enabling
        !> the reorganization of the bc_sections initialized by
        !> the buffer layer
        !
        !>@param bc_sections
        !> buffer layer whose bc_sections are investigated
        !
        !>@param gen_coords
        !> coordinates expressed in the general frame identifying the
        !> SW and NE corners of the grid-points extracted
        !
        !>@param grdpts_id
        !> configuration of the grid-points
        !--------------------------------------------------------------
        subroutine reinitialize_overlaps(bc_sections)

          implicit none

          integer(ikind), dimension(:,:), intent(inout) :: bc_sections

          integer :: k


          do k=1, size(bc_sections,2)

             select case(bc_sections(1,k))
               case(N_edge_type,S_edge_type,
     $              E_edge_type,W_edge_type)

                 bc_sections(5,k) = no_overlap

               case(NE_corner_type,NW_corner_type,
     $              SE_corner_type,SW_corner_type,
     $              NE_edge_type,NW_edge_type,
     $              SE_edge_type,SW_edge_type)
               
                 bc_sections(4,k) = no_overlap
                 bc_sections(5,k) = no_overlap

              case default

                 print '(''mainlayer_interface_time'')'
                 print '(''reinitialize_overlaps'')'
                 print '(''bc_section not recognized'')'
                 stop ''

            end select

          end do

        end subroutine reinitialize_overlaps
        
      end module mainlayer_interface_time_class
