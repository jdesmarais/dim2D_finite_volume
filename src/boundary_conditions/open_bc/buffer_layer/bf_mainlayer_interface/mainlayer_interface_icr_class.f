      !> @file
      !> mainlayer_interface_grdpts_id_update augmented
      !> with procedures to include which bc_interior_pt
      !> should be updated
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> mainlayer_interface_grdpts_id_update augmented
      !> with procedures to include which bc_interior_pt
      !> should be updated
      !
      !> @date
      ! 22_04_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module mainlayer_interface_icr_class
      
        use bf_layer_bc_sections_merge_module, only :
     $       test_grdpts_id_config,
     $       get_edge_test_param,
     $       get_extent_bc_section_edge

        use bf_layer_class, only :
     $       bf_layer

        use mainlayer_interface_grdpts_id_update_class, only :
     $       mainlayer_interface_grdpts_id_update

        use parameters_constant, only :
     $       left,right

        use parameters_kind, only :
     $       ikind        

        implicit none

        private
        public :: mainlayer_interface_icr

        
        !>@class mainlayer_interface_icr
        !> mainlayer_interface_grdpts_id_update augmented
        !> with procedures to include which bc_interior_pt
        !> should be updated
        !
        !> @param analyze_bc_section_edge
        !> check whether there are anti-corners at the sides of the
        !> edge to know whether the bc_interior_pt of the edge should
        !> be all updated to prevent a crenel structure inside the
        !> computational domain
        !--------------------------------------------------------------
        type, extends(mainlayer_interface_grdpts_id_update) :: mainlayer_interface_icr

          contains

          procedure, pass :: analyze_bc_section_edge

        end type mainlayer_interface_icr


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> check whether there are anti-corners at the sides of the
        !> edge to know whether the bc_interior_pt of the edge should
        !> be all updated to prevent a crenel structure inside the
        !> computational domain
        !
        !> @date
        !> 22_04_2015 - initial version - J.L. Desmarais
        !
        !> @param[in] this
        !> mainlayer_interface_grdpts_id_update augmented
        !> with procedures to include which bc_interior_pt
        !> should be updated
        !
        !> @param[in] bc_section
        !> edge bc_section
        !
        !> @param[in] bf_layer_used
        !> buffer layer analyzed to which the bc_section belongs. If this
        !> is an interior bc_section, the pointer is not associated and 
        !> the edge is not checked
        !
        !> @return update_entire_bc_section
        !> if the configuration of the edge matches the requirements, all
        !> its bc_interior_pt should be updated
        !--------------------------------------------------------------
        function analyze_bc_section_edge(
     $       this,
     $       bc_section,
     $       bf_layer_used)
     $       result(update_entire_bc_section)

          implicit none

          class(mainlayer_interface_icr), intent(in) :: this
          integer(ikind), dimension(5)  , intent(in) :: bc_section
          class(bf_layer)               , intent(in) :: bf_layer_used
          logical                                    :: update_entire_bc_section


          integer(ikind)       , parameter :: edge_extent_limit = 2
          logical, dimension(2), parameter :: side = [left,right]

          integer(ikind) :: extent
          integer        :: j

          
          ! parameters to test the anti-corners at
          ! the sides of the edge
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


          ! parameters to extract the grdpts to check
          ! whether the anti-corners exist at the sides
          ! of the edge
          integer(ikind), dimension(2)   :: match_table
          integer(ikind), dimension(2,2) :: gen_coords
          integer       , dimension(4,4) :: grdptsid_side
          integer                        :: size_x
          integer                        :: size_y
          logical                        :: modify_edge

          
          ! check whether the extent of the edge is large
          ! enough to test for anti-corners on both side
          extent = get_extent_bc_section_edge(bc_section)


          if((extent.gt.0).and.(extent.le.edge_extent_limit)) then


             update_entire_bc_section = .true.


             ! loop over the sides of the edges
             do j=1,2
             
                ! get the parameters to test whether there
                ! are anti-corners overlap with corners at
                ! each sides of the edges
                call get_edge_test_param(
     $               bc_section,
     $               side(j),
     $               
     $               grdpts_ex_borders,
     $               
     $               test_merge_loc_borders,
     $               test_merge_array,
     $            
     $               test_over_loc_borders,
     $               test_over_array,
     $               
     $               merge_anticorner_type,
     $               merge_anticorner_position,
     $               
     $               over_corner_type,
     $               over_corner_position,
     $               over_corner_overlap,
     $               
     $               edge_new_position)
             
             
                ! determine the borders of the grdpts_id
                ! to be extracted to check for anti-corners
                match_table = bf_layer_used%get_general_to_local_coord_tab()
                gen_coords(1,1) = grdpts_ex_borders(1,1) + match_table(1)
                gen_coords(1,2) = grdpts_ex_borders(1,2) + match_table(1)
                gen_coords(2,1) = grdpts_ex_borders(2,1) + match_table(2)
                gen_coords(2,2) = grdpts_ex_borders(2,2) + match_table(2)
                
                size_x = gen_coords(1,2)-gen_coords(1,1)+1
                size_y = gen_coords(2,2)-gen_coords(2,1)+1

             
                ! check whether both sides have anti-corners
                call this%extract_grdpts_id_for_merge(
     $               bf_layer_used,
     $               gen_coords,
     $               grdptsid_side(1:size_x,1:size_y))
                
             
                ! test whether the grdpts_id on the side
                ! have the configuration of an anti-corner
                ! overlap by a corner
                !
                ! ex:
                !
                !     corner       corner
                !        |   edge    |
                !       / \  __|__  / \
                !     _____ /     \ _____
                !     3 3 3|_ _ _ _|3 3 3
                !     2 2 3 3 ... 3 3 2 2
                !      |2 2 2 ... 2 2 2|
                !-----------------------------------------
                modify_edge = test_grdpts_id_config(
     $            grdptsid_side(
     $               test_over_loc_borders(1,1):test_over_loc_borders(1,2),
     $               test_over_loc_borders(2,1):test_over_loc_borders(2,2)),
     $            test_over_array(
     $               test_over_loc_borders(1,1):test_over_loc_borders(1,2),
     $               test_over_loc_borders(2,1):test_over_loc_borders(2,2)))
                
                update_entire_bc_section =
     $               update_entire_bc_section.and.modify_edge


             end do

          else
          
             update_entire_bc_section = .false.
          
          end if

        end function analyze_bc_section_edge


      end module mainlayer_interface_icr_class
