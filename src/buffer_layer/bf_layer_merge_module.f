      module bf_layer_merge_module
      
        use bf_layer_allocate_module, only : get_additional_blocks_N,
     $                                       get_additional_blocks_S,
     $                                       get_additional_blocks_E,
     $                                       get_additional_blocks_W,
     $                                       get_match_interior

        use parameters_bf_layer     , only : no_pt,
     $                                       align_N, align_S, align_E, align_W
        use parameters_constant     , only : x_direction, y_direction
        use parameters_input        , only : debug, nx, ny, ne, bc_size
        use parameters_kind         , only : ikind, rkind

        private
        public :: merge_bf_layers_N,
     $            merge_bf_layers_S,
     $            merge_bf_layers_E,
     $            merge_bf_layers_W,
     $            get_new_size



        contains

        !< merge northern buffer layers
        subroutine merge_bf_layers_N(
     $       nodes1, nodes2, interior_nodes,
     $       grdpts_id1, grdpts_id2,
     $       alignment1, alignment2, final_alignment_i)

          implicit none

          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: nodes1
          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: nodes2
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer       , dimension(:,:)     , allocatable, intent(inout) :: grdpts_id1
          integer       , dimension(:,:)     , allocatable, intent(inout) :: grdpts_id2
          integer(ikind), dimension(2,2)                  , intent(inout) :: alignment1
          integer(ikind), dimension(2,2)                  , intent(in)    :: alignment2
          integer(ikind), dimension(2,2)     , optional   , intent(in)    :: final_alignment_i


          integer(ikind), dimension(2,2)                :: final_alignment
          integer(ikind), dimension(2,2)                :: bf_alignment
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
          integer(ikind), dimension(2)                  :: new_size
          integer(ikind)                                :: outside_i_max1
          integer(ikind)                                :: outside_i_max2
          integer(ikind)                                :: interior_i_max1
          integer(ikind)                                :: interior_i_max2
          integer(ikind)                                :: interior_i_max3
          integer(ikind)                                :: i_min1
          integer(ikind)                                :: i_min3
          integer(ikind)                                :: i_min4
          integer(ikind)                                :: i_min5
          integer(ikind)                                :: i_min6
          integer(ikind)                                :: i_min8
          integer(ikind)                                :: j_min1
          integer(ikind)                                :: j_min2
          integer(ikind)                                :: interior_i_max11
          integer(ikind)                                :: interior_i_max13
          integer(ikind)                                :: interior_i_max21
          integer(ikind)                                :: interior_i_max23
          integer(ikind)                                :: interior_i_max31
          integer(ikind)                                :: interior_i_max33
          integer(ikind)                                :: i_min11, i_min13
          integer(ikind)                                :: i_min21, i_min23
          integer(ikind)                                :: i_min31, i_min33
          


          !get the new size of the tables
          !and the indices to match easily during the copy of the tables
          if(present(final_alignment_i)) then

             final_alignment      = final_alignment_i
             final_alignment(2,1) = align_N

             new_size = get_new_size(alignment1,
     $                               alignment2,
     $                               final_alignment)

             call get_match(
     $            y_direction,
     $            outside_i_max1, outside_i_max2,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $            j_min1, j_min2,
     $            alignment1, alignment2, final_alignment)

             bf_alignment(1,1) = min(alignment1(1,1),
     $                               alignment2(1,1),
     $                               final_alignment(1,1))

             bf_alignment(2,1) = align_N

             bf_alignment(1,2) = max(alignment1(1,2),
     $                               alignment2(1,2),
     $                               final_alignment(1,2))

             bf_alignment(2,2) = max(alignment1(2,2),
     $                               alignment2(2,2),
     $                               final_alignment(2,2))

          else

             new_size = get_new_size(alignment1, alignment2)

             call get_match(
     $            y_direction,
     $            outside_i_max1, outside_i_max2,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $            j_min1, j_min2,
     $            alignment1, alignment2)
             
             bf_alignment(1,1) = min(alignment1(1,1), alignment2(1,1))
             bf_alignment(2,1) = align_N
             bf_alignment(1,2) = max(alignment1(1,2), alignment2(1,2))
             bf_alignment(2,2) = max(alignment1(2,2), alignment2(2,2))

          end if


          !get the subdivision of the interior blocks
          call get_match_interior_layers(
     $         x_direction, nx,
     $         bf_alignment,
     $         i_min1, i_min3, interior_i_max1,
     $         i_min4, i_min5, interior_i_max2,
     $         i_min6, i_min8, interior_i_max3,
     $         interior_i_max11, interior_i_max13,
     $         interior_i_max21, interior_i_max23,
     $         interior_i_max31, interior_i_max33,
     $         i_min11, i_min13,
     $         i_min21, i_min23,
     $         i_min31, i_min33)


          !allocate the nodes and copy the tables
          allocate(new_nodes(new_size(1), new_size(2), ne))
          call merge_nodes_N(
     $         new_nodes, nodes1, nodes2, interior_nodes,
     $         alignment1, alignment2, bf_alignment,
     $         i_min1, i_min3, i_min4, i_min5, i_min6,
     $         interior_i_max1, interior_i_max2, interior_i_max3,
     $         j_min1, j_min2)
          deallocate(nodes2)
          call MOVE_ALLOC(new_nodes,nodes1)


          !allocate the gridpts_id and copy the tables
          allocate(new_grdpts_id(new_size(1), new_size(2)))
          call merge_grdpts_id_N(
     $         new_grdpts_id,
     $         grdpts_id1, grdpts_id2,
     $         alignment1, alignment2, bf_alignment,
     $         outside_i_max1, outside_i_max2,
     $         interior_i_max1, interior_i_max2, interior_i_max3,
     $         i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $         interior_i_max11, interior_i_max13,
     $         interior_i_max21, interior_i_max23,
     $         interior_i_max31, interior_i_max33,
     $         i_min11, i_min13, i_min21, i_min23, i_min31, i_min33,
     $         j_min1, j_min2)
          deallocate(grdpts_id2)
          call MOVE_ALLOC(new_grdpts_id,grdpts_id1)


          !update the alignment
          alignment1 = bf_alignment          

        end subroutine merge_bf_layers_N


        !< merge southern buffer layers
        subroutine merge_bf_layers_S(
     $       nodes1, nodes2, interior_nodes,
     $       grdpts_id1, grdpts_id2,
     $       alignment1, alignment2, final_alignment_i)

          implicit none

          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes1
          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes2
          real(rkind)   , dimension(nx,ny,ne)             , intent(in) :: interior_nodes
          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id1
          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id2
          integer(ikind), dimension(2,2)               , intent(inout) :: alignment1
          integer(ikind), dimension(2,2)               , intent(in)    :: alignment2
          integer(ikind), dimension(2,2)  , optional   , intent(in)    :: final_alignment_i

          integer(ikind), dimension(2,2)                :: final_alignment
          integer(ikind), dimension(2,2)                :: bf_alignment
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
          integer(ikind), dimension(2)                  :: new_size
          integer(ikind)                                :: outside_i_max1
          integer(ikind)                                :: outside_i_max2
          integer(ikind)                                :: interior_i_max1
          integer(ikind)                                :: interior_i_max2
          integer(ikind)                                :: interior_i_max3
          integer(ikind)                                :: i_min1
          integer(ikind)                                :: i_min3
          integer(ikind)                                :: i_min4
          integer(ikind)                                :: i_min5
          integer(ikind)                                :: i_min6
          integer(ikind)                                :: i_min8
          integer(ikind)                                :: j_min1
          integer(ikind)                                :: j_min2
          integer(ikind)                                :: interior_i_max11, interior_i_max13
          integer(ikind)                                :: interior_i_max21, interior_i_max23
          integer(ikind)                                :: interior_i_max31, interior_i_max33
          integer(ikind)                                :: i_min11, i_min13
          integer(ikind)                                :: i_min21, i_min23
          integer(ikind)                                :: i_min31, i_min33


          !get the new size of the tables
          !and the indices to match easily during the copy of the tables
          if(present(final_alignment_i)) then

             final_alignment      = final_alignment_i
             final_alignment(2,2) = align_S
             
             new_size = get_new_size(alignment1,
     $                               alignment2,
     $                               final_alignment)

             call get_match(
     $            y_direction,
     $            outside_i_max1, outside_i_max2,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $            j_min1, j_min2,
     $            alignment1, alignment2, final_alignment)

             bf_alignment(1,1) = min(alignment1(1,1),
     $                               alignment2(1,1),
     $                               final_alignment(1,1))

             bf_alignment(2,1) = min(alignment1(2,1),
     $                               alignment2(2,1),
     $                               final_alignment(2,1))

             bf_alignment(1,2) = max(alignment1(1,2),
     $                               alignment2(1,2),
     $                               final_alignment(1,2))

             bf_alignment(2,2) = align_S

          else

             new_size = get_new_size(alignment1, alignment2)

             call get_match(
     $            y_direction,
     $            outside_i_max1, outside_i_max2,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $            j_min1, j_min2,
     $            alignment1, alignment2)

             bf_alignment(1,1) = min(alignment1(1,1), alignment2(1,1))
             bf_alignment(2,1) = min(alignment1(2,1), alignment2(2,1))
             bf_alignment(1,2) = max(alignment1(1,2), alignment2(1,2))
             bf_alignment(2,2) = align_S           

          end if


          !get the subdivision of the interior blocks
          call get_match_interior_layers(
     $         x_direction, nx,
     $         bf_alignment,
     $         i_min1, i_min3, interior_i_max1,
     $         i_min4, i_min5, interior_i_max2,
     $         i_min6, i_min8, interior_i_max3,
     $         interior_i_max11, interior_i_max13,
     $         interior_i_max21, interior_i_max23,
     $         interior_i_max31, interior_i_max33,
     $         i_min11, i_min13,
     $         i_min21, i_min23,
     $         i_min31, i_min33)


          !allocate the nodes and copy the tables
          allocate(new_nodes(new_size(1), new_size(2), ne))
          call merge_nodes_S(
     $         new_nodes, nodes1, nodes2, interior_nodes,
     $         alignment1, alignment2, bf_alignment,
     $         i_min1, i_min3, i_min4, i_min5, i_min6,
     $         interior_i_max1, interior_i_max2, interior_i_max3)
          deallocate(nodes2)
          call MOVE_ALLOC(new_nodes,nodes1)


          !allocate the gridpts_id and copy the tables
          allocate(new_grdpts_id(new_size(1), new_size(2)))
          call merge_grdpts_id_S(
     $         new_grdpts_id,
     $         grdpts_id1, grdpts_id2,
     $         alignment1, alignment2, bf_alignment,
     $         outside_i_max1, outside_i_max2,
     $         interior_i_max1, interior_i_max2, interior_i_max3,
     $         i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $         interior_i_max11, interior_i_max13,
     $         interior_i_max21, interior_i_max23,
     $         interior_i_max31, interior_i_max33,
     $         i_min11, i_min13, i_min21, i_min23, i_min31, i_min33,
     $         j_min1, j_min2)
          deallocate(grdpts_id2)
          call MOVE_ALLOC(new_grdpts_id,grdpts_id1)

          !update the alignment
          alignment1 = bf_alignment

        end subroutine merge_bf_layers_S


        !< merge eastern buffer layers
        subroutine merge_bf_layers_E(
     $       nodes1, nodes2, interior_nodes,
     $       grdpts_id1, grdpts_id2,
     $       alignment1, alignment2, final_alignment_i)

          implicit none

          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: nodes1
          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: nodes2
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer       , dimension(:,:)     , allocatable, intent(inout) :: grdpts_id1
          integer       , dimension(:,:)     , allocatable, intent(inout) :: grdpts_id2
          integer(ikind), dimension(2,2)                  , intent(inout) :: alignment1
          integer(ikind), dimension(2,2)                  , intent(in)    :: alignment2
          integer(ikind), dimension(2,2)     , optional   , intent(in)    :: final_alignment_i


          integer(ikind), dimension(2,2)                :: final_alignment
          integer(ikind), dimension(2,2)                :: bf_alignment
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
          integer(ikind), dimension(2)                  :: new_size
          integer(ikind)                                :: outside_j_max1
          integer(ikind)                                :: outside_j_max2
          integer(ikind)                                :: interior_j_max1
          integer(ikind)                                :: interior_j_max2
          integer(ikind)                                :: interior_j_max3
          integer(ikind)                                :: j_min1
          integer(ikind)                                :: j_min3
          integer(ikind)                                :: j_min4
          integer(ikind)                                :: j_min5
          integer(ikind)                                :: j_min6
          integer(ikind)                                :: j_min8
          integer(ikind)                                :: i_min1
          integer(ikind)                                :: i_min2
          integer(ikind)                                :: interior_j_max11, interior_j_max13
          integer(ikind)                                :: interior_j_max21, interior_j_max23
          integer(ikind)                                :: interior_j_max31, interior_j_max33
          integer(ikind)                                :: j_min11, j_min13
          integer(ikind)                                :: j_min21, j_min23
          integer(ikind)                                :: j_min31, j_min33


          !get the new size of the tables
          !and the indices to match easily during the copy of the tables
          if(present(final_alignment_i)) then

             final_alignment      = final_alignment_i
             final_alignment(1,1) = align_E

             new_size = get_new_size(alignment1,
     $                               alignment2,
     $                               final_alignment)

             call get_match(
     $            x_direction,
     $            outside_j_max1, outside_j_max2,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            j_min1, j_min3, j_min4, j_min5, j_min6, j_min8,
     $            i_min1, i_min2,
     $            alignment1, alignment2, final_alignment)

             bf_alignment(1,1) = align_E

             bf_alignment(2,1) = min(alignment1(2,1),
     $                               alignment2(2,1),
     $                               final_alignment(2,1))

             bf_alignment(1,2) = max(alignment1(1,2),
     $                               alignment2(1,2),
     $                               final_alignment(1,2))

             bf_alignment(2,2) = max(alignment1(2,2),
     $                               alignment2(2,2),
     $                               final_alignment(2,2))

          else

             new_size = get_new_size(alignment1, alignment2)

             call get_match(
     $            x_direction,
     $            outside_j_max1, outside_j_max2,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            j_min1, j_min3, j_min4, j_min5, j_min6, j_min8,
     $            i_min1, i_min2,
     $            alignment1, alignment2)

             bf_alignment(1,1) = align_E
             bf_alignment(2,1) = min(alignment1(2,1), alignment2(1,2))
             bf_alignment(1,2) = max(alignment1(1,2), alignment2(1,2))
             bf_alignment(2,2) = max(alignment1(2,2), alignment2(2,2))

          end if


          !get the subdivision of the interior layers blocks
          call get_match_interior_layers(
     $            y_direction, ny,
     $            bf_alignment,
     $            j_min1, j_min3, interior_j_max1,
     $            j_min4, j_min5, interior_j_max2,
     $            j_min6, j_min8, interior_j_max3,
     $            interior_j_max11, interior_j_max13,
     $            interior_j_max21, interior_j_max23,
     $            interior_j_max31, interior_j_max33,
     $            j_min11, j_min13,
     $            j_min21, j_min23,
     $            j_min31, j_min33)


          !allocate the nodes and copy the tables
          allocate(new_nodes(new_size(1), new_size(2), ne))
          call merge_nodes_E(
     $         new_nodes, nodes1, nodes2, interior_nodes,
     $         alignment1, alignment2, bf_alignment,
     $         j_min1, j_min3, j_min4, j_min5, j_min6,
     $         interior_j_max1, interior_j_max2, interior_j_max3)
          deallocate(nodes2)
          call MOVE_ALLOC(new_nodes,nodes1)


          !allocate the gridpts_id and copy the tables
          allocate(new_grdpts_id(new_size(1), new_size(2)))
          call merge_grdpts_id_E(
     $         new_grdpts_id,
     $         grdpts_id1, grdpts_id2,
     $         alignment1, alignment2, bf_alignment,
     $         outside_j_max1, outside_j_max2,
     $         interior_j_max1, interior_j_max2, interior_j_max3,
     $         j_min1, j_min3, j_min4, j_min5, j_min6, j_min8,
     $         interior_j_max11, interior_j_max13,
     $         interior_j_max21, interior_j_max23,
     $         interior_j_max31, interior_j_max33,
     $         j_min11, j_min13, j_min21, j_min23, j_min31, j_min33)
          deallocate(grdpts_id2)
          call MOVE_ALLOC(new_grdpts_id,grdpts_id1)

          !update the alignment
          alignment1 = bf_alignment

        end subroutine merge_bf_layers_E


        !< merge western buffer layers
        subroutine merge_bf_layers_W(
     $       nodes1, nodes2, interior_nodes,
     $       grdpts_id1, grdpts_id2,
     $       alignment1, alignment2, final_alignment_i)

          implicit none

          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes1
          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes2
          real(rkind)   , dimension(nx,ny,ne)          , intent(in)    :: interior_nodes
          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id1
          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id2
          integer(ikind), dimension(2,2)               , intent(inout) :: alignment1
          integer(ikind), dimension(2,2)               , intent(in)    :: alignment2
          integer(ikind), dimension(2,2)  , optional   , intent(in)    :: final_alignment_i


          integer(ikind), dimension(2,2)                :: final_alignment
          integer(ikind), dimension(2,2)                :: bf_alignment
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
          integer(ikind), dimension(2)                  :: new_size
          integer(ikind)                                :: outside_j_max1
          integer(ikind)                                :: outside_j_max2
          integer(ikind)                                :: interior_j_max1
          integer(ikind)                                :: interior_j_max2
          integer(ikind)                                :: interior_j_max3
          integer(ikind)                                :: j_min1
          integer(ikind)                                :: j_min3
          integer(ikind)                                :: j_min4
          integer(ikind)                                :: j_min5
          integer(ikind)                                :: j_min6
          integer(ikind)                                :: j_min8
          integer(ikind)                                :: i_min1
          integer(ikind)                                :: i_min2
          integer(ikind)                                :: interior_j_max11
          integer(ikind)                                :: interior_j_max13
          integer(ikind)                                :: interior_j_max21
          integer(ikind)                                :: interior_j_max23
          integer(ikind)                                :: interior_j_max31
          integer(ikind)                                :: interior_j_max33
          integer(ikind)                                :: j_min11, j_min13
          integer(ikind)                                :: j_min21, j_min23
          integer(ikind)                                :: j_min31, j_min33



          !get the new size of the tables
          !and the indices to match easily during the copy of the tables
          if(present(final_alignment_i)) then

             final_alignment      = final_alignment_i
             final_alignment(1,2) = align_W

             new_size = get_new_size(alignment1,
     $                               alignment2,
     $                               final_alignment)

             call get_match(
     $            x_direction,
     $            outside_j_max1, outside_j_max2,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            j_min1, j_min3, j_min4, j_min5, j_min6, j_min8,
     $            i_min1, i_min2,
     $            alignment1, alignment2, final_alignment)

             bf_alignment(1,1) = min(alignment1(1,1),
     $                               alignment2(1,1),
     $                               final_alignment(1,1))

             bf_alignment(2,1) = min(alignment1(2,1),
     $                               alignment2(2,1),
     $                               final_alignment(2,1))

             bf_alignment(1,2) = align_W

             bf_alignment(2,2) = max(alignment1(2,2),
     $                               alignment2(2,2),
     $                               final_alignment(2,2))

          else

             new_size = get_new_size(alignment1, alignment2)

             call get_match(
     $            x_direction,
     $            outside_j_max1, outside_j_max2,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            j_min1, j_min3, j_min4, j_min5, j_min6, j_min8,
     $            i_min1, i_min2,
     $            alignment1, alignment2)

             bf_alignment(1,1) = min(alignment1(1,1), alignment2(1,1))
             bf_alignment(2,1) = min(alignment1(2,1), alignment2(2,1))
             bf_alignment(1,2) = align_W
             bf_alignment(2,2) = max(alignment1(2,2), alignment2(2,2))

          end if


          !get the subdivision of the interior layers blocks
          call get_match_interior_layers(
     $            y_direction, ny,
     $            bf_alignment,
     $            j_min1, j_min3, interior_j_max1,
     $            j_min4, j_min5, interior_j_max2,
     $            j_min6, j_min8, interior_j_max3,
     $            interior_j_max11, interior_j_max13,
     $            interior_j_max21, interior_j_max23,
     $            interior_j_max31, interior_j_max33,
     $            j_min11, j_min13,
     $            j_min21, j_min23,
     $            j_min31, j_min33)


          !allocate the nodes and copy the tables
          allocate(new_nodes(new_size(1), new_size(2), ne))
          call merge_nodes_W(
     $         new_nodes, nodes1, nodes2, interior_nodes,
     $         alignment1, alignment2, bf_alignment,
     $         j_min1, j_min3, j_min4, j_min5, j_min6,
     $         interior_j_max1, interior_j_max2, interior_j_max3)
          deallocate(nodes2)
          call MOVE_ALLOC(new_nodes,nodes1)


          !allocate the gridpts_id and copy the tables
          allocate(new_grdpts_id(new_size(1), new_size(2)))
          call merge_grdpts_id_W(
     $         new_grdpts_id,
     $         grdpts_id1, grdpts_id2,
     $         alignment1, alignment2, bf_alignment,
     $         outside_j_max1, outside_j_max2,
     $         interior_j_max1, interior_j_max2, interior_j_max3,
     $         j_min1, j_min3, j_min4, j_min5, j_min6, j_min8,
     $         interior_j_max11, interior_j_max13,
     $         interior_j_max21, interior_j_max23,
     $         interior_j_max31, interior_j_max33,
     $         j_min11, j_min13, j_min21, j_min23, j_min31, j_min33)
          deallocate(grdpts_id2)
          call MOVE_ALLOC(new_grdpts_id,grdpts_id1)


          !update the alignment
          alignment1 = bf_alignment

        end subroutine merge_bf_layers_W


        !< get the new size of the nodes and gridpts_id tables
        !> after merging northern and southern buffer layers
        function get_new_size(alignment1, alignment2, final_alignment)
     $     result(new_size)

          implicit none

          integer(ikind), dimension(2,2)          , intent(in)  :: alignment1
          integer(ikind), dimension(2,2)          , intent(in)  :: alignment2
          integer(ikind), dimension(2,2), optional, intent(in)  :: final_alignment
          integer(ikind), dimension(2)                          :: new_size

          if(present(final_alignment)) then
             new_size(1) = max(alignment1(1,2),alignment2(1,2),final_alignment(1,2))-
     $                     min(alignment1(1,1),alignment2(1,1),final_alignment(1,1))+
     $                     2*bc_size+1
             new_size(2) = max(alignment1(2,2),alignment2(2,2),final_alignment(2,2))-
     $                     min(alignment1(2,1),alignment2(2,1),final_alignment(2,1))+
     $                     2*bc_size+1
          else
             new_size(1) = max(alignment1(1,2),alignment2(1,2))-
     $                     min(alignment1(1,1),alignment2(1,1))+
     $                     2*bc_size+1
             new_size(2) = max(alignment1(2,2),alignment2(2,2))-
     $                     min(alignment1(2,1),alignment2(2,1))+
     $                     2*bc_size+1
          end if


        end function get_new_size


        !< get the indices needed to match the tables copied
        !> when merging southern and northerm buffer layers
        subroutine get_match(
     $     dir,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $     j_min1, j_min2,
     $     alignment1, alignment2, final_alignment)

          implicit none

          integer                                 , intent(in)  :: dir
          integer(ikind)                          , intent(out) :: outside_i_max1
          integer(ikind)                          , intent(out) :: outside_i_max2
          integer(ikind)                          , intent(out) :: interior_i_max1
          integer(ikind)                          , intent(out) :: interior_i_max2
          integer(ikind)                          , intent(out) :: interior_i_max3
          integer(ikind)                          , intent(out) :: i_min1
          integer(ikind)                          , intent(out) :: i_min3
          integer(ikind)                          , intent(out) :: i_min4
          integer(ikind)                          , intent(out) :: i_min5
          integer(ikind)                          , intent(out) :: i_min6
          integer(ikind)                          , intent(out) :: i_min8
          integer(ikind)                          , intent(out) :: j_min1
          integer(ikind)                          , intent(out) :: j_min2
          integer(ikind), dimension(2,2)          , intent(in)  :: alignment1
          integer(ikind), dimension(2,2)          , intent(in)  :: alignment2
          integer(ikind), dimension(2,2), optional, intent(in)  :: final_alignment


          integer(ikind)                 :: size1, size2, ndir
          integer                        :: dir1, dir2
          integer(ikind), dimension(2,2) :: border_changes

          select case(dir)
            case(y_direction)
               dir1 = x_direction
               dir2 = y_direction
               ndir = nx
            case(x_direction)
               dir1 = y_direction
               dir2 = x_direction
               ndir = ny
            case default
               print '(''bf_layer_merge_module'')'
               print '(''get_match'')'
               print '(''dir not recognized'')'
               print '(''dir: '',I2)', dir
               stop 'either use x_direction or y_direction'
          end select


          if(present(final_alignment)) then

             !compute the border changes due to the final alignment
             border_changes(1,1) = final_alignment(1,1) - min(alignment1(1,1),alignment2(1,1))
             border_changes(2,1) = final_alignment(2,1) - min(alignment1(2,1),alignment2(2,1))
             border_changes(1,2) = final_alignment(1,2) - max(alignment1(1,2),alignment2(1,2))
             border_changes(2,2) = final_alignment(2,2) - max(alignment1(2,2),alignment2(2,2))


             !define the length of the block 1
             !the length is restricted by the length of the
             !interior domain
             if(min(alignment1(dir1,1),alignment2(dir1,1)).le.(bc_size+1)) then
                outside_i_max1 = max(0, -border_changes(dir1,1))
             else
                if(final_alignment(dir1,1).le.(bc_size+1)) then
                   outside_i_max1 = bc_size + 1 - final_alignment(dir1,1)
                else
                   outside_i_max1 = 0
                end if
             end if


             !define the length of the blocks 9
             !the length is restricted by the length of the
             !interior domain
             if(max(alignment1(dir1,2),alignment2(dir1,2)).ge.(ndir-bc_size)) then
                outside_i_max2 = max(0, border_changes(dir1,2))
             else
                if(final_alignment(dir1,2).ge.(ndir-bc_size)) then
                   outside_i_max2 = final_alignment(dir1,2)-(ndir-bc_size)
                else
                   outside_i_max2 = 0
                end if
             end if


             !define the length of the blocks 2+3
             !the length is restricted by the length of the
             !interior domain
             if(min(alignment1(dir1,1),alignment2(dir1,1)).ge.(bc_size+1)) then
                if(final_alignment(dir1,1).ge.(bc_size+1)) then
                   interior_i_max1 = max(0,-border_changes(dir1,1))
                else
                   interior_i_max1 = min(alignment1(dir1,1),alignment2(dir1,1))-(bc_size+1)
                end if
             else
                interior_i_max1 = 0
             end if


             !define the length of the blocks 7+8
             !the length is restricted by the length of the
             !interior domain
             if(max(alignment1(dir1,2),alignment2(dir1,2)).le.(ndir-bc_size)) then
                if(final_alignment(dir1,2).le.(ndir-bc_size)) then
                   interior_i_max3 = max(0, border_changes(dir1,2))
                else
                   interior_i_max3 = (ndir-bc_size) - max(alignment1(dir1,2),alignment2(dir1,2))
                end if
             else
                interior_i_max3 = 0
             end if

          else
             outside_i_max1  = 0
             outside_i_max2  = 0
             interior_i_max1 = 0
             interior_i_max3 = 0
          end if             

          !define the length of the blocks 5, 12, and 17
          interior_i_max2 = max(alignment1(dir1,1), alignment2(dir1,1)) -
     $                      min(alignment1(dir1,2), alignment2(dir1,2)) -
     $                      (2*bc_size+1)

          if(debug) then
             if(interior_i_max2.lt.0) then
                print '(''bf_layer_merge_module'')'
                print '(''get_match_NS'')'
                print '(''the two tables are superposed'')'
                stop 'check the alignment of the two tables'
             end if
          end if
          
          if(alignment1(dir1,1).lt.alignment2(dir1,1)) then
             size1 = alignment1(dir1,2)-alignment1(dir1,1) + 2*bc_size + 1
             size2 = alignment2(dir1,2)-alignment2(dir1,1) + 2*bc_size + 1
          else
             size1 = alignment2(dir1,2)-alignment2(dir1,1) + 2*bc_size + 1
             size2 = alignment1(dir1,2)-alignment1(dir1,1) + 2*bc_size + 1
          end if

          i_min1 = outside_i_max1
          i_min3 = i_min1 + interior_i_max1
          i_min4 = i_min3 + size1
          i_min5 = i_min4 + interior_i_max2
          i_min6 = i_min5 + size2
          i_min8 = i_min6 + interior_i_max3

          j_min1 = min((alignment1(dir2,2)-alignment1(dir2,1)+2*bc_size+1),
     $                 (alignment2(dir2,2)-alignment2(dir2,1)+2*bc_size+1))
          j_min2 = max((alignment1(dir2,2)-alignment1(dir2,1)+2*bc_size+1),
     $                 (alignment2(dir2,2)-alignment2(dir2,1)+2*bc_size+1))

        end subroutine get_match


        subroutine get_match_interior_layers(
     $     dir, ndir,
     $     bf_alignment,
     $     i_min1, i_min3, interior_i_max1,
     $     i_min4, i_min5, interior_i_max2,
     $     i_min6, i_min8, interior_i_max3,
     $     interior_i_max11, interior_i_max13,
     $     interior_i_max21, interior_i_max23,
     $     interior_i_max31, interior_i_max33,
     $     i_min11, i_min13,
     $     i_min21, i_min23,
     $     i_min31, i_min33)

          implicit none

          integer                       , intent(in)  :: dir
          integer(ikind)                , intent(in)  :: ndir
          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
          integer(ikind)                , intent(in)  :: i_min1
          integer(ikind)                , intent(in)  :: i_min3
          integer(ikind)                , intent(in)  :: interior_i_max1
          integer(ikind)                , intent(in)  :: i_min4
          integer(ikind)                , intent(in)  :: i_min5
          integer(ikind)                , intent(in)  :: interior_i_max2
          integer(ikind)                , intent(in)  :: i_min6
          integer(ikind)                , intent(in)  :: i_min8
          integer(ikind)                , intent(out) :: interior_i_max3
          integer(ikind)                , intent(out) :: interior_i_max11
          integer(ikind)                , intent(out) :: interior_i_max13
          integer(ikind)                , intent(out) :: interior_i_max21
          integer(ikind)                , intent(out) :: interior_i_max23
          integer(ikind)                , intent(out) :: interior_i_max31
          integer(ikind)                , intent(out) :: interior_i_max33
          integer(ikind)                , intent(out) :: i_min11
          integer(ikind)                , intent(out) :: i_min13
          integer(ikind)                , intent(out) :: i_min21
          integer(ikind)                , intent(out) :: i_min23
          integer(ikind)                , intent(out) :: i_min31
          integer(ikind)                , intent(out) :: i_min33


          !define the lengths of the blocks 2-3-4
          call get_match_interior(
     $         dir, ndir,
     $         bf_alignment, i_min1, interior_i_max1,
     $         interior_i_max11, interior_i_max13)

          !update the block indices
          i_min11 = i_min1 + interior_i_max11
          i_min13 = i_min3 - interior_i_max13


          !define the lengths of the blocks 6-7-8
          call get_match_interior(
     $         dir, ndir,
     $         bf_alignment, i_min4, interior_i_max2,
     $         interior_i_max21, interior_i_max23)

          !update the block indices
          i_min21 = i_min4 + interior_i_max21
          i_min23 = i_min5 - interior_i_max23


          !define the lengths of the blocks 6-7-8
          call get_match_interior(
     $         dir, ndir,
     $         bf_alignment, i_min6, interior_i_max3,
     $         interior_i_max31, interior_i_max33)

          !update the block indices
          i_min31 = i_min6 + interior_i_max31
          i_min33 = i_min8 - interior_i_max33

        end subroutine get_match_interior_layers

        
        !> merge the nodes for northern buffer layers
        subroutine merge_nodes_N(
     $     new_nodes,
     $     nodes1, nodes2, interior_nodes,
     $     alignment1, alignment2, bf_alignment,
     $     i_min1, i_min3, i_min4, i_min5, i_min6,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     j_min1, j_min2)

          implicit none

          real(rkind), dimension(:,:,:)   , intent(out):: new_nodes
          real(rkind), dimension(:,:,:)   , intent(in) :: nodes1
          real(rkind), dimension(:,:,:)   , intent(in) :: nodes2
          real(rkind), dimension(nx,ny,ne), intent(in) :: interior_nodes
          integer(ikind), dimension(2,2)  , intent(in) :: alignment1
          integer(ikind), dimension(2,2)  , intent(in) :: alignment2
          integer(ikind), dimension(2,2)  , intent(in) :: bf_alignment
          integer(ikind)                  , intent(in) :: i_min1
          integer(ikind)                  , intent(in) :: i_min3
          integer(ikind)                  , intent(in) :: i_min4
          integer(ikind)                  , intent(in) :: i_min5
          integer(ikind)                  , intent(in) :: i_min6
          integer(ikind)                  , intent(in) :: interior_i_max1
          integer(ikind)                  , intent(in) :: interior_i_max2
          integer(ikind)                  , intent(in) :: interior_i_max3
          integer(ikind)                  , intent(in) :: j_min1
          integer(ikind)                  , intent(in) :: j_min2

          integer :: k
          integer(ikind) :: j_match1, j_match2, j_matchI

          j_match1 = 0
          j_match2 = 0
          j_matchI = ny-(2*bc_size)

          !nodes1 - nodes2
          if(alignment1(1,1).lt.alignment2(1,1)) then

             if(alignment1(2,2).gt.alignment2(2,2)) then

                do k=1,ne
                   
                   call add_nodes_blocks_2_to_12_NS(
     $                  new_nodes,
     $                  nodes1, nodes2, interior_nodes,
     $                  bf_alignment,
     $                  k, 1, 2*bc_size,
     $                  j_match1, j_match2, j_matchI,
     $                  interior_i_max1,interior_i_max2,interior_i_max3,
     $                  i_min1, i_min3, i_min4, i_min5, i_min6)

                   call add_nodes_blocks_15_to_17_NS(
     $                  new_nodes,
     $                  nodes1, nodes2,
     $                  k, 2*bc_size+1, j_min1,
     $                  j_match1, j_match2,
     $                  i_min3, i_min5)
                   
                   call add_nodes_blocks_20_to_22_NS(
     $                  new_nodes,
     $                  nodes1,
     $                  k, j_min1+1, j_min2,
     $                  j_match1,
     $                  i_min3)

                end do

             else
                do k=1, ne
                   call add_nodes_blocks_2_to_12_NS(
     $                  new_nodes,
     $                  nodes1, nodes2, interior_nodes,
     $                  bf_alignment,
     $                  k, 1, 2*bc_size,
     $                  j_match1, j_match2, j_matchI,
     $                  interior_i_max1,interior_i_max2,interior_i_max3,
     $                  i_min1, i_min3, i_min4, i_min5, i_min6)

                   call add_nodes_blocks_15_to_17_NS(
     $                  new_nodes,
     $                  nodes1, nodes2,
     $                  k, 2*bc_size+1, j_min1,
     $                  j_match1, j_match2,
     $                  i_min3, i_min5)
                   
                   call add_nodes_blocks_20_to_22_NS(
     $                  new_nodes,
     $                  nodes2,
     $                  k, j_min1+1, j_min2,
     $                  j_match2,
     $                  i_min5)
                end do

             end if

          !nodes2 - nodes1
          else
             if(alignment1(2,2).gt.alignment2(2,2)) then
                do k=1, ne

                   call add_nodes_blocks_2_to_12_NS(
     $                  new_nodes,
     $                  nodes2, nodes1, interior_nodes,
     $                  bf_alignment,
     $                  k, 1, 2*bc_size,
     $                  j_match2, j_match1, j_matchI,
     $                  interior_i_max1,interior_i_max2,interior_i_max3,
     $                  i_min1, i_min3, i_min4, i_min5, i_min6)

                   call add_nodes_blocks_15_to_17_NS(
     $                  new_nodes,
     $                  nodes2, nodes1,
     $                  k, 2*bc_size+1, j_min1,
     $                  j_match2, j_match1,
     $                  i_min3, i_min5)
                   
                   call add_nodes_blocks_20_to_22_NS(
     $                  new_nodes,
     $                  nodes1,
     $                  k, j_min1+1, j_min2,
     $                  j_match1,
     $                  i_min5)

                end do
             else
                do k=1, ne
                   call add_nodes_blocks_2_to_12_NS(
     $                  new_nodes,
     $                  nodes2, nodes1, interior_nodes,
     $                  bf_alignment,
     $                  k, 1, 2*bc_size,
     $                  j_match2, j_match1, j_matchI,
     $                  interior_i_max1,interior_i_max2,interior_i_max3,
     $                  i_min1, i_min3, i_min4, i_min5, i_min6)

                   call add_nodes_blocks_15_to_17_NS(
     $                  new_nodes,
     $                  nodes2, nodes1,
     $                  k, 2*bc_size+1, j_min1,
     $                  j_match2, j_match1,
     $                  i_min3, i_min5)
                   
                   call add_nodes_blocks_20_to_22_NS(
     $                  new_nodes,
     $                  nodes2,
     $                  k, j_min1+1, j_min2,
     $                  j_match2,
     $                  i_min3)
                end do
             end if
          end if

        end subroutine merge_nodes_N



        !> merge the nodes for southern buffer layers
        subroutine merge_nodes_S(
     $     new_nodes,
     $     nodes1, nodes2, interior_nodes,
     $     alignment1, alignment2, bf_alignment,
     $     i_min1, i_min3, i_min4, i_min5, i_min6,
     $     interior_i_max1, interior_i_max2, interior_i_max3)

          implicit none

          real(rkind), dimension(:,:,:)   , intent(out):: new_nodes
          real(rkind), dimension(:,:,:)   , intent(in) :: nodes1
          real(rkind), dimension(:,:,:)   , intent(in) :: nodes2
          real(rkind), dimension(nx,ny,ne), intent(in) :: interior_nodes
          integer(ikind), dimension(2,2)  , intent(in) :: alignment1
          integer(ikind), dimension(2,2)  , intent(in) :: alignment2
          integer(ikind), dimension(2,2)  , intent(in) :: bf_alignment
          integer(ikind)                  , intent(in) :: i_min1
          integer(ikind)                  , intent(in) :: i_min3
          integer(ikind)                  , intent(in) :: i_min4
          integer(ikind)                  , intent(in) :: i_min5
          integer(ikind)                  , intent(in) :: i_min6
          integer(ikind)                  , intent(in) :: interior_i_max1
          integer(ikind)                  , intent(in) :: interior_i_max2
          integer(ikind)                  , intent(in) :: interior_i_max3

          integer        :: k
          integer(ikind) :: j_match1, j_match2, j_matchI

          j_match1 = - size(new_nodes,2)+size(nodes1,2)
          j_match2 = - size(new_nodes,2)+size(nodes2,2)
          j_matchI = - size(new_nodes,2)+(2*bc_size)
          
          !nodes1 - nodes2
          if(alignment1(1,1).lt.alignment2(1,1)) then

             if(alignment1(2,1).lt.alignment2(2,1)) then

                do k=1,ne

                   call add_nodes_blocks_20_to_22_NS(
     $                  new_nodes,
     $                  nodes1,
     $                  k,
     $                  size(new_nodes,2)-size(nodes1,2)+1,
     $                  size(new_nodes,2)-size(nodes2,2),
     $                  j_match1,
     $                  i_min3)

                   call add_nodes_blocks_15_to_17_NS(
     $                  new_nodes,
     $                  nodes1, nodes2,
     $                  k,
     $                  size(new_nodes,2)-size(nodes2,2)+1,
     $                  size(new_nodes,2)-(2*bc_size),
     $                  j_match1, j_match2,
     $                  i_min3, i_min5)
                   
                   call add_nodes_blocks_2_to_12_NS(
     $                  new_nodes,
     $                  nodes1, nodes2, interior_nodes,
     $                  bf_alignment,
     $                  k,
     $                  size(new_nodes,2)-(2*bc_size)+1,
     $                  size(new_nodes,2),
     $                  j_match1, j_match2, j_matchI,
     $                  interior_i_max1,interior_i_max2,interior_i_max3,
     $                  i_min1, i_min3, i_min4, i_min5, i_min6)                   

                end do

             else                

                do k=1, ne
                   call add_nodes_blocks_20_to_22_NS(
     $                  new_nodes,
     $                  nodes2,
     $                  k,
     $                  size(new_nodes,2)-size(nodes2,2)+1,
     $                  size(new_nodes,2)-size(nodes1,2),
     $                  j_match2,
     $                  i_min5)

                   call add_nodes_blocks_15_to_17_NS(
     $                  new_nodes,
     $                  nodes1, nodes2,
     $                  k,
     $                  size(new_nodes,2)-size(nodes1,2)+1,
     $                  size(new_nodes,2)-(2*bc_size),
     $                  j_match1, j_match2,
     $                  i_min3, i_min5)

                   call add_nodes_blocks_2_to_12_NS(
     $                  new_nodes,
     $                  nodes1, nodes2, interior_nodes,
     $                  bf_alignment,
     $                  k,
     $                  size(new_nodes,2)-(2*bc_size)+1,
     $                  size(new_nodes,2),
     $                  j_match1, j_match2, j_matchI,
     $                  interior_i_max1,interior_i_max2,interior_i_max3,
     $                  i_min1, i_min3, i_min4, i_min5, i_min6)                   
                   
                end do

             end if

          !nodes2 - nodes1
          else
             if(alignment1(2,1).lt.alignment2(2,1)) then
                do k=1, ne

                   call add_nodes_blocks_20_to_22_NS(
     $                  new_nodes,
     $                  nodes1,
     $                  k,
     $                  size(new_nodes,2)-size(nodes1,2)+1,
     $                  size(new_nodes,2)-size(nodes2,2),
     $                  j_match1,
     $                  i_min5)

                   call add_nodes_blocks_15_to_17_NS(
     $                  new_nodes,
     $                  nodes2, nodes1,
     $                  k,
     $                  size(new_nodes,2)-size(nodes2,2)+1,
     $                  size(new_nodes,2)-(2*bc_size),
     $                  j_match2, j_match1,
     $                  i_min3, i_min5)

                   call add_nodes_blocks_2_to_12_NS(
     $                  new_nodes,
     $                  nodes2, nodes1, interior_nodes,
     $                  bf_alignment,
     $                  k,
     $                  size(new_nodes,2)-(2*bc_size)+1,
     $                  size(new_nodes,2),
     $                  j_match2, j_match1, j_matchI,
     $                  interior_i_max1,interior_i_max2,interior_i_max3,
     $                  i_min1, i_min3, i_min4, i_min5, i_min6)                   

                end do
             else
                do k=1, ne

                   call add_nodes_blocks_20_to_22_NS(
     $                  new_nodes,
     $                  nodes2,
     $                  k,
     $                  size(new_nodes,2)-size(nodes2,2)+1,
     $                  size(new_nodes,2)-size(nodes1,2),
     $                  j_match2,
     $                  i_min3)

                   call add_nodes_blocks_15_to_17_NS(
     $                  new_nodes,
     $                  nodes2, nodes1,
     $                  k,
     $                  size(new_nodes,2)-size(nodes1,2)+1,
     $                  size(new_nodes,2)-(2*bc_size),
     $                  j_match2, j_match1,
     $                  i_min3, i_min5)

                   call add_nodes_blocks_2_to_12_NS(
     $                  new_nodes,
     $                  nodes2, nodes1, interior_nodes,
     $                  bf_alignment,
     $                  k,
     $                  size(new_nodes,2)-(2*bc_size)+1,
     $                  size(new_nodes,2),
     $                  j_match2, j_match1, j_matchI,
     $                  interior_i_max1,interior_i_max2,interior_i_max3,
     $                  i_min1, i_min3, i_min4, i_min5, i_min6)                   
                   
                end do
             end if
          end if

        end subroutine merge_nodes_S


        !> merge the nodes for southern buffer layers
        subroutine merge_nodes_E(
     $     new_nodes,
     $     nodes1, nodes2, interior_nodes,
     $     alignment1, alignment2, bf_alignment,
     $     j_min1, j_min3, j_min4, j_min5, j_min6,
     $     interior_j_max1, interior_j_max2, interior_j_max3)

          implicit none

          real(rkind), dimension(:,:,:)   , intent(out):: new_nodes
          real(rkind), dimension(:,:,:)   , intent(in) :: nodes1
          real(rkind), dimension(:,:,:)   , intent(in) :: nodes2
          real(rkind), dimension(nx,ny,ne), intent(in) :: interior_nodes
          integer(ikind), dimension(2,2)  , intent(in) :: alignment1
          integer(ikind), dimension(2,2)  , intent(in) :: alignment2
          integer(ikind), dimension(2,2)  , intent(in) :: bf_alignment
          integer(ikind)                  , intent(in) :: j_min1
          integer(ikind)                  , intent(in) :: j_min3
          integer(ikind)                  , intent(in) :: j_min4
          integer(ikind)                  , intent(in) :: j_min5
          integer(ikind)                  , intent(in) :: j_min6
          integer(ikind)                  , intent(in) :: interior_j_max1
          integer(ikind)                  , intent(in) :: interior_j_max2
          integer(ikind)                  , intent(in) :: interior_j_max3

          integer        :: k
          integer(ikind) :: i_matchN, i_matchI

          i_matchN = 0
          i_matchI = nx-(2*bc_size)
          
          !nodes1 - nodes2
          if(alignment1(2,1).lt.alignment2(2,1)) then
             
             do k=1, ne
                call add_nodes_interior_blocks_EW(
     $               new_nodes, interior_nodes,
     $               bf_alignment,
     $               k, j_min1+1, j_min1+interior_j_max1,
     $               i_matchN, i_matchI)
                
                call add_nodes_sublayer_block_EW(
     $               new_nodes, nodes1,
     $               k, j_min3+1, j_min3+size(nodes1,2),
     $               i_matchN)
                
                call add_nodes_interior_blocks_EW(
     $               new_nodes, interior_nodes,
     $               bf_alignment,
     $               k, j_min4+1, j_min4+interior_j_max2,
     $               i_matchN, i_matchI)
                
                call add_nodes_sublayer_block_EW(
     $               new_nodes, nodes2,
     $               k, j_min5+1, j_min5+size(nodes2,2),
     $               i_matchN)
                
                call add_nodes_interior_blocks_EW(
     $               new_nodes, interior_nodes,
     $               bf_alignment,
     $               k, j_min6+1, j_min6+interior_j_max3,
     $               i_matchN, i_matchI)

             end do
                
          !nodes2 - nodes1
          else
             
             do k=1, ne

                call add_nodes_interior_blocks_EW(
     $               new_nodes, interior_nodes,
     $               bf_alignment,
     $               k, j_min1+1, j_min1+interior_j_max1,
     $               i_matchN, i_matchI)
                
                call add_nodes_sublayer_block_EW(
     $               new_nodes, nodes2,
     $               k, j_min3+1, j_min3+size(nodes2,2),
     $               i_matchN)

                call add_nodes_interior_blocks_EW(
     $               new_nodes, interior_nodes,
     $               bf_alignment,
     $               k, j_min4+1, j_min4+interior_j_max2,
     $               i_matchN, i_matchI)

                call add_nodes_sublayer_block_EW(
     $               new_nodes, nodes1,
     $               k, j_min5+1, j_min5+size(nodes1,2),
     $               i_matchN)

                call add_nodes_interior_blocks_EW(
     $               new_nodes, interior_nodes,
     $               bf_alignment,
     $               k, j_min6+1, j_min6+interior_j_max3,
     $               i_matchN, i_matchI)

             end do

          end if

        end subroutine merge_nodes_E


        !> merge the nodes for western buffer layers
        subroutine merge_nodes_W(
     $     new_nodes,
     $     nodes1, nodes2, interior_nodes,
     $     alignment1, alignment2, bf_alignment,
     $     j_min1, j_min3, j_min4, j_min5, j_min6,
     $     interior_j_max1, interior_j_max2, interior_j_max3)

          implicit none

          real(rkind), dimension(:,:,:)   , intent(out):: new_nodes
          real(rkind), dimension(:,:,:)   , intent(in) :: nodes1
          real(rkind), dimension(:,:,:)   , intent(in) :: nodes2
          real(rkind), dimension(nx,ny,ne), intent(in) :: interior_nodes
          integer(ikind), dimension(2,2)  , intent(in) :: alignment1
          integer(ikind), dimension(2,2)  , intent(in) :: alignment2
          integer(ikind), dimension(2,2)  , intent(in) :: bf_alignment
          integer(ikind)                  , intent(in) :: j_min1
          integer(ikind)                  , intent(in) :: j_min3
          integer(ikind)                  , intent(in) :: j_min4
          integer(ikind)                  , intent(in) :: j_min5
          integer(ikind)                  , intent(in) :: j_min6
          integer(ikind)                  , intent(in) :: interior_j_max1
          integer(ikind)                  , intent(in) :: interior_j_max2
          integer(ikind)                  , intent(in) :: interior_j_max3

          integer        :: k
          integer(ikind) :: i_match1, i_match2, i_matchN, i_matchI

          i_match1 = size(new_nodes,1)-size(nodes1,1)
          i_match2 = size(new_nodes,1)-size(nodes2,1)
          i_matchN = size(new_nodes,1)-(2*bc_size)
          i_matchI = 0
          
          !nodes1 - nodes2
          if(alignment1(2,1).lt.alignment2(2,1)) then
             
             do k=1, ne
                call add_nodes_interior_blocks_EW(
     $               new_nodes, interior_nodes,
     $               bf_alignment,
     $               k, j_min1+1, j_min1+interior_j_max1,
     $               i_matchN, i_matchI)
                
                call add_nodes_sublayer_block_EW(
     $               new_nodes, nodes1,
     $               k, j_min3+1, j_min3+size(nodes1,2),
     $               i_match1)
                
                call add_nodes_interior_blocks_EW(
     $               new_nodes, interior_nodes,
     $               bf_alignment,
     $               k, j_min4+1, j_min4+interior_j_max2,
     $               i_matchN, i_matchI)
                
                call add_nodes_sublayer_block_EW(
     $               new_nodes, nodes2,
     $               k, j_min5+1, j_min5+size(nodes2,2),
     $               i_match2)
                
                call add_nodes_interior_blocks_EW(
     $               new_nodes, interior_nodes,
     $               bf_alignment,
     $               k, j_min6+1, j_min6+interior_j_max3,
     $               i_matchN, i_matchI)

             end do
                
          !nodes2 - nodes1
          else
             
             do k=1, ne

                call add_nodes_interior_blocks_EW(
     $               new_nodes, interior_nodes,
     $               bf_alignment,
     $               k, j_min1+1, j_min1+interior_j_max1,
     $               i_matchN, i_matchI)
                
                call add_nodes_sublayer_block_EW(
     $               new_nodes, nodes2,
     $               k, j_min3+1, j_min3+size(nodes2,2),
     $               i_match2)

                call add_nodes_interior_blocks_EW(
     $               new_nodes, interior_nodes,
     $               bf_alignment,
     $               k, j_min4+1, j_min4+interior_j_max2,
     $               i_matchN, i_matchI)

                call add_nodes_sublayer_block_EW(
     $               new_nodes, nodes1,
     $               k, j_min5+1, j_min5+size(nodes1,2),
     $               i_match1)

                call add_nodes_interior_blocks_EW(
     $               new_nodes, interior_nodes,
     $               bf_alignment,
     $               k, j_min6+1, j_min6+interior_j_max3,
     $               i_matchN, i_matchI)

             end do

          end if

        end subroutine merge_nodes_W


        !> merge the nodes for northern buffer layers
        subroutine merge_grdpts_id_N(
     $     new_grdpts_id,
     $     grdpts_id1, grdpts_id2,
     $     alignment1, alignment2, bf_alignment,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $     interior_i_max11, interior_i_max13,
     $     interior_i_max21, interior_i_max23,
     $     interior_i_max31, interior_i_max33,
     $     i_min11, i_min13, i_min21, i_min23, i_min31, i_min33,
     $     j_min1, j_min2)

          implicit none

          integer       , dimension(:,:), intent(out):: new_grdpts_id
          integer       , dimension(:,:), intent(in) :: grdpts_id1
          integer       , dimension(:,:), intent(in) :: grdpts_id2
          integer(ikind), dimension(2,2), intent(in) :: alignment1
          integer(ikind), dimension(2,2), intent(in) :: alignment2          
          integer(ikind), dimension(2,2), intent(in) :: bf_alignment
          integer(ikind)                , intent(in) :: outside_i_max1
          integer(ikind)                , intent(in) :: outside_i_max2
          integer(ikind)                , intent(in) :: interior_i_max1
          integer(ikind)                , intent(in) :: interior_i_max2
          integer(ikind)                , intent(in) :: interior_i_max3
          integer(ikind)                , intent(in) :: i_min1
          integer(ikind)                , intent(in) :: i_min3
          integer(ikind)                , intent(in) :: i_min4
          integer(ikind)                , intent(in) :: i_min5
          integer(ikind)                , intent(in) :: i_min6
          integer(ikind)                , intent(in) :: i_min8
          integer(ikind)                , intent(in) :: interior_i_max11
          integer(ikind)                , intent(in) :: interior_i_max13
          integer(ikind)                , intent(in) :: interior_i_max21
          integer(ikind)                , intent(in) :: interior_i_max23
          integer(ikind)                , intent(in) :: interior_i_max31
          integer(ikind)                , intent(in) :: interior_i_max33
          integer(ikind)                , intent(in) :: i_min11, i_min13
          integer(ikind)                , intent(in) :: i_min21, i_min23
          integer(ikind)                , intent(in) :: i_min31, i_min33
          integer(ikind)                , intent(in) :: j_min1
          integer(ikind)                , intent(in) :: j_min2


          integer, dimension(bc_size, 2*bc_size) :: border_W1, border_E1
          integer, dimension(bc_size, 2*bc_size) :: border_W2, border_E2
          integer, dimension(bc_size, 2*bc_size) :: border_W3, border_E3
          integer, dimension(2*bc_size)          :: interior_profile
          integer(ikind)                         :: j_match1, j_match2

          j_match1 = 0
          j_match2 = 0
          
          !get the additional blocks
          call get_additional_blocks_N(
     $         bf_alignment(1,1)-bc_size+i_min1,
     $         bf_alignment(1,1)-bc_size+i_min13,
     $         border_W1, border_E1, interior_profile)

          call get_additional_blocks_N(
     $         bf_alignment(1,1)-bc_size+i_min4,
     $         bf_alignment(1,1)-bc_size+i_min23,
     $         border_W2, border_E2, interior_profile)

          call get_additional_blocks_N(
     $         bf_alignment(1,1)-bc_size+i_min6,
     $         bf_alignment(1,1)-bc_size+i_min33,
     $         border_W3, border_E3, interior_profile)
          

          !nodes1 - nodes2
          if(alignment1(1,1).lt.alignment2(1,1)) then

             if(alignment1(2,2).gt.alignment2(2,2)) then
                
                call add_grdpts_id_blocks_1_to_13_NS(
     $               new_grdpts_id,
     $               grdpts_id1, grdpts_id2,
     $               border_W1, border_E1,
     $               border_W2, border_E2,
     $               border_W3, border_E3,
     $               interior_profile,
     $               1, 2*bc_size,
     $               j_match1, j_match2,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $               interior_i_max11, interior_i_max13,
     $               interior_i_max21, interior_i_max23,
     $               interior_i_max31, interior_i_max33,
     $               i_min11, i_min13, i_min21, i_min23, i_min31, i_min33)
                
                call add_grdpts_id_blocks_14_to_18_NS(
     $               new_grdpts_id,
     $               grdpts_id1, grdpts_id2,
     $               2*bc_size+1, j_min1,
     $               j_match1, j_match2,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4, i_min5, i_min6)

                call add_grdpts_id_blocks_19_to_23_with_20_taller_NS(
     $               new_grdpts_id, 
     $               grdpts_id1, grdpts_id2,
     $               j_min1+1, j_min2,
     $               j_match1,
     $               outside_i_max1,  outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4)
                
                call add_grdpts_id_block_24_NS(
     $               new_grdpts_id,
     $               j_min2+1, size(new_grdpts_id,2))
                
             else

                call add_grdpts_id_blocks_1_to_13_NS(
     $               new_grdpts_id,
     $               grdpts_id1, grdpts_id2,
     $               border_W1, border_E1,
     $               border_W2, border_E2,
     $               border_W3, border_E3,
     $               interior_profile,
     $               1, 2*bc_size,
     $               j_match1, j_match2,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $               interior_i_max11, interior_i_max13,
     $               interior_i_max21, interior_i_max23,
     $               interior_i_max31, interior_i_max33,
     $               i_min11, i_min13, i_min21, i_min23, i_min31, i_min33)

                call add_grdpts_id_blocks_14_to_18_NS(
     $               new_grdpts_id,
     $               grdpts_id1, grdpts_id2,
     $               2*bc_size+1, j_min1,
     $               j_match1, j_match2,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4, i_min5, i_min6)

                call add_grdpts_id_blocks_19_to_23_with_22_taller_NS(
     $               new_grdpts_id, 
     $               grdpts_id1, grdpts_id2,
     $               j_min1+1, j_min2,
     $               j_match2,
     $               outside_i_max1,  outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min5, i_min6)
                
                call add_grdpts_id_block_24_NS(
     $               new_grdpts_id,
     $               j_min2+1, size(new_grdpts_id,2))

             end if
          else
             
             if(alignment1(2,2).gt.alignment2(2,2)) then

                call add_grdpts_id_blocks_1_to_13_NS(
     $               new_grdpts_id,
     $               grdpts_id2, grdpts_id1,
     $               border_W1, border_E1,
     $               border_W2, border_E2,
     $               border_W3, border_E3,
     $               interior_profile,
     $               1, 2*bc_size,
     $               j_match2, j_match1,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $               interior_i_max11, interior_i_max13,
     $               interior_i_max21, interior_i_max23,
     $               interior_i_max31, interior_i_max33,
     $               i_min11, i_min13, i_min21, i_min23, i_min31, i_min33)

                call add_grdpts_id_blocks_14_to_18_NS(
     $               new_grdpts_id,
     $               grdpts_id2, grdpts_id1,
     $               2*bc_size+1, j_min1,
     $               j_match2, j_match1,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4, i_min5, i_min6)

                call add_grdpts_id_blocks_19_to_23_with_20_taller_NS(
     $               new_grdpts_id, 
     $               grdpts_id2, grdpts_id1,
     $               j_min1+1, j_min2,
     $               j_match2,
     $               outside_i_max1,  outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4)
                
                call add_grdpts_id_block_24_NS(
     $               new_grdpts_id,
     $               j_min2+1, size(new_grdpts_id,2))
                
             else

                call add_grdpts_id_blocks_1_to_13_NS(
     $               new_grdpts_id,
     $               grdpts_id2, grdpts_id1,
     $               border_W1, border_E1,
     $               border_W2, border_E2,
     $               border_W3, border_E3,
     $               interior_profile,
     $               1, 2*bc_size,
     $               j_match2, j_match1,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $               interior_i_max11, interior_i_max13,
     $               interior_i_max21, interior_i_max23,
     $               interior_i_max31, interior_i_max33,
     $               i_min11, i_min13, i_min21, i_min23, i_min31, i_min33)

                call add_grdpts_id_blocks_14_to_18_NS(
     $               new_grdpts_id,
     $               grdpts_id2, grdpts_id1,
     $               2*bc_size+1, j_min1,
     $               j_match2, j_match1,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4, i_min5, i_min6)

                call add_grdpts_id_blocks_19_to_23_with_22_taller_NS(
     $               new_grdpts_id, 
     $               grdpts_id2, grdpts_id1,
     $               j_min1+1, j_min2,
     $               j_match1,
     $               outside_i_max1,  outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min5, i_min6)
                
                call add_grdpts_id_block_24_NS(
     $               new_grdpts_id,
     $               j_min2+1, size(new_grdpts_id,2))
             end if
          end if

        end subroutine merge_grdpts_id_N


        !> merge the nodes for southern buffer layers
        subroutine merge_grdpts_id_S(
     $     new_grdpts_id,
     $     grdpts_id1, grdpts_id2,
     $     alignment1, alignment2, bf_alignment,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $     interior_i_max11, interior_i_max13,
     $     interior_i_max21, interior_i_max23,
     $     interior_i_max31, interior_i_max33,
     $     i_min11, i_min13, i_min21, i_min23, i_min31, i_min33,
     $     j_min1, j_min2)

          implicit none

          integer       , dimension(:,:), intent(out):: new_grdpts_id
          integer       , dimension(:,:), intent(in) :: grdpts_id1
          integer       , dimension(:,:), intent(in) :: grdpts_id2
          integer(ikind), dimension(2,2), intent(in) :: alignment1
          integer(ikind), dimension(2,2), intent(in) :: alignment2          
          integer(ikind), dimension(2,2), intent(in) :: bf_alignment
          integer(ikind)                , intent(in) :: outside_i_max1
          integer(ikind)                , intent(in) :: outside_i_max2
          integer(ikind)                , intent(in) :: interior_i_max1
          integer(ikind)                , intent(in) :: interior_i_max2
          integer(ikind)                , intent(in) :: interior_i_max3
          integer(ikind)                , intent(in) :: i_min1
          integer(ikind)                , intent(in) :: i_min3
          integer(ikind)                , intent(in) :: i_min4
          integer(ikind)                , intent(in) :: i_min5
          integer(ikind)                , intent(in) :: i_min6
          integer(ikind)                , intent(in) :: i_min8
          integer(ikind)                , intent(in) :: interior_i_max11
          integer(ikind)                , intent(in) :: interior_i_max13
          integer(ikind)                , intent(in) :: interior_i_max21
          integer(ikind)                , intent(in) :: interior_i_max23
          integer(ikind)                , intent(in) :: interior_i_max31
          integer(ikind)                , intent(in) :: interior_i_max33
          integer(ikind)                , intent(in) :: i_min11, i_min13
          integer(ikind)                , intent(in) :: i_min21, i_min23
          integer(ikind)                , intent(in) :: i_min31, i_min33
          integer(ikind)                , intent(in) :: j_min1
          integer(ikind)                , intent(in) :: j_min2


          integer, dimension(bc_size, 2*bc_size) :: border_W1, border_E1
          integer, dimension(bc_size, 2*bc_size) :: border_W2, border_E2
          integer, dimension(bc_size, 2*bc_size) :: border_W3, border_E3
          integer, dimension(2*bc_size)          :: interior_profile
          integer(ikind)                         :: j_match1, j_match2

          j_match1 = - size(new_grdpts_id,2)+size(grdpts_id1,2)
          j_match2 = - size(new_grdpts_id,2)+size(grdpts_id2,2)

          
          !get the additional blocks
          call get_additional_blocks_S(
     $         bf_alignment(1,1)-bc_size+i_min1,
     $         bf_alignment(1,1)-bc_size+i_min13,
     $         border_W1, border_E1, interior_profile)

          call get_additional_blocks_S(
     $         bf_alignment(1,1)-bc_size+i_min4,
     $         bf_alignment(1,1)-bc_size+i_min23,
     $         border_W2, border_E2, interior_profile)

          call get_additional_blocks_S(
     $         bf_alignment(1,1)-bc_size+i_min6,
     $         bf_alignment(1,1)-bc_size+i_min33,
     $         border_W3, border_E3, interior_profile)
          

          !nodes1 - nodes2
          if(alignment1(1,1).lt.alignment2(1,1)) then

             if(alignment1(2,1).lt.alignment2(2,1)) then
                
                call add_grdpts_id_block_24_NS(
     $               new_grdpts_id,
     $               1, size(new_grdpts_id,2)-j_min2)

                call add_grdpts_id_blocks_19_to_23_with_20_taller_NS(
     $               new_grdpts_id, 
     $               grdpts_id1, grdpts_id2,
     $               size(new_grdpts_id,2)-j_min2+1,
     $               size(new_grdpts_id,2)-j_min1,
     $               j_match1,
     $               outside_i_max1,  outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4)

                call add_grdpts_id_blocks_14_to_18_NS(
     $               new_grdpts_id,
     $               grdpts_id1, grdpts_id2,
     $               size(new_grdpts_id,2)-j_min1+1,
     $               size(new_grdpts_id,2)-(2*bc_size),
     $               j_match1, j_match2,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4, i_min5, i_min6)

                call add_grdpts_id_blocks_1_to_13_NS(
     $               new_grdpts_id,
     $               grdpts_id1, grdpts_id2,
     $               border_W1, border_E1,
     $               border_W2, border_E2,
     $               border_W3, border_E3,
     $               interior_profile,
     $               size(new_grdpts_id,2)-(2*bc_size)+1,
     $               size(new_grdpts_id,2),
     $               j_match1, j_match2,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $               interior_i_max11, interior_i_max13,
     $               interior_i_max21, interior_i_max23,
     $               interior_i_max31, interior_i_max33,
     $               i_min11, i_min13, i_min21, i_min23, i_min31, i_min33)
                
             else

                call add_grdpts_id_block_24_NS(
     $               new_grdpts_id,
     $               1, size(new_grdpts_id,2)-j_min2)

                call add_grdpts_id_blocks_19_to_23_with_22_taller_NS(
     $               new_grdpts_id, 
     $               grdpts_id1, grdpts_id2,
     $               size(new_grdpts_id,2)-j_min2+1, size(new_grdpts_id,2)-j_min1,
     $               j_match2,
     $               outside_i_max1,  outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min5, i_min6)

                call add_grdpts_id_blocks_14_to_18_NS(
     $               new_grdpts_id,
     $               grdpts_id1, grdpts_id2,
     $               size(new_grdpts_id,2)-j_min1+1,
     $               size(new_grdpts_id,2)-(2*bc_size),
     $               j_match1, j_match2,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4, i_min5, i_min6)

                call add_grdpts_id_blocks_1_to_13_NS(
     $               new_grdpts_id,
     $               grdpts_id1, grdpts_id2,
     $               border_W1, border_E1,
     $               border_W2, border_E2,
     $               border_W3, border_E3,
     $               interior_profile,
     $               size(new_grdpts_id,2)-(2*bc_size)+1,
     $               size(new_grdpts_id,2),
     $               j_match1, j_match2,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $               interior_i_max11, interior_i_max13,
     $               interior_i_max21, interior_i_max23,
     $               interior_i_max31, interior_i_max33,
     $               i_min11, i_min13, i_min21, i_min23, i_min31, i_min33)

             end if
          else
             
             if(alignment1(2,2).gt.alignment2(2,2)) then

                call add_grdpts_id_block_24_NS(
     $               new_grdpts_id,
     $               1, size(new_grdpts_id,2)-j_min2)

                call add_grdpts_id_blocks_19_to_23_with_20_taller_NS(
     $               new_grdpts_id, 
     $               grdpts_id2, grdpts_id1,
     $               size(new_grdpts_id,2)-j_min2+1,
     $               size(new_grdpts_id,2)-j_min1,
     $               j_match2,
     $               outside_i_max1,  outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4)

                call add_grdpts_id_blocks_14_to_18_NS(
     $               new_grdpts_id,
     $               grdpts_id2, grdpts_id1,
     $               size(new_grdpts_id,2)-j_min1+1,
     $               size(new_grdpts_id,2)-(2*bc_size),
     $               j_match2, j_match1,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4, i_min5, i_min6)

                call add_grdpts_id_blocks_1_to_13_NS(
     $               new_grdpts_id,
     $               grdpts_id2, grdpts_id1,
     $               border_W1, border_E1,
     $               border_W2, border_E2,
     $               border_W3, border_E3,
     $               interior_profile,
     $               size(new_grdpts_id,2)-(2*bc_size)+1,
     $               size(new_grdpts_id,2),
     $               j_match2, j_match1,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $               interior_i_max11, interior_i_max13,
     $               interior_i_max21, interior_i_max23,
     $               interior_i_max31, interior_i_max33,
     $               i_min11, i_min13, i_min21, i_min23, i_min31, i_min33)
                
             else

                call add_grdpts_id_block_24_NS(
     $               new_grdpts_id,
     $               1, size(new_grdpts_id,2)-j_min2)

                call add_grdpts_id_blocks_19_to_23_with_22_taller_NS(
     $               new_grdpts_id, 
     $               grdpts_id2, grdpts_id1,
     $               size(new_grdpts_id,2)-j_min2+1,
     $               size(new_grdpts_id,2)-j_min1,
     $               j_match1,
     $               outside_i_max1,  outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min5, i_min6)

                call add_grdpts_id_blocks_14_to_18_NS(
     $               new_grdpts_id,
     $               grdpts_id2, grdpts_id1,
     $               size(new_grdpts_id,2)-j_min1+1,
     $               size(new_grdpts_id,2)-(2*bc_size),
     $               j_match2, j_match1,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4, i_min5, i_min6)

                call add_grdpts_id_blocks_1_to_13_NS(
     $               new_grdpts_id,
     $               grdpts_id2, grdpts_id1,
     $               border_W1, border_E1,
     $               border_W2, border_E2,
     $               border_W3, border_E3,
     $               interior_profile,
     $               size(new_grdpts_id,2)-(2*bc_size)+1,
     $               size(new_grdpts_id,2),
     $               j_match2, j_match1,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $               interior_i_max11, interior_i_max13,
     $               interior_i_max21, interior_i_max23,
     $               interior_i_max31, interior_i_max33,
     $               i_min11, i_min13, i_min21, i_min23, i_min31, i_min33)
                
             end if
          end if

        end subroutine merge_grdpts_id_S


        !> merge the nodes for eastern buffer layers
        subroutine merge_grdpts_id_E(
     $     new_grdpts_id,
     $     grdpts_id1, grdpts_id2,
     $     alignment1, alignment2, bf_alignment,
     $     outside_j_max1, outside_j_max2,
     $     interior_j_max1, interior_j_max2, interior_j_max3,
     $     j_min1, j_min3, j_min4, j_min5, j_min6, j_min8,
     $     interior_j_max11, interior_j_max13,
     $     interior_j_max21, interior_j_max23,
     $     interior_j_max31, interior_j_max33,
     $     j_min11, j_min13, j_min21, j_min23, j_min31, j_min33)

          implicit none

          integer       , dimension(:,:), intent(out):: new_grdpts_id
          integer       , dimension(:,:), intent(in) :: grdpts_id1
          integer       , dimension(:,:), intent(in) :: grdpts_id2
          integer(ikind), dimension(2,2), intent(in) :: alignment1
          integer(ikind), dimension(2,2), intent(in) :: alignment2          
          integer(ikind), dimension(2,2), intent(in) :: bf_alignment
          integer(ikind)                , intent(in) :: outside_j_max1
          integer(ikind)                , intent(in) :: outside_j_max2
          integer(ikind)                , intent(in) :: interior_j_max1
          integer(ikind)                , intent(in) :: interior_j_max2
          integer(ikind)                , intent(in) :: interior_j_max3
          integer(ikind)                , intent(in) :: j_min1
          integer(ikind)                , intent(in) :: j_min3
          integer(ikind)                , intent(in) :: j_min4
          integer(ikind)                , intent(in) :: j_min5
          integer(ikind)                , intent(in) :: j_min6
          integer(ikind)                , intent(in) :: j_min8
          integer(ikind)                , intent(in) :: interior_j_max11
          integer(ikind)                , intent(in) :: interior_j_max13
          integer(ikind)                , intent(in) :: interior_j_max21
          integer(ikind)                , intent(in) :: interior_j_max23
          integer(ikind)                , intent(in) :: interior_j_max31
          integer(ikind)                , intent(in) :: interior_j_max33
          integer(ikind)                , intent(in) :: j_min11, j_min13
          integer(ikind)                , intent(in) :: j_min21, j_min23
          integer(ikind)                , intent(in) :: j_min31, j_min33


          integer, dimension(2*bc_size, bc_size) :: border_S1, border_N1
          integer, dimension(2*bc_size, bc_size) :: border_S2, border_N2
          integer, dimension(2*bc_size, bc_size) :: border_S3, border_N3
          integer, dimension(2*bc_size)          :: interior_profile
          

          !get the additional blocks
          call get_additional_blocks_E(
     $         bf_alignment(2,1)-bc_size+j_min1,
     $         bf_alignment(2,1)-bc_size+j_min13,
     $         border_S1, border_N1, interior_profile)

          call get_additional_blocks_E(
     $         bf_alignment(2,1)-bc_size+j_min4,
     $         bf_alignment(2,1)-bc_size+j_min23,
     $         border_S2, border_N2, interior_profile)

          call get_additional_blocks_E(
     $         bf_alignment(2,1)-bc_size+j_min6,
     $         bf_alignment(2,1)-bc_size+j_min33,
     $         border_S3, border_N3, interior_profile)
          

          !nodes1 - nodes2
          if(alignment1(2,1).lt.alignment2(2,1)) then

             call add_grdpts_id_outside_blocks_EW(
     $            new_grdpts_id,
     $            1, outside_j_max1)

             call add_grdpts_id_edge_block_E(
     $            new_grdpts_id,
     $            border_S1,
     $            j_min1+1,
     $            j_min1+interior_j_max11)

             call add_grdpts_id_interior_block_E(
     $            new_grdpts_id,
     $            interior_profile,
     $            j_min11+1,
     $            j_min11+interior_j_max1-
     $            (interior_j_max11+interior_j_max13))

             call add_grdpts_id_edge_block_E(
     $            new_grdpts_id,
     $            border_N1,
     $            j_min13+1,
     $            j_min13+interior_j_max13)

             call add_grdpts_id_sublayer_block_E(
     $            new_grdpts_id, grdpts_id1,
     $            j_min3+1, j_min3+size(grdpts_id1,2))

             call add_grdpts_id_edge_block_E(
     $            new_grdpts_id,
     $            border_S2,
     $            j_min4+1,
     $            j_min4+interior_j_max21)

             call add_grdpts_id_interior_block_E(
     $            new_grdpts_id,
     $            interior_profile,
     $            j_min21+1,
     $            j_min21+interior_j_max2-
     $            (interior_j_max21+interior_j_max23))

             call add_grdpts_id_edge_block_E(
     $            new_grdpts_id,
     $            border_N2,
     $            j_min23+1,
     $            j_min23+interior_j_max23)

             call add_grdpts_id_sublayer_block_E(
     $            new_grdpts_id, grdpts_id2,
     $            j_min5+1, j_min5+size(grdpts_id2,2))

             call add_grdpts_id_edge_block_E(
     $            new_grdpts_id,
     $            border_S3,
     $            j_min6+1,
     $            j_min6+interior_j_max31)

             call add_grdpts_id_interior_block_E(
     $            new_grdpts_id,
     $            interior_profile,
     $            j_min31+1,
     $            j_min31+interior_j_max3-
     $            (interior_j_max31+interior_j_max33))

             call add_grdpts_id_edge_block_E(
     $            new_grdpts_id,
     $            border_N3,
     $            j_min33+1,
     $            j_min33+interior_j_max33)

             call add_grdpts_id_outside_blocks_EW(
     $            new_grdpts_id,
     $            j_min8+1, j_min8+outside_j_max2)

          else

             call add_grdpts_id_outside_blocks_EW(
     $            new_grdpts_id,
     $            1, outside_j_max1)

             call add_grdpts_id_edge_block_E(
     $            new_grdpts_id,
     $            border_S1,
     $            j_min1+1,
     $            j_min1+interior_j_max11)

             call add_grdpts_id_interior_block_E(
     $            new_grdpts_id,
     $            interior_profile,
     $            j_min11+1,
     $            j_min11+interior_j_max1-
     $            (interior_j_max11+interior_j_max13))

             call add_grdpts_id_edge_block_E(
     $            new_grdpts_id,
     $            border_N1,
     $            j_min13+1,
     $            j_min13+interior_j_max13)

             call add_grdpts_id_sublayer_block_E(
     $            new_grdpts_id, grdpts_id2,
     $            j_min3+1, j_min3+size(grdpts_id2,2))

             call add_grdpts_id_edge_block_E(
     $            new_grdpts_id,
     $            border_S2,
     $            j_min4+1,
     $            j_min4+interior_j_max21)

             call add_grdpts_id_interior_block_E(
     $            new_grdpts_id,
     $            interior_profile,
     $            j_min21+1,
     $            j_min21+interior_j_max2-
     $            (interior_j_max21+interior_j_max23))

             call add_grdpts_id_edge_block_E(
     $            new_grdpts_id,
     $            border_N2,
     $            j_min23+1,
     $            j_min23+interior_j_max23)

             call add_grdpts_id_sublayer_block_E(
     $            new_grdpts_id, grdpts_id1,
     $            j_min5+1, j_min5+size(grdpts_id1,2))

             call add_grdpts_id_edge_block_E(
     $            new_grdpts_id,
     $            border_S3,
     $            j_min6+1,
     $            j_min6+interior_j_max31)

             call add_grdpts_id_interior_block_E(
     $            new_grdpts_id,
     $            interior_profile,
     $            j_min31+1,
     $            j_min31+interior_j_max3-
     $            (interior_j_max31+interior_j_max33))

             call add_grdpts_id_edge_block_E(
     $            new_grdpts_id,
     $            border_N3,
     $            j_min33+1,
     $            j_min33+interior_j_max33)

             call add_grdpts_id_outside_blocks_EW(
     $            new_grdpts_id,
     $            j_min8+1, size(new_grdpts_id,2))

          end if

        end subroutine merge_grdpts_id_E


        !> merge the nodes for eastern buffer layers
        subroutine merge_grdpts_id_W(
     $     new_grdpts_id,
     $     grdpts_id1, grdpts_id2,
     $     alignment1, alignment2, bf_alignment,
     $     outside_j_max1, outside_j_max2,
     $     interior_j_max1, interior_j_max2, interior_j_max3,
     $     j_min1, j_min3, j_min4, j_min5, j_min6, j_min8,
     $     interior_j_max11, interior_j_max13,
     $     interior_j_max21, interior_j_max23,
     $     interior_j_max31, interior_j_max33,
     $     j_min11, j_min13, j_min21, j_min23, j_min31, j_min33)

          implicit none

          integer       , dimension(:,:), intent(out):: new_grdpts_id
          integer       , dimension(:,:), intent(in) :: grdpts_id1
          integer       , dimension(:,:), intent(in) :: grdpts_id2
          integer(ikind), dimension(2,2), intent(in) :: alignment1
          integer(ikind), dimension(2,2), intent(in) :: alignment2          
          integer(ikind), dimension(2,2), intent(in) :: bf_alignment
          integer(ikind)                , intent(in) :: outside_j_max1
          integer(ikind)                , intent(in) :: outside_j_max2
          integer(ikind)                , intent(in) :: interior_j_max1
          integer(ikind)                , intent(in) :: interior_j_max2
          integer(ikind)                , intent(in) :: interior_j_max3
          integer(ikind)                , intent(in) :: j_min1
          integer(ikind)                , intent(in) :: j_min3
          integer(ikind)                , intent(in) :: j_min4
          integer(ikind)                , intent(in) :: j_min5
          integer(ikind)                , intent(in) :: j_min6
          integer(ikind)                , intent(in) :: j_min8
          integer(ikind)                , intent(in) :: interior_j_max11
          integer(ikind)                , intent(in) :: interior_j_max13
          integer(ikind)                , intent(in) :: interior_j_max21
          integer(ikind)                , intent(in) :: interior_j_max23
          integer(ikind)                , intent(in) :: interior_j_max31
          integer(ikind)                , intent(in) :: interior_j_max33
          integer(ikind)                , intent(in) :: j_min11, j_min13
          integer(ikind)                , intent(in) :: j_min21, j_min23
          integer(ikind)                , intent(in) :: j_min31, j_min33


          integer, dimension(2*bc_size, bc_size) :: border_S1, border_N1
          integer, dimension(2*bc_size, bc_size) :: border_S2, border_N2
          integer, dimension(2*bc_size, bc_size) :: border_S3, border_N3
          integer, dimension(2*bc_size)          :: interior_profile
          

          !get the additional blocks
          call get_additional_blocks_W(
     $         bf_alignment(2,1)-bc_size+j_min1,
     $         bf_alignment(2,1)-bc_size+j_min13,
     $         border_S1, border_N1, interior_profile)

          call get_additional_blocks_W(
     $         bf_alignment(2,1)-bc_size+j_min4,
     $         bf_alignment(2,1)-bc_size+j_min23,
     $         border_S2, border_N2, interior_profile)

          call get_additional_blocks_W(
     $         bf_alignment(2,1)-bc_size+j_min6,
     $         bf_alignment(2,1)-bc_size+j_min33,
     $         border_S3, border_N3, interior_profile)
          

          !nodes1 - nodes2
          if(alignment1(2,1).lt.alignment2(2,1)) then

             call add_grdpts_id_outside_blocks_EW(
     $            new_grdpts_id,
     $            1, outside_j_max1)

             call add_grdpts_id_edge_block_W(
     $            new_grdpts_id,
     $            border_S1,
     $            j_min1+1,
     $            j_min1+interior_j_max11)

             call add_grdpts_id_interior_block_W(
     $            new_grdpts_id,
     $            interior_profile,
     $            j_min11+1,
     $            j_min11+interior_j_max1-
     $            (interior_j_max11+interior_j_max13))

             call add_grdpts_id_edge_block_W(
     $            new_grdpts_id,
     $            border_N1,
     $            j_min13+1,
     $            j_min13+interior_j_max13)

             call add_grdpts_id_sublayer_block_W(
     $            new_grdpts_id, grdpts_id1,
     $            j_min3+1, j_min3+size(grdpts_id1,2))

             call add_grdpts_id_edge_block_W(
     $            new_grdpts_id,
     $            border_S2,
     $            j_min4+1,
     $            j_min4+interior_j_max21)

             call add_grdpts_id_interior_block_W(
     $            new_grdpts_id,
     $            interior_profile,
     $            j_min21+1,
     $            j_min21+interior_j_max2-
     $            (interior_j_max21+interior_j_max23))

             call add_grdpts_id_edge_block_W(
     $            new_grdpts_id,
     $            border_N2,
     $            j_min23+1,
     $            j_min23+interior_j_max23)

             call add_grdpts_id_sublayer_block_W(
     $            new_grdpts_id, grdpts_id2,
     $            j_min5+1, j_min5+size(grdpts_id2,2))

             call add_grdpts_id_edge_block_W(
     $            new_grdpts_id,
     $            border_S3,
     $            j_min6+1,
     $            j_min6+interior_j_max31)

             call add_grdpts_id_interior_block_W(
     $            new_grdpts_id,
     $            interior_profile,
     $            j_min31+1,
     $            j_min31+interior_j_max3-
     $            (interior_j_max31+interior_j_max33))

             call add_grdpts_id_edge_block_W(
     $            new_grdpts_id,
     $            border_N3,
     $            j_min33+1,
     $            j_min33+interior_j_max33)

             call add_grdpts_id_outside_blocks_EW(
     $            new_grdpts_id,
     $            j_min8+1, j_min8+outside_j_max2)

          else

             call add_grdpts_id_outside_blocks_EW(
     $            new_grdpts_id,
     $            1, outside_j_max1)

             call add_grdpts_id_edge_block_W(
     $            new_grdpts_id,
     $            border_S1,
     $            j_min1+1,
     $            j_min1+interior_j_max11)

             call add_grdpts_id_interior_block_W(
     $            new_grdpts_id,
     $            interior_profile,
     $            j_min11+1,
     $            j_min11+interior_j_max1-
     $            (interior_j_max11+interior_j_max13))

             call add_grdpts_id_edge_block_W(
     $            new_grdpts_id,
     $            border_N1,
     $            j_min13+1,
     $            j_min13+interior_j_max13)

             call add_grdpts_id_sublayer_block_W(
     $            new_grdpts_id, grdpts_id2,
     $            j_min3+1, j_min3+size(grdpts_id2,2))

             call add_grdpts_id_edge_block_W(
     $            new_grdpts_id,
     $            border_S2,
     $            j_min4+1,
     $            j_min4+interior_j_max21)

             call add_grdpts_id_interior_block_W(
     $            new_grdpts_id,
     $            interior_profile,
     $            j_min21+1,
     $            j_min21+interior_j_max2-
     $            (interior_j_max21+interior_j_max23))

             call add_grdpts_id_edge_block_W(
     $            new_grdpts_id,
     $            border_N2,
     $            j_min23+1,
     $            j_min23+interior_j_max23)

             call add_grdpts_id_sublayer_block_W(
     $            new_grdpts_id, grdpts_id1,
     $            j_min5+1, j_min5+size(grdpts_id1,2))

             call add_grdpts_id_edge_block_W(
     $            new_grdpts_id,
     $            border_S3,
     $            j_min6+1,
     $            j_min6+interior_j_max31)

             call add_grdpts_id_interior_block_W(
     $            new_grdpts_id,
     $            interior_profile,
     $            j_min31+1,
     $            j_min31+interior_j_max3-
     $            (interior_j_max31+interior_j_max33))

             call add_grdpts_id_edge_block_W(
     $            new_grdpts_id,
     $            border_N3,
     $            j_min33+1,
     $            j_min33+interior_j_max33)

             call add_grdpts_id_outside_blocks_EW(
     $            new_grdpts_id,
     $            j_min8+1, size(new_grdpts_id,2))

          end if

        end subroutine merge_grdpts_id_W


        subroutine add_nodes_blocks_2_to_12_NS(
     $     new_nodes,
     $     nodes_block4, nodes_block6, interior_nodes,
     $     bf_alignment,
     $     k, j_min, j_max,
     $     j_match_block4, j_match_block6, j_match_interior,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min1, i_min3, i_min4, i_min5, i_min6)

          implicit none

          real(rkind), dimension(:,:,:), intent(out):: new_nodes
          real(rkind), dimension(:,:,:), intent(in) :: nodes_block4
          real(rkind), dimension(:,:,:), intent(in) :: nodes_block6
          real(rkind), dimension(nx,ny,ne), intent(in) :: interior_nodes
          integer(ikind), dimension(2,2), intent(in) :: bf_alignment
          integer       , intent(in) :: k
          integer(ikind), intent(in) :: j_min
          integer(ikind), intent(in) :: j_max
          integer(ikind), intent(in) :: j_match_block4
          integer(ikind), intent(in) :: j_match_block6
          integer(ikind), intent(in) :: j_match_interior
          integer(ikind), intent(in) :: interior_i_max1
          integer(ikind), intent(in) :: interior_i_max2
          integer(ikind), intent(in) :: interior_i_max3
          integer(ikind), intent(in) :: i_min1
          integer(ikind), intent(in) :: i_min3
          integer(ikind), intent(in) :: i_min4
          integer(ikind), intent(in) :: i_min5
          integer(ikind), intent(in) :: i_min6
          
          
          integer(ikind) :: i,j

          do j=j_min, j_max
             do i=1, interior_i_max1
                new_nodes(i_min1+i,j,k) = interior_nodes(
     $               bf_alignment(1,1)-(bc_size+1)+i_min1+i,
     $               j_match_interior+j,
     $               k)
             end do

             do i=1, size(nodes_block4,1)
                new_nodes(i_min3+i,j,k) = nodes_block4(
     $               i,
     $               j_match_block4+j,
     $               k)
             end do
             
             do i=1, interior_i_max2
                new_nodes(i_min4+i,j,k) = interior_nodes(
     $               bf_alignment(1,1)-(bc_size+1)+i_min4+i,
     $               j_match_interior+j,
     $               k)
             end do

             do i=1, size(nodes_block6,1)
                new_nodes(i_min5+i,j,k) = nodes_block6(
     $               i,
     $               j_match_block6+j,
     $               k)
             end do

             do i=1, interior_i_max3
                new_nodes(i_min6+i,j,k) = interior_nodes(
     $               bf_alignment(1,1)-(bc_size+1)+i_min6+i,
     $               j_match_interior+j,
     $               k)
             end do
          end do

        end subroutine add_nodes_blocks_2_to_12_NS


        subroutine add_nodes_blocks_15_to_17_NS(
     $     new_nodes,
     $     nodes_block11, nodes_block13,
     $     k, j_min, j_max,
     $     j_match_block11, j_match_block13,
     $     i_min3, i_min5)

          implicit none

          real(rkind), dimension(:,:,:), intent(out):: new_nodes
          real(rkind), dimension(:,:,:), intent(in) :: nodes_block11
          real(rkind), dimension(:,:,:), intent(in) :: nodes_block13
          integer       , intent(in) :: k
          integer(ikind), intent(in) :: j_min
          integer(ikind), intent(in) :: j_max
          integer(ikind), intent(in) :: j_match_block11
          integer(ikind), intent(in) :: j_match_block13
          integer(ikind), intent(in) :: i_min3
          integer(ikind), intent(in) :: i_min5


          integer(ikind) :: i,j


          do j=j_min, j_max
             do i=1, size(nodes_block11,1)
                new_nodes(i_min3+i,j,k) = nodes_block11(
     $               i,
     $               j_match_block11+j,
     $               k)
             end do
             
             do i=1, size(nodes_block13,1)
                new_nodes(i_min5+i,j,k) = nodes_block13(
     $               i,
     $               j_match_block13+j,
     $               k)
             end do
          end do

        end subroutine add_nodes_blocks_15_to_17_NS


        subroutine add_nodes_blocks_20_to_22_NS(
     $     new_nodes,
     $     nodes_block16,
     $     k, j_min, j_max,
     $     j_match_block16,
     $     i_min3)

          implicit none

          real(rkind), dimension(:,:,:), intent(out):: new_nodes
          real(rkind), dimension(:,:,:), intent(in) :: nodes_block16
          integer       , intent(in) :: k
          integer(ikind), intent(in) :: j_min
          integer(ikind), intent(in) :: j_max
          integer(ikind), intent(in) :: j_match_block16
          integer(ikind), intent(in) :: i_min3


          integer(ikind) :: i,j


          do j=j_min, j_max
             do i=1, size(nodes_block16,1)
                new_nodes(i_min3+i,j,k) = nodes_block16(
     $               i,
     $               j_match_block16+j,
     $               k)
             end do
          end do

        end subroutine add_nodes_blocks_20_to_22_NS


        subroutine add_grdpts_id_blocks_1_to_13_NS(
     $     new_grdpts_id,
     $     grdpts_block4, grdpts_block6,
     $     border_W1, border_E1,
     $     border_W2, border_E2,
     $     border_W3, border_E3,
     $     interior_profile,
     $     j_min, j_max,
     $     j_match_block4, j_match_block6,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min1, i_min3, i_min4, i_min5, i_min6, i_min8,
     $     interior_i_max11, interior_i_max13,
     $     interior_i_max21, interior_i_max23,
     $     interior_i_max31, interior_i_max33,
     $     i_min11, i_min13, i_min21, i_min23, i_min31, i_min33)

          implicit none

          integer, dimension(:,:)              , intent(out):: new_grdpts_id
          integer, dimension(:,:)              , intent(in) :: grdpts_block4
          integer, dimension(:,:)              , intent(in) :: grdpts_block6
          integer, dimension(bc_size,2*bc_size), intent(in) :: border_W1
          integer, dimension(bc_size,2*bc_size), intent(in) :: border_E1
          integer, dimension(bc_size,2*bc_size), intent(in) :: border_W2
          integer, dimension(bc_size,2*bc_size), intent(in) :: border_E2
          integer, dimension(bc_size,2*bc_size), intent(in) :: border_W3
          integer, dimension(bc_size,2*bc_size), intent(in) :: border_E3
          integer, dimension(2*bc_size)        , intent(in) :: interior_profile
          integer(ikind)                       , intent(in) :: j_min
          integer(ikind)                       , intent(in) :: j_max
          integer(ikind)                       , intent(in) :: j_match_block4
          integer(ikind)                       , intent(in) :: j_match_block6
          integer(ikind)                       , intent(in) :: outside_i_max1
          integer(ikind)                       , intent(in) :: outside_i_max2
          integer(ikind)                       , intent(in) :: interior_i_max1
          integer(ikind)                       , intent(in) :: interior_i_max2
          integer(ikind)                       , intent(in) :: interior_i_max3
          integer(ikind)                       , intent(in) :: i_min1
          integer(ikind)                       , intent(in) :: i_min3
          integer(ikind)                       , intent(in) :: i_min4
          integer(ikind)                       , intent(in) :: i_min5
          integer(ikind)                       , intent(in) :: i_min6
          integer(ikind)                       , intent(in) :: i_min8
          integer(ikind)                       , intent(in) :: interior_i_max11
          integer(ikind)                       , intent(in) :: interior_i_max13
          integer(ikind)                       , intent(in) :: interior_i_max21
          integer(ikind)                       , intent(in) :: interior_i_max23
          integer(ikind)                       , intent(in) :: interior_i_max31
          integer(ikind)                       , intent(in) :: interior_i_max33
          integer(ikind)                       , intent(in) :: i_min11
          integer(ikind)                       , intent(in) :: i_min13
          integer(ikind)                       , intent(in) :: i_min21
          integer(ikind)                       , intent(in) :: i_min23
          integer(ikind)                       , intent(in) :: i_min31
          integer(ikind)                       , intent(in) :: i_min33

          integer(ikind) :: i,j


          do j=j_min, j_max

             !block 1
             do i=1, outside_i_max1
                new_grdpts_id(i,j) = no_pt
             end do


             !block 2
             do i=1, interior_i_max11
                new_grdpts_id(i_min1+i,j) = border_W1(i,j-(j_min-1))
             end do

             !block 3
             do i=1, interior_i_max1-(interior_i_max11+interior_i_max13)
                new_grdpts_id(i_min11+i,j) = interior_profile(j-(j_min-1))
             end do

             !block 4
             do i=1, interior_i_max13
                new_grdpts_id(i_min13+i,j) = border_E1(i,j-(j_min-1))
             end do


             !block 5
             do i=1, size(grdpts_block4,1)
                new_grdpts_id(i_min3+i,j) = grdpts_block4(i,j_match_block4+j)
             end do


             !block 6
             do i=1, interior_i_max21
                new_grdpts_id(i_min4+i,j) = border_W2(i,j-(j_min-1))
             end do

             !block 7
             do i=1, interior_i_max2-(interior_i_max21+interior_i_max23)
                new_grdpts_id(i_min21+i,j) = interior_profile(j-(j_min-1))
             end do

             !block 8
             do i=1, interior_i_max23
                new_grdpts_id(i_min23+i,j) = border_E2(i,j-(j_min-1))
             end do


             !block 9
             do i=1, size(grdpts_block6,1)
                new_grdpts_id(i_min5+i,j) = grdpts_block6(i,j_match_block6+j)
             end do


             !block 10
             do i=1, interior_i_max31
                new_grdpts_id(i_min6+i,j) = border_W3(i,j-(j_min-1))
             end do

             !block 11
             do i=1, interior_i_max3-(interior_i_max31+interior_i_max33)
                new_grdpts_id(i_min31+i,j) = interior_profile(j-(j_min-1))
             end do

             !block 12
             do i=1, interior_i_max33
                new_grdpts_id(i_min33+i,j) = border_E3(i,j-(j_min-1))
             end do


             !block 13
             do i=1, outside_i_max2
                new_grdpts_id(i_min8+i,j) = no_pt
             end do
          end do

        end subroutine add_grdpts_id_blocks_1_to_13_NS
        

        subroutine add_grdpts_id_blocks_14_to_18_NS(
     $     new_grdpts_id,
     $     grdpts_block11, grdpts_block13,
     $     j_min, j_max,
     $     j_match_block11, j_match_block13,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min3, i_min4, i_min5, i_min6)

          implicit none

          integer, dimension(:,:), intent(out):: new_grdpts_id
          integer, dimension(:,:), intent(in) :: grdpts_block11
          integer, dimension(:,:), intent(in) :: grdpts_block13
          integer(ikind), intent(in) :: j_min
          integer(ikind), intent(in) :: j_max
          integer(ikind), intent(in) :: j_match_block11
          integer(ikind), intent(in) :: j_match_block13
          integer(ikind), intent(in) :: outside_i_max1
          integer(ikind), intent(in) :: outside_i_max2
          integer(ikind), intent(in) :: interior_i_max1
          integer(ikind), intent(in) :: interior_i_max2
          integer(ikind), intent(in) :: interior_i_max3
          integer(ikind), intent(in) :: i_min3
          integer(ikind), intent(in) :: i_min4
          integer(ikind), intent(in) :: i_min5
          integer(ikind), intent(in) :: i_min6


          integer(ikind) :: i,j
                    

          do j=j_min, j_max
            !block 10
             do i=1, outside_i_max1+interior_i_max1
                new_grdpts_id(i,j) = no_pt
             end do
             
             !block 11
             do i=1, size(grdpts_block11,1)
                new_grdpts_id(i_min3+i,j) = grdpts_block11(i,j_match_block11+j)
             end do

             !block 12
             do i=1, interior_i_max2
                new_grdpts_id(i_min4+i,j) = no_pt
             end do

             !block 13
             do i=1, size(grdpts_block13,1)
                new_grdpts_id(i_min5+i,j) = grdpts_block13(i,j_match_block13+j)
             end do

             !block 14
             do i=1, interior_i_max3+outside_i_max2
                new_grdpts_id(i_min6+i,j) = no_pt
             end do
          end do

        end subroutine add_grdpts_id_blocks_14_to_18_NS


        subroutine add_grdpts_id_blocks_19_to_23_with_20_taller_NS(
     $     new_grdpts_id,
     $     grdpts_block16, grdpts_block18,
     $     j_min, j_max,
     $     j_match_block16,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min3, i_min4)


          implicit none

          integer, dimension(:,:), intent(out) :: new_grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_block16
          integer, dimension(:,:), intent(in)  :: grdpts_block18
          integer(ikind), intent(in) :: j_min
          integer(ikind), intent(in) :: j_max
          integer(ikind), intent(in) :: j_match_block16
          integer(ikind), intent(in) :: outside_i_max1
          integer(ikind), intent(in) :: outside_i_max2
          integer(ikind), intent(in) :: interior_i_max1
          integer(ikind), intent(in) :: interior_i_max2
          integer(ikind), intent(in) :: interior_i_max3
          integer(ikind), intent(in) :: i_min3
          integer(ikind), intent(in) :: i_min4
          

          integer(ikind) :: i,j


          do j=j_min, j_max
            !block 15
             do i=1, outside_i_max1+interior_i_max1
                new_grdpts_id(i,j) = no_pt
             end do
             
             !block 16
             do i=1, size(grdpts_block16,1)
                new_grdpts_id(i_min3+i,j) = grdpts_block16(i,j_match_block16+j)
             end do
             
             !block 17+18+19
             do i=1, interior_i_max2+size(grdpts_block18,1)+interior_i_max3+outside_i_max2
                new_grdpts_id(i_min4+i,j) = no_pt
             end do
          end do

        end subroutine add_grdpts_id_blocks_19_to_23_with_20_taller_NS


        subroutine add_grdpts_id_blocks_19_to_23_with_22_taller_NS(
     $     new_grdpts_id,
     $     grdpts_block16, grdpts_block18,
     $     j_min, j_max,
     $     j_match_block18,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min5, i_min6)


          implicit none

          integer, dimension(:,:), intent(out) :: new_grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_block16
          integer, dimension(:,:), intent(in)  :: grdpts_block18
          integer(ikind), intent(in) :: j_min
          integer(ikind), intent(in) :: j_max
          integer(ikind), intent(in) :: j_match_block18
          integer(ikind), intent(in) :: outside_i_max1
          integer(ikind), intent(in) :: outside_i_max2
          integer(ikind), intent(in) :: interior_i_max1
          integer(ikind), intent(in) :: interior_i_max2
          integer(ikind), intent(in) :: interior_i_max3
          integer(ikind), intent(in) :: i_min5
          integer(ikind), intent(in) :: i_min6
          

          integer(ikind) :: i,j
                    

          do j=j_min, j_max
             !block 15+16+17
             do i=1, outside_i_max1+interior_i_max1+size(grdpts_block16,1)+interior_i_max2
                new_grdpts_id(i,j) = no_pt
             end do
             
             !block 18
             do i=1, size(grdpts_block18,1)
                new_grdpts_id(i_min5+i,j) = grdpts_block18(i,j_match_block18+j)
             end do

             !block 19
             do i=1, interior_i_max3+outside_i_max2
                new_grdpts_id(i_min6+i,j) = no_pt
             end do
          end do

        end subroutine add_grdpts_id_blocks_19_to_23_with_22_taller_NS

      
        subroutine add_grdpts_id_block_24_NS(
     $     new_grdpts_id,
     $     j_min, j_max)

          implicit none

          integer, dimension(:,:), intent(out) :: new_grdpts_id
          integer(ikind)         , intent(in)  :: j_min
          integer(ikind)         , intent(in)  :: j_max


          integer(ikind) :: i,j
          

          do j=j_min, j_max
             do i=1, size(new_grdpts_id,1)
                new_grdpts_id(i,j) = no_pt
             end do
          end do

        end subroutine add_grdpts_id_block_24_NS


        subroutine add_nodes_interior_blocks_EW(
     $     new_nodes, interior_nodes,
     $     bf_alignment,
     $     k, j_min, j_max,
     $     i_matchN, i_matchI)

          implicit none

          real(rkind)   , dimension(:,:,:)   , intent(out) :: new_nodes
          real(rkind)   , dimension(nx,ny,ne), intent(in)  :: interior_nodes
          integer(ikind), dimension(2,2)     , intent(in)  :: bf_alignment
          integer(ikind)                     , intent(in)  :: k
          integer(ikind)                     , intent(in)  :: j_min
          integer(ikind)                     , intent(in)  :: j_max
          integer(ikind)                     , intent(in)  :: i_matchN
          integer(ikind)                     , intent(in)  :: i_matchI          

          integer(ikind) :: i,j

          do j=j_min, j_max
             do i=1, 2*bc_size
                new_nodes(i_matchN+i,j,k) = interior_nodes(
     $               i_matchI+i,
     $               bf_alignment(2,1)-(bc_size+1)+j,
     $               k)
             end do
          end do

        end subroutine add_nodes_interior_blocks_EW


        subroutine add_nodes_sublayer_block_EW(
     $     new_nodes, sublayer_nodes,
     $     k, j_min, j_max,
     $     i_matchN)

          implicit none

          real(rkind)   , dimension(:,:,:), intent(out) :: new_nodes
          real(rkind)   , dimension(:,:,:), intent(in)  :: sublayer_nodes
          integer(ikind)                  , intent(in)  :: k
          integer(ikind)                  , intent(in)  :: j_min
          integer(ikind)                  , intent(in)  :: j_max
          integer(ikind)                  , intent(in)  :: i_matchN

          integer(ikind) :: i,j  


          do j=j_min, j_max
             do i=1, size(sublayer_nodes,1)
                new_nodes(i_matchN+i,j,k) = sublayer_nodes(i,j-(j_min-1),k)
             end do
          end do

       end subroutine add_nodes_sublayer_block_EW


       subroutine add_grdpts_id_outside_blocks_EW(
     $     new_grdpts_id,
     $     j_min, j_max)

         implicit none

         integer, dimension(:,:), intent(out) :: new_grdpts_id
         integer(ikind)         , intent(in)  :: j_min
         integer(ikind)         , intent(in)  :: j_max

         integer(ikind) :: i,j

         do j=j_min, j_max
            do i=1, size(new_grdpts_id,1)
               new_grdpts_id(i,j) = no_pt
            end do
         end do

       end subroutine add_grdpts_id_outside_blocks_EW


        subroutine add_grdpts_id_edge_block_E(
     $     new_grdpts_id,
     $     border, j_min, j_max)

          implicit none

          integer, dimension(:,:)              , intent(out) :: new_grdpts_id
          integer, dimension(2*bc_size,bc_size), intent(in)  :: border
          integer(ikind)                       , intent(in)  :: j_min
          integer(ikind)                       , intent(in)  :: j_max

          integer(ikind) :: i,j
          
          do j=j_min, j_max
             do i=1, 2*bc_size
                new_grdpts_id(i,j) = border(i,j-(j_min-1))
             end do
             do i=2*bc_size+1, size(new_grdpts_id,1)
                new_grdpts_id(i,j) = no_pt
             end do
          end do

       end subroutine add_grdpts_id_edge_block_E


       subroutine add_grdpts_id_edge_block_W(
     $     new_grdpts_id,
     $     border, j_min, j_max)

          implicit none

          integer, dimension(:,:)              , intent(out) :: new_grdpts_id
          integer, dimension(2*bc_size,bc_size), intent(in)  :: border
          integer(ikind)                       , intent(in)  :: j_min
          integer(ikind)                       , intent(in)  :: j_max

          integer(ikind) :: i,j
          
          do j=j_min, j_max
             do i=1, size(new_grdpts_id,1)-(2*bc_size)
                new_grdpts_id(i,j) = no_pt
             end do
             do i=size(new_grdpts_id,1)-(2*bc_size)+1, size(new_grdpts_id,1)
                new_grdpts_id(i,j) = border(
     $               i-(size(new_grdpts_id,1)-(2*bc_size)),
     $               j-(j_min-1))
             end do             
          end do
          
       end subroutine add_grdpts_id_edge_block_W


       subroutine add_grdpts_id_sublayer_block_E(
     $     new_grdpts_id,
     $     sublayer_grdpts,
     $     j_min, j_max)

         implicit none

         integer       , dimension(:,:), intent(out) :: new_grdpts_id
         integer       , dimension(:,:), intent(in)  :: sublayer_grdpts
         integer(ikind)                , intent(in)  :: j_min
         integer(ikind)                , intent(in)  :: j_max

         integer(ikind) :: i,j

         do j=j_min, j_max
            do i=1, size(sublayer_grdpts,1)
               new_grdpts_id(i,j) = sublayer_grdpts(i,j-(j_min-1))
            end do
            do i=size(sublayer_grdpts,1)+1, size(new_grdpts_id,1)
               new_grdpts_id(i,j) = no_pt
            end do
         end do

       end subroutine add_grdpts_id_sublayer_block_E

      
       subroutine add_grdpts_id_sublayer_block_W(
     $     new_grdpts_id,
     $     sublayer_grdpts,
     $     j_min, j_max)

         implicit none

         integer       , dimension(:,:), intent(out) :: new_grdpts_id
         integer       , dimension(:,:), intent(in)  :: sublayer_grdpts
         integer(ikind)                , intent(in)  :: j_min
         integer(ikind)                , intent(in)  :: j_max

         integer(ikind) :: i,j

         do j=j_min, j_max
            do i=1, size(new_grdpts_id,1)-size(sublayer_grdpts,1)
               new_grdpts_id(i,j) = no_pt
            end do
            do i=size(new_grdpts_id,1)-size(sublayer_grdpts,1)+1, size(new_grdpts_id,1)
               new_grdpts_id(i,j) = sublayer_grdpts(
     $              i-(size(new_grdpts_id,1)-size(sublayer_grdpts,1)),
     $              j-(j_min-1))
            end do
            
         end do

       end subroutine add_grdpts_id_sublayer_block_W


       subroutine add_grdpts_id_interior_block_E(
     $     new_grdpts_id,
     $     interior_profile,
     $     j_min, j_max)

         implicit none

         integer, dimension(:,:)      , intent(out) :: new_grdpts_id
         integer, dimension(2*bc_size), intent(in)  :: interior_profile
         integer(ikind)               , intent(in)  :: j_min
         integer(ikind)               , intent(in)  :: j_max

         integer(ikind) :: i,j

         do j=j_min, j_max
            do i=1, 2*bc_size
               new_grdpts_id(i,j) = interior_profile(i)
            end do
            do i=2*bc_size+1, size(new_grdpts_id,1)
               new_grdpts_id(i,j) = no_pt
            end do
         end do

       end subroutine add_grdpts_id_interior_block_E


       subroutine add_grdpts_id_interior_block_W(
     $     new_grdpts_id,
     $     interior_profile,
     $     j_min, j_max)

         implicit none

         integer, dimension(:,:)      , intent(out) :: new_grdpts_id
         integer, dimension(2*bc_size), intent(in)  :: interior_profile
         integer(ikind)               , intent(in)  :: j_min
         integer(ikind)               , intent(in)  :: j_max

         integer(ikind) :: i,j

         do j=j_min, j_max
            do i=1, size(new_grdpts_id,1)-(2*bc_size)
               new_grdpts_id(i,j) = no_pt
            end do
            do i=size(new_grdpts_id,1)-(2*bc_size)+1, size(new_grdpts_id,1)
               new_grdpts_id(i,j) = interior_profile(i-(size(new_grdpts_id,1)-(2*bc_size)))
            end do            
         end do

       end subroutine add_grdpts_id_interior_block_W


       subroutine add_grdpts_id_edge_N_blocks_E(
     $     new_grdpts_id,
     $     border_N, interior_profile,
     $     j_min6, j_min7, j_match_borderN,
     $     interior_j_max3)

         implicit none

         integer, dimension(:,:)              , intent(out) :: new_grdpts_id
         integer, dimension(2*bc_size,bc_size), intent(in)  :: border_N
         integer, dimension(2*bc_size)        , intent(in)  :: interior_profile
         integer(ikind)                       , intent(in)  :: j_min6
         integer(ikind)                       , intent(in)  :: j_min7
         integer(ikind)                       , intent(in)  :: j_match_borderN
         integer(ikind)                       , intent(in)  :: interior_j_max3

         integer(ikind) :: i,j

         do j=1, interior_j_max3-bc_size
            do i=1, 2*bc_size
               new_grdpts_id(i,j_min6+j) = interior_profile(i)
            end do
            do i=2*bc_size+1, size(new_grdpts_id,1)
               new_grdpts_id(i,j_min6+j) = no_pt
            end do            
         end do

         do j=1, min(interior_j_max3,bc_size)
            do i=1, 2*bc_size
               new_grdpts_id(i,j_min7+j) = border_N(i,j_match_borderN+j)
            end do
            do i=2*bc_size+1, size(new_grdpts_id,1)
               new_grdpts_id(i,j_min7+j) = no_pt
            end do
         end do         

       end subroutine add_grdpts_id_edge_N_blocks_E


      subroutine add_grdpts_id_edge_N_blocks_W(
     $     new_grdpts_id,
     $     border_N, interior_profile,
     $     j_min6, j_min7, j_match_borderN,
     $     interior_j_max3)

         implicit none

         integer, dimension(:,:)              , intent(out) :: new_grdpts_id
         integer, dimension(2*bc_size,bc_size), intent(in)  :: border_N
         integer, dimension(2*bc_size)        , intent(in)  :: interior_profile
         integer(ikind)                       , intent(in)  :: j_min6
         integer(ikind)                       , intent(in)  :: j_min7
         integer(ikind)                       , intent(in)  :: j_match_borderN
         integer(ikind)                       , intent(in)  :: interior_j_max3

         integer(ikind) :: i,j

         do j=1, interior_j_max3-bc_size
            do i=1, size(new_grdpts_id,1)-2*bc_size
               new_grdpts_id(i,j_min6+j) = no_pt
            end do
            do i=size(new_grdpts_id,1)-2*bc_size+1, size(new_grdpts_id,1)
               new_grdpts_id(i,j_min6+j) = interior_profile(i-(size(new_grdpts_id,1)-2*bc_size))
            end do            
         end do

         do j=1, min(interior_j_max3,bc_size)
            do i=1, size(new_grdpts_id,1)-2*bc_size
               new_grdpts_id(i,j_min7+j) = no_pt
            end do
            do i=size(new_grdpts_id,1)-2*bc_size+1, size(new_grdpts_id,1)
               new_grdpts_id(i,j_min7+j) = border_N(i-(size(new_grdpts_id,1)-2*bc_size),j_match_borderN+j)
            end do
         end do         

       end subroutine add_grdpts_id_edge_N_blocks_W

      end module bf_layer_merge_module
