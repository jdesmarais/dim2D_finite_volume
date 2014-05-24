      module bf_layer_merge_module
      
        use bf_layer_reallocate_module, only : get_additional_blocks_N
        use parameters_bf_layer       , only : no_pt
        use parameters_constant       , only : x_direction, y_direction
        use parameters_input          , only : debug, nx, ny, ne, bc_size
        use parameters_kind           , only : ikind, rkind

        private
        public :: merge_bf_layers_N,
     $            get_new_size
c$$$     $            merge_bf_layers_S,
c$$$     $            merge_bf_layers_E,
c$$$     $            merge_bf_layers_W,


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
          integer                                       :: outside_i_max1
          integer                                       :: outside_i_max2
          integer                                       :: interior_i_max1
          integer                                       :: interior_i_max2
          integer                                       :: interior_i_max3
          integer                                       :: i_min1
          integer                                       :: i_min3
          integer                                       :: i_min4
          integer                                       :: i_min5
          integer                                       :: i_min6
          integer                                       :: i_min7
          integer                                       :: i_min8
          integer                                       :: j_min1
          integer                                       :: j_min2
          integer                                       :: i_match_borderE


          !get the new size of the tables
          !and the indices to match easily during the copy of the tables
          if(present(final_alignment_i)) then

             final_alignment      = final_alignment_i
             final_alignment(2,1) = ny+1

             new_size = get_new_size(alignment1,
     $                               alignment2,
     $                               final_alignment)

             call get_match(
     $            y_direction,
     $            outside_i_max1, outside_i_max2,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_min1, i_min3, i_min4, i_min5, i_min6, i_min7, i_min8,
     $            i_match_borderE,
     $            j_min1, j_min2,
     $            alignment1, alignment2, final_alignment)

             bf_alignment(1,1) = min(alignment1(1,1),
     $                               alignment2(1,1),
     $                               final_alignment(1,1))

             bf_alignment(2,1) = ny+1

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
     $            i_min1, i_min3, i_min4, i_min5, i_min6, i_min7, i_min8,
     $            i_match_borderE,
     $            j_min1, j_min2,
     $            alignment1, alignment2)
             
             bf_alignment(1,1) = min(alignment1(1,1), alignment2(1,1))
             bf_alignment(2,1) = ny+1
             bf_alignment(1,2) = max(alignment1(1,2), alignment2(1,2))
             bf_alignment(2,2) = max(alignment1(2,2), alignment2(2,2))

          end if


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
     $         i_min1, i_min3, i_min4, i_min5, i_min6, i_min7, i_min8,
     $         j_min1, j_min2,
     $         i_match_borderE)
          deallocate(grdpts_id2)
          call MOVE_ALLOC(new_grdpts_id,grdpts_id1)


          !update the alignment
          alignment1 = bf_alignment          

        end subroutine merge_bf_layers_N


c$$$        !< merge southern buffer layers
c$$$        subroutine merge_bf_layers_S(
c$$$     $       nodes1, nodes2,
c$$$     $       grdpts_id1, grdpts_id2,
c$$$     $       alignment1, alignment2, final_alignment_i)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes1
c$$$          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes2
c$$$          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id1
c$$$          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id2
c$$$          integer(ikind), dimension(2,2)               , intent(inout) :: alignment1
c$$$          integer(ikind), dimension(2,2)               , intent(in)    :: alignment2
c$$$          integer(ikind), dimension(2,2)  , optional   , intent(in)    :: final_alignment_i
c$$$
c$$$          integer(ikind), dimension(2,2)                :: final_alignment
c$$$          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
c$$$          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
c$$$          integer(ikind), dimension(2) :: new_size
c$$$          integer                      :: interior_i_max1
c$$$          integer                      :: interior_i_max2
c$$$          integer                      :: interior_i_max3
c$$$          integer                      :: i_min1
c$$$          integer                      :: i_min2
c$$$          integer                      :: i_min3
c$$$          integer                      :: i_min4
c$$$          integer                      :: j_min1
c$$$          integer                      :: j_min2
c$$$          integer                      :: j_max
c$$$
c$$$
c$$$          !get the new size of the tables
c$$$          !and the indices to match easily during the copy of the tables
c$$$          if(present(final_alignment_i)) then
c$$$
c$$$             final_alignment      = final_alignment_i
c$$$             final_alignment(2,2) = 0
c$$$             
c$$$             new_size = get_new_size(alignment1, alignment2, final_alignment)
c$$$
c$$$             call get_match(
c$$$     $            y_direction,
c$$$     $            interior_i_max1, interior_i_max2, interior_i_max3,
c$$$     $            i_min1, i_min2, i_min3, i_min4,
c$$$     $            j_min1, j_min2,
c$$$     $            alignment1, alignment2, final_alignment)
c$$$
c$$$          else
c$$$
c$$$             new_size = get_new_size(alignment1, alignment2)
c$$$
c$$$             call get_match(
c$$$     $            y_direction,
c$$$     $            interior_i_max1, interior_i_max2, interior_i_max3,
c$$$     $            i_min1, i_min2, i_min3, i_min4,
c$$$     $            j_min1, j_min2,
c$$$     $            alignment1, alignment2)
c$$$          end if
c$$$          j_max = new_size(2)
c$$$
c$$$
c$$$          !allocate the nodes and copy the tables
c$$$          allocate(new_nodes(new_size(1), new_size(2), ne))
c$$$          call merge_nodes_S(
c$$$     $         new_nodes,
c$$$     $         nodes1, nodes2,
c$$$     $         alignment1, alignment2,
c$$$     $         i_min1, i_min3,
c$$$     $         j_min1, j_min2, j_max)
c$$$          deallocate(nodes2)
c$$$          call MOVE_ALLOC(new_nodes,nodes1)
c$$$
c$$$
c$$$          !allocate the gridpts_id and copy the tables
c$$$          allocate(new_grdpts_id(new_size(1), new_size(2)))
c$$$          call merge_grdpts_id_S(
c$$$     $         new_grdpts_id,
c$$$     $         grdpts_id1, grdpts_id2,
c$$$     $         alignment1, alignment2,
c$$$     $         interior_i_max1, interior_i_max2, interior_i_max3,
c$$$     $         i_min1, i_min2, i_min3, i_min4,
c$$$     $         j_min1, j_min2, j_max)
c$$$          deallocate(grdpts_id2)
c$$$          call MOVE_ALLOC(new_grdpts_id,grdpts_id1)
c$$$
c$$$          !update the alignment
c$$$          if(present(final_alignment_i)) then
c$$$             alignment1(1,1) = min(alignment1(1,1), alignment2(1,1), final_alignment(1,1))
c$$$             alignment1(2,1) = min(alignment1(2,1), alignment2(2,1), final_alignment(2,1))
c$$$             alignment1(1,2) = max(alignment1(1,2), alignment2(1,2), final_alignment(1,2))
c$$$             alignment1(2,2) = 0
c$$$          else
c$$$             alignment1(1,1) = min(alignment1(1,1), alignment2(1,1))
c$$$             alignment1(2,1) = min(alignment1(2,1), alignment2(2,1))
c$$$             alignment1(1,2) = max(alignment1(1,2), alignment2(1,2))
c$$$             alignment1(2,2) = 0
c$$$          end if
c$$$
c$$$        end subroutine merge_bf_layers_S
c$$$
c$$$
c$$$        !< merge eastern buffer layers
c$$$        subroutine merge_bf_layers_E(
c$$$     $       nodes1, nodes2,
c$$$     $       grdpts_id1, grdpts_id2,
c$$$     $       alignment1, alignment2, final_alignment_i)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes1
c$$$          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes2
c$$$          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id1
c$$$          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id2
c$$$          integer(ikind), dimension(2,2)               , intent(inout) :: alignment1
c$$$          integer(ikind), dimension(2,2)               , intent(in)    :: alignment2
c$$$          integer(ikind), dimension(2,2)  , optional   , intent(in)    :: final_alignment_i
c$$$
c$$$
c$$$          integer(ikind), dimension(2,2)                :: final_alignment
c$$$          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
c$$$          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
c$$$          integer(ikind), dimension(2) :: new_size
c$$$          integer                      :: interior_j_max1
c$$$          integer                      :: interior_j_max2
c$$$          integer                      :: interior_j_max3
c$$$          integer                      :: j_min1
c$$$          integer                      :: j_min2
c$$$          integer                      :: j_min3
c$$$          integer                      :: j_min4
c$$$          integer                      :: i_min1
c$$$          integer                      :: i_min2
c$$$          integer                      :: i_max
c$$$
c$$$
c$$$          !get the new size of the tables
c$$$          !and the indices to match easily during the copy of the tables
c$$$          if(present(final_alignment_i)) then
c$$$
c$$$             final_alignment      = final_alignment_i
c$$$             final_alignment(1,1) = nx+1
c$$$
c$$$             new_size = get_new_size(alignment1, alignment2, final_alignment)
c$$$
c$$$             call get_match(
c$$$     $            x_direction,
c$$$     $            interior_j_max1, interior_j_max2, interior_j_max3,
c$$$     $            j_min1, j_min2, j_min3, j_min4,
c$$$     $            i_min1, i_min2,
c$$$     $            alignment1, alignment2, final_alignment)
c$$$
c$$$          else
c$$$
c$$$             new_size = get_new_size(alignment1, alignment2)
c$$$
c$$$             call get_match(
c$$$     $            x_direction,
c$$$     $            interior_j_max1, interior_j_max2, interior_j_max3,
c$$$     $            j_min1, j_min2, j_min3, j_min4,
c$$$     $            i_min1, i_min2,
c$$$     $            alignment1, alignment2)
c$$$          end if
c$$$          i_max = new_size(2)
c$$$
c$$$          !allocate the nodes and copy the tables
c$$$          allocate(new_nodes(new_size(1), new_size(2), ne))
c$$$          call merge_nodes_E(
c$$$     $         new_nodes,
c$$$     $         nodes1, nodes2,
c$$$     $         alignment1, alignment2,
c$$$     $         j_min1, j_min3)
c$$$          deallocate(nodes2)
c$$$          call MOVE_ALLOC(new_nodes,nodes1)
c$$$
c$$$
c$$$          !allocate the gridpts_id and copy the tables
c$$$          allocate(new_grdpts_id(new_size(1), new_size(2)))
c$$$          call merge_grdpts_id_E(
c$$$     $         new_grdpts_id,
c$$$     $         grdpts_id1, grdpts_id2,
c$$$     $         alignment1, alignment2,
c$$$     $         interior_j_max1, interior_j_max2, interior_j_max3,
c$$$     $         j_min1, j_min2, j_min3, j_min4)
c$$$          deallocate(grdpts_id2)
c$$$          call MOVE_ALLOC(new_grdpts_id,grdpts_id1)
c$$$
c$$$          !update the alignment
c$$$          if(present(final_alignment_i)) then
c$$$             alignment1(1,1) = nx+1
c$$$             alignment1(2,1) = min(alignment1(2,1), alignment2(2,1), final_alignment(2,1))
c$$$             alignment1(1,2) = max(alignment1(1,2), alignment2(1,2), final_alignment(1,2))
c$$$             alignment1(2,2) = max(alignment1(2,2), alignment2(2,2), final_alignment(2,2))
c$$$          else
c$$$             alignment1(1,1) = nx+1
c$$$             alignment1(2,1) = min(alignment1(2,1), alignment2(2,1))
c$$$             alignment1(1,2) = max(alignment1(1,2), alignment2(1,2))
c$$$             alignment1(2,2) = max(alignment1(2,2), alignment2(2,2))
c$$$          end if
c$$$
c$$$        end subroutine merge_bf_layers_E
c$$$
c$$$
c$$$        !< merge western buffer layers
c$$$        subroutine merge_bf_layers_W(
c$$$     $       nodes1, nodes2,
c$$$     $       grdpts_id1, grdpts_id2,
c$$$     $       alignment1, alignment2, final_alignment_i)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes1
c$$$          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes2
c$$$          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id1
c$$$          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id2
c$$$          integer(ikind), dimension(2,2)               , intent(inout) :: alignment1
c$$$          integer(ikind), dimension(2,2)               , intent(in)    :: alignment2
c$$$          integer(ikind), dimension(2,2)  , optional   , intent(in)    :: final_alignment_i
c$$$
c$$$
c$$$          integer(ikind), dimension(2,2)                :: final_alignment
c$$$          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
c$$$          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
c$$$          integer(ikind), dimension(2) :: new_size
c$$$          integer                      :: interior_j_max1
c$$$          integer                      :: interior_j_max2
c$$$          integer                      :: interior_j_max3
c$$$          integer                      :: j_min1
c$$$          integer                      :: j_min2
c$$$          integer                      :: j_min3
c$$$          integer                      :: j_min4
c$$$          integer                      :: i_min1
c$$$          integer                      :: i_min2
c$$$
c$$$
c$$$          !get the new size of the tables
c$$$          !and the indices to match easily during the copy of the tables
c$$$          if(present(final_alignment_i)) then
c$$$
c$$$             final_alignment      = final_alignment_i
c$$$             final_alignment(1,2) = 0
c$$$
c$$$             new_size = get_new_size(alignment1, alignment2, final_alignment)
c$$$
c$$$             call get_match(
c$$$     $            x_direction,
c$$$     $            interior_j_max1, interior_j_max2, interior_j_max3,
c$$$     $            j_min1, j_min2, j_min3, j_min4,
c$$$     $            i_min1, i_min2,
c$$$     $            alignment1, alignment2, final_alignment)
c$$$
c$$$          else
c$$$
c$$$             new_size = get_new_size(alignment1, alignment2)
c$$$
c$$$             call get_match(
c$$$     $            x_direction,
c$$$     $            interior_j_max1, interior_j_max2, interior_j_max3,
c$$$     $            j_min1, j_min2, j_min3, j_min4,
c$$$     $            i_min1, i_min2,
c$$$     $            alignment1, alignment2)
c$$$          end if
c$$$
c$$$          !allocate the nodes and copy the tables
c$$$          allocate(new_nodes(new_size(1), new_size(2), ne))
c$$$          call merge_nodes_W(
c$$$     $         new_nodes,
c$$$     $         nodes1, nodes2,
c$$$     $         alignment1, alignment2,
c$$$     $         j_min1, j_min3)
c$$$          deallocate(nodes2)
c$$$          call MOVE_ALLOC(new_nodes,nodes1)
c$$$
c$$$
c$$$          !allocate the gridpts_id and copy the tables
c$$$          allocate(new_grdpts_id(new_size(1), new_size(2)))
c$$$          call merge_grdpts_id_W(
c$$$     $         new_grdpts_id,
c$$$     $         grdpts_id1, grdpts_id2,
c$$$     $         alignment1, alignment2,
c$$$     $         interior_j_max1, interior_j_max2, interior_j_max3,
c$$$     $         j_min1, j_min2, j_min3, j_min4)
c$$$          deallocate(grdpts_id2)
c$$$          call MOVE_ALLOC(new_grdpts_id,grdpts_id1)
c$$$
c$$$          !update the alignment
c$$$          if(present(final_alignment_i)) then
c$$$             alignment1(1,1) = min(alignment1(1,1), alignment2(1,1), final_alignment(1,1))
c$$$             alignment1(2,1) = min(alignment1(2,1), alignment2(2,1), final_alignment(2,1))
c$$$             alignment1(1,2) = 0
c$$$             alignment1(2,2) = max(alignment1(2,2), alignment2(2,2), final_alignment(2,2))
c$$$          else
c$$$             alignment1(1,1) = min(alignment1(1,1), alignment2(1,1), final_alignment(1,1))
c$$$             alignment1(2,1) = min(alignment1(2,1), alignment2(2,1))
c$$$             alignment1(1,2) = 0
c$$$             alignment1(2,2) = max(alignment1(2,2), alignment2(2,2))
c$$$          end if
c$$$
c$$$        end subroutine merge_bf_layers_W


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
     $     i_min1, i_min3, i_min4, i_min5, i_min6, i_min7, i_min8,
     $     i_match_borderE,
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
          integer(ikind)                          , intent(out) :: i_min7
          integer(ikind)                          , intent(out) :: i_min8
          integer(ikind)                          , intent(out) :: i_match_borderE
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
          i_min7 = i_min6 + max(0, interior_i_max3-bc_size)
          i_min8 = i_min7 + interior_i_max3

          j_min1 = min((alignment1(dir2,2)-alignment1(dir2,1)+2*bc_size+1),
     $                 (alignment2(dir2,2)-alignment2(dir2,1)+2*bc_size+1))
          j_min2 = max((alignment1(dir2,2)-alignment1(dir2,1)+2*bc_size+1),
     $                 (alignment2(dir2,2)-alignment2(dir2,1)+2*bc_size+1))

          if(present(final_alignment)) then
             i_match_borderE = max(0,
     $            min(alignment1(dir1,1),
     $                alignment2(dir1,1),
     $                final_alignment(dir1,1))-
     $            (bc_size+1)+i_min7-(ndir-bc_size))
          else
             i_match_borderE = max(0,
     $            min(alignment1(dir1,1),
     $                alignment2(dir1,1))-
     $            (bc_size+1)+i_min7-(ndir-bc_size))
          end if
        end subroutine get_match

        
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

          !nodes1 - nodes2
          if(alignment1(1,1).lt.alignment2(1,1)) then

             if(alignment1(2,2).gt.alignment2(2,2)) then

                do k=1,ne
                   
                   call add_nodes_blocks_2_to_8_NS(
     $                  new_nodes,
     $                  nodes1, nodes2, interior_nodes,
     $                  bf_alignment,
     $                  k, 1, 2*bc_size,
     $                  interior_i_max1,interior_i_max2,interior_i_max3,
     $                  i_min1, i_min3, i_min4, i_min5, i_min6)

                   call add_nodes_blocks_11_to_13_NS(
     $                  new_nodes,
     $                  nodes1, nodes2,
     $                  k, 2*bc_size+1, j_min1,
     $                  i_min3, i_min5)
                   
                   call add_nodes_blocks_16_to_18_NS(
     $                  new_nodes,
     $                  nodes1,
     $                  k, j_min1+1, j_min2,
     $                  i_min3)

                end do

             else
                do k=1, ne
                   call add_nodes_blocks_2_to_8_NS(
     $                  new_nodes,
     $                  nodes1, nodes2, interior_nodes,
     $                  bf_alignment,
     $                  k, 1, 2*bc_size,
     $                  interior_i_max1,interior_i_max2,interior_i_max3,
     $                  i_min1, i_min3, i_min4, i_min5, i_min6)

                   call add_nodes_blocks_11_to_13_NS(
     $                  new_nodes,
     $                  nodes1, nodes2,
     $                  k, 2*bc_size+1, j_min1,
     $                  i_min3, i_min5)
                   
                   call add_nodes_blocks_16_to_18_NS(
     $                  new_nodes,
     $                  nodes2,
     $                  k, j_min1+1, j_min2,
     $                  i_min5)
                end do

             end if

          !nodes2 - nodes1
          else
             if(alignment1(2,2).gt.alignment2(2,2)) then
                do k=1, ne

                   call add_nodes_blocks_2_to_8_NS(
     $                  new_nodes,
     $                  nodes2, nodes1, interior_nodes,
     $                  bf_alignment,
     $                  k, 1, 2*bc_size,
     $                  interior_i_max1,interior_i_max2,interior_i_max3,
     $                  i_min1, i_min3, i_min4, i_min5, i_min6)

                   call add_nodes_blocks_11_to_13_NS(
     $                  new_nodes,
     $                  nodes2, nodes1,
     $                  k, 2*bc_size+1, j_min1,
     $                  i_min3, i_min5)
                   
                   call add_nodes_blocks_16_to_18_NS(
     $                  new_nodes,
     $                  nodes1,
     $                  k, j_min1+1, j_min2,
     $                  i_min5)

                end do
             else
                do k=1, ne
                   call add_nodes_blocks_2_to_8_NS(
     $                  new_nodes,
     $                  nodes2, nodes1, interior_nodes,
     $                  bf_alignment,
     $                  k, 1, 2*bc_size,
     $                  interior_i_max1,interior_i_max2,interior_i_max3,
     $                  i_min1, i_min3, i_min4, i_min5, i_min6)

                   call add_nodes_blocks_11_to_13_NS(
     $                  new_nodes,
     $                  nodes2, nodes1,
     $                  k, 2*bc_size+1, j_min1,
     $                  i_min3, i_min5)
                   
                   call add_nodes_blocks_16_to_18_NS(
     $                  new_nodes,
     $                  nodes2,
     $                  k, j_min1+1, j_min2,
     $                  i_min3)
                end do
             end if
          end if

        end subroutine merge_nodes_N


c$$$        !> merge the nodes for southern buffer layers
c$$$        subroutine merge_nodes_S(
c$$$     $     new_nodes,
c$$$     $     nodes1, nodes2,
c$$$     $     alignment1, alignment2,
c$$$     $     i_min1, i_min3,
c$$$     $     j_min1, j_min2, j_max)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind), dimension(:,:,:) , intent(out):: new_nodes
c$$$          real(rkind), dimension(:,:,:) , intent(in) :: nodes1
c$$$          real(rkind), dimension(:,:,:) , intent(in) :: nodes2
c$$$          integer(ikind), dimension(2,2), intent(in) :: alignment1
c$$$          integer(ikind), dimension(2,2), intent(in) :: alignment2
c$$$          integer(ikind)                , intent(in) :: i_min1
c$$$          integer(ikind)                , intent(in) :: i_min3
c$$$          integer(ikind)                , intent(in) :: j_min1
c$$$          integer(ikind)                , intent(in) :: j_min2
c$$$          integer(ikind)                , intent(in) :: j_max
c$$$
c$$$          integer(ikind) :: i,j
c$$$          integer        :: k
c$$$
c$$$          !nodes1 - nodes2
c$$$          if(alignment1(1,1).lt.alignment2(1,1)) then
c$$$
c$$$             if(alignment1(2,1).lt.alignment2(2,1)) then
c$$$                do k=1, ne
c$$$                   do j=j_max-j_min2+1, j_max-j_min1
c$$$                      do i=1, size(nodes1,1)
c$$$                         new_nodes(i_min1+i,j,k) = nodes1(i,j-(j_max-j_min2),k)
c$$$                      end do
c$$$                   end do
c$$$
c$$$                   do j=j_max-j_min1+1, j_max
c$$$                      do i=1, size(nodes1,1)
c$$$                         new_nodes(i_min1+i,j,k) = nodes1(i,j-(j_max-j_min2),k)
c$$$                      end do
c$$$                      
c$$$                      do i=1, size(nodes2,1)
c$$$                         new_nodes(i_min3+i,j,k) = nodes2(i,j-(j_max-j_min1),k)
c$$$                      end do
c$$$                   end do
c$$$                end do                
c$$$             else
c$$$                do k=1, ne
c$$$                   do j=j_max-j_min2+1, j_max-j_min1
c$$$                      do i=1, size(nodes2,1)
c$$$                         new_nodes(i_min3+i,j,k) = nodes2(i,j-(j_max-j_min2),k)
c$$$                      end do
c$$$                   end do
c$$$
c$$$                   do j=j_max-j_min1+1, j_max
c$$$                      do i=1, size(nodes1,1)
c$$$                         new_nodes(i_min1+i,j,k) = nodes1(i,j-(j_max-j_min1),k)
c$$$                      end do
c$$$                      
c$$$                      do i=1, size(nodes2,1)
c$$$                         new_nodes(i_min3+i,j,k) = nodes2(i,j-(j_max-j_min2),k)
c$$$                      end do
c$$$                   end do
c$$$                end do
c$$$             end if            
c$$$
c$$$          !nodes2 - nodes1
c$$$          else
c$$$
c$$$             if(alignment1(2,1).lt.alignment2(2,1)) then
c$$$                !print *, 'nodes:', 'x1>x2 : y1>y2'
c$$$                do k=1, ne
c$$$                   do j=j_max-j_min2+1, j_max-j_min1
c$$$                      do i=1, size(nodes1,1)
c$$$                         new_nodes(i_min3+i,j,k) = nodes1(i,j-(j_max-j_min2),k)
c$$$                      end do
c$$$                   end do
c$$$
c$$$                   do j=j_max-j_min1+1, j_max
c$$$                      do i=1, size(nodes2,1)
c$$$                         new_nodes(i_min1+i,j,k) = nodes2(i,j-(j_max-j_min1),k)
c$$$                      end do
c$$$
c$$$                      do i=1, size(nodes1,1)
c$$$                         new_nodes(i_min3+i,j,k) = nodes1(i,j-(j_max-j_min2),k)
c$$$                      end do
c$$$                   end do
c$$$                end do
c$$$             else
c$$$                do k=1, ne
c$$$                   do j=j_max-j_min2, j_max-j_min1
c$$$                      do i=1, size(nodes2,1)
c$$$                         new_nodes(i_min1+i,j,k) = nodes2(i,j-(j_max-j_min2),k)
c$$$                      end do
c$$$                   end do
c$$$
c$$$                   do j=j_max-j_min1+1, j_max
c$$$                      do i=1, size(nodes2,1)
c$$$                         new_nodes(i_min1+i,j,k) = nodes2(i,j-(j_max-j_min2),k)
c$$$                      end do
c$$$
c$$$                      do i=1, size(nodes1,1)
c$$$                         new_nodes(i_min3+i,j,k) = nodes1(i,j-(j_max-j_min1),k)
c$$$                      end do
c$$$                   end do
c$$$                end do
c$$$             end if
c$$$          end if
c$$$
c$$$        end subroutine merge_nodes_S
c$$$
c$$$
c$$$        !> merge the nodes for eastern buffer layers
c$$$        subroutine merge_nodes_E(
c$$$     $     new_nodes,
c$$$     $     nodes1, nodes2,
c$$$     $     alignment1, alignment2,
c$$$     $     j_min1, j_min3)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind), dimension(:,:,:) , intent(out):: new_nodes
c$$$          real(rkind), dimension(:,:,:) , intent(in) :: nodes1
c$$$          real(rkind), dimension(:,:,:) , intent(in) :: nodes2
c$$$          integer(ikind), dimension(2,2), intent(in) :: alignment1
c$$$          integer(ikind), dimension(2,2), intent(in) :: alignment2
c$$$          integer(ikind)                , intent(in) :: j_min1
c$$$          integer(ikind)                , intent(in) :: j_min3
c$$$
c$$$          integer(ikind) :: i,j
c$$$          integer        :: k
c$$$
c$$$          !nodes1 - nodes2
c$$$          if(alignment1(2,1).lt.alignment2(2,1)) then
c$$$             print *, 'y1<y2'
c$$$
c$$$             do k=1, ne
c$$$
c$$$                do j=1, size(nodes1,2)
c$$$                   do i=1, size(nodes1,1)
c$$$                      new_nodes(i,j_min1+j,k) = nodes1(i,j,k)
c$$$                   end do
c$$$                end do
c$$$                   
c$$$                do j=1, size(nodes2,1)
c$$$                   do i=1, size(nodes2,1)
c$$$                      new_nodes(i,j_min3+j,k) = nodes2(i,j,k)
c$$$                   end do
c$$$                end do
c$$$
c$$$             end do
c$$$
c$$$          else
c$$$             print *, 'y2<y1'
c$$$
c$$$             do k=1, ne
c$$$
c$$$                do j=1, size(nodes2,2)
c$$$                   do i=1, size(nodes2,1)
c$$$                      new_nodes(i,j_min1+j,k) = nodes2(i,j,k)
c$$$                   end do
c$$$                end do
c$$$
c$$$                do j=1, size(nodes1,2)
c$$$                   do i=1, size(nodes1,1)
c$$$                      new_nodes(i,j_min3+j,k) = nodes1(i,j,k)
c$$$                   end do
c$$$                end do
c$$$
c$$$             end do
c$$$
c$$$          end if
c$$$
c$$$        end subroutine merge_nodes_E
c$$$
c$$$
c$$$        !> merge the nodes for western buffer layers
c$$$        subroutine merge_nodes_W(
c$$$     $     new_nodes,
c$$$     $     nodes1, nodes2,
c$$$     $     alignment1, alignment2,
c$$$     $     j_min1, j_min3)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind), dimension(:,:,:) , intent(out):: new_nodes
c$$$          real(rkind), dimension(:,:,:) , intent(in) :: nodes1
c$$$          real(rkind), dimension(:,:,:) , intent(in) :: nodes2
c$$$          integer(ikind), dimension(2,2), intent(in) :: alignment1
c$$$          integer(ikind), dimension(2,2), intent(in) :: alignment2
c$$$          integer(ikind)                , intent(in) :: j_min1
c$$$          integer(ikind)                , intent(in) :: j_min3
c$$$
c$$$          integer(ikind) :: i,j, i_min1, i_min2
c$$$          integer        :: k
c$$$
c$$$          !nodes1 - nodes2
c$$$          if(alignment1(2,1).lt.alignment2(2,1)) then
c$$$             !print *, 'y1<y2'
c$$$
c$$$             i_min1 = size(new_nodes,1)-size(nodes1,1)
c$$$             i_min2 = size(new_nodes,1)-size(nodes2,1)
c$$$
c$$$             do k=1, ne
c$$$                
c$$$                do j=1, size(nodes1,2)
c$$$                   do i=i_min1+1, size(new_nodes,1)
c$$$                      new_nodes(i,j_min1+j,k) = nodes1(i-i_min1,j,k)
c$$$                   end do
c$$$                end do
c$$$                   
c$$$                do j=1, size(nodes2,1)
c$$$                   do i=i_min2+1, size(new_nodes,1)
c$$$                      new_nodes(i,j_min3+j,k) = nodes2(i-i_min2,j,k)
c$$$                   end do
c$$$                end do
c$$$
c$$$             end do
c$$$
c$$$          else
c$$$             !print *, 'y2<y1'
c$$$
c$$$             i_min1 = size(new_nodes,1)-size(nodes2,1)
c$$$             i_min2 = size(new_nodes,1)-size(nodes1,1)
c$$$
c$$$             do k=1, ne
c$$$
c$$$                do j=1, size(nodes2,2)
c$$$                   do i=i_min1+1, size(new_nodes,1)
c$$$                      new_nodes(i,j_min1+j,k) = nodes2(i-i_min1,j,k)
c$$$                   end do
c$$$                end do
c$$$
c$$$                do j=1, size(nodes1,2)
c$$$                   do i=i_min2+1, size(new_nodes,1)
c$$$                      new_nodes(i,j_min3+j,k) = nodes1(i-i_min2,j,k)
c$$$                   end do
c$$$                end do
c$$$
c$$$             end do
c$$$
c$$$          end if
c$$$
c$$$        end subroutine merge_nodes_W
c$$$
c$$$
        !> merge the nodes for northern buffer layers
        subroutine merge_grdpts_id_N(
     $     new_grdpts_id,
     $     grdpts_id1, grdpts_id2,
     $     alignment1, alignment2, bf_alignment,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min1, i_min3, i_min4, i_min5, i_min6, i_min7, i_min8,
     $     j_min1, j_min2,
     $     i_match_borderE)

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
          integer(ikind)                , intent(in) :: i_min7
          integer(ikind)                , intent(in) :: i_min8
          integer(ikind)                , intent(in) :: j_min1
          integer(ikind)                , intent(in) :: j_min2
          integer(ikind)                , intent(in) :: i_match_borderE


          integer, dimension(bc_size, 2*bc_size) :: border_W
          integer, dimension(bc_size, 2*bc_size) :: border_E
          integer, dimension(2*bc_size)          :: interior_profile

          
          !get the additional blocks
          call get_additional_blocks_N(
     $         bf_alignment,
     $         border_W, border_E, interior_profile)
          

          !nodes1 - nodes2
          if(alignment1(1,1).lt.alignment2(1,1)) then

             if(alignment1(2,2).gt.alignment2(2,2)) then
                
                call add_grdpts_id_blocks_1_to_9_NS(
     $               new_grdpts_id,
     $               grdpts_id1, grdpts_id2,
     $               border_W, border_E, interior_profile,
     $               1, 2*bc_size,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min1, i_min3, i_min4, i_min5, i_min6, i_min7, i_min8,
     $               i_match_borderE)

                call add_grdpts_id_blocks_10_to_14_NS(
     $               new_grdpts_id,
     $               grdpts_id1, grdpts_id2,
     $               2*bc_size+1, j_min1,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4, i_min5, i_min6)

                call add_grdpts_id_blocks_15_to_19_with_16_taller_NS(
     $               new_grdpts_id, 
     $               grdpts_id1, grdpts_id2,
     $               j_min1+1, j_min2,
     $               outside_i_max1,  outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4)
                
                call add_grdpts_id_block_20_NS(
     $               new_grdpts_id,
     $               j_min2+1, size(new_grdpts_id,2))
                
             else

                call add_grdpts_id_blocks_1_to_9_NS(
     $               new_grdpts_id,
     $               grdpts_id1, grdpts_id2,
     $               border_W, border_E, interior_profile,
     $               1, 2*bc_size,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min1, i_min3, i_min4, i_min5, i_min6, i_min7, i_min8,
     $               i_match_borderE)

                call add_grdpts_id_blocks_10_to_14_NS(
     $               new_grdpts_id,
     $               grdpts_id1, grdpts_id2,
     $               2*bc_size+1, j_min1,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4, i_min5, i_min6)

                call add_grdpts_id_blocks_15_to_19_with_18_taller_NS(
     $               new_grdpts_id, 
     $               grdpts_id1, grdpts_id2,
     $               j_min1+1, j_min2,
     $               outside_i_max1,  outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min5, i_min6)
                
                call add_grdpts_id_block_20_NS(
     $               new_grdpts_id,
     $               j_min2+1, size(new_grdpts_id,2))

             end if
          else
             
             if(alignment1(2,2).gt.alignment2(2,2)) then

                call add_grdpts_id_blocks_1_to_9_NS(
     $               new_grdpts_id,
     $               grdpts_id2, grdpts_id1,
     $               border_W, border_E, interior_profile,
     $               1, 2*bc_size,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min1, i_min3, i_min4, i_min5, i_min6, i_min7, i_min8,
     $               i_match_borderE)

                call add_grdpts_id_blocks_10_to_14_NS(
     $               new_grdpts_id,
     $               grdpts_id2, grdpts_id1,
     $               2*bc_size+1, j_min1,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4, i_min5, i_min6)

                call add_grdpts_id_blocks_15_to_19_with_16_taller_NS(
     $               new_grdpts_id, 
     $               grdpts_id2, grdpts_id1,
     $               j_min1+1, j_min2,
     $               outside_i_max1,  outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4)
                
                call add_grdpts_id_block_20_NS(
     $               new_grdpts_id,
     $               j_min2+1, size(new_grdpts_id,2))
                
             else

                call add_grdpts_id_blocks_1_to_9_NS(
     $               new_grdpts_id,
     $               grdpts_id2, grdpts_id1,
     $               border_W, border_E, interior_profile,
     $               1, 2*bc_size,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min1, i_min3, i_min4, i_min5, i_min6, i_min7, i_min8,
     $               i_match_borderE)

                call add_grdpts_id_blocks_10_to_14_NS(
     $               new_grdpts_id,
     $               grdpts_id2, grdpts_id1,
     $               2*bc_size+1, j_min1,
     $               outside_i_max1, outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min3, i_min4, i_min5, i_min6)

                call add_grdpts_id_blocks_15_to_19_with_18_taller_NS(
     $               new_grdpts_id, 
     $               grdpts_id2, grdpts_id1,
     $               j_min1+1, j_min2,
     $               outside_i_max1,  outside_i_max2,
     $               interior_i_max1, interior_i_max2, interior_i_max3,
     $               i_min5, i_min6)
                
                call add_grdpts_id_block_20_NS(
     $               new_grdpts_id,
     $               j_min2+1, size(new_grdpts_id,2))
             end if
          end if

        end subroutine merge_grdpts_id_N


c$$$        !> merge the nodes for southern buffer layers
c$$$        subroutine merge_grdpts_id_S(
c$$$     $     new_grdpts_id,
c$$$     $     grdpts_id1, grdpts_id2,
c$$$     $     alignment1, alignment2,
c$$$     $     interior_i_max1, interior_i_max2, interior_i_max3,
c$$$     $     i_min1, i_min2, i_min3, i_min4,
c$$$     $     j_min1, j_min2, j_max)
c$$$
c$$$          implicit none
c$$$
c$$$          integer       , dimension(:,:), intent(out):: new_grdpts_id
c$$$          integer       , dimension(:,:), intent(in) :: grdpts_id1
c$$$          integer       , dimension(:,:), intent(in) :: grdpts_id2
c$$$          integer(ikind), dimension(2,2), intent(in) :: alignment1
c$$$          integer(ikind), dimension(2,2), intent(in) :: alignment2
c$$$          integer(ikind)                , intent(in) :: interior_i_max1
c$$$          integer(ikind)                , intent(in) :: interior_i_max2
c$$$          integer(ikind)                , intent(in) :: interior_i_max3
c$$$          integer(ikind)                , intent(in) :: i_min1
c$$$          integer(ikind)                , intent(in) :: i_min2
c$$$          integer(ikind)                , intent(in) :: i_min3
c$$$          integer(ikind)                , intent(in) :: i_min4
c$$$          integer(ikind)                , intent(in) :: j_min1
c$$$          integer(ikind)                , intent(in) :: j_min2
c$$$          integer(ikind)                , intent(in) :: j_max
c$$$
c$$$          integer(ikind) :: i,j
c$$$
c$$$          !grdpts_id1 - grdpts_id2
c$$$          if(alignment1(1,1).lt.alignment2(1,1)) then
c$$$
c$$$             if(alignment1(2,1).lt.alignment2(2,1)) then
c$$$
c$$$                print *, 'x1<x2 : y1>y2'
c$$$                do j=1, j_max-j_min2
c$$$                   do i=1, size(new_grdpts_id,1)
c$$$                      new_grdpts_id(i,j) = no_pt
c$$$                   end do
c$$$                end do
c$$$
c$$$                do j=j_max-j_min2+1, j_max-j_min1
c$$$                   do i=1, interior_i_max1
c$$$                      new_grdpts_id(i,j) = no_pt
c$$$                   end do
c$$$                   
c$$$                   do i=1, size(grdpts_id1,1)
c$$$                      new_grdpts_id(i_min1+i,j) = grdpts_id1(i,j-(j_max-j_min2))
c$$$                   end do
c$$$                   
c$$$                   do i=1, interior_i_max2+size(grdpts_id2,1)+interior_i_max3
c$$$                      new_grdpts_id(i_min2+i,j) = no_pt
c$$$                   end do
c$$$                end do
c$$$                
c$$$                do j=j_max-j_min1+1, j_max
c$$$                   do i=1, interior_i_max1
c$$$                      new_grdpts_id(i,j) = no_pt
c$$$                   end do
c$$$
c$$$                   do i=1, size(grdpts_id1,1)
c$$$                      new_grdpts_id(i_min1+i,j) = grdpts_id1(i,j-(j_max-j_min2))
c$$$                   end do
c$$$
c$$$                   do i=1, interior_i_max2
c$$$                      new_grdpts_id(i_min2+i,j) = no_pt
c$$$                   end do
c$$$                   
c$$$                   do i=1, size(grdpts_id2,1)
c$$$                      new_grdpts_id(i_min3+i,j) = grdpts_id2(i,j-(j_max-j_min1))
c$$$                   end do
c$$$
c$$$                   do i=1, interior_i_max3
c$$$                      new_grdpts_id(i_min4+i,j) = no_pt
c$$$                   end do
c$$$                end do
c$$$             else
c$$$                print *, 'x1<x2 : y1<y2'
c$$$                do j=1, j_max-j_min2
c$$$                   do i=1, size(new_grdpts_id,1)
c$$$                      new_grdpts_id(i,j) = no_pt
c$$$                   end do
c$$$                end do
c$$$
c$$$                do j=j_max-j_min2+1, j_max-j_min1
c$$$                   do i=1, interior_i_max1+size(grdpts_id1,1)+interior_i_max2
c$$$                      new_grdpts_id(i,j) = no_pt
c$$$                   end do
c$$$
c$$$                   do i=1, size(grdpts_id2,1)
c$$$                      new_grdpts_id(i_min3+i,j) = grdpts_id2(i,j-(j_max-j_min2))
c$$$                   end do
c$$$
c$$$                   do i=1, interior_i_max3
c$$$                      new_grdpts_id(i_min4+i,j) = no_pt
c$$$                   end do
c$$$                end do
c$$$
c$$$                do j=j_max-j_min1+1, j_max
c$$$                   do i=1, interior_i_max1
c$$$                      new_grdpts_id(i,j) = no_pt
c$$$                   end do
c$$$
c$$$                   do i=1, size(grdpts_id1,1)
c$$$                      new_grdpts_id(i_min1+i,j) = grdpts_id1(i,j-(j_max-j_min1))
c$$$                   end do
c$$$
c$$$                   do i=1, interior_i_max2
c$$$                      new_grdpts_id(i_min2+i,j) = no_pt
c$$$                   end do
c$$$                   
c$$$                   do i=1, size(grdpts_id2,1)
c$$$                      new_grdpts_id(i_min3+i,j) = grdpts_id2(i,j-(j_max-j_min2))
c$$$                   end do
c$$$
c$$$                   do i=1, interior_i_max3
c$$$                      new_grdpts_id(i_min4+i,j) = no_pt
c$$$                   end do
c$$$                end do
c$$$             end if            
c$$$
c$$$          !grdpts_id2 - grdpts_id1
c$$$          else
c$$$
c$$$             if(alignment1(2,1).lt.alignment2(2,1)) then
c$$$                print *, 'x1>x2 : y1>y2'
c$$$                do j=1, j_max-j_min2
c$$$                   do i=1, size(new_grdpts_id,1)
c$$$                      new_grdpts_id(i,j) = no_pt
c$$$                   end do
c$$$                end do
c$$$
c$$$                do j=j_max-j_min2+1, j_max-j_min1
c$$$                   do i=1, interior_i_max1+size(grdpts_id2,1)+interior_i_max2
c$$$                      new_grdpts_id(i,j) = no_pt
c$$$                   end do
c$$$
c$$$                   do i=1, size(grdpts_id1,1)
c$$$                      new_grdpts_id(i_min3+i,j) = grdpts_id1(i,j-(j_max-j_min2))
c$$$                   end do
c$$$
c$$$                   do i=1, interior_i_max3
c$$$                      new_grdpts_id(i_min4+i,j) = no_pt
c$$$                   end do
c$$$
c$$$                end do
c$$$
c$$$                do j=j_max-j_min1+1, j_max
c$$$                   do i=1, interior_i_max1
c$$$                      new_grdpts_id(i,j) = no_pt
c$$$                   end do
c$$$
c$$$                   do i=1, size(grdpts_id2,1)
c$$$                      new_grdpts_id(i_min1+i,j) = grdpts_id2(i,j-(j_max-j_min1))
c$$$                   end do
c$$$
c$$$                   do i=1, interior_i_max2
c$$$                      new_grdpts_id(i_min2+i,j) = no_pt
c$$$                   end do
c$$$
c$$$                   do i=1, size(grdpts_id1,1)
c$$$                      new_grdpts_id(i_min3+i,j) = grdpts_id1(i,j-(j_max-j_min2))
c$$$                   end do
c$$$
c$$$                   do i=1, interior_i_max3
c$$$                      new_grdpts_id(i_min4+i,j) = no_pt
c$$$                   end do
c$$$
c$$$                end do
c$$$
c$$$             else
c$$$
c$$$                print *, 'x1>x2 : y1<y2'
c$$$                do j=1, j_max-j_min2
c$$$                   do i=1, size(new_grdpts_id,1)
c$$$                      new_grdpts_id(i,j) = no_pt
c$$$                   end do
c$$$                end do
c$$$
c$$$                do j=j_max-j_min2+1, j_max-j_min1
c$$$                   do i=1, interior_i_max1
c$$$                      new_grdpts_id(i,j) = no_pt
c$$$                   end do
c$$$
c$$$                   do i=1, size(grdpts_id2,1)
c$$$                      new_grdpts_id(i_min1+i,j) = grdpts_id2(i,j-(j_max-j_min2))
c$$$                   end do
c$$$
c$$$                   do i=1, interior_i_max2+size(grdpts_id1,1)+interior_i_max3
c$$$                      new_grdpts_id(i_min2+i,j) = no_pt
c$$$                   end do
c$$$                end do
c$$$
c$$$                do j=j_max-j_min1+1, j_max
c$$$                   do i=1, interior_i_max1
c$$$                      new_grdpts_id(i,j) = no_pt
c$$$                   end do
c$$$
c$$$                   do i=1, size(grdpts_id2,1)
c$$$                      new_grdpts_id(i_min1+i,j) = grdpts_id2(i,j-(j_max-j_min2))
c$$$                   end do
c$$$
c$$$                   do i=1, interior_i_max2
c$$$                      new_grdpts_id(i_min2+i,j) = no_pt
c$$$                   end do
c$$$                   
c$$$                   do i=1, size(grdpts_id1,1)
c$$$                      new_grdpts_id(i_min3+i,j) = grdpts_id1(i,j-(j_max-j_min1))
c$$$                   end do
c$$$
c$$$                   do i=1, interior_i_max3
c$$$                      new_grdpts_id(i_min4+i,j) = no_pt
c$$$                   end do
c$$$                end do
c$$$             end if
c$$$
c$$$          end if
c$$$
c$$$        end subroutine merge_grdpts_id_S
c$$$
c$$$        !> merge the nodes for eastern buffer layers
c$$$        subroutine merge_grdpts_id_E(
c$$$     $     new_grdpts_id,
c$$$     $     grdpts_id1, grdpts_id2,
c$$$     $     alignment1, alignment2,
c$$$     $     interior_j_max1, interior_j_max2, interior_j_max3,
c$$$     $     j_min1, j_min2, j_min3, j_min4)
c$$$
c$$$          implicit none
c$$$
c$$$          integer       , dimension(:,:), intent(out):: new_grdpts_id
c$$$          integer       , dimension(:,:), intent(in) :: grdpts_id1
c$$$          integer       , dimension(:,:), intent(in) :: grdpts_id2
c$$$          integer(ikind), dimension(2,2), intent(in) :: alignment1
c$$$          integer(ikind), dimension(2,2), intent(in) :: alignment2
c$$$          integer(ikind)                , intent(in) :: interior_j_max1
c$$$          integer(ikind)                , intent(in) :: interior_j_max2
c$$$          integer(ikind)                , intent(in) :: interior_j_max3
c$$$          integer(ikind)                , intent(in) :: j_min1
c$$$          integer(ikind)                , intent(in) :: j_min2
c$$$          integer(ikind)                , intent(in) :: j_min3
c$$$          integer(ikind)                , intent(in) :: j_min4
c$$$
c$$$          integer(ikind) :: i,j
c$$$
c$$$          !nodes1 - nodes2
c$$$          if(alignment1(2,1).lt.alignment2(2,1)) then
c$$$             print *, 'y1<y2'
c$$$
c$$$             do j=1, interior_j_max1
c$$$                do i=1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, size(grdpts_id1,2)
c$$$                do i=1, size(grdpts_id1,1)
c$$$                   new_grdpts_id(i,j_min1+j) = grdpts_id1(i,j)
c$$$                end do
c$$$
c$$$                do i=size(grdpts_id1,1)+1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min1+j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, interior_j_max2
c$$$                do i=1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min2+j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, size(grdpts_id2,2)
c$$$                do i=1, size(grdpts_id2,1)
c$$$                   new_grdpts_id(i,j_min3+j) = grdpts_id2(i,j)
c$$$                end do
c$$$
c$$$                do i=size(grdpts_id2,1)+1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min3+j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, interior_j_max3
c$$$                do i=1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min4+j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$          else
c$$$             print *, 'y2<y1'
c$$$
c$$$             do j=1, interior_j_max1
c$$$                do i=1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, size(grdpts_id2,2)
c$$$                do i=1, size(grdpts_id2,1)
c$$$                   new_grdpts_id(i,j_min1+j) = grdpts_id2(i,j)
c$$$                end do
c$$$
c$$$                do i=size(grdpts_id2,1)+1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min1+j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, interior_j_max2
c$$$                do i=1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min2+j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, size(grdpts_id1,2)
c$$$                do i=1, size(grdpts_id1,1)
c$$$                   new_grdpts_id(i,j_min3+j) = grdpts_id1(i,j)
c$$$                end do
c$$$
c$$$                do i=size(grdpts_id1,1)+1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min3+j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, interior_j_max3
c$$$                do i=1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min4+j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$          end if
c$$$
c$$$        end subroutine merge_grdpts_id_E
c$$$
c$$$
c$$$        !> merge the grdpts_id for western buffer layers
c$$$        subroutine merge_grdpts_id_W(
c$$$     $     new_grdpts_id,
c$$$     $     grdpts_id1, grdpts_id2,
c$$$     $     alignment1, alignment2,
c$$$     $     interior_j_max1, interior_j_max2, interior_j_max3,
c$$$     $     j_min1, j_min2, j_min3, j_min4)
c$$$
c$$$          implicit none
c$$$
c$$$          integer       , dimension(:,:), intent(out):: new_grdpts_id
c$$$          integer       , dimension(:,:), intent(in) :: grdpts_id1
c$$$          integer       , dimension(:,:), intent(in) :: grdpts_id2
c$$$          integer(ikind), dimension(2,2), intent(in) :: alignment1
c$$$          integer(ikind), dimension(2,2), intent(in) :: alignment2
c$$$          integer(ikind)                , intent(in) :: interior_j_max1
c$$$          integer(ikind)                , intent(in) :: interior_j_max2
c$$$          integer(ikind)                , intent(in) :: interior_j_max3
c$$$          integer(ikind)                , intent(in) :: j_min1
c$$$          integer(ikind)                , intent(in) :: j_min2
c$$$          integer(ikind)                , intent(in) :: j_min3
c$$$          integer(ikind)                , intent(in) :: j_min4
c$$$
c$$$          integer(ikind) :: i,j, i_min1, i_min2
c$$$
c$$$          !grdpts_id1 - grdpts_id2
c$$$          if(alignment1(2,1).lt.alignment2(2,1)) then
c$$$             !print *, 'y1<y2'
c$$$
c$$$             i_min1 = size(new_grdpts_id,1)-size(grdpts_id1,1)
c$$$             i_min2 = size(new_grdpts_id,1)-size(grdpts_id2,1)
c$$$
c$$$             do j=1, interior_j_max1
c$$$                do i=1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, size(grdpts_id1,2)
c$$$                do i=1, i_min1
c$$$                   new_grdpts_id(i,j_min1+j) = no_pt
c$$$                end do
c$$$
c$$$                do i=i_min1+1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min1+j) = grdpts_id1(i-i_min1,j)
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, interior_j_max2
c$$$                do i=1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min2+j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, size(grdpts_id2,2)
c$$$                do i=1, i_min2
c$$$                   new_grdpts_id(i,j_min3+j) = no_pt
c$$$                end do
c$$$
c$$$                do i=i_min2+1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min3+j) = grdpts_id2(i-i_min2,j)
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, interior_j_max3
c$$$                do i=1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min4+j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$          else
c$$$             !print *, 'y2<y1'
c$$$
c$$$             i_min1 = size(new_grdpts_id,1)-size(grdpts_id2,1)
c$$$             i_min2 = size(new_grdpts_id,1)-size(grdpts_id1,1)
c$$$
c$$$             do j=1, interior_j_max1
c$$$                do i=1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, size(grdpts_id2,2)
c$$$                do i=1, i_min1
c$$$                   new_grdpts_id(i,j_min1+j) = no_pt
c$$$                end do
c$$$
c$$$                do i=i_min1+1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min1+j) = grdpts_id2(i-i_min1,j)
c$$$                end do
c$$$             end do
c$$$             
c$$$             do j=1, interior_j_max2
c$$$                do i=1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min2+j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, size(grdpts_id1,2)
c$$$                do i=1, i_min2
c$$$                   new_grdpts_id(i,j_min3+j) = no_pt
c$$$                end do
c$$$
c$$$                do i=i_min2+1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min3+j) = grdpts_id1(i-i_min2,j)
c$$$                end do
c$$$             end do
c$$$
c$$$             do j=1, interior_j_max3
c$$$                do i=1, size(new_grdpts_id,1)
c$$$                   new_grdpts_id(i,j_min4+j) = no_pt
c$$$                end do
c$$$             end do
c$$$
c$$$          end if
c$$$
c$$$        end subroutine merge_grdpts_id_W


        subroutine add_nodes_blocks_2_to_8_NS(
     $     new_nodes,
     $     nodes_block4, nodes_block6, interior_nodes,
     $     bf_alignment,
     $     k, j_min, j_max,
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
     $               ny-(2*bc_size)+j,
     $               k)
             end do

             do i=1, size(nodes_block4,1)
                new_nodes(i_min3+i,j,k) = nodes_block4(i,j,k)
             end do
             
             do i=1, interior_i_max2
                new_nodes(i_min4+i,j,k) = interior_nodes(
     $               bf_alignment(1,1)-(bc_size+1)+i_min4+i,
     $               ny-(2*bc_size)+j,
     $               k)
             end do

             do i=1, size(nodes_block6,1)
                new_nodes(i_min5+i,j,k) = nodes_block6(i,j,k)
             end do

             do i=1, interior_i_max3
                new_nodes(i_min6+i,j,k) = interior_nodes(
     $               bf_alignment(1,1)-(bc_size+1)+i_min6+i,
     $               ny-(2*bc_size)+j,
     $               k)
             end do
          end do

        end subroutine add_nodes_blocks_2_to_8_NS


        subroutine add_nodes_blocks_11_to_13_NS(
     $     new_nodes,
     $     nodes_block11, nodes_block13,
     $     k, j_min, j_max,
     $     i_min3, i_min5)

          implicit none

          real(rkind), dimension(:,:,:), intent(out):: new_nodes
          real(rkind), dimension(:,:,:), intent(in) :: nodes_block11
          real(rkind), dimension(:,:,:), intent(in) :: nodes_block13
          integer       , intent(in) :: k
          integer(ikind), intent(in) :: j_min
          integer(ikind), intent(in) :: j_max
          integer(ikind), intent(in) :: i_min3
          integer(ikind), intent(in) :: i_min5


          integer(ikind) :: i,j


          do j=j_min, j_max
             do i=1, size(nodes_block11,1)
                new_nodes(i_min3+i,j,k) = nodes_block11(i,j,k)
             end do
             
             do i=1, size(nodes_block13,1)
                new_nodes(i_min5+i,j,k) = nodes_block13(i,j,k)
             end do
          end do

        end subroutine add_nodes_blocks_11_to_13_NS


        subroutine add_nodes_blocks_16_to_18_NS(
     $     new_nodes,
     $     nodes_block16,
     $     k, j_min, j_max,
     $     i_min3)

          implicit none

          real(rkind), dimension(:,:,:), intent(out):: new_nodes
          real(rkind), dimension(:,:,:), intent(in) :: nodes_block16
          integer       , intent(in) :: k
          integer(ikind), intent(in) :: j_min
          integer(ikind), intent(in) :: j_max
          integer(ikind), intent(in) :: i_min3


          integer(ikind) :: i,j


          do j=j_min, j_max
             do i=1, size(nodes_block16,1)
                new_nodes(i_min3+i,j,k) = nodes_block16(i,j,k)
             end do
          end do

        end subroutine add_nodes_blocks_16_to_18_NS


        subroutine add_grdpts_id_blocks_1_to_9_NS(
     $     new_grdpts_id,
     $     grdpts_block4, grdpts_block6,
     $     border_W, border_E, interior_profile,
     $     j_min, j_max,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min1, i_min3, i_min4, i_min5, i_min6, i_min7, i_min8,
     $     i_match_borderE)

          implicit none

          integer, dimension(:,:), intent(out):: new_grdpts_id
          integer, dimension(:,:), intent(in) :: grdpts_block4
          integer, dimension(:,:), intent(in) :: grdpts_block6
          integer, dimension(bc_size,2*bc_size), intent(in):: border_W
          integer, dimension(bc_size,2*bc_size), intent(in):: border_E
          integer, dimension(2*bc_size), intent(in):: interior_profile
          integer(ikind), intent(in) :: j_min
          integer(ikind), intent(in) :: j_max
          integer(ikind), intent(in) :: outside_i_max1
          integer(ikind), intent(in) :: outside_i_max2
          integer(ikind), intent(in) :: interior_i_max1
          integer(ikind), intent(in) :: interior_i_max2
          integer(ikind), intent(in) :: interior_i_max3
          integer(ikind), intent(in) :: i_min1
          integer(ikind), intent(in) :: i_min3
          integer(ikind), intent(in) :: i_min4
          integer(ikind), intent(in) :: i_min5
          integer(ikind), intent(in) :: i_min6
          integer(ikind), intent(in) :: i_min7
          integer(ikind), intent(in) :: i_min8
          integer(ikind), intent(in) :: i_match_borderE


          integer(ikind) :: i,j


          do j=j_min, j_max

             !block 1
             do i=1, outside_i_max1
                new_grdpts_id(i,j) = no_pt
             end do

             !block 2
             do i=1, min(bc_size, interior_i_max1)
                new_grdpts_id(i_min1+i,j) = border_W(i,j)
             end do

             !block 3
             do i=bc_size+1, interior_i_max1
                new_grdpts_id(i_min1+i,j) = interior_profile(j)
             end do

             !block 4
             do i=1, size(grdpts_block4,1)
                new_grdpts_id(i_min3+i,j) = grdpts_block4(i,j)
             end do

             !block 5
             do i=1, interior_i_max2
                new_grdpts_id(i_min4+i,j) = interior_profile(j)
             end do

             !block 6
             do i=1, size(grdpts_block6,1)
                new_grdpts_id(i_min5+i,j) = grdpts_block6(i,j)
             end do

             !block 7
             do i=1, interior_i_max3-bc_size
                new_grdpts_id(i_min6+i,j) = interior_profile(j)
             end do

             !block 8
             do i=1, min(interior_i_max3,bc_size)
                new_grdpts_id(i_min7+i,j) = border_E(
     $               i_match_borderE+i,j)
             end do

             !block 9
             do i=1, outside_i_max2
                new_grdpts_id(i_min8+i,j) = no_pt
             end do
          end do

        end subroutine add_grdpts_id_blocks_1_to_9_NS
        


        subroutine add_grdpts_id_blocks_10_to_14_NS(
     $     new_grdpts_id,
     $     grdpts_block11, grdpts_block13,
     $     j_min, j_max,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min3, i_min4, i_min5, i_min6)

          implicit none

          integer, dimension(:,:), intent(out):: new_grdpts_id
          integer, dimension(:,:), intent(in) :: grdpts_block11
          integer, dimension(:,:), intent(in) :: grdpts_block13
          integer(ikind), intent(in) :: j_min
          integer(ikind), intent(in) :: j_max
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
                new_grdpts_id(i_min3+i,j) = grdpts_block11(i,j)
             end do

             !block 12
             do i=1, interior_i_max2
                new_grdpts_id(i_min4+i,j) = no_pt
             end do

             !block 13
             do i=1, size(grdpts_block13,1)
                new_grdpts_id(i_min5+i,j) = grdpts_block13(i,j)
             end do

             !block 14
             do i=1, interior_i_max3+outside_i_max2
                new_grdpts_id(i_min6+i,j) = no_pt
             end do
          end do

        end subroutine add_grdpts_id_blocks_10_to_14_NS


        subroutine add_grdpts_id_blocks_15_to_19_with_16_taller_NS(
     $     new_grdpts_id,
     $     grdpts_block16, grdpts_block18,
     $     j_min, j_max,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min3, i_min4)


          implicit none

          integer, dimension(:,:), intent(out) :: new_grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_block16
          integer, dimension(:,:), intent(in)  :: grdpts_block18
          integer(ikind), intent(in) :: j_min
          integer(ikind), intent(in) :: j_max
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
                new_grdpts_id(i_min3+i,j) = grdpts_block16(i,j)
             end do
             
             !block 17+18+19
             do i=1, interior_i_max2+size(grdpts_block18,1)+interior_i_max3+outside_i_max2
                new_grdpts_id(i_min4+i,j) = no_pt
             end do
          end do

        end subroutine add_grdpts_id_blocks_15_to_19_with_16_taller_NS


        subroutine add_grdpts_id_blocks_15_to_19_with_18_taller_NS(
     $     new_grdpts_id,
     $     grdpts_block16, grdpts_block18,
     $     j_min, j_max,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min5, i_min6)


          implicit none

          integer, dimension(:,:), intent(out) :: new_grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_block16
          integer, dimension(:,:), intent(in)  :: grdpts_block18
          integer(ikind), intent(in) :: j_min
          integer(ikind), intent(in) :: j_max
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
                new_grdpts_id(i_min5+i,j) = grdpts_block18(i,j)
             end do

             !block 19
             do i=1, interior_i_max3+outside_i_max2
                new_grdpts_id(i_min6+i,j) = no_pt
             end do
          end do

        end subroutine add_grdpts_id_blocks_15_to_19_with_18_taller_NS

      
        subroutine add_grdpts_id_block_20_NS(
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

        end subroutine add_grdpts_id_block_20_NS

      end module bf_layer_merge_module
