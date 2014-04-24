      module bf_layer_merge_module

        use parameters_bf_layer, only : no_pt
        use parameters_constant, only : x_direction, y_direction
        use parameters_input   , only : debug, nx, ny, ne, bc_size
        use parameters_kind    , only : ikind, rkind

        private
        public :: merge_bf_layers_N,
     $            merge_bf_layers_S,
     $            merge_bf_layers_E,
     $            merge_bf_layers_W

        contains

        !< merge northern buffer layers
        subroutine merge_bf_layers_N(
     $       nodes1, nodes2,
     $       grdpts_id1, grdpts_id2,
     $       alignment1, alignment2, final_alignment_i)

          implicit none

          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes1
          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes2
          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id1
          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id2
          integer(ikind), dimension(2,2)               , intent(inout) :: alignment1
          integer(ikind), dimension(2,2)               , intent(in)    :: alignment2
          integer(ikind), dimension(2,2)  , optional   , intent(in)    :: final_alignment_i

          integer(ikind), dimension(2,2)                :: final_alignment
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
          integer(ikind), dimension(2) :: new_size
          integer                      :: interior_i_max1
          integer                      :: interior_i_max2
          integer                      :: interior_i_max3
          integer                      :: i_min1
          integer                      :: i_min2
          integer                      :: i_min3
          integer                      :: i_min4
          integer                      :: j_min1
          integer                      :: j_min2
          integer                      :: j_max


          !get the new size of the tables
          !and the indices to match easily during the copy of the tables
          if(present(final_alignment_i)) then

             final_alignment      = final_alignment_i
             final_alignment(2,1) = ny+1

             new_size = get_new_size(alignment1, alignment2, final_alignment)

             call get_match_indices(
     $            y_direction,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_min1, i_min2, i_min3, i_min4,
     $            j_min1, j_min2,
     $            alignment1, alignment2, final_alignment)

          else

             new_size = get_new_size(alignment1, alignment2)

             call get_match_indices(
     $            y_direction,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_min1, i_min2, i_min3, i_min4,
     $            j_min1, j_min2,
     $            alignment1, alignment2)
          end if
          j_max = new_size(2)

          !allocate the nodes and copy the tables
          allocate(new_nodes(new_size(1), new_size(2), ne))
          call merge_nodes_N(
     $         new_nodes,
     $         nodes1, nodes2,
     $         alignment1, alignment2,
     $         i_min1, i_min3,
     $         j_min1)
          deallocate(nodes2)
          call MOVE_ALLOC(new_nodes,nodes1)


          !allocate the gridpts_id and copy the tables
          allocate(new_grdpts_id(new_size(1), new_size(2)))
          call merge_grdpts_id_N(
     $         new_grdpts_id,
     $         grdpts_id1, grdpts_id2,
     $         alignment1, alignment2,
     $         interior_i_max1, interior_i_max2, interior_i_max3,
     $         i_min1, i_min2, i_min3, i_min4,
     $         j_min1, j_max)
          deallocate(grdpts_id2)
          call MOVE_ALLOC(new_grdpts_id,grdpts_id1)

          !update the alignment
          if(present(final_alignment_i)) then
             alignment1(1,1) = min(alignment1(1,1), alignment2(1,1), final_alignment(1,1))
             alignment1(2,1) = ny+1
             alignment1(1,2) = max(alignment1(1,2), alignment2(1,2), final_alignment(1,2))
             alignment1(2,2) = max(alignment1(2,2), alignment2(2,2), final_alignment(2,2))
          else
             alignment1(1,1) = min(alignment1(1,1), alignment2(1,1))
             alignment1(2,1) = ny+1
             alignment1(1,2) = max(alignment1(1,2), alignment2(1,2))
             alignment1(2,2) = max(alignment1(2,2), alignment2(2,2))
          end if

        end subroutine merge_bf_layers_N


        !< merge southern buffer layers
        subroutine merge_bf_layers_S(
     $       nodes1, nodes2,
     $       grdpts_id1, grdpts_id2,
     $       alignment1, alignment2, final_alignment_i)

          implicit none

          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes1
          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes2
          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id1
          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id2
          integer(ikind), dimension(2,2)               , intent(inout) :: alignment1
          integer(ikind), dimension(2,2)               , intent(in)    :: alignment2
          integer(ikind), dimension(2,2)  , optional   , intent(in)    :: final_alignment_i

          integer(ikind), dimension(2,2)                :: final_alignment
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
          integer(ikind), dimension(2) :: new_size
          integer                      :: interior_i_max1
          integer                      :: interior_i_max2
          integer                      :: interior_i_max3
          integer                      :: i_min1
          integer                      :: i_min2
          integer                      :: i_min3
          integer                      :: i_min4
          integer                      :: j_min1
          integer                      :: j_min2
          integer                      :: j_max


          !get the new size of the tables
          !and the indices to match easily during the copy of the tables
          if(present(final_alignment_i)) then

             final_alignment      = final_alignment_i
             final_alignment(2,2) = 0
             
             new_size = get_new_size(alignment1, alignment2, final_alignment)

             call get_match_indices(
     $            y_direction,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_min1, i_min2, i_min3, i_min4,
     $            j_min1, j_min2,
     $            alignment1, alignment2, final_alignment)

          else

             new_size = get_new_size(alignment1, alignment2)

             call get_match_indices(
     $            y_direction,
     $            interior_i_max1, interior_i_max2, interior_i_max3,
     $            i_min1, i_min2, i_min3, i_min4,
     $            j_min1, j_min2,
     $            alignment1, alignment2)
          end if
          j_max = new_size(2)


          !allocate the nodes and copy the tables
          allocate(new_nodes(new_size(1), new_size(2), ne))
          call merge_nodes_S(
     $         new_nodes,
     $         nodes1, nodes2,
     $         alignment1, alignment2,
     $         i_min1, i_min3,
     $         j_min1, j_min2, j_max)
          deallocate(nodes2)
          call MOVE_ALLOC(new_nodes,nodes1)


          !allocate the gridpts_id and copy the tables
          allocate(new_grdpts_id(new_size(1), new_size(2)))
          call merge_grdpts_id_S(
     $         new_grdpts_id,
     $         grdpts_id1, grdpts_id2,
     $         alignment1, alignment2,
     $         interior_i_max1, interior_i_max2, interior_i_max3,
     $         i_min1, i_min2, i_min3, i_min4,
     $         j_min1, j_min2, j_max)
          deallocate(grdpts_id2)
          call MOVE_ALLOC(new_grdpts_id,grdpts_id1)

          !update the alignment
          if(present(final_alignment_i)) then
             alignment1(1,1) = min(alignment1(1,1), alignment2(1,1), final_alignment(1,1))
             alignment1(2,1) = min(alignment1(2,1), alignment2(2,1), final_alignment(2,1))
             alignment1(1,2) = max(alignment1(1,2), alignment2(1,2), final_alignment(1,2))
             alignment1(2,2) = 0
          else
             alignment1(1,1) = min(alignment1(1,1), alignment2(1,1))
             alignment1(2,1) = min(alignment1(2,1), alignment2(2,1))
             alignment1(1,2) = max(alignment1(1,2), alignment2(1,2))
             alignment1(2,2) = 0
          end if

        end subroutine merge_bf_layers_S


        !< merge eastern buffer layers
        subroutine merge_bf_layers_E(
     $       nodes1, nodes2,
     $       grdpts_id1, grdpts_id2,
     $       alignment1, alignment2, final_alignment_i)

          implicit none

          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes1
          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes2
          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id1
          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id2
          integer(ikind), dimension(2,2)               , intent(inout) :: alignment1
          integer(ikind), dimension(2,2)               , intent(in)    :: alignment2
          integer(ikind), dimension(2,2)  , optional   , intent(in)    :: final_alignment_i


          integer(ikind), dimension(2,2)                :: final_alignment
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
          integer(ikind), dimension(2) :: new_size
          integer                      :: interior_j_max1
          integer                      :: interior_j_max2
          integer                      :: interior_j_max3
          integer                      :: j_min1
          integer                      :: j_min2
          integer                      :: j_min3
          integer                      :: j_min4
          integer                      :: i_min1
          integer                      :: i_min2
          integer                      :: i_max


          !get the new size of the tables
          !and the indices to match easily during the copy of the tables
          if(present(final_alignment_i)) then

             final_alignment      = final_alignment_i
             final_alignment(1,1) = nx+1

             new_size = get_new_size(alignment1, alignment2, final_alignment)

             call get_match_indices(
     $            x_direction,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            j_min1, j_min2, j_min3, j_min4,
     $            i_min1, i_min2,
     $            alignment1, alignment2, final_alignment)

          else

             new_size = get_new_size(alignment1, alignment2)

             call get_match_indices(
     $            x_direction,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            j_min1, j_min2, j_min3, j_min4,
     $            i_min1, i_min2,
     $            alignment1, alignment2)
          end if
          i_max = new_size(2)

          !allocate the nodes and copy the tables
          allocate(new_nodes(new_size(1), new_size(2), ne))
          call merge_nodes_E(
     $         new_nodes,
     $         nodes1, nodes2,
     $         alignment1, alignment2,
     $         j_min1, j_min3)
          deallocate(nodes2)
          call MOVE_ALLOC(new_nodes,nodes1)


          !allocate the gridpts_id and copy the tables
          allocate(new_grdpts_id(new_size(1), new_size(2)))
          call merge_grdpts_id_E(
     $         new_grdpts_id,
     $         grdpts_id1, grdpts_id2,
     $         alignment1, alignment2,
     $         interior_j_max1, interior_j_max2, interior_j_max3,
     $         j_min1, j_min2, j_min3, j_min4)
          deallocate(grdpts_id2)
          call MOVE_ALLOC(new_grdpts_id,grdpts_id1)

          !update the alignment
          if(present(final_alignment_i)) then
             alignment1(1,1) = nx+1
             alignment1(2,1) = min(alignment1(2,1), alignment2(2,1), final_alignment(2,1))
             alignment1(1,2) = max(alignment1(1,2), alignment2(1,2), final_alignment(1,2))
             alignment1(2,2) = max(alignment1(2,2), alignment2(2,2), final_alignment(2,2))
          else
             alignment1(1,1) = nx+1
             alignment1(2,1) = min(alignment1(2,1), alignment2(2,1))
             alignment1(1,2) = max(alignment1(1,2), alignment2(1,2))
             alignment1(2,2) = max(alignment1(2,2), alignment2(2,2))
          end if

        end subroutine merge_bf_layers_E


        !< merge western buffer layers
        subroutine merge_bf_layers_W(
     $       nodes1, nodes2,
     $       grdpts_id1, grdpts_id2,
     $       alignment1, alignment2, final_alignment_i)

          implicit none

          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes1
          real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: nodes2
          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id1
          integer       , dimension(:,:)  , allocatable, intent(inout) :: grdpts_id2
          integer(ikind), dimension(2,2)               , intent(inout) :: alignment1
          integer(ikind), dimension(2,2)               , intent(in)    :: alignment2
          integer(ikind), dimension(2,2)  , optional   , intent(in)    :: final_alignment_i


          integer(ikind), dimension(2,2)                :: final_alignment
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
          integer(ikind), dimension(2) :: new_size
          integer                      :: interior_j_max1
          integer                      :: interior_j_max2
          integer                      :: interior_j_max3
          integer                      :: j_min1
          integer                      :: j_min2
          integer                      :: j_min3
          integer                      :: j_min4
          integer                      :: i_min1
          integer                      :: i_min2


          !get the new size of the tables
          !and the indices to match easily during the copy of the tables
          if(present(final_alignment_i)) then

             final_alignment      = final_alignment_i
             final_alignment(1,2) = 0

             new_size = get_new_size(alignment1, alignment2, final_alignment)

             call get_match_indices(
     $            x_direction,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            j_min1, j_min2, j_min3, j_min4,
     $            i_min1, i_min2,
     $            alignment1, alignment2, final_alignment)

          else

             new_size = get_new_size(alignment1, alignment2)

             call get_match_indices(
     $            x_direction,
     $            interior_j_max1, interior_j_max2, interior_j_max3,
     $            j_min1, j_min2, j_min3, j_min4,
     $            i_min1, i_min2,
     $            alignment1, alignment2)
          end if

          !allocate the nodes and copy the tables
          allocate(new_nodes(new_size(1), new_size(2), ne))
          call merge_nodes_W(
     $         new_nodes,
     $         nodes1, nodes2,
     $         alignment1, alignment2,
     $         j_min1, j_min3)
          deallocate(nodes2)
          call MOVE_ALLOC(new_nodes,nodes1)


          !allocate the gridpts_id and copy the tables
          allocate(new_grdpts_id(new_size(1), new_size(2)))
          call merge_grdpts_id_W(
     $         new_grdpts_id,
     $         grdpts_id1, grdpts_id2,
     $         alignment1, alignment2,
     $         interior_j_max1, interior_j_max2, interior_j_max3,
     $         j_min1, j_min2, j_min3, j_min4)
          deallocate(grdpts_id2)
          call MOVE_ALLOC(new_grdpts_id,grdpts_id1)

          !update the alignment
          if(present(final_alignment_i)) then
             alignment1(1,1) = min(alignment1(1,1), alignment2(1,1), final_alignment(1,1))
             alignment1(2,1) = min(alignment1(2,1), alignment2(2,1), final_alignment(2,1))
             alignment1(1,2) = 0
             alignment1(2,2) = max(alignment1(2,2), alignment2(2,2), final_alignment(2,2))
          else
             alignment1(1,1) = min(alignment1(1,1), alignment2(1,1), final_alignment(1,1))
             alignment1(2,1) = min(alignment1(2,1), alignment2(2,1))
             alignment1(1,2) = 0
             alignment1(2,2) = max(alignment1(2,2), alignment2(2,2))
          end if

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
        subroutine get_match_indices_NS(
     $       interior_i_max1, interior_i_max2, interior_i_max3,
     $       i_min1, i_min2, i_min3, i_min4,
     $       j_min1, j_min2,
     $       alignment1, alignment2, final_alignment)

          implicit none

          integer                                 , intent(out) :: interior_i_max1
          integer                                 , intent(out) :: interior_i_max2
          integer                                 , intent(out) :: interior_i_max3
          integer                                 , intent(out) :: i_min1
          integer                                 , intent(out) :: i_min2
          integer                                 , intent(out) :: i_min3
          integer                                 , intent(out) :: i_min4
          integer                                 , intent(out) :: j_min1
          integer                                 , intent(out) :: j_min2
          integer(ikind), dimension(2,2)          , intent(in)  :: alignment1
          integer(ikind), dimension(2,2)          , intent(in)  :: alignment2
          integer(ikind), dimension(2,2), optional, intent(in)  :: final_alignment

          integer :: size1, size2


          if(present(final_alignment)) then
             interior_i_max1 = max(min(alignment1(1,1), alignment2(1,1)) -
     $                         final_alignment(1,1),0)

             interior_i_max3 = max(final_alignment(1,2) -
     $                         max(alignment1(1,2),alignment2(1,2)),0)
          else
             interior_i_max1 = 0
             interior_i_max3 = 0
          end if             

          interior_i_max2 = max(alignment1(1,1), alignment2(1,1)) -
     $                      min(alignment1(1,2), alignment2(1,2)) -
     $                      (2*bc_size+1)

          if(debug) then
             if(interior_i_max2.lt.0) then
                print '(''bf_layer_merge_module'')'
                print '(''get_match_indices_NS'')'
                print '(''the two tables are superposed'')'
                stop 'check the alignment of the two tables'
             end if
          end if
          
          if(alignment1(1,1).lt.alignment2(1,1)) then
             size1 = alignment1(1,2)-alignment1(1,1) + 2*bc_size + 1
             size2 = alignment2(1,2)-alignment2(1,1) + 2*bc_size + 1
          else
             size1 = alignment2(1,2)-alignment2(1,1) + 2*bc_size + 1
             size2 = alignment1(1,2)-alignment1(1,1) + 2*bc_size + 1
          end if

          i_min1 = interior_i_max1
          i_min2 = i_min1 + size1
          i_min3 = i_min2 + interior_i_max2
          i_min4 = i_min3 + size2

          j_min1 = min((alignment1(2,2)-alignment1(2,1)+2*bc_size+1),
     $                 (alignment2(2,2)-alignment2(2,1)+2*bc_size+1))
          j_min2 = max((alignment1(2,2)-alignment1(2,1)+2*bc_size+1),
     $                 (alignment2(2,2)-alignment2(2,1)+2*bc_size+1))

        end subroutine get_match_indices_NS


        !< get the indices needed to match the tables copied
        !> when merging southern and northerm buffer layers
        subroutine get_match_indices(
     $     direction,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min1, i_min2, i_min3, i_min4,
     $     j_min1, j_min2,
     $     alignment1, alignment2, final_alignment)

          implicit none

          integer                                 , intent(in)  :: direction
          integer                                 , intent(out) :: interior_i_max1
          integer                                 , intent(out) :: interior_i_max2
          integer                                 , intent(out) :: interior_i_max3
          integer                                 , intent(out) :: i_min1
          integer                                 , intent(out) :: i_min2
          integer                                 , intent(out) :: i_min3
          integer                                 , intent(out) :: i_min4
          integer                                 , intent(out) :: j_min1
          integer                                 , intent(out) :: j_min2
          integer(ikind), dimension(2,2)          , intent(in)  :: alignment1
          integer(ikind), dimension(2,2)          , intent(in)  :: alignment2
          integer(ikind), dimension(2,2), optional, intent(in)  :: final_alignment


          integer :: size1, size2
          integer :: dir1, dir2

          select case(direction)
            case(y_direction)
               dir1 = x_direction
               dir2 = y_direction
            case(x_direction)
               dir1 = y_direction
               dir2 = x_direction
            case default
               print '(''bf_layer_merge_module'')'
               print '(''get_match_indices'')'
               print '(''direction not recognized'')'
               print '(''direction: '',I2)', direction
               stop 'either use x_direction or y_direction'
          end select

          if(present(final_alignment)) then
             interior_i_max1 = max(min(alignment1(dir1,1), alignment2(dir1,1)) -
     $                             final_alignment(dir1,1),0)

             interior_i_max3 = max(final_alignment(dir1,2) -
     $                             max(alignment1(dir1,2),alignment2(dir1,2)),0)
          else
             interior_i_max1 = 0
             interior_i_max3 = 0
          end if             

          interior_i_max2 = max(alignment1(dir1,1), alignment2(dir1,1)) -
     $                      min(alignment1(dir1,2), alignment2(dir1,2)) -
     $                      (2*bc_size+1)

          if(debug) then
             if(interior_i_max2.lt.0) then
                print '(''bf_layer_merge_module'')'
                print '(''get_match_indices_NS'')'
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

          i_min1 = interior_i_max1
          i_min2 = i_min1 + size1
          i_min3 = i_min2 + interior_i_max2
          i_min4 = i_min3 + size2

          j_min1 = min((alignment1(dir2,2)-alignment1(dir2,1)+2*bc_size+1),
     $                 (alignment2(dir2,2)-alignment2(dir2,1)+2*bc_size+1))
          j_min2 = max((alignment1(dir2,2)-alignment1(dir2,1)+2*bc_size+1),
     $                 (alignment2(dir2,2)-alignment2(dir2,1)+2*bc_size+1))

        end subroutine get_match_indices

        
        !> merge the nodes for northern buffer layers
        subroutine merge_nodes_N(
     $     new_nodes,
     $     nodes1, nodes2,
     $     alignment1, alignment2,
     $     i_min1, i_min3,
     $     j_min1)

          implicit none

          real(rkind), dimension(:,:,:) , intent(out):: new_nodes
          real(rkind), dimension(:,:,:) , intent(in) :: nodes1
          real(rkind), dimension(:,:,:) , intent(in) :: nodes2
          integer(ikind), dimension(2,2), intent(in) :: alignment1
          integer(ikind), dimension(2,2), intent(in) :: alignment2
          integer(ikind)                , intent(in) :: i_min1
          integer(ikind)                , intent(in) :: i_min3
          integer(ikind)                , intent(in) :: j_min1

          integer(ikind) :: i,j
          integer        :: k

          !nodes1 - nodes2
          if(alignment1(1,1).lt.alignment2(1,1)) then

             if(alignment1(2,2).gt.alignment2(2,2)) then
                do k=1, ne
                   do j=1, j_min1
                      do i=1, size(nodes1,1)
                         new_nodes(i_min1+i,j,k) = nodes1(i,j,k)
                      end do
                      
                      do i=1, size(nodes2,1)
                         new_nodes(i_min3+i,j,k) = nodes2(i,j,k)
                      end do
                   end do

                   do j=j_min1+1, size(nodes1,2)
                      do i=1, size(nodes1,1)
                         new_nodes(i_min1+i,j,k) = nodes1(i,j,k)
                      end do
                   end do
                end do
             else
                do k=1, ne
                   do j=1, j_min1
                      do i=1, size(nodes1,1)
                         new_nodes(i_min1+i,j,k) = nodes1(i,j,k)
                      end do
                      
                      do i=1, size(nodes2,1)
                         new_nodes(i_min3+i,j,k) = nodes2(i,j,k)
                      end do
                   end do

                   do j=j_min1+1, size(nodes1,2)
                      do i=1, size(nodes1,1)
                         new_nodes(i_min1+i,j,k) = nodes1(i,j,k)
                      end do
                   end do
                end do
             end if

          !nodes2 - nodes1
          else
             if(alignment1(2,2).gt.alignment2(2,2)) then
                do k=1, ne
                   do j=1, j_min1
                      do i=1, size(nodes2,1)
                         new_nodes(i_min1+i,j,k) = nodes2(i,j,k)
                      end do
                      
                      do i=1, size(nodes1,1)
                         new_nodes(i_min3+i,j,k) = nodes1(i,j,k)
                      end do
                   end do

                   do j=j_min1+1, size(nodes1,2)
                      do i=1, size(nodes1,1)
                         new_nodes(i_min3+i,j,k) = nodes1(i,j,k)
                      end do
                   end do
                end do
             else
                do k=1, ne
                   do j=1, j_min1
                      do i=1, size(nodes2,1)
                         new_nodes(i_min1+i,j,k) = nodes2(i,j,k)
                      end do
                      
                      do i=1, size(nodes1,1)
                         new_nodes(i_min3+i,j,k) = nodes1(i,j,k)
                      end do
                   end do

                   do j=j_min1+1, size(nodes2,2)
                      do i=1, size(nodes2,1)
                         new_nodes(i_min1+i,j,k) = nodes2(i,j,k)
                      end do
                   end do
                end do
             end if

          end if

        end subroutine merge_nodes_N


        !> merge the nodes for southern buffer layers
        subroutine merge_nodes_S(
     $     new_nodes,
     $     nodes1, nodes2,
     $     alignment1, alignment2,
     $     i_min1, i_min3,
     $     j_min1, j_min2, j_max)

          implicit none

          real(rkind), dimension(:,:,:) , intent(out):: new_nodes
          real(rkind), dimension(:,:,:) , intent(in) :: nodes1
          real(rkind), dimension(:,:,:) , intent(in) :: nodes2
          integer(ikind), dimension(2,2), intent(in) :: alignment1
          integer(ikind), dimension(2,2), intent(in) :: alignment2
          integer(ikind)                , intent(in) :: i_min1
          integer(ikind)                , intent(in) :: i_min3
          integer(ikind)                , intent(in) :: j_min1
          integer(ikind)                , intent(in) :: j_min2
          integer(ikind)                , intent(in) :: j_max

          integer(ikind) :: i,j
          integer        :: k

          !nodes1 - nodes2
          if(alignment1(1,1).lt.alignment2(1,1)) then

             if(alignment1(2,1).lt.alignment2(2,1)) then
                do k=1, ne
                   do j=j_max-j_min2+1, j_max-j_min1
                      do i=1, size(nodes1,1)
                         new_nodes(i_min1+i,j,k) = nodes1(i,j-(j_max-j_min2),k)
                      end do
                   end do

                   do j=j_max-j_min1+1, j_max
                      do i=1, size(nodes1,1)
                         new_nodes(i_min1+i,j,k) = nodes1(i,j-(j_max-j_min2),k)
                      end do
                      
                      do i=1, size(nodes2,1)
                         new_nodes(i_min3+i,j,k) = nodes2(i,j-(j_max-j_min1),k)
                      end do
                   end do
                end do                
             else
                do k=1, ne
                   do j=j_max-j_min2+1, j_max-j_min1
                      do i=1, size(nodes2,1)
                         new_nodes(i_min3+i,j,k) = nodes2(i,j-(j_max-j_min2),k)
                      end do
                   end do

                   do j=j_max-j_min1+1, j_max
                      do i=1, size(nodes1,1)
                         new_nodes(i_min1+i,j,k) = nodes1(i,j-(j_max-j_min1),k)
                      end do
                      
                      do i=1, size(nodes2,1)
                         new_nodes(i_min3+i,j,k) = nodes2(i,j-(j_max-j_min2),k)
                      end do
                   end do
                end do
             end if            

          !nodes2 - nodes1
          else

             if(alignment1(2,1).lt.alignment2(2,1)) then
                !print *, 'nodes:', 'x1>x2 : y1>y2'
                do k=1, ne
                   do j=j_max-j_min2+1, j_max-j_min1
                      do i=1, size(nodes1,1)
                         new_nodes(i_min3+i,j,k) = nodes1(i,j-(j_max-j_min2),k)
                      end do
                   end do

                   do j=j_max-j_min1+1, j_max
                      do i=1, size(nodes2,1)
                         new_nodes(i_min1+i,j,k) = nodes2(i,j-(j_max-j_min1),k)
                      end do

                      do i=1, size(nodes1,1)
                         new_nodes(i_min3+i,j,k) = nodes1(i,j-(j_max-j_min2),k)
                      end do
                   end do
                end do
             else
                do k=1, ne
                   do j=j_max-j_min2, j_max-j_min1
                      do i=1, size(nodes2,1)
                         new_nodes(i_min1+i,j,k) = nodes2(i,j-(j_max-j_min2),k)
                      end do
                   end do

                   do j=j_max-j_min1+1, j_max
                      do i=1, size(nodes2,1)
                         new_nodes(i_min1+i,j,k) = nodes2(i,j-(j_max-j_min2),k)
                      end do

                      do i=1, size(nodes1,1)
                         new_nodes(i_min3+i,j,k) = nodes1(i,j-(j_max-j_min1),k)
                      end do
                   end do
                end do
             end if
          end if

        end subroutine merge_nodes_S


        !> merge the nodes for eastern buffer layers
        subroutine merge_nodes_E(
     $     new_nodes,
     $     nodes1, nodes2,
     $     alignment1, alignment2,
     $     j_min1, j_min3)

          implicit none

          real(rkind), dimension(:,:,:) , intent(out):: new_nodes
          real(rkind), dimension(:,:,:) , intent(in) :: nodes1
          real(rkind), dimension(:,:,:) , intent(in) :: nodes2
          integer(ikind), dimension(2,2), intent(in) :: alignment1
          integer(ikind), dimension(2,2), intent(in) :: alignment2
          integer(ikind)                , intent(in) :: j_min1
          integer(ikind)                , intent(in) :: j_min3

          integer(ikind) :: i,j
          integer        :: k

          !nodes1 - nodes2
          if(alignment1(2,1).lt.alignment2(2,1)) then
             print *, 'y1<y2'

             do k=1, ne

                do j=1, size(nodes1,2)
                   do i=1, size(nodes1,1)
                      new_nodes(i,j_min1+j,k) = nodes1(i,j,k)
                   end do
                end do
                   
                do j=1, size(nodes2,1)
                   do i=1, size(nodes2,1)
                      new_nodes(i,j_min3+j,k) = nodes2(i,j,k)
                   end do
                end do

             end do

          else
             print *, 'y2<y1'

             do k=1, ne

                do j=1, size(nodes2,2)
                   do i=1, size(nodes2,1)
                      new_nodes(i,j_min1+j,k) = nodes2(i,j,k)
                   end do
                end do

                do j=1, size(nodes1,2)
                   do i=1, size(nodes1,1)
                      new_nodes(i,j_min3+j,k) = nodes1(i,j,k)
                   end do
                end do

             end do

          end if

        end subroutine merge_nodes_E


        !> merge the nodes for western buffer layers
        subroutine merge_nodes_W(
     $     new_nodes,
     $     nodes1, nodes2,
     $     alignment1, alignment2,
     $     j_min1, j_min3)

          implicit none

          real(rkind), dimension(:,:,:) , intent(out):: new_nodes
          real(rkind), dimension(:,:,:) , intent(in) :: nodes1
          real(rkind), dimension(:,:,:) , intent(in) :: nodes2
          integer(ikind), dimension(2,2), intent(in) :: alignment1
          integer(ikind), dimension(2,2), intent(in) :: alignment2
          integer(ikind)                , intent(in) :: j_min1
          integer(ikind)                , intent(in) :: j_min3

          integer(ikind) :: i,j, i_min1, i_min2
          integer        :: k

          !nodes1 - nodes2
          if(alignment1(2,1).lt.alignment2(2,1)) then
             !print *, 'y1<y2'

             i_min1 = size(new_nodes,1)-size(nodes1,1)
             i_min2 = size(new_nodes,1)-size(nodes2,1)

             do k=1, ne
                
                do j=1, size(nodes1,2)
                   do i=i_min1+1, size(new_nodes,1)
                      new_nodes(i,j_min1+j,k) = nodes1(i-i_min1,j,k)
                   end do
                end do
                   
                do j=1, size(nodes2,1)
                   do i=i_min2+1, size(new_nodes,1)
                      new_nodes(i,j_min3+j,k) = nodes2(i-i_min2,j,k)
                   end do
                end do

             end do

          else
             !print *, 'y2<y1'

             i_min1 = size(new_nodes,1)-size(nodes2,1)
             i_min2 = size(new_nodes,1)-size(nodes1,1)

             do k=1, ne

                do j=1, size(nodes2,2)
                   do i=i_min1+1, size(new_nodes,1)
                      new_nodes(i,j_min1+j,k) = nodes2(i-i_min1,j,k)
                   end do
                end do

                do j=1, size(nodes1,2)
                   do i=i_min2+1, size(new_nodes,1)
                      new_nodes(i,j_min3+j,k) = nodes1(i-i_min2,j,k)
                   end do
                end do

             end do

          end if

        end subroutine merge_nodes_W


        !> merge the nodes for northern buffer layers
        subroutine merge_grdpts_id_N(
     $     new_grdpts_id,
     $     grdpts_id1, grdpts_id2,
     $     alignment1, alignment2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min1, i_min2, i_min3, i_min4,
     $     j_min1, j_max)

          implicit none

          integer       , dimension(:,:), intent(out):: new_grdpts_id
          integer       , dimension(:,:), intent(in) :: grdpts_id1
          integer       , dimension(:,:), intent(in) :: grdpts_id2
          integer(ikind), dimension(2,2), intent(in) :: alignment1
          integer(ikind), dimension(2,2), intent(in) :: alignment2          
          integer(ikind)                , intent(in) :: interior_i_max1
          integer(ikind)                , intent(in) :: interior_i_max2
          integer(ikind)                , intent(in) :: interior_i_max3
          integer(ikind)                , intent(in) :: i_min1
          integer(ikind)                , intent(in) :: i_min2
          integer(ikind)                , intent(in) :: i_min3
          integer(ikind)                , intent(in) :: i_min4
          integer(ikind)                , intent(in) :: j_min1
          integer(ikind)                , intent(in) :: j_max

          integer(ikind) :: i,j

          !nodes1 - nodes2
          if(alignment1(1,1).lt.alignment2(1,1)) then

             do j=1, j_min1
                do i=1, interior_i_max1
                   new_grdpts_id(i,j) = no_pt
                end do

                do i=1, size(grdpts_id1,1)
                   new_grdpts_id(i_min1+i,j) = grdpts_id1(i,j)
                end do

                do i=1, interior_i_max2
                   new_grdpts_id(i_min2+i,j) = no_pt
                end do

                do i=1, size(grdpts_id2,1)
                   new_grdpts_id(i_min3+i,j) = grdpts_id2(i,j)
                end do

                do i=1, interior_i_max3
                   new_grdpts_id(i_min4+i,j) = no_pt
                end do
             end do

             if(alignment1(2,2).gt.alignment2(2,2)) then
                do j=j_min1+1, size(grdpts_id1,2)
                   do i=1, interior_i_max1
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=1, size(grdpts_id1,1)
                      new_grdpts_id(i_min1+i,j) = grdpts_id1(i,j)
                   end do

                   do i=1, interior_i_max2+size(grdpts_id2,1)+interior_i_max3
                      new_grdpts_id(i_min2+i,j) = no_pt
                   end do
                end do

                do j=size(grdpts_id1,2)+1, j_max
                   do i=1, size(new_grdpts_id,1)
                      new_grdpts_id(i,j)=no_pt
                   end do
                end do
             else
                do j=j_min1+1, size(grdpts_id2,2)
                   do i=1, interior_i_max1+size(grdpts_id1,1)+interior_i_max2
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=1, size(grdpts_id2,1)
                      new_grdpts_id(i_min3+i,j) = grdpts_id2(i,j)
                   end do

                   do i=1, interior_i_max3
                      new_grdpts_id(i_min4+i,j) = no_pt
                   end do

                end do

                do j=size(grdpts_id2,2)+1, j_max
                   do i=1, size(new_grdpts_id,1)
                      new_grdpts_id(i,j)=no_pt
                   end do
                end do
             end if

          !nodes2 - nodes1
          else

             do j=1, j_min1
                do i=1, interior_i_max1
                   new_grdpts_id(i,j) = no_pt
                end do

                do i=1, size(grdpts_id2,1)
                   new_grdpts_id(i_min1+i,j) = grdpts_id2(i,j)
                end do

                do i=1, interior_i_max2
                   new_grdpts_id(i_min2+i,j) = no_pt
                end do

                do i=1, size(grdpts_id1,1)
                   new_grdpts_id(i_min3+i,j) = grdpts_id1(i,j)
                end do

                do i=1, interior_i_max3
                   new_grdpts_id(i_min4+i,j) = no_pt
                end do                
                
             end do

             if(alignment1(2,2).gt.alignment2(2,2)) then
                do j=j_min1+1, size(grdpts_id1,2)
                   do i=1, interior_i_max1+size(grdpts_id2,1)+interior_i_max2
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=1, size(grdpts_id1,1)
                      new_grdpts_id(i_min3+i,j) = grdpts_id1(i,j)
                   end do
                   
                   do i=1, interior_i_max3
                      new_grdpts_id(i_min4+i,j) = no_pt
                   end do
                end do

                do j=size(grdpts_id1,2)+1, j_max
                   do i=1, size(new_grdpts_id,1)
                      new_grdpts_id(i,j)=no_pt
                   end do
                end do

             else
                do j=j_min1+1, size(grdpts_id2,2)
                   do i=1, interior_i_max1
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=1, size(grdpts_id2,1)
                      new_grdpts_id(i_min1+i,j) = grdpts_id2(i,j)
                   end do

                   do i=1, interior_i_max2+size(grdpts_id1,1)+interior_i_max3
                      new_grdpts_id(i_min2+i,j) = no_pt
                   end do
                end do

                do j=size(grdpts_id2,2)+1, j_max
                   do i=1, size(new_grdpts_id,1)
                      new_grdpts_id(i,j)=no_pt
                   end do
                end do
             end if

          end if

        end subroutine merge_grdpts_id_N

        
        !> merge the nodes for southern buffer layers
        subroutine merge_grdpts_id_S(
     $     new_grdpts_id,
     $     grdpts_id1, grdpts_id2,
     $     alignment1, alignment2,
     $     interior_i_max1, interior_i_max2, interior_i_max3,
     $     i_min1, i_min2, i_min3, i_min4,
     $     j_min1, j_min2, j_max)

          implicit none

          integer       , dimension(:,:), intent(out):: new_grdpts_id
          integer       , dimension(:,:), intent(in) :: grdpts_id1
          integer       , dimension(:,:), intent(in) :: grdpts_id2
          integer(ikind), dimension(2,2), intent(in) :: alignment1
          integer(ikind), dimension(2,2), intent(in) :: alignment2
          integer(ikind)                , intent(in) :: interior_i_max1
          integer(ikind)                , intent(in) :: interior_i_max2
          integer(ikind)                , intent(in) :: interior_i_max3
          integer(ikind)                , intent(in) :: i_min1
          integer(ikind)                , intent(in) :: i_min2
          integer(ikind)                , intent(in) :: i_min3
          integer(ikind)                , intent(in) :: i_min4
          integer(ikind)                , intent(in) :: j_min1
          integer(ikind)                , intent(in) :: j_min2
          integer(ikind)                , intent(in) :: j_max

          integer(ikind) :: i,j

          !grdpts_id1 - grdpts_id2
          if(alignment1(1,1).lt.alignment2(1,1)) then

             if(alignment1(2,1).lt.alignment2(2,1)) then

                print *, 'x1<x2 : y1>y2'
                do j=1, j_max-j_min2
                   do i=1, size(new_grdpts_id,1)
                      new_grdpts_id(i,j) = no_pt
                   end do
                end do

                do j=j_max-j_min2+1, j_max-j_min1
                   do i=1, interior_i_max1
                      new_grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=1, size(grdpts_id1,1)
                      new_grdpts_id(i_min1+i,j) = grdpts_id1(i,j-(j_max-j_min2))
                   end do
                   
                   do i=1, interior_i_max2+size(grdpts_id2,1)+interior_i_max3
                      new_grdpts_id(i_min2+i,j) = no_pt
                   end do
                end do
                
                do j=j_max-j_min1+1, j_max
                   do i=1, interior_i_max1
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=1, size(grdpts_id1,1)
                      new_grdpts_id(i_min1+i,j) = grdpts_id1(i,j-(j_max-j_min2))
                   end do

                   do i=1, interior_i_max2
                      new_grdpts_id(i_min2+i,j) = no_pt
                   end do
                   
                   do i=1, size(grdpts_id2,1)
                      new_grdpts_id(i_min3+i,j) = grdpts_id2(i,j-(j_max-j_min1))
                   end do

                   do i=1, interior_i_max3
                      new_grdpts_id(i_min4+i,j) = no_pt
                   end do
                end do
             else
                print *, 'x1<x2 : y1<y2'
                do j=1, j_max-j_min2
                   do i=1, size(new_grdpts_id,1)
                      new_grdpts_id(i,j) = no_pt
                   end do
                end do

                do j=j_max-j_min2+1, j_max-j_min1
                   do i=1, interior_i_max1+size(grdpts_id1,1)+interior_i_max2
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=1, size(grdpts_id2,1)
                      new_grdpts_id(i_min3+i,j) = grdpts_id2(i,j-(j_max-j_min2))
                   end do

                   do i=1, interior_i_max3
                      new_grdpts_id(i_min4+i,j) = no_pt
                   end do
                end do

                do j=j_max-j_min1+1, j_max
                   do i=1, interior_i_max1
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=1, size(grdpts_id1,1)
                      new_grdpts_id(i_min1+i,j) = grdpts_id1(i,j-(j_max-j_min1))
                   end do

                   do i=1, interior_i_max2
                      new_grdpts_id(i_min2+i,j) = no_pt
                   end do
                   
                   do i=1, size(grdpts_id2,1)
                      new_grdpts_id(i_min3+i,j) = grdpts_id2(i,j-(j_max-j_min2))
                   end do

                   do i=1, interior_i_max3
                      new_grdpts_id(i_min4+i,j) = no_pt
                   end do
                end do
             end if            

          !grdpts_id2 - grdpts_id1
          else

             if(alignment1(2,1).lt.alignment2(2,1)) then
                print *, 'x1>x2 : y1>y2'
                do j=1, j_max-j_min2
                   do i=1, size(new_grdpts_id,1)
                      new_grdpts_id(i,j) = no_pt
                   end do
                end do

                do j=j_max-j_min2+1, j_max-j_min1
                   do i=1, interior_i_max1+size(grdpts_id2,1)+interior_i_max2
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=1, size(grdpts_id1,1)
                      new_grdpts_id(i_min3+i,j) = grdpts_id1(i,j-(j_max-j_min2))
                   end do

                   do i=1, interior_i_max3
                      new_grdpts_id(i_min4+i,j) = no_pt
                   end do

                end do

                do j=j_max-j_min1+1, j_max
                   do i=1, interior_i_max1
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=1, size(grdpts_id2,1)
                      new_grdpts_id(i_min1+i,j) = grdpts_id2(i,j-(j_max-j_min1))
                   end do

                   do i=1, interior_i_max2
                      new_grdpts_id(i_min2+i,j) = no_pt
                   end do

                   do i=1, size(grdpts_id1,1)
                      new_grdpts_id(i_min3+i,j) = grdpts_id1(i,j-(j_max-j_min2))
                   end do

                   do i=1, interior_i_max3
                      new_grdpts_id(i_min4+i,j) = no_pt
                   end do

                end do

             else

                print *, 'x1>x2 : y1<y2'
                do j=1, j_max-j_min2
                   do i=1, size(new_grdpts_id,1)
                      new_grdpts_id(i,j) = no_pt
                   end do
                end do

                do j=j_max-j_min2+1, j_max-j_min1
                   do i=1, interior_i_max1
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=1, size(grdpts_id2,1)
                      new_grdpts_id(i_min1+i,j) = grdpts_id2(i,j-(j_max-j_min2))
                   end do

                   do i=1, interior_i_max2+size(grdpts_id1,1)+interior_i_max3
                      new_grdpts_id(i_min2+i,j) = no_pt
                   end do
                end do

                do j=j_max-j_min1+1, j_max
                   do i=1, interior_i_max1
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=1, size(grdpts_id2,1)
                      new_grdpts_id(i_min1+i,j) = grdpts_id2(i,j-(j_max-j_min2))
                   end do

                   do i=1, interior_i_max2
                      new_grdpts_id(i_min2+i,j) = no_pt
                   end do
                   
                   do i=1, size(grdpts_id1,1)
                      new_grdpts_id(i_min3+i,j) = grdpts_id1(i,j-(j_max-j_min1))
                   end do

                   do i=1, interior_i_max3
                      new_grdpts_id(i_min4+i,j) = no_pt
                   end do
                end do
             end if

          end if

        end subroutine merge_grdpts_id_S

        !> merge the nodes for eastern buffer layers
        subroutine merge_grdpts_id_E(
     $     new_grdpts_id,
     $     grdpts_id1, grdpts_id2,
     $     alignment1, alignment2,
     $     interior_j_max1, interior_j_max2, interior_j_max3,
     $     j_min1, j_min2, j_min3, j_min4)

          implicit none

          integer       , dimension(:,:), intent(out):: new_grdpts_id
          integer       , dimension(:,:), intent(in) :: grdpts_id1
          integer       , dimension(:,:), intent(in) :: grdpts_id2
          integer(ikind), dimension(2,2), intent(in) :: alignment1
          integer(ikind), dimension(2,2), intent(in) :: alignment2
          integer(ikind)                , intent(in) :: interior_j_max1
          integer(ikind)                , intent(in) :: interior_j_max2
          integer(ikind)                , intent(in) :: interior_j_max3
          integer(ikind)                , intent(in) :: j_min1
          integer(ikind)                , intent(in) :: j_min2
          integer(ikind)                , intent(in) :: j_min3
          integer(ikind)                , intent(in) :: j_min4

          integer(ikind) :: i,j

          !nodes1 - nodes2
          if(alignment1(2,1).lt.alignment2(2,1)) then
             print *, 'y1<y2'

             do j=1, interior_j_max1
                do i=1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j) = no_pt
                end do
             end do

             do j=1, size(grdpts_id1,2)
                do i=1, size(grdpts_id1,1)
                   new_grdpts_id(i,j_min1+j) = grdpts_id1(i,j)
                end do

                do i=size(grdpts_id1,1)+1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min1+j) = no_pt
                end do
             end do

             do j=1, interior_j_max2
                do i=1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min2+j) = no_pt
                end do
             end do

             do j=1, size(grdpts_id2,2)
                do i=1, size(grdpts_id2,1)
                   new_grdpts_id(i,j_min3+j) = grdpts_id2(i,j)
                end do

                do i=size(grdpts_id2,1)+1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min3+j) = no_pt
                end do
             end do

             do j=1, interior_j_max3
                do i=1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min4+j) = no_pt
                end do
             end do

          else
             print *, 'y2<y1'

             do j=1, interior_j_max1
                do i=1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j) = no_pt
                end do
             end do

             do j=1, size(grdpts_id2,2)
                do i=1, size(grdpts_id2,1)
                   new_grdpts_id(i,j_min1+j) = grdpts_id2(i,j)
                end do

                do i=size(grdpts_id2,1)+1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min1+j) = no_pt
                end do
             end do

             do j=1, interior_j_max2
                do i=1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min2+j) = no_pt
                end do
             end do

             do j=1, size(grdpts_id1,2)
                do i=1, size(grdpts_id1,1)
                   new_grdpts_id(i,j_min3+j) = grdpts_id1(i,j)
                end do

                do i=size(grdpts_id1,1)+1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min3+j) = no_pt
                end do
             end do

             do j=1, interior_j_max3
                do i=1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min4+j) = no_pt
                end do
             end do

          end if

        end subroutine merge_grdpts_id_E


        !> merge the grdpts_id for western buffer layers
        subroutine merge_grdpts_id_W(
     $     new_grdpts_id,
     $     grdpts_id1, grdpts_id2,
     $     alignment1, alignment2,
     $     interior_j_max1, interior_j_max2, interior_j_max3,
     $     j_min1, j_min2, j_min3, j_min4)

          implicit none

          integer       , dimension(:,:), intent(out):: new_grdpts_id
          integer       , dimension(:,:), intent(in) :: grdpts_id1
          integer       , dimension(:,:), intent(in) :: grdpts_id2
          integer(ikind), dimension(2,2), intent(in) :: alignment1
          integer(ikind), dimension(2,2), intent(in) :: alignment2
          integer(ikind)                , intent(in) :: interior_j_max1
          integer(ikind)                , intent(in) :: interior_j_max2
          integer(ikind)                , intent(in) :: interior_j_max3
          integer(ikind)                , intent(in) :: j_min1
          integer(ikind)                , intent(in) :: j_min2
          integer(ikind)                , intent(in) :: j_min3
          integer(ikind)                , intent(in) :: j_min4

          integer(ikind) :: i,j, i_min1, i_min2

          !grdpts_id1 - grdpts_id2
          if(alignment1(2,1).lt.alignment2(2,1)) then
             !print *, 'y1<y2'

             i_min1 = size(new_grdpts_id,1)-size(grdpts_id1,1)
             i_min2 = size(new_grdpts_id,1)-size(grdpts_id2,1)

             do j=1, interior_j_max1
                do i=1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j) = no_pt
                end do
             end do

             do j=1, size(grdpts_id1,2)
                do i=1, i_min1
                   new_grdpts_id(i,j_min1+j) = no_pt
                end do

                do i=i_min1+1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min1+j) = grdpts_id1(i-i_min1,j)
                end do
             end do

             do j=1, interior_j_max2
                do i=1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min2+j) = no_pt
                end do
             end do

             do j=1, size(grdpts_id2,2)
                do i=1, i_min2
                   new_grdpts_id(i,j_min3+j) = no_pt
                end do

                do i=i_min2+1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min3+j) = grdpts_id2(i-i_min2,j)
                end do
             end do

             do j=1, interior_j_max3
                do i=1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min4+j) = no_pt
                end do
             end do

          else
             !print *, 'y2<y1'

             i_min1 = size(new_grdpts_id,1)-size(grdpts_id2,1)
             i_min2 = size(new_grdpts_id,1)-size(grdpts_id1,1)

             do j=1, interior_j_max1
                do i=1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j) = no_pt
                end do
             end do

             do j=1, size(grdpts_id2,2)
                do i=1, i_min1
                   new_grdpts_id(i,j_min1+j) = no_pt
                end do

                do i=i_min1+1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min1+j) = grdpts_id2(i-i_min1,j)
                end do
             end do
             
             do j=1, interior_j_max2
                do i=1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min2+j) = no_pt
                end do
             end do

             do j=1, size(grdpts_id1,2)
                do i=1, i_min2
                   new_grdpts_id(i,j_min3+j) = no_pt
                end do

                do i=i_min2+1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min3+j) = grdpts_id1(i-i_min2,j)
                end do
             end do

             do j=1, interior_j_max3
                do i=1, size(new_grdpts_id,1)
                   new_grdpts_id(i,j_min4+j) = no_pt
                end do
             end do

          end if

        end subroutine merge_grdpts_id_W

      end module bf_layer_merge_module
