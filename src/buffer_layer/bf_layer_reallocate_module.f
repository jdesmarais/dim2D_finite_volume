      module bf_layer_reallocate_module

        use bf_layer_allocate_module, only : get_additional_blocks_N,
     $                                       get_additional_blocks_S,
     $                                       get_additional_blocks_E,
     $                                       get_additional_blocks_W

        use parameters_bf_layer     , only : interior_pt, bc_interior_pt,
     $                                       bc_pt, no_pt
        use parameters_constant     , only : x_direction, y_direction
        use parameters_input        , only : nx,ny,ne,bc_size
        use parameters_kind         , only : ikind, rkind


        implicit none


        private
        public :: reallocate_bf_layer_N,
     $            reallocate_bf_layer_S,
     $            reallocate_bf_layer_E,
     $            reallocate_bf_layer_W


        contains

        
        subroutine reallocate_bf_layer_N(
     $       bf_nodes, interior_nodes,
     $       bf_grdpts_id,
     $       bf_alignment, final_alignment)

          implicit none
          
          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer       , dimension(:,:)     , allocatable, intent(inout) :: bf_grdpts_id
          integer(ikind), dimension(2,2)                  , intent(inout) :: bf_alignment
          integer(ikind), dimension(2,2)                  , intent(in)    :: final_alignment


          integer(ikind), dimension(2,2)                :: border_changes
          integer(ikind)                                :: i_match, i_max, j_max
          integer(ikind)                                :: i_min1, i_min3, i_min4, i_min5, i_min6
          integer(ikind)                                :: outside_i_max1, outside_i_max2
          integer(ikind)                                :: interior_i_max1, interior_i_max2
          integer(ikind)                                :: i_match_border_E
          integer(ikind), dimension(2)                  :: new_sizes


          call get_border_changes_and_new_alignment_N(
     $         border_changes, final_alignment, bf_alignment)

          call get_match(
     $         x_direction,
     $         bf_alignment, border_changes,
     $         size(bf_nodes,1), size(bf_nodes,2),
     $         i_min1, i_min3, i_min4, i_min5, i_min6,
     $         i_match, i_max, j_max,
     $         outside_i_max1, outside_i_max2,
     $         interior_i_max1, interior_i_max2,
     $         i_match_border_E)

          new_sizes = get_new_sizes(border_changes, size(bf_nodes,1), size(bf_nodes,2))

          call reallocate_nodes_N(
     $         bf_nodes, interior_nodes,
     $         i_min1, i_min3, i_min4,
     $         i_match, i_max, j_max,
     $         interior_i_max1, interior_i_max2,
     $         bf_alignment, new_sizes)

          call reallocate_grdpts_id_N(
     $         bf_grdpts_id,
     $         i_min1, i_min3, i_min4, i_min5, i_min6,
     $         i_match, i_max, j_max,
     $         outside_i_max1, outside_i_max2,
     $         interior_i_max1, interior_i_max2,
     $         bf_alignment ,new_sizes,
     $         i_match_border_E)

        end subroutine reallocate_bf_layer_N


        subroutine reallocate_bf_layer_S(
     $       bf_nodes, interior_nodes,
     $       bf_grdpts_id,
     $       bf_alignment, final_alignment)

          implicit none
          
          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer       , dimension(:,:)     , allocatable, intent(inout) :: bf_grdpts_id
          integer(ikind), dimension(2,2)                  , intent(inout) :: bf_alignment
          integer(ikind), dimension(2,2)                  , intent(in)    :: final_alignment


          integer(ikind), dimension(2,2)                :: border_changes
          integer(ikind)                                :: i_match, i_max, j_max
          integer(ikind)                                :: i_min1, i_min3, i_min4, i_min5, i_min6
          integer(ikind)                                :: outside_i_max1, outside_i_max2
          integer(ikind)                                :: interior_i_max1, interior_i_max2
          integer(ikind)                                :: i_match_border_E
          integer(ikind), dimension(2)                  :: new_sizes


          call get_border_changes_and_new_alignment_S(
     $         border_changes, final_alignment, bf_alignment)

          call get_match(
     $         x_direction,
     $         bf_alignment, border_changes,
     $         size(bf_nodes,1), size(bf_nodes,2),
     $         i_min1, i_min3, i_min4, i_min5, i_min6,
     $         i_match, i_max, j_max,
     $         outside_i_max1, outside_i_max2,
     $         interior_i_max1, interior_i_max2,
     $         i_match_border_E)

          new_sizes = get_new_sizes(border_changes, size(bf_nodes,1), size(bf_nodes,2))

          call reallocate_nodes_S(
     $         bf_nodes, interior_nodes,
     $         i_min1, i_min3, i_min4,
     $         i_match, i_max,
     $         interior_i_max1, interior_i_max2,
     $         bf_alignment, new_sizes)

          call reallocate_grdpts_id_S(
     $         bf_grdpts_id,
     $         i_min1, i_min3, i_min4, i_min5, i_min6,
     $         i_match, i_max, j_max,
     $         outside_i_max1, outside_i_max2,
     $         interior_i_max1, interior_i_max2,
     $         bf_alignment ,new_sizes,
     $         i_match_border_E)

        end subroutine reallocate_bf_layer_S


        subroutine reallocate_bf_layer_E(
     $       bf_nodes, interior_nodes,
     $       bf_grdpts_id,
     $       bf_alignment, final_alignment)

          implicit none
          
          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer       , dimension(:,:)     , allocatable, intent(inout) :: bf_grdpts_id
          integer(ikind), dimension(2,2)                  , intent(inout) :: bf_alignment
          integer(ikind), dimension(2,2)                  , intent(in)    :: final_alignment


          integer(ikind), dimension(2,2) :: border_changes
          integer(ikind)                 :: j_match, i_max, j_max
          integer(ikind)                 :: j_min1, j_min3, j_min4, j_min5, j_min6
          integer(ikind)                 :: outside_j_max1, outside_j_max2
          integer(ikind)                 :: interior_j_max1, interior_j_max2
          integer(ikind)                 :: j_match_border_N
          integer(ikind), dimension(2)   :: new_sizes


          call get_border_changes_and_new_alignment_E(
     $         border_changes, final_alignment, bf_alignment)

          call get_match(
     $         y_direction,
     $         bf_alignment, border_changes,
     $         size(bf_nodes,1), size(bf_nodes,2),
     $         j_min1, j_min3, j_min4, j_min5, j_min6,
     $         j_match, i_max, j_max,
     $         outside_j_max1, outside_j_max2,
     $         interior_j_max1, interior_j_max2,
     $         j_match_border_N)

          new_sizes = get_new_sizes(border_changes, size(bf_nodes,1), size(bf_nodes,2))

          call reallocate_nodes_E(
     $         bf_nodes, interior_nodes,
     $         j_min1, j_min3, j_min4,
     $         j_match, i_max, j_max,
     $         interior_j_max1, interior_j_max2,
     $         bf_alignment, new_sizes)

          call reallocate_grdpts_id_E(
     $         bf_grdpts_id,
     $         j_min1, j_min3, j_min4, j_min5, j_min6,
     $         j_match, i_max, j_max,
     $         outside_j_max1, outside_j_max2,
     $         interior_j_max1, interior_j_max2,
     $         bf_alignment ,new_sizes,
     $         j_match_border_N)

        end subroutine reallocate_bf_layer_E


        !> reallocate the western buffer layer
        subroutine reallocate_bf_layer_W(
     $       bf_nodes, interior_nodes,
     $       bf_grdpts_id,
     $       bf_alignment, final_alignment)

          implicit none
          
          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer       , dimension(:,:)     , allocatable, intent(inout) :: bf_grdpts_id
          integer(ikind), dimension(2,2)                  , intent(inout) :: bf_alignment
          integer(ikind), dimension(2,2)                  , intent(in)    :: final_alignment


          integer(ikind), dimension(2,2) :: border_changes
          integer(ikind)                 :: j_match, i_max, j_max
          integer(ikind)                 :: j_min1, j_min3, j_min4, j_min5, j_min6
          integer(ikind)                 :: outside_j_max1, outside_j_max2
          integer(ikind)                 :: interior_j_max1, interior_j_max2
          integer(ikind)                 :: j_match_border_N
          integer(ikind), dimension(2)   :: new_sizes


          call get_border_changes_and_new_alignment_W(
     $         border_changes, final_alignment, bf_alignment)

          call get_match(
     $         y_direction,
     $         bf_alignment, border_changes,
     $         size(bf_nodes,1), size(bf_nodes,2),
     $         j_min1, j_min3, j_min4, j_min5, j_min6,
     $         j_match, i_max, j_max,
     $         outside_j_max1, outside_j_max2,
     $         interior_j_max1, interior_j_max2,
     $         j_match_border_N)

          new_sizes = get_new_sizes(border_changes, size(bf_nodes,1), size(bf_nodes,2))

          call reallocate_nodes_W(
     $         bf_nodes, interior_nodes,
     $         j_min1, j_min3, j_min4,
     $         j_match, j_max,
     $         interior_j_max1, interior_j_max2,
     $         bf_alignment, new_sizes)

          call reallocate_grdpts_id_W(
     $         bf_grdpts_id,
     $         j_min1, j_min3, j_min4, j_min5, j_min6,
     $         j_match, j_max,
     $         outside_j_max1, outside_j_max2,
     $         interior_j_max1, interior_j_max2,
     $         bf_alignment ,new_sizes,
     $         j_match_border_N)

        end subroutine reallocate_bf_layer_W


        !< reallocate nodes for northern buffer layers
        subroutine reallocate_nodes_N(
     $     bf_nodes, interior_nodes,
     $     i_min1, i_min3, i_min4,
     $     i_match, i_max, j_max,
     $     interior_i_max1, interior_i_max2,
     $     bf_alignment, new_sizes)

          implicit none

          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer(ikind)                                  , intent(in)    :: i_min1
          integer(ikind)                                  , intent(in)    :: i_min3
          integer(ikind)                                  , intent(in)    :: i_min4
          integer(ikind)                                  , intent(in)    :: i_match
          integer(ikind)                                  , intent(in)    :: i_max
          integer(ikind)                                  , intent(in)    :: j_max
          integer(ikind)                                  , intent(in)    :: interior_i_max1
          integer(ikind)                                  , intent(in)    :: interior_i_max2
          integer(ikind), dimension(2,2)                  , intent(in)    :: bf_alignment
          integer(ikind), dimension(2)                    , intent(in)    :: new_sizes

          real(rkind), dimension(:,:,:), allocatable :: new_nodes
          integer(ikind) :: i,j
          integer        :: k

          allocate(new_nodes(new_sizes(1),new_sizes(2),ne))

          do k=1, ne
             do j=1, 2*bc_size
                do i=1, interior_i_max1
                   new_nodes(i_min1+i,j,k) = interior_nodes(
     $                  bf_alignment(1,1)-(bc_size+1)+i_min1+i,
     $                  ny-(2*bc_size)+j,
     $                  k)
                end do

                do i=1, i_max
                   new_nodes(i_min3+i,j,k) = bf_nodes(i_match+i,j,k)
                end do

                do i=1, interior_i_max2
                   new_nodes(i_min4+i,j,k) = interior_nodes(
     $                  bf_alignment(1,1)-(bc_size+1)+i_min4+i,
     $                  ny-(2*bc_size)+j,
     $                  k)
                end do
             end do
          
             do j=2*bc_size+1, j_max
                 do i=1, i_max
                   new_nodes(i_min3+i,j,k) = bf_nodes(i_match+i,j,k)
                end do
             end do
          end do

          call MOVE_ALLOC(new_nodes,bf_nodes)

        end subroutine reallocate_nodes_N


        !< reallocate gridpts_id for northern buffer layers
        subroutine reallocate_grdpts_id_N(
     $     bf_grdpts_id,
     $     i_min1, i_min3, i_min4, i_min5, i_min6,
     $     i_match, i_max, j_max,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2,
     $     bf_alignment, new_sizes,
     $     i_match_border_E)

          implicit none

          integer       , dimension(:,:)     , allocatable, intent(inout) :: bf_grdpts_id
          integer(ikind)                                  , intent(in)    :: i_min1
          integer(ikind)                                  , intent(in)    :: i_min3
          integer(ikind)                                  , intent(in)    :: i_min4
          integer(ikind)                                  , intent(in)    :: i_min5
          integer(ikind)                                  , intent(in)    :: i_min6
          integer(ikind)                                  , intent(in)    :: i_match
          integer(ikind)                                  , intent(in)    :: i_max
          integer(ikind)                                  , intent(in)    :: j_max
          integer(ikind)                                  , intent(in)    :: outside_i_max1
          integer(ikind)                                  , intent(in)    :: outside_i_max2
          integer(ikind)                                  , intent(in)    :: interior_i_max1
          integer(ikind)                                  , intent(in)    :: interior_i_max2
          integer(ikind), dimension(2,2)                  , intent(in)    :: bf_alignment
          integer(ikind), dimension(2)                    , intent(in)    :: new_sizes
          integer(ikind)                                  , intent(in)    :: i_match_border_E


          integer, dimension(:,:), allocatable :: new_grdpts_id
          
          integer(ikind) :: i,j

          integer, dimension(bc_size, 2*bc_size) :: border_W
          integer, dimension(bc_size, 2*bc_size) :: border_E
          integer, dimension(2*bc_size)          :: interior_profile

          allocate(new_grdpts_id(new_sizes(1),new_sizes(2)))
          
          !get the special border for the left and right
          !as well as the interior profile
          call get_additional_blocks_N(
     $         bf_alignment,
     $         border_W, border_E, interior_profile)


          !fill the blocks according to Fig.1
          do j=1, 2*bc_size
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
             do i=1, i_max
                new_grdpts_id(i_min3+i,j) = bf_grdpts_id(i_match+i,j)
             end do

             !block 5
             do i=1, interior_i_max2-bc_size
                new_grdpts_id(i_min4+i,j) = interior_profile(j)
             end do

             !block 6
             do i=1, min(interior_i_max2,bc_size)
                new_grdpts_id(i_min5+i,j) = border_E(i_match_border_E+i,j)
             end do

             !block 7
             do i=1, outside_i_max2
                new_grdpts_id(i_min6+i,j) = no_pt
             end do
          end do
          
          do j=2*bc_size+1, j_max
             !block 8
             do i=1, outside_i_max1+interior_i_max1
                new_grdpts_id(i,j) = no_pt
             end do

             !block 9
             do i=1, i_max
                new_grdpts_id(i_min3+i,j) = bf_grdpts_id(i_match+i,j)
             end do

             !block 10
             do i=1, interior_i_max2+outside_i_max2
                new_grdpts_id(i_min4+i,j) = no_pt
             end do
          end do

          !block 11
          do j=j_max+1, size(new_grdpts_id,2)
             do i=1, size(new_grdpts_id,1)
                new_grdpts_id(i,j) = no_pt
             end do
          end do

          call MOVE_ALLOC(new_grdpts_id,bf_grdpts_id)

        end subroutine reallocate_grdpts_id_N


        !< reallocate nodes for southern buffer layers
        subroutine reallocate_nodes_S(
     $     bf_nodes, interior_nodes,
     $     i_min1, i_min3, i_min4,
     $     i_match, i_max,
     $     interior_i_max1, interior_i_max2,
     $     bf_alignment, new_sizes)

          implicit none

          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer(ikind)                                  , intent(in)    :: i_min1
          integer(ikind)                                  , intent(in)    :: i_min3
          integer(ikind)                                  , intent(in)    :: i_min4
          integer(ikind)                                  , intent(in)    :: i_match
          integer(ikind)                                  , intent(in)    :: i_max
          integer(ikind)                                  , intent(in)    :: interior_i_max1
          integer(ikind)                                  , intent(in)    :: interior_i_max2
          integer(ikind), dimension(2,2)                  , intent(in)    :: bf_alignment
          integer(ikind), dimension(2)                    , intent(in)    :: new_sizes

          real(rkind), dimension(:,:,:), allocatable :: new_nodes
          integer(ikind) :: i,j,j_min,j_start
          integer        :: k

          allocate(new_nodes(new_sizes(1),new_sizes(2),ne))

          j_min = size(new_nodes,2)-size(bf_nodes,2)
          if(j_min.lt.0) then
             j_start = 0
          else
             j_start = j_min
          end if
          
          do k=1, ne

             do j=j_start+1, size(new_nodes,2)-(2*bc_size)
                 do i=1, i_max
                   new_nodes(i_min3+i,j,k) = bf_nodes(i_match+i,j-j_min,k)
                end do
             end do

             do j=size(new_nodes,2)-(2*bc_size)+1, size(new_nodes,2)
                do i=1, interior_i_max1
                   new_nodes(i_min1+i,j,k) = interior_nodes(
     $                  bf_alignment(1,1)-(bc_size+1)+i_min1+i,
     $                  j-(size(new_nodes,2)-2*bc_size),
     $                  k)
                end do

                do i=1, i_max
                   new_nodes(i_min3+i,j,k) = bf_nodes(i_match+i,j-j_min,k)
                end do

                do i=1, interior_i_max2
                   new_nodes(i_min4+i,j,k) = interior_nodes(
     $                  bf_alignment(1,1)-(bc_size+1)+i_min4+i,
     $                  j-(size(new_nodes,2)-2*bc_size),
     $                  k)
                end do
             end do
             
          end do

          call MOVE_ALLOC(new_nodes,bf_nodes)

        end subroutine reallocate_nodes_S


        !< reallocate gridpts_id for northern buffer layers
        subroutine reallocate_grdpts_id_S(
     $     bf_grdpts_id,
     $     i_min1, i_min3, i_min4, i_min5, i_min6,
     $     i_match, i_max, j_max,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2,
     $     bf_alignment, new_sizes,
     $     i_match_border_E)

          implicit none

          integer       , dimension(:,:)     , allocatable, intent(inout) :: bf_grdpts_id
          integer(ikind)                                  , intent(in)    :: i_min1
          integer(ikind)                                  , intent(in)    :: i_min3
          integer(ikind)                                  , intent(in)    :: i_min4
          integer(ikind)                                  , intent(in)    :: i_min5
          integer(ikind)                                  , intent(in)    :: i_min6
          integer(ikind)                                  , intent(in)    :: i_match
          integer(ikind)                                  , intent(in)    :: i_max
          integer(ikind)                                  , intent(in)    :: j_max
          integer(ikind)                                  , intent(in)    :: outside_i_max1
          integer(ikind)                                  , intent(in)    :: outside_i_max2
          integer(ikind)                                  , intent(in)    :: interior_i_max1
          integer(ikind)                                  , intent(in)    :: interior_i_max2
          integer(ikind), dimension(2,2)                  , intent(in)    :: bf_alignment
          integer(ikind), dimension(2)                    , intent(in)    :: new_sizes
          integer(ikind)                                  , intent(in)    :: i_match_border_E


          integer, dimension(:,:), allocatable :: new_grdpts_id
          
          integer(ikind) :: i,j,j_min,j_match_profile,j_start

          integer, dimension(bc_size, 2*bc_size) :: border_W
          integer, dimension(bc_size, 2*bc_size) :: border_E
          integer, dimension(2*bc_size)          :: interior_profile


          allocate(new_grdpts_id(new_sizes(1),new_sizes(2)))

          j_min           = size(new_grdpts_id,2)-j_max
          j_match_profile = size(new_grdpts_id,2)-(2*bc_size)
          if(j_min.lt.0) then
             j_start = 0
          else
             j_start = j_min
          end if

          !get the special border for the left and right
          !as well as the interior profile
          call get_additional_blocks_S(
     $         bf_alignment,
     $         border_W, border_E, interior_profile)


          !fill the blocks according to Fig.2

          !block 11
          do j=1, j_min
             do i=1, size(new_grdpts_id,1)
                new_grdpts_id(i,j) = no_pt
             end do
          end do


          do j=j_start+1, size(new_grdpts_id,2)-(2*bc_size)
             !block 8
             do i=1, outside_i_max1+interior_i_max1
                new_grdpts_id(i,j) = no_pt
             end do

             !block 9
             do i=1, i_max
                new_grdpts_id(i_min3+i,j) = bf_grdpts_id(i_match+i,j-j_min)
             end do

             !block 10
             do i=1, interior_i_max2+outside_i_max2
                new_grdpts_id(i_min4+i,j) = no_pt
             end do
          end do


          do j=size(new_grdpts_id,2)-(2*bc_size)+1, size(new_grdpts_id,2)
             !block 1
             do i=1, outside_i_max1
                new_grdpts_id(i,j) = no_pt
             end do

             !block 2
             do i=1, min(bc_size, interior_i_max1)
                new_grdpts_id(i_min1+i,j) = border_W(i,j-j_match_profile)
             end do

             !block 3
             do i=bc_size+1, interior_i_max1
                new_grdpts_id(i_min1+i,j) = interior_profile(j-j_match_profile)
             end do

             !block 4
             do i=1, i_max
                new_grdpts_id(i_min3+i,j) = bf_grdpts_id(i_match+i,j-j_min)
             end do

             !block 5
             do i=1, interior_i_max2-bc_size
                new_grdpts_id(i_min4+i,j) = interior_profile(j-j_match_profile)
             end do

             !block 6
             do i=1, min(interior_i_max2,bc_size)
                new_grdpts_id(i_min5+i,j) = border_E(i_match_border_E+i,j-j_match_profile)
             end do

             !block 7
             do i=1, outside_i_max2
                new_grdpts_id(i_min6+i,j) = no_pt
             end do
          end do          

          call MOVE_ALLOC(new_grdpts_id,bf_grdpts_id)

        end subroutine reallocate_grdpts_id_S


        !< reallocate nodes for eastern buffer layers
        subroutine reallocate_nodes_E(
     $     bf_nodes, interior_nodes,
     $     j_min1, j_min3, j_min4,
     $     j_match, i_max, j_max,
     $     interior_j_max1, interior_j_max2,
     $     bf_alignment, new_sizes)

          implicit none

          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer(ikind)                                  , intent(in)    :: j_min1
          integer(ikind)                                  , intent(in)    :: j_min3
          integer(ikind)                                  , intent(in)    :: j_min4
          integer(ikind)                                  , intent(in)    :: j_match
          integer(ikind)                                  , intent(in)    :: i_max
          integer(ikind)                                  , intent(in)    :: j_max
          integer(ikind)                                  , intent(in)    :: interior_j_max1
          integer(ikind)                                  , intent(in)    :: interior_j_max2
          integer(ikind), dimension(2,2)                  , intent(in)    :: bf_alignment
          integer(ikind), dimension(2)                    , intent(in)    :: new_sizes

          real(rkind), dimension(:,:,:), allocatable :: new_nodes
          integer(ikind) :: i,j
          integer        :: k

          allocate(new_nodes(new_sizes(1),new_sizes(2),ne))

          do k=1, ne
             do j=1, interior_j_max1
                do i=1, 2*bc_size
                   new_nodes(i,j_min1+j,k) = interior_nodes(
     $                  nx-(2*bc_size)+i,
     $                  bf_alignment(2,1)-(bc_size+1)+j_min1+j,
     $                  k)
                end do
             end do

             do j=1, j_max
                do i=1, i_max
                   new_nodes(i,j_min3+j,k) = bf_nodes(i,j_match+j,k)
                end do
             end do
             
             do j=1, interior_j_max2
                do i=1, 2*bc_size
                   new_nodes(i,j_min4+j,k) = interior_nodes(
     $                  nx-(2*bc_size)+i,
     $                  bf_alignment(2,1)-(bc_size+1)+j_min4+j,
     $                  k)
                end do
             end do
          end do

          call MOVE_ALLOC(new_nodes,bf_nodes)

        end subroutine reallocate_nodes_E


        !< reallocate gridpts_id for eastern buffer layers
        subroutine reallocate_grdpts_id_E(
     $     bf_grdpts_id,
     $     j_min1, j_min3, j_min4, j_min5, j_min6,
     $     j_match, i_max, j_max,
     $     outside_j_max1, outside_j_max2,
     $     interior_j_max1, interior_j_max2,
     $     bf_alignment, new_sizes,
     $     j_match_border_N)

          implicit none

          integer       , dimension(:,:), allocatable, intent(inout) :: bf_grdpts_id
          integer(ikind)                             , intent(in)    :: j_min1
          integer(ikind)                             , intent(in)    :: j_min3
          integer(ikind)                             , intent(in)    :: j_min4
          integer(ikind)                             , intent(in)    :: j_min5
          integer(ikind)                             , intent(in)    :: j_min6
          integer(ikind)                             , intent(in)    :: j_match
          integer(ikind)                             , intent(in)    :: i_max
          integer(ikind)                             , intent(in)    :: j_max
          integer(ikind)                             , intent(in)    :: outside_j_max1
          integer(ikind)                             , intent(in)    :: outside_j_max2
          integer(ikind)                             , intent(in)    :: interior_j_max1
          integer(ikind)                             , intent(in)    :: interior_j_max2
          integer(ikind), dimension(2,2)             , intent(in)    :: bf_alignment
          integer(ikind), dimension(2)               , intent(in)    :: new_sizes
          integer(ikind)                             , intent(in)    :: j_match_border_N


          integer, dimension(:,:), allocatable :: new_grdpts_id
          
          integer(ikind) :: i,j

          integer, dimension(2*bc_size,bc_size) :: border_N
          integer, dimension(2*bc_size,bc_size) :: border_S
          integer, dimension(2*bc_size)         :: interior_profile


          allocate(new_grdpts_id(new_sizes(1),new_sizes(2)))

          
          !get the special border for the left and right
          !as well as the interior profile
          call get_additional_blocks_E(
     $         bf_alignment,
     $         border_S, border_N, interior_profile)


          !fill the blocks 1 and partially 8 and 11
          do j=1, outside_j_max1
             do i=1, size(new_grdpts_id,1)
                new_grdpts_id(i,j) = no_pt
             end do
          end do

          
          !fill the blocks 2 and partially 8 and 11
          do j=1, min(interior_j_max1, bc_size)
             do i=1, 2*bc_size
                new_grdpts_id(i,j_min1+j) = border_S(i,j)
             end do             
             do i=2*bc_size+1, size(new_grdpts_id,1)
                new_grdpts_id(i,j_min1+j) = no_pt
             end do
          end do


          !fill the blocks 3 and partially 8 and 11
          do j=bc_size+1, interior_j_max1
             do i=1, 2*bc_size
                new_grdpts_id(i,j_min1+j) = interior_profile(i)
             end do
             do i=2*bc_size+1, size(new_grdpts_id,1)
                new_grdpts_id(i,j_min1+j) = no_pt
             end do
          end do


          !fill the blocks 4, 9 and partially 11
          do j=1, j_max
             do i=1, i_max
                new_grdpts_id(i,j_min3+j) = bf_grdpts_id(i,j_match+j)
             end do
             do i=i_max+1, size(new_grdpts_id,1)
                new_grdpts_id(i,j_min3+j) = no_pt
             end do
          end do


          !fill the blocks 5 and partially 10 and 11
          do j=1, interior_j_max2-bc_size
             do i=1, 2*bc_size
                new_grdpts_id(i,j_min4+j) = interior_profile(i)
             end do
             do i=2*bc_size+1, size(new_grdpts_id,1)
                new_grdpts_id(i,j_min4+j) = no_pt
             end do
          end do
          

          !fill the blocks 6 and partially 10 and 11
          do j=1, min(interior_j_max2,bc_size)
             do i=1, 2*bc_size
                new_grdpts_id(i,j_min5+j) = border_N(i,j_match_border_N+j)
             end do
             do i=2*bc_size+1, size(new_grdpts_id,1)
                new_grdpts_id(i,j_min5+j) = no_pt
             end do
          end do


          !fill the blocks 7 and partially 10 and 11
          do j=1, outside_j_max2
             do i=1, size(new_grdpts_id,1)
                new_grdpts_id(i,j_min6+j) = no_pt
             end do
          end do


          !replace the content of bf_grdpts_id with new_grdpts_id
          call MOVE_ALLOC(new_grdpts_id,bf_grdpts_id)

        end subroutine reallocate_grdpts_id_E


        !< reallocate nodes for western buffer layers
        subroutine reallocate_nodes_W(
     $     bf_nodes, interior_nodes,
     $     j_min1, j_min3, j_min4,
     $     j_match, j_max,
     $     interior_j_max1, interior_j_max2,
     $     bf_alignment, new_sizes)

          implicit none

          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer(ikind)                                  , intent(in)    :: j_min1
          integer(ikind)                                  , intent(in)    :: j_min3
          integer(ikind)                                  , intent(in)    :: j_min4
          integer(ikind)                                  , intent(in)    :: j_match
          integer(ikind)                                  , intent(in)    :: j_max
          integer(ikind)                                  , intent(in)    :: interior_j_max1
          integer(ikind)                                  , intent(in)    :: interior_j_max2
          integer(ikind), dimension(2,2)                  , intent(in)    :: bf_alignment
          integer(ikind), dimension(2)                    , intent(in)    :: new_sizes

          real(rkind), dimension(:,:,:), allocatable :: new_nodes
          integer(ikind) :: i,j,i_min,i_start
          integer        :: k

          allocate(new_nodes(new_sizes(1),new_sizes(2),ne))

          i_min = size(new_nodes,1)-size(bf_nodes,1)
          if(i_min.lt.0) then
             i_start = 0
          else
             i_start = i_min
          end if

          do k=1, ne
             do j=1, interior_j_max1
                do i=size(new_nodes,1)-(2*bc_size)+1, size(new_nodes,1)
                   new_nodes(i,j_min1+j,k) = interior_nodes(
     $                  i-(size(new_nodes,1)-2*bc_size),
     $                  bf_alignment(2,1)-(bc_size+1)+j_min1+j,
     $                  k)
                end do
             end do

             do j=1, j_max
                do i=i_start+1, size(new_nodes,1)
                   new_nodes(i,j_min3+j,k) = bf_nodes(i-i_min,j_match+j,k)
                end do
             end do
             
             do j=1, interior_j_max2
                do i=size(new_nodes,1)-(2*bc_size)+1, size(new_nodes,1)
                   new_nodes(i,j_min4+j,k) = interior_nodes(
     $                  i-(size(new_nodes,1)-2*bc_size),
     $                  bf_alignment(2,1)-(bc_size+1)+j_min4+j,
     $                  k)
                end do
             end do
          end do

          call MOVE_ALLOC(new_nodes,bf_nodes)

        end subroutine reallocate_nodes_W


        !< reallocate gridpts_id for eastern buffer layers
        subroutine reallocate_grdpts_id_W(
     $     bf_grdpts_id,
     $     j_min1, j_min3, j_min4, j_min5, j_min6,
     $     j_match, j_max,
     $     outside_j_max1, outside_j_max2,
     $     interior_j_max1, interior_j_max2,
     $     bf_alignment, new_sizes,
     $     j_match_border_N)

          implicit none

          integer       , dimension(:,:), allocatable, intent(inout) :: bf_grdpts_id
          integer(ikind)                             , intent(in)    :: j_min1
          integer(ikind)                             , intent(in)    :: j_min3
          integer(ikind)                             , intent(in)    :: j_min4
          integer(ikind)                             , intent(in)    :: j_min5
          integer(ikind)                             , intent(in)    :: j_min6
          integer(ikind)                             , intent(in)    :: j_match
          integer(ikind)                             , intent(in)    :: j_max
          integer(ikind)                             , intent(in)    :: outside_j_max1
          integer(ikind)                             , intent(in)    :: outside_j_max2
          integer(ikind)                             , intent(in)    :: interior_j_max1
          integer(ikind)                             , intent(in)    :: interior_j_max2
          integer(ikind), dimension(2,2)             , intent(in)    :: bf_alignment
          integer(ikind), dimension(2)               , intent(in)    :: new_sizes
          integer(ikind)                             , intent(in)    :: j_match_border_N


          integer, dimension(:,:), allocatable :: new_grdpts_id
          
          integer(ikind) :: i,j,i_min,i_match_profile,i_start

          integer, dimension(2*bc_size,bc_size) :: border_N
          integer, dimension(2*bc_size,bc_size) :: border_S
          integer, dimension(2*bc_size)         :: interior_profile


          allocate(new_grdpts_id(new_sizes(1),new_sizes(2)))

          i_min           = size(new_grdpts_id,1)-size(bf_grdpts_id,1)
          i_match_profile = size(new_grdpts_id,1)-(2*bc_size)
          if(i_min.lt.0) then
             i_start = 0
          else
             i_start = i_min
          end if

          !get the special border for the left and right
          !as well as the interior profile
          call get_additional_blocks_W(
     $         bf_alignment,
     $         border_S, border_N, interior_profile)


          !fill the blocks 1 and partially 8 and 11
          do j=1, outside_j_max1
             do i=1, size(new_grdpts_id,1)
                new_grdpts_id(i,j) = no_pt
             end do
          end do

          
          !fill the blocks 2 and partially 8 and 11
          do j=1, min(interior_j_max1, bc_size)
             do i=1, i_match_profile
                new_grdpts_id(i,j_min1+j) = no_pt
             end do
             do i=i_match_profile+1, size(new_grdpts_id,1)
                new_grdpts_id(i,j_min1+j)= border_S(i-i_match_profile,j)
             end do                          
          end do


          !fill the blocks 3 and partially 8 and 11
          do j=bc_size+1, interior_j_max1
             do i=1, i_match_profile
                new_grdpts_id(i,j_min1+j) = no_pt
             end do
             do i=i_match_profile+1, size(new_grdpts_id,1)
                new_grdpts_id(i,j_min1+j) = interior_profile(i-i_match_profile)
             end do             
          end do


          !fill the blocks 4, 9 and partially 11
          do j=1, j_max
             do i=1, i_min
                new_grdpts_id(i,j_min3+j) = no_pt
             end do
             do i=i_start+1, size(new_grdpts_id,1)
                new_grdpts_id(i,j_min3+j) = bf_grdpts_id(i-i_min,j_match+j)
             end do
          end do


          !fill the blocks 5 and partially 10 and 11
          do j=1, interior_j_max2-bc_size
             do i=1, i_match_profile
                new_grdpts_id(i,j_min4+j) = no_pt
             end do
             do i=i_match_profile+1, size(new_grdpts_id,1)
                new_grdpts_id(i,j_min4+j) = interior_profile(i-i_match_profile)
             end do             
          end do
          

          !fill the blocks 6 and partially 10 and 11
          do j=1, min(interior_j_max2,bc_size)
             do i=1, i_match_profile
                new_grdpts_id(i,j_min5+j) = no_pt
             end do
             do i=i_match_profile+1, size(new_grdpts_id,1)
                new_grdpts_id(i,j_min5+j) = border_N(i-i_match_profile,
     $                                               j+j_match_border_N)
             end do             
          end do


          !fill the blocks 7 and partially 10 and 11
          do j=1, outside_j_max2
             do i=1, size(new_grdpts_id,1)
                new_grdpts_id(i,j_min6+j) = no_pt
             end do
          end do


          !replace the content of bf_grdpts_id with new_grdpts_id
          call MOVE_ALLOC(new_grdpts_id,bf_grdpts_id)

        end subroutine reallocate_grdpts_id_W


        subroutine get_match(
     $     dir,
     $     bf_alignment, border_changes,
     $     size_x, size_y,
     $     i_min1, i_min3, i_min4, i_min5, i_min6,
     $     i_match, i_max, j_max,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1, interior_i_max2,
     $     i_match_border_W)
        
          implicit none

          integer                       , intent(in)  :: dir
          integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
          integer(ikind), dimension(2,2), intent(in)  :: border_changes
          integer(ikind)                , intent(in)  :: size_x
          integer(ikind)                , intent(in)  :: size_y
          integer(ikind)                , intent(out) :: i_min1
          integer(ikind)                , intent(out) :: i_min3
          integer(ikind)                , intent(out) :: i_min4
          integer(ikind)                , intent(out) :: i_min5
          integer(ikind)                , intent(out) :: i_min6
          integer(ikind)                , intent(out) :: i_match
          integer(ikind)                , intent(out) :: i_max
          integer(ikind)                , intent(out) :: j_max
          integer(ikind)                , intent(out) :: outside_i_max1
          integer(ikind)                , intent(out) :: outside_i_max2
          integer(ikind)                , intent(out) :: interior_i_max1
          integer(ikind)                , intent(out) :: interior_i_max2
          integer(ikind)                , intent(out) :: i_match_border_W

          integer(ikind), dimension(2) :: old_bf_alignment
          integer(ikind) :: ndir

          select case(dir)
            case(x_direction)
               ndir = nx
            case(y_direction)
               ndir = ny
            case default
               print '(''bf_layer_reallocate_module'')'
               print '(''get_match'')'
               stop 'dir not recognized'
          end select


          !compute the previosu alignment
          old_bf_alignment(1) = bf_alignment(dir,1)-border_changes(dir,1)
          old_bf_alignment(2) = bf_alignment(dir,2)-border_changes(dir,2)


          !define the length of the blocks 1 and 8
          !the length is restricted by the length of the
          !interior domain
          if(old_bf_alignment(1).le.(bc_size+1)) then
             outside_i_max1 = max(0, -border_changes(dir,1))
          else
             if(bf_alignment(dir,1).le.(bc_size+1)) then
                outside_i_max1 = bc_size + 1 - bf_alignment(dir,1)
             else
                outside_i_max1 = 0
             end if
          end if


          !define the length of the blocks 7 and 12
          !the length is restricted by the length of the
          !interior domain
          if(old_bf_alignment(2).ge.(ndir-bc_size)) then
             outside_i_max2 = max(0, border_changes(dir,2))
          else
             if(bf_alignment(dir,2).ge.(ndir-bc_size)) then
                outside_i_max2 = bf_alignment(dir,2)-(ndir-bc_size)
             else
                outside_i_max2 = 0
             end if
          end if


          !define the length of the blocks 2 and 3
          !the length is restricted by the length of the
          !interior domain
          if(old_bf_alignment(1).ge.(bc_size+1)) then
             if(bf_alignment(dir,1).ge.(bc_size+1)) then
                interior_i_max1 = max(0,-border_changes(dir,1))
             else
                interior_i_max1 = old_bf_alignment(1)-(bc_size+1)
             end if
          else
             interior_i_max1 = 0
          end if


          !define the length of the blocks 5 and 6
          !the length is restricted by the length of the
          !interior domain
          if(old_bf_alignment(2).le.(ndir-bc_size)) then
             if(bf_alignment(dir,2).le.(ndir-bc_size)) then
                interior_i_max2 = max(0, border_changes(dir,2))
             else
                interior_i_max2 = (ndir-bc_size) - old_bf_alignment(2)
             end if
          else
             interior_i_max2 = 0
          end if


          !define the length of the blocks 4 and 9
          !in the case the reallocated table in smaller than the
          !original one, we need to change the borders identifying
          !which part of the original table is copied
          i_match = 0 + max(0,border_changes(dir,1))
          i_max   = size_x -
     $         max(0,border_changes(1,1)) -
     $         max(0,-border_changes(1,2))
          j_max   = size_y -
     $         max(0,border_changes(2,1)) -
     $         max(0,-border_changes(2,2))

          i_min1 = outside_i_max1
          i_min3 = i_min1 + interior_i_max1
          if(dir.eq.x_direction) then
             i_min4 = i_min3 + i_max
          else
             i_min4 = i_min3 + j_max
          end if
          i_min5 = i_min4 + max(0,interior_i_max2-bc_size)
          i_min6 = i_min4 + interior_i_max2

          i_match_border_W = max(0,bf_alignment(dir,1)-
     $         (bc_size+1)+i_min5-(ndir-bc_size))

        end subroutine get_match


c$$$        subroutine get_additional_blocks_N(
c$$$     $     bf_alignment,
c$$$     $     border_W, border_E, interior_profile)
c$$$
c$$$          implicit none
c$$$
c$$$          integer(ikind), dimension(2,2)              , intent(in)  :: bf_alignment
c$$$          integer       , dimension(bc_size,2*bc_size), intent(out) :: border_W
c$$$          integer       , dimension(bc_size,2*bc_size), intent(out) :: border_E
c$$$          integer       , dimension(2*bc_size)        , intent(out) :: interior_profile
c$$$
c$$$          integer :: i,j
c$$$
c$$$          !interior profile
c$$$          do i=1, bc_size
c$$$             interior_profile(i) = interior_pt
c$$$          end do
c$$$          interior_profile(bc_size+1) = bc_interior_pt
c$$$          do i=bc_size+2, 2*bc_size
c$$$             interior_profile(i) = bc_pt
c$$$          end do
c$$$
c$$$          !border_W
c$$$          if(bf_alignment(1,1).gt.(bc_size+2)) then
c$$$             do j=1, bc_size
c$$$                do i=1, bc_size
c$$$                   border_W(i,j) = interior_pt
c$$$                end do
c$$$             end do
c$$$             j=bc_size+1
c$$$             do i=1, bc_size
c$$$                border_W(i,j) = bc_interior_pt
c$$$             end do
c$$$             do j=bc_size+2, 2*bc_size
c$$$                do i=1, bc_size
c$$$                   border_W(i,j) = bc_pt
c$$$                end do
c$$$             end do
c$$$          else
c$$$             if(bf_alignment(1,1).eq.bc_size+2) then
c$$$                do j=1, 2*bc_size-2
c$$$                   border_W(1,j) = bc_interior_pt
c$$$                   border_W(2,j) = interior_pt
c$$$               end do
c$$$
c$$$               j=2*bc_size-1
c$$$               border_W(1,j) = bc_interior_pt
c$$$               border_W(2,j) = bc_interior_pt
c$$$
c$$$               j=2*bc_size
c$$$               border_W(1,j) = bc_pt
c$$$               border_W(2,j) = bc_pt
c$$$
c$$$             else
c$$$               do j=1, 2*bc_size-1
c$$$                  border_W(1,j) = bc_pt
c$$$                  border_W(2,j) = bc_interior_pt
c$$$               end do
c$$$               j=2*bc_size
c$$$               border_W(1,j) = bc_pt
c$$$               border_W(2,j) = bc_pt
c$$$             end if
c$$$          end if
c$$$
c$$$          !border_E
c$$$          if(bf_alignment(1,2).lt.(nx-bc_size-1)) then
c$$$             do j=1, bc_size
c$$$                do i=1, bc_size
c$$$                   border_E(i,j) = interior_pt
c$$$                end do
c$$$             end do
c$$$             j=bc_size+1
c$$$             do i=1, bc_size
c$$$                border_E(i,j) = bc_interior_pt
c$$$             end do
c$$$             do j=bc_size+2, 2*bc_size
c$$$                do i=1, bc_size
c$$$                   border_E(i,j) = bc_pt
c$$$                end do
c$$$             end do
c$$$          else
c$$$             if(bf_alignment(1,2).eq.(nx-bc_size-1)) then
c$$$                do j=1, 2*bc_size-2
c$$$                   border_E(1,j) = interior_pt
c$$$                   border_E(2,j) = bc_interior_pt
c$$$               end do
c$$$
c$$$               j=2*bc_size-1
c$$$               border_E(1,j) = bc_interior_pt
c$$$               border_E(2,j) = bc_interior_pt
c$$$
c$$$               j=2*bc_size
c$$$               border_E(1,j) = bc_pt
c$$$               border_E(2,j) = bc_pt
c$$$
c$$$             else
c$$$               do j=1, 2*bc_size-1
c$$$                  border_E(1,j) = bc_interior_pt
c$$$                  border_E(2,j) = bc_pt
c$$$               end do
c$$$               j=2*bc_size
c$$$               border_E(1,j) = bc_pt
c$$$               border_E(2,j) = bc_pt
c$$$             end if
c$$$          end if
c$$$
c$$$        end subroutine get_additional_blocks_N
c$$$
c$$$
c$$$        subroutine get_additional_blocks_S(
c$$$     $     bf_alignment,
c$$$     $     border_W, border_E, interior_profile)
c$$$
c$$$          implicit none
c$$$
c$$$          integer(ikind), dimension(2,2)              , intent(in)  :: bf_alignment
c$$$          integer       , dimension(bc_size,2*bc_size), intent(out) :: border_W
c$$$          integer       , dimension(bc_size,2*bc_size), intent(out) :: border_E
c$$$          integer       , dimension(2*bc_size)        , intent(out) :: interior_profile
c$$$
c$$$          integer :: i,j
c$$$
c$$$          !interior profile
c$$$          do i=1, bc_size-1
c$$$             interior_profile(i) = bc_pt
c$$$          end do
c$$$          interior_profile(bc_size) = bc_interior_pt
c$$$          do i=bc_size+1, 2*bc_size
c$$$             interior_profile(i) = interior_pt
c$$$          end do          
c$$$
c$$$          !border_W
c$$$          if(bf_alignment(1,1).gt.(bc_size+2)) then
c$$$             do j=1, bc_size-1
c$$$                do i=1, bc_size
c$$$                   border_W(i,j) = bc_pt
c$$$                end do
c$$$             end do
c$$$             j=bc_size
c$$$             do i=1, bc_size
c$$$                border_W(i,j) = bc_interior_pt
c$$$             end do
c$$$             do j=bc_size+1, 2*bc_size
c$$$                do i=1, bc_size
c$$$                   border_W(i,j) = interior_pt
c$$$                end do
c$$$             end do
c$$$             
c$$$          else
c$$$             if(bf_alignment(1,1).eq.bc_size+2) then
c$$$                j=1
c$$$                border_W(1,j) = bc_pt
c$$$                border_W(2,j) = bc_pt
c$$$               
c$$$                j=2
c$$$                border_W(1,j) = bc_interior_pt
c$$$                border_W(2,j) = bc_interior_pt
c$$$
c$$$                do j=3, 2*bc_size
c$$$                   border_W(1,j) = bc_interior_pt
c$$$                   border_W(2,j) = interior_pt
c$$$                end do               
c$$$
c$$$             else
c$$$                j=1
c$$$                border_W(1,j) = bc_pt
c$$$                border_W(2,j) = bc_pt
c$$$
c$$$                do j=2, 2*bc_size
c$$$                   border_W(1,j) = bc_pt
c$$$                   border_W(2,j) = bc_interior_pt
c$$$                end do
c$$$               
c$$$             end if
c$$$          end if
c$$$
c$$$          !border_E
c$$$          if(bf_alignment(1,2).lt.(nx-bc_size-1)) then
c$$$             do j=1, bc_size-1
c$$$                do i=1, bc_size
c$$$                   border_E(i,j) = bc_pt
c$$$                end do
c$$$             end do
c$$$             j=bc_size
c$$$             do i=1, bc_size
c$$$                border_E(i,j) = bc_interior_pt
c$$$             end do
c$$$             do j=bc_size+1, 2*bc_size
c$$$                do i=1, bc_size
c$$$                   border_E(i,j) = interior_pt
c$$$                end do
c$$$             end do
c$$$             
c$$$          else
c$$$             if(bf_alignment(1,2).eq.(nx-bc_size-1)) then
c$$$                j=1
c$$$                border_E(1,j) = bc_pt
c$$$                border_E(2,j) = bc_pt
c$$$
c$$$                j=2
c$$$                border_E(1,j) = bc_interior_pt
c$$$                border_E(2,j) = bc_interior_pt
c$$$                
c$$$                do j=3, 2*bc_size
c$$$                   border_E(1,j) = interior_pt
c$$$                   border_E(2,j) = bc_interior_pt
c$$$                end do
c$$$
c$$$             else
c$$$               j=1
c$$$               border_E(1,j) = bc_pt
c$$$               border_E(2,j) = bc_pt
c$$$
c$$$               do j=2, 2*bc_size
c$$$                  border_E(1,j) = bc_interior_pt
c$$$                  border_E(2,j) = bc_pt
c$$$               end do
c$$$             end if
c$$$          end if
c$$$
c$$$        end subroutine get_additional_blocks_S
c$$$
c$$$
c$$$        subroutine get_additional_blocks_E(
c$$$     $     bf_alignment,
c$$$     $     border_S, border_N, interior_profile)
c$$$
c$$$          implicit none
c$$$
c$$$          integer(ikind), dimension(2,2)              , intent(in)  :: bf_alignment
c$$$          integer       , dimension(2*bc_size,bc_size), intent(out) :: border_S
c$$$          integer       , dimension(2*bc_size,bc_size), intent(out) :: border_N
c$$$          integer       , dimension(2*bc_size)        , intent(out) :: interior_profile
c$$$
c$$$          integer :: i,j
c$$$
c$$$          !interior profile
c$$$          do i=1, bc_size
c$$$             interior_profile(i) = interior_pt
c$$$          end do
c$$$          interior_profile(bc_size+1) = bc_interior_pt
c$$$          do i=bc_size+2, 2*bc_size
c$$$             interior_profile(i) = bc_pt
c$$$          end do
c$$$
c$$$          !border_S
c$$$          if(bf_alignment(2,1).gt.(bc_size+2)) then
c$$$             do j=1, bc_size
c$$$                do i=1, bc_size
c$$$                   border_S(i,j) = interior_pt
c$$$                end do
c$$$
c$$$                i = bc_size+1
c$$$                border_S(i,j) = bc_interior_pt
c$$$
c$$$                do i=bc_size+2, 2*bc_size
c$$$                   border_S(i,j) = bc_pt
c$$$                end do
c$$$             end do
c$$$
c$$$          else
c$$$             if(bf_alignment(2,1).eq.bc_size+2) then
c$$$                j=1
c$$$                do i=1, bc_size+1
c$$$                   border_S(i,j) = bc_interior_pt
c$$$                end do
c$$$                do i=bc_size+2,2*bc_size
c$$$                   border_S(i,j) = bc_pt
c$$$                end do
c$$$                
c$$$                do j=2, bc_size
c$$$                   do i=1, bc_size
c$$$                      border_S(i,j) = interior_pt
c$$$                   end do
c$$$                   i=bc_size+1
c$$$                   border_S(i,j) = bc_interior_pt
c$$$                   do i=bc_size+2, 2*bc_size
c$$$                      border_S(i,j) = bc_pt
c$$$                   end do
c$$$                end do
c$$$
c$$$             else
c$$$                j=1
c$$$                do i=1, 2*bc_size
c$$$                   border_S(i,j) = bc_pt
c$$$                end do
c$$$                do j=2, bc_size
c$$$                   do i=1, bc_size+1
c$$$                      border_S(i,j) = bc_interior_pt
c$$$                   end do
c$$$                   do i=bc_size+2, 2*bc_size
c$$$                      border_S(i,j) = bc_pt
c$$$                   end do
c$$$                end do
c$$$             end if
c$$$          end if
c$$$
c$$$          !border_N
c$$$          if(bf_alignment(2,2).lt.(ny-bc_size-1)) then
c$$$             do j=1, bc_size
c$$$                do i=1, bc_size
c$$$                   border_N(i,j) = interior_pt
c$$$                end do
c$$$                i=bc_size+1
c$$$                border_N(i,j) = bc_interior_pt
c$$$                do i=bc_size+2, 2*bc_size
c$$$                   border_N(i,j) = bc_pt
c$$$                end do
c$$$             end do                
c$$$
c$$$          else
c$$$             if(bf_alignment(2,2).eq.(ny-bc_size-1)) then
c$$$                j=1
c$$$                do i=1, bc_size
c$$$                   border_N(i,j) = interior_pt
c$$$                end do
c$$$                i=bc_size+1
c$$$                border_N(i,j) = bc_interior_pt
c$$$                do i=bc_size+2, 2*bc_size
c$$$                   border_N(i,j) = bc_pt
c$$$                end do
c$$$
c$$$                do j=2, bc_size
c$$$                   do i=1, bc_size+1
c$$$                      border_N(i,j) = bc_interior_pt
c$$$                   end do
c$$$                   do i=bc_size+2, 2*bc_size
c$$$                      border_N(i,j) = bc_pt
c$$$                   end do
c$$$                end do
c$$$
c$$$             else
c$$$                 j=1
c$$$                 do i=1, bc_size+1
c$$$                    border_N(i,j) = bc_interior_pt
c$$$                 end do
c$$$                 do i=bc_size+2, 2*bc_size
c$$$                    border_N(i,j) = bc_pt
c$$$                 end do
c$$$                 do j=2, bc_size
c$$$                    do i=1, 2*bc_size
c$$$                       border_N(i,j) = bc_pt
c$$$                    end do
c$$$                 end do
c$$$             end if
c$$$          end if
c$$$
c$$$        end subroutine get_additional_blocks_E
c$$$
c$$$
c$$$        !< compute the additional block W for Fig.4
c$$$        subroutine get_additional_blocks_W(
c$$$     $     bf_alignment,
c$$$     $     border_S, border_N, interior_profile)
c$$$
c$$$          implicit none
c$$$
c$$$          integer(ikind), dimension(2,2)              , intent(in)  :: bf_alignment
c$$$          integer       , dimension(2*bc_size,bc_size), intent(out) :: border_S
c$$$          integer       , dimension(2*bc_size,bc_size), intent(out) :: border_N
c$$$          integer       , dimension(2*bc_size)        , intent(out) :: interior_profile
c$$$
c$$$          integer :: i,j
c$$$
c$$$          !interior profile
c$$$          do i=1, bc_size-1
c$$$             interior_profile(i) = bc_pt
c$$$          end do
c$$$          i=bc_size
c$$$          interior_profile(bc_size) = bc_interior_pt
c$$$          do i=bc_size+1, 2*bc_size
c$$$             interior_profile(i) = interior_pt
c$$$          end do
c$$$
c$$$          !border_S
c$$$          if(bf_alignment(2,1).gt.(bc_size+2)) then
c$$$             do j=1, bc_size
c$$$                do i=1, bc_size
c$$$                   border_S(i,j) = bc_pt
c$$$                end do
c$$$
c$$$                i = bc_size
c$$$                border_S(i,j) = bc_interior_pt
c$$$
c$$$                do i=bc_size+1, 2*bc_size
c$$$                   border_S(i,j) = interior_pt
c$$$                end do
c$$$             end do
c$$$
c$$$          else
c$$$             if(bf_alignment(2,1).eq.bc_size+2) then
c$$$                j=1
c$$$                do i=1,bc_size-1
c$$$                   border_S(i,j) = bc_pt
c$$$                end do
c$$$                do i=bc_size, 2*bc_size
c$$$                   border_S(i,j) = bc_interior_pt
c$$$                end do
c$$$
c$$$                do j=2, bc_size
c$$$                   do i=1,bc_size-1
c$$$                      border_S(i,j) = bc_pt
c$$$                   end do
c$$$                   i=bc_size
c$$$                   border_S(i,j) = bc_interior_pt
c$$$                   do i=bc_size+1, 2*bc_size
c$$$                      border_S(i,j) = interior_pt
c$$$                   end do
c$$$                end do                
c$$$
c$$$             else
c$$$                do j=1, bc_size-1
c$$$                   do i=1,2*bc_size
c$$$                      border_S(i,j) = bc_pt
c$$$                   end do
c$$$                end do
c$$$
c$$$                j=bc_size
c$$$                do i=1,bc_size-1
c$$$                   border_S(i,j) = bc_pt
c$$$                end do
c$$$                do i=bc_size, 2*bc_size
c$$$                   border_S(i,j) = bc_interior_pt
c$$$                end do
c$$$                
c$$$             end if
c$$$          end if
c$$$
c$$$          !border_N
c$$$          if(bf_alignment(2,2).lt.(ny-bc_size-1)) then
c$$$             do j=1, bc_size
c$$$                do i=1, bc_size-1
c$$$                   border_N(i,j) = bc_pt
c$$$                end do
c$$$                i=bc_size
c$$$                border_N(i,j) = bc_interior_pt
c$$$                do i=bc_size+1, 2*bc_size
c$$$                   border_N(i,j) = interior_pt
c$$$                end do
c$$$             end do
c$$$
c$$$          else
c$$$             if(bf_alignment(2,2).eq.(ny-bc_size-1)) then
c$$$                j=1
c$$$                do i=1, bc_size-1
c$$$                   border_N(i,j) = bc_pt
c$$$                end do
c$$$                i=bc_size
c$$$                border_N(i,j) = bc_interior_pt
c$$$                do i=bc_size+1, 2*bc_size
c$$$                   border_N(i,j) = interior_pt
c$$$                end do
c$$$
c$$$                do j=2, bc_size
c$$$                   do i=1, bc_size-1
c$$$                      border_N(i,j) = bc_pt
c$$$                   end do
c$$$                   do i=bc_size, 2*bc_size
c$$$                      border_N(i,j) = bc_interior_pt
c$$$                   end do
c$$$                end do
c$$$
c$$$             else
c$$$                j=1
c$$$                do i=1, bc_size-1
c$$$                   border_N(i,j) = bc_pt
c$$$                end do
c$$$                do i=bc_size, 2*bc_size
c$$$                   border_N(i,j) = bc_interior_pt
c$$$                end do
c$$$
c$$$                do j=2, bc_size
c$$$                   do i=1, 2*bc_size
c$$$                      border_N(i,j) = bc_pt
c$$$                   end do
c$$$                end do
c$$$
c$$$             end if
c$$$          end if
c$$$
c$$$        end subroutine get_additional_blocks_W


        !< new size of tables for the buffer layer
        function get_new_sizes(border_changes, size_x, size_y)
     $     result(new_sizes)

          implicit none

          integer(ikind), dimension(2,2), intent(in) :: border_changes
          integer(ikind)                , intent(in) :: size_x
          integer(ikind)                , intent(in) :: size_y
          integer(ikind), dimension(2)               :: new_sizes

          new_sizes(1) = size_x - border_changes(1,1) + border_changes(1,2)
          new_sizes(2) = size_y - border_changes(2,1) + border_changes(2,2)

        end function get_new_sizes


        subroutine get_border_changes_and_new_alignment_N(
     $     border_changes, final_alignment, bf_alignment)

          implicit none

          integer(ikind), dimension(2,2), intent(out)  :: border_changes
          integer(ikind), dimension(2,2), intent(in)   :: final_alignment
          integer(ikind), dimension(2,2), intent(inout):: bf_alignment

          border_changes(1,1) = final_alignment(1,1) - bf_alignment(1,1)
          border_changes(2,1) = 0
          border_changes(1,2) = final_alignment(1,2) - bf_alignment(1,2)
          border_changes(2,2) = final_alignment(2,2) - bf_alignment(2,2)

          bf_alignment(1,1) = final_alignment(1,1)
          !bf_alignment(2,1) = final_alignment(2,1)
          bf_alignment(1,2) = final_alignment(1,2)
          bf_alignment(2,2) = final_alignment(2,2)

        end subroutine get_border_changes_and_new_alignment_N


        subroutine get_border_changes_and_new_alignment_S(
     $     border_changes, final_alignment, bf_alignment)

          implicit none

          integer(ikind), dimension(2,2), intent(out)  :: border_changes
          integer(ikind), dimension(2,2), intent(in)   :: final_alignment
          integer(ikind), dimension(2,2), intent(inout):: bf_alignment
          
          
          border_changes(1,1) = final_alignment(1,1) - bf_alignment(1,1)
          border_changes(2,1) = final_alignment(2,1) - bf_alignment(2,1)
          border_changes(1,2) = final_alignment(1,2) - bf_alignment(1,2)
          border_changes(2,2) = 0

          bf_alignment(1,1) = final_alignment(1,1)
          bf_alignment(2,1) = final_alignment(2,1)
          bf_alignment(1,2) = final_alignment(1,2)
          !bf_alignment(2,2) = final_alignment(2,2)

        end subroutine get_border_changes_and_new_alignment_S


        subroutine get_border_changes_and_new_alignment_E(
     $     border_changes, final_alignment, bf_alignment)

          implicit none

          integer(ikind), dimension(2,2), intent(out)  :: border_changes
          integer(ikind), dimension(2,2), intent(in)   :: final_alignment
          integer(ikind), dimension(2,2), intent(inout):: bf_alignment
          
          border_changes(1,1) = 0
          border_changes(2,1) = final_alignment(2,1) - bf_alignment(2,1)
          border_changes(1,2) = final_alignment(1,2) - bf_alignment(1,2)
          border_changes(2,2) = final_alignment(2,2) - bf_alignment(2,2)

          !bf_alignment(1,1) = final_alignment(1,1)
          bf_alignment(2,1) = final_alignment(2,1)
          bf_alignment(1,2) = final_alignment(1,2)
          bf_alignment(2,2) = final_alignment(2,2)

        end subroutine get_border_changes_and_new_alignment_E


        subroutine get_border_changes_and_new_alignment_W(
     $     border_changes, final_alignment, bf_alignment)

          implicit none

          integer(ikind), dimension(2,2), intent(out)  :: border_changes
          integer(ikind), dimension(2,2), intent(in)   :: final_alignment
          integer(ikind), dimension(2,2), intent(inout):: bf_alignment
          
          border_changes(1,1) = final_alignment(1,1) - bf_alignment(1,1)
          border_changes(2,1) = final_alignment(2,1) - bf_alignment(2,1)
          border_changes(1,2) = 0
          border_changes(2,2) = final_alignment(2,2) - bf_alignment(2,2)

          bf_alignment(1,1) = final_alignment(1,1)
          bf_alignment(2,1) = final_alignment(2,1)
          !bf_alignment(1,2) = final_alignment(1,2)
          bf_alignment(2,2) = final_alignment(2,2)

        end subroutine get_border_changes_and_new_alignment_W



      end module bf_layer_reallocate_module
