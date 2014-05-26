      module bf_layer_allocate_module

        use parameters_bf_layer, only : no_pt, interior_pt,
     $                                  bc_interior_pt, bc_pt
        use parameters_constant, only : N,S,E,W, x_direction, y_direction
        use parameters_input   , only : nx,ny,ne,bc_size
        use parameters_kind    , only : ikind, rkind     

        implicit none


        private
        public :: allocate_bf_layer_N,
     $            allocate_bf_layer_S,
     $            allocate_bf_layer_E,
     $            allocate_bf_layer_W,
     $            get_additional_blocks_N,
     $            get_additional_blocks_S,
     $            get_additional_blocks_E,
     $            get_additional_blocks_W


        contains


        subroutine allocate_bf_layer_N(
     $       bf_nodes, interior_nodes,
     $       bf_grdpts_id,
     $       bf_alignment, final_alignment)

          implicit none

          real(rkind), dimension(:,:,:), allocatable, intent(inout) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in)    :: interior_nodes
          integer    , dimension(:,:)  , allocatable, intent(inout) :: bf_grdpts_id
          integer(ikind), dimension(2,2)            , intent(out)   :: bf_alignment
          integer(ikind), dimension(2,2)            , intent(in)    :: final_alignment


          integer(ikind), dimension(2) :: new_sizes
          integer(ikind)               :: i_min1
          integer(ikind)               :: i_min3
          integer(ikind)               :: i_min4
          integer(ikind)               :: outside_i_max1
          integer(ikind)               :: outside_i_max2
          integer(ikind)               :: interior_i_max1


          !compute the new sizes
          call get_new_sizes_and_alignment_N(
     $         final_alignment, new_sizes, bf_alignment)

          !get match indices to copy the tables
          call get_match(
     $         x_direction,
     $         bf_alignment,
     $         i_min1, i_min3, i_min4,
     $         outside_i_max1, outside_i_max2,
     $         interior_i_max1)

          !allocate the nodes
          call allocate_nodes_N(
     $         bf_nodes, interior_nodes,
     $         i_min1,
     $         interior_i_max1,
     $         bf_alignment, new_sizes)

          !allocate the grdpts_id
          call allocate_grdpts_id_N(
     $         bf_grdpts_id,
     $         outside_i_max1, outside_i_max2,
     $         interior_i_max1,
     $         i_min1, i_min3, i_min4,
     $         bf_alignment, new_sizes)

        end subroutine allocate_bf_layer_N


        subroutine allocate_bf_layer_S(
     $       bf_nodes, interior_nodes,
     $       bf_grdpts_id,
     $       bf_alignment, final_alignment)

          implicit none

          real(rkind), dimension(:,:,:), allocatable, intent(inout) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in)    :: interior_nodes
          integer    , dimension(:,:)  , allocatable, intent(inout) :: bf_grdpts_id
          integer(ikind), dimension(2,2)            , intent(out)   :: bf_alignment
          integer(ikind), dimension(2,2)            , intent(in)    :: final_alignment


          integer(ikind), dimension(2) :: new_sizes
          integer(ikind)               :: i_min1
          integer(ikind)               :: i_min3
          integer(ikind)               :: i_min4
          integer(ikind)               :: outside_i_max1
          integer(ikind)               :: outside_i_max2
          integer(ikind)               :: interior_i_max1

          !compute the new sizes
          call get_new_sizes_and_alignment_S(
     $         final_alignment, new_sizes, bf_alignment)

          !get match indices to copy the tables
          call get_match(
     $         x_direction,
     $         bf_alignment,
     $         i_min1, i_min3, i_min4,
     $         outside_i_max1, outside_i_max2,
     $         interior_i_max1)

          !allocate the nodes
          call allocate_nodes_S(
     $         bf_nodes, interior_nodes,
     $         i_min1,
     $         interior_i_max1,
     $         bf_alignment, new_sizes)

          !allocate the grdpts_id
          call allocate_grdpts_id_S(
     $         bf_grdpts_id,
     $         outside_i_max1, outside_i_max2,
     $         interior_i_max1,
     $         i_min1, i_min3, i_min4,
     $         bf_alignment, new_sizes)

        end subroutine allocate_bf_layer_S


        subroutine allocate_bf_layer_E(
     $       bf_nodes, interior_nodes,
     $       bf_grdpts_id,
     $       bf_alignment, final_alignment)

          implicit none

          real(rkind), dimension(:,:,:), allocatable, intent(inout) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in)    :: interior_nodes
          integer    , dimension(:,:)  , allocatable, intent(inout) :: bf_grdpts_id
          integer(ikind), dimension(2,2)            , intent(out)   :: bf_alignment
          integer(ikind), dimension(2,2)            , intent(in)    :: final_alignment


          integer(ikind), dimension(2) :: new_sizes
          integer(ikind)               :: j_min1
          integer(ikind)               :: j_min3
          integer(ikind)               :: j_min4
          integer(ikind)               :: outside_j_max1
          integer(ikind)               :: outside_j_max2
          integer(ikind)               :: interior_j_max1

          !compute the new sizes
          call get_new_sizes_and_alignment_E(
     $         final_alignment, new_sizes, bf_alignment)

          !get match indices to copy the tables
          call get_match(
     $         y_direction,
     $         bf_alignment,
     $         j_min1, j_min3, j_min4,
     $         outside_j_max1, outside_j_max2,
     $         interior_j_max1)

          !allocate the nodes
          call allocate_nodes_E(
     $         bf_nodes, interior_nodes,
     $         j_min1,
     $         interior_j_max1,
     $         bf_alignment, new_sizes)

          !allocate the grdpts_id
          call allocate_grdpts_id_E(
     $         bf_grdpts_id,
     $         outside_j_max1, outside_j_max2,
     $         interior_j_max1,
     $         j_min1, j_min3, j_min4,
     $         bf_alignment, new_sizes)

        end subroutine allocate_bf_layer_E


        subroutine allocate_bf_layer_W(
     $       bf_nodes, interior_nodes,
     $       bf_grdpts_id,
     $       bf_alignment, final_alignment)

          implicit none

          real(rkind), dimension(:,:,:), allocatable, intent(inout) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in)    :: interior_nodes
          integer    , dimension(:,:)  , allocatable, intent(inout) :: bf_grdpts_id
          integer(ikind), dimension(2,2)            , intent(out)   :: bf_alignment
          integer(ikind), dimension(2,2)            , intent(in)    :: final_alignment


          integer(ikind), dimension(2) :: new_sizes
          integer(ikind)               :: j_min1
          integer(ikind)               :: j_min3
          integer(ikind)               :: j_min4
          integer(ikind)               :: outside_j_max1
          integer(ikind)               :: outside_j_max2
          integer(ikind)               :: interior_j_max1

          !compute the new sizes
          call get_new_sizes_and_alignment_W(
     $         final_alignment, new_sizes, bf_alignment)

          !get match indices to copy the tables
          call get_match(
     $         y_direction,
     $         bf_alignment,
     $         j_min1, j_min3, j_min4,
     $         outside_j_max1, outside_j_max2,
     $         interior_j_max1)

          !allocate the nodes
          call allocate_nodes_W(
     $         bf_nodes, interior_nodes,
     $         j_min1,
     $         interior_j_max1,
     $         bf_alignment, new_sizes)

          !allocate the grdpts_id
          call allocate_grdpts_id_W(
     $         bf_grdpts_id,
     $         outside_j_max1, outside_j_max2,
     $         interior_j_max1,
     $         j_min1, j_min3, j_min4,
     $         bf_alignment, new_sizes)
          
        end subroutine allocate_bf_layer_W


        subroutine allocate_nodes_N(
     $     bf_nodes, interior_nodes,
     $     i_min1,
     $     interior_i_max1,
     $     bf_alignment, new_sizes)

          implicit none

          real(rkind), dimension(:,:,:), allocatable, intent(out) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in)  :: interior_nodes
          integer(ikind)                            , intent(in)  :: i_min1
          integer(ikind)                            , intent(in)  :: interior_i_max1
          integer(ikind), dimension(2,2)            , intent(in)  :: bf_alignment
          integer(ikind), dimension(2)              , intent(in)  :: new_sizes



          integer(ikind) :: i_match, j_match
          integer(ikind) :: i,j
          integer        :: k

          allocate(bf_nodes(new_sizes(1), new_sizes(2), ne))

          !copy of grid points from the interior
          i_match = bf_alignment(1,1)-(bc_size+1)+i_min1
          j_match = ny-2*bc_size
          do k=1, ne
             do j=1, 2*bc_size
                do i=1, interior_i_max1
                   bf_nodes(i_min1+i,j,k) = interior_nodes(
     $                  i_match+i,j_match+j,k)
                end do
             end do
          end do

        end subroutine allocate_nodes_N


        subroutine allocate_nodes_S(
     $     bf_nodes, interior_nodes,
     $     i_min1,
     $     interior_i_max1,
     $     bf_alignment, new_sizes)

          implicit none

          real(rkind), dimension(:,:,:), allocatable, intent(out) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in)  :: interior_nodes
          integer(ikind)                            , intent(in)  :: i_min1
          integer(ikind)                            , intent(in)  :: interior_i_max1
          integer(ikind), dimension(2,2)            , intent(in)  :: bf_alignment
          integer(ikind), dimension(2)              , intent(in)  :: new_sizes

          integer(ikind) :: i_match, j_match
          integer(ikind) :: i,j
          integer        :: k

          allocate(bf_nodes(new_sizes(1), new_sizes(2), ne))

          !copy of grid points from the interior
          i_match = bf_alignment(1,1)-(bc_size+1)+i_min1
          j_match = 0
          do k=1, ne
             do j=1, 2*bc_size
                do i=1, interior_i_max1
                   bf_nodes(i_min1+i,j+1,k) = interior_nodes(i_match+i,j_match+j,k)
                end do
             end do
          end do

        end subroutine allocate_nodes_S


        subroutine allocate_nodes_E(
     $     bf_nodes, interior_nodes,
     $     j_min1,
     $     interior_j_max1,
     $     bf_alignment, new_sizes)

          implicit none

          real(rkind), dimension(:,:,:), allocatable, intent(out) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in)  :: interior_nodes
          integer(ikind)                            , intent(in)  :: j_min1
          integer(ikind)                            , intent(in)  :: interior_j_max1
          integer(ikind), dimension(2,2)            , intent(in)  :: bf_alignment
          integer(ikind), dimension(2)              , intent(in)  :: new_sizes

          integer(ikind) :: i_match, j_match
          integer(ikind) :: i,j
          integer        :: k

          allocate(bf_nodes(new_sizes(1), new_sizes(2), ne))

          !copy of grid points from the interior
          i_match = nx-2*bc_size
          j_match = bf_alignment(2,1)-(bc_size+1)+j_min1
          do k=1, ne
             do j=1, interior_j_max1
                do i=1, 2*bc_size
                   bf_nodes(i,j_min1+j,k) = interior_nodes(i_match+i,j_match+j,k)
                end do
             end do
          end do

        end subroutine allocate_nodes_E


        subroutine allocate_nodes_W(
     $     bf_nodes, interior_nodes,
     $     j_min1,
     $     interior_j_max1,
     $     bf_alignment, new_sizes)

          implicit none

          real(rkind), dimension(:,:,:), allocatable, intent(out) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in)  :: interior_nodes
          integer(ikind)                            , intent(in)  :: j_min1
          integer(ikind)                            , intent(in)  :: interior_j_max1
          integer(ikind), dimension(2,2)            , intent(in)  :: bf_alignment
          integer(ikind), dimension(2)              , intent(in)  :: new_sizes

          integer(ikind) :: i_match, j_match
          integer(ikind) :: i,j
          integer        :: k

          allocate(bf_nodes(new_sizes(1), new_sizes(2), ne))

          !copy of grid points from the interior
          i_match = 0
          j_match = bf_alignment(2,1)-(bc_size+1)+j_min1
          do k=1, ne
             do j=1, interior_j_max1
                do i=1, 2*bc_size
                   bf_nodes(i+1,j_min1+j,k) = interior_nodes(i_match+i,j_match+j,k)
                end do
             end do
          end do

        end subroutine allocate_nodes_W

      
        subroutine allocate_grdpts_id_N(
     $     bf_grdpts_id,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1,
     $     i_min1, i_min3, i_min4,
     $     bf_alignment, new_sizes)

          implicit none

          integer       , dimension(:,:), allocatable, intent(out):: bf_grdpts_id
          integer(ikind)                             , intent(in) :: outside_i_max1
          integer(ikind)                             , intent(in) :: outside_i_max2
          integer(ikind)                             , intent(in) :: interior_i_max1
          integer(ikind)                             , intent(in) :: i_min1
          integer(ikind)                             , intent(in) :: i_min3
          integer(ikind)                             , intent(in) :: i_min4
          integer(ikind), dimension(2,2)             , intent(in) :: bf_alignment
          integer(ikind), dimension(2)               , intent(in) :: new_sizes


          integer, dimension(bc_size,2*bc_size) :: border_W
          integer, dimension(bc_size,2*bc_size) :: border_E
          integer, dimension(2*bc_size)         :: interior_profile


          allocate(bf_grdpts_id(new_sizes(1), new_sizes(2)))

          !get the special border for the left and right
          !as well as the interior profile
          call get_additional_blocks_N(
     $         bf_alignment(1,1)-bc_size+i_min1,
     $         bf_alignment(1,1)-bc_size+i_min3,
     $         border_W, border_E, interior_profile)

          !blocks 1-5
          call add_grdpts_id_blocks_1_to_5_NS(
     $         bf_grdpts_id,
     $         border_W, border_E, interior_profile,
     $         outside_i_max1, outside_i_max2,
     $         interior_i_max1,
     $         i_min1, i_min3, i_min4,
     $         0)

          !block 6
          call add_grdpts_id_block_6_NS(
     $         bf_grdpts_id,
     $         2*bc_size+1, size(bf_grdpts_id,2))

        end subroutine allocate_grdpts_id_N


        subroutine allocate_grdpts_id_S(
     $     bf_grdpts_id,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1,
     $     i_min1, i_min3, i_min4,
     $     bf_alignment, new_sizes)

          implicit none

          integer       , dimension(:,:), allocatable, intent(out):: bf_grdpts_id
          integer(ikind)                             , intent(in) :: outside_i_max1
          integer(ikind)                             , intent(in) :: outside_i_max2
          integer(ikind)                             , intent(in) :: interior_i_max1
          integer(ikind)                             , intent(in) :: i_min1
          integer(ikind)                             , intent(in) :: i_min3
          integer(ikind)                             , intent(in) :: i_min4
          integer(ikind), dimension(2,2)             , intent(in) :: bf_alignment
          integer(ikind), dimension(2)               , intent(in) :: new_sizes


          integer, dimension(bc_size,2*bc_size) :: border_W
          integer, dimension(bc_size,2*bc_size) :: border_E
          integer, dimension(2*bc_size)         :: interior_profile


          allocate(bf_grdpts_id(new_sizes(1), new_sizes(2)))

          !get the special border for the left and right
          !as well as the interior profile
          call get_additional_blocks_S(
     $         bf_alignment(1,1)-bc_size+i_min1,
     $         bf_alignment(1,1)-bc_size+i_min3,
     $         border_W, border_E, interior_profile)

          !block 6
          call add_grdpts_id_block_6_NS(
     $         bf_grdpts_id,
     $         1, size(bf_grdpts_id,2)-2*bc_size)

          !blocks 1-5
          call add_grdpts_id_blocks_1_to_5_NS(
     $         bf_grdpts_id,
     $         border_W, border_E, interior_profile,
     $         outside_i_max1, outside_i_max2,
     $         interior_i_max1,
     $         i_min1, i_min3, i_min4,
     $         1)          

        end subroutine allocate_grdpts_id_S


        subroutine allocate_grdpts_id_E(
     $     bf_grdpts_id,
     $     outside_j_max1, outside_j_max2,
     $     interior_j_max1,
     $     j_min1, j_min3, j_min4,
     $     bf_alignment, new_sizes)

          implicit none

          integer       , dimension(:,:), allocatable, intent(out):: bf_grdpts_id
          integer(ikind)                             , intent(in) :: outside_j_max1
          integer(ikind)                             , intent(in) :: outside_j_max2
          integer(ikind)                             , intent(in) :: interior_j_max1
          integer(ikind)                             , intent(in) :: j_min1
          integer(ikind)                             , intent(in) :: j_min3
          integer(ikind)                             , intent(in) :: j_min4
          integer(ikind), dimension(2,2)             , intent(in) :: bf_alignment
          integer(ikind), dimension(2)               , intent(in) :: new_sizes


          integer, dimension(2*bc_size,bc_size) :: border_S
          integer, dimension(2*bc_size,bc_size) :: border_N
          integer, dimension(2*bc_size)         :: interior_profile


          allocate(bf_grdpts_id(new_sizes(1), new_sizes(2)))

          !get the special border for the left and right
          !as well as the interior profile
          call get_additional_blocks_E(
     $         bf_alignment(2,1)-bc_size+j_min1,
     $         bf_alignment(2,1)-bc_size+j_min3,
     $         border_S, border_N, interior_profile)

          !blocks 1-5
          call add_grdpts_id_blocks_1_to_5_E(
     $         bf_grdpts_id,
     $         border_S, border_N, interior_profile,
     $         outside_j_max1, outside_j_max2,
     $         interior_j_max1,
     $         j_min1, j_min3, j_min4)

        end subroutine allocate_grdpts_id_E


        subroutine allocate_grdpts_id_W(
     $     bf_grdpts_id,
     $     outside_j_max1, outside_j_max2,
     $     interior_j_max1,
     $     j_min1, j_min3, j_min4,
     $     bf_alignment, new_sizes)

          implicit none

          integer       , dimension(:,:), allocatable, intent(out):: bf_grdpts_id
          integer(ikind)                             , intent(in) :: outside_j_max1
          integer(ikind)                             , intent(in) :: outside_j_max2
          integer(ikind)                             , intent(in) :: interior_j_max1
          integer(ikind)                             , intent(in) :: j_min1
          integer(ikind)                             , intent(in) :: j_min3
          integer(ikind)                             , intent(in) :: j_min4
          integer(ikind), dimension(2,2)             , intent(in) :: bf_alignment
          integer(ikind), dimension(2)               , intent(in) :: new_sizes


          integer, dimension(2*bc_size,bc_size) :: border_S
          integer, dimension(2*bc_size,bc_size) :: border_N
          integer, dimension(2*bc_size)         :: interior_profile


          allocate(bf_grdpts_id(new_sizes(1), new_sizes(2)))

          !get the special border for the left and right
          !as well as the interior profile
          call get_additional_blocks_W(
     $         bf_alignment(2,1)-bc_size+j_min1,
     $         bf_alignment(2,1)-bc_size+j_min3,
     $         border_S, border_N, interior_profile)

          !blocks 1-5
          call add_grdpts_id_blocks_1_to_5_W(
     $         bf_grdpts_id,
     $         border_S, border_N, interior_profile,
     $         outside_j_max1, outside_j_max2,
     $         interior_j_max1,
     $         j_min1, j_min3, j_min4)

        end subroutine allocate_grdpts_id_W


        subroutine get_new_sizes_and_alignment_N(
     $     final_alignment,
     $     new_sizes, bf_alignment)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: final_alignment
          integer(ikind), dimension(2)  , intent(out) :: new_sizes
          integer(ikind), dimension(2,2), intent(out) :: bf_alignment
          
          !compute the new_sizes
          new_sizes(1) = final_alignment(1,2) - final_alignment(1,1) +
     $                   2*bc_size + 1
          new_sizes(2) = final_alignment(2,2) - (ny+1) +
     $                   2*bc_size + 1

          !compute the new alignment
          bf_alignment(1,1) = final_alignment(1,1)
          bf_alignment(2,1) = ny+1
          bf_alignment(1,2) = final_alignment(1,2)
          bf_alignment(2,2) = final_alignment(2,2)

        end subroutine get_new_sizes_and_alignment_N


        subroutine get_new_sizes_and_alignment_S(
     $     final_alignment,
     $     new_sizes, bf_alignment)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: final_alignment
          integer(ikind), dimension(2)  , intent(out) :: new_sizes
          integer(ikind), dimension(2,2), intent(out) :: bf_alignment
          
          !compute the new_sizes
          new_sizes(1) = final_alignment(1,2) - final_alignment(1,1) +
     $                   2*bc_size + 1
          new_sizes(2) = 0 - final_alignment(2,1) +
     $                   2*bc_size + 1

          !compute the new alignment
          bf_alignment(1,1) = final_alignment(1,1)
          bf_alignment(2,1) = final_alignment(2,1)
          bf_alignment(1,2) = final_alignment(1,2)
          bf_alignment(2,2) = 0

        end subroutine get_new_sizes_and_alignment_S


        subroutine get_new_sizes_and_alignment_E(
     $     final_alignment,
     $     new_sizes, bf_alignment)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: final_alignment
          integer(ikind), dimension(2)  , intent(out) :: new_sizes
          integer(ikind), dimension(2,2), intent(out) :: bf_alignment
          
          !compute the new_sizes
          new_sizes(1) = final_alignment(1,2) - (nx+1) +
     $                   2*bc_size + 1
          new_sizes(2) = final_alignment(2,2) - final_alignment(2,1) +
     $                   2*bc_size + 1

          !compute the new alignment
          bf_alignment(1,1) = nx+1
          bf_alignment(2,1) = final_alignment(2,1)
          bf_alignment(1,2) = final_alignment(1,2)
          bf_alignment(2,2) = final_alignment(2,2)

        end subroutine get_new_sizes_and_alignment_E


        subroutine get_new_sizes_and_alignment_W(
     $     final_alignment,
     $     new_sizes, bf_alignment)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: final_alignment
          integer(ikind), dimension(2)  , intent(out) :: new_sizes
          integer(ikind), dimension(2,2), intent(out) :: bf_alignment
          
          !compute the new_sizes
          new_sizes(1) = 0 - final_alignment(1,1) +
     $                   2*bc_size + 1
          new_sizes(2) = final_alignment(2,2) - final_alignment(2,1) +
     $                   2*bc_size + 1

          !compute the new alignment
          bf_alignment(1,1) = final_alignment(1,1)
          bf_alignment(2,1) = final_alignment(2,1)
          bf_alignment(1,2) = 0
          bf_alignment(2,2) = final_alignment(2,2)

        end subroutine get_new_sizes_and_alignment_W

      
         subroutine get_match(
     $     dir,
     $     bf_alignment,
     $     i_min1, i_min3, i_min4,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1)

           implicit none

           integer                       , intent(in)  :: dir
           integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
           integer(ikind)                , intent(out) :: i_min1
           integer(ikind)                , intent(out) :: i_min3
           integer(ikind)                , intent(out) :: i_min4
           integer(ikind)                , intent(out) :: outside_i_max1
           integer(ikind)                , intent(out) :: outside_i_max2
           integer(ikind)                , intent(out) :: interior_i_max1

           integer(ikind) :: ndir
           integer(ikind) :: min_border, max_border


           !detect the direction
           select case(dir)
             case(x_direction)
                ndir = nx
             case(y_direction)
                ndir = ny
             case default
                print '(''bf_layer_allocate_module'')'
                print '(''get_match'')'
                stop 'dir not recognized'
           end select


           !define the length of the block 1
           if(bf_alignment(dir,2).le.0) then
              outside_i_max1 = 2*bc_size + 1 +
     $                         bf_alignment(dir,2) - bf_alignment(dir,1) -
     $                         max(0, bf_alignment(dir,2)+bc_size)
           else
              outside_i_max1 = max(0,(bc_size+1) - bf_alignment(dir,1))
           end if


           !define the length of the block 2+3+4
           if(bf_alignment(dir,1).lt.(bc_size+1)) then
              min_border = 1
           else
              min_border = bf_alignment(dir,1)-bc_size
           end if
           if(bf_alignment(dir,2).gt.(ndir-bc_size)) then
              max_border = ndir
           else
              max_border = bf_alignment(dir,2)+bc_size
           end if
           interior_i_max1 = max_border-min_border+1

           if(bf_alignment(dir,2).le.0) then
              interior_i_max1 = max(0, bf_alignment(dir,2)+bc_size)
           end if
           if(bf_alignment(dir,1).ge.(ndir+1)) then
              interior_i_max1 = max(0, ndir+bc_size-bf_alignment(dir,1)+1)
           end if


           !define the length of the block 5
           if(bf_alignment(dir,1).ge.(ndir+1)) then
              outside_i_max2 = 2*bc_size + 1 + 
     $                         bf_alignment(dir,2) - bf_alignment(dir,1) -
     $                         interior_i_max1
           else
              outside_i_max2 = max(0, bf_alignment(dir,2)-(ndir-bc_size))
           end if


           !define the border indices
           i_min1 = outside_i_max1
           i_min3 = i_min1 +
     $              min(bc_size, interior_i_max1) +
     $              max(0, interior_i_max1 - 2*bc_size)
           i_min4 = i_min1 + interior_i_max1

         end subroutine get_match


         subroutine add_grdpts_id_blocks_1_to_5_NS(
     $     bf_grdpts_id,
     $     border_W, border_E, interior_profile,
     $     outside_i_max1, outside_i_max2,
     $     interior_i_max1,
     $     i_min1, i_min3, i_min4,
     $     j_match)

           implicit none

           integer      , dimension(:,:)              , intent(out):: bf_grdpts_id
           integer      , dimension(bc_size,2*bc_size), intent(in) :: border_W
           integer      , dimension(bc_size,2*bc_size), intent(in) :: border_E
           integer      , dimension(2*bc_size)        , intent(in) :: interior_profile
           integer(ikind)                             , intent(in) :: outside_i_max1
           integer(ikind)                             , intent(in) :: outside_i_max2
           integer(ikind)                             , intent(in) :: interior_i_max1
           integer(ikind)                             , intent(in) :: i_min1
           integer(ikind)                             , intent(in) :: i_min3
           integer(ikind)                             , intent(in) :: i_min4
           integer(ikind)                             , intent(in) :: j_match

           integer(ikind) :: i_start,i,j

           if((interior_i_max1-bc_size).gt.0) then
              i_start = max(bc_size, interior_i_max1-bc_size)
           else
              i_start = interior_i_max1
           end if


           do j=1, 2*bc_size

              !block 1
              do i=1, outside_i_max1
                 bf_grdpts_id(i,j_match+j) = no_pt
              end do

              !block 2
              do i=1, min(bc_size, interior_i_max1)
                 bf_grdpts_id(i_min1+i,j_match+j) = border_W(i,j)
              end do
                 
              !block 3
              do i=min(bc_size, interior_i_max1)+1, interior_i_max1-bc_size
                 bf_grdpts_id(i_min1+i,j_match+j) = interior_profile(j)
              end do

              !block 4
c$$$              print *, 'interior_i_max1', interior_i_max1
c$$$              print *, 'interior_i_max1-i_start', interior_i_max1-i_start
              do i=1, interior_i_max1-i_start
                 bf_grdpts_id(i_min3+i,j_match+j) = border_E(i,j)
              end do

              !block 5
              do i=1, outside_i_max2
                 bf_grdpts_id(i_min4+i,j_match+j) = no_pt
              end do
           end do

         end subroutine add_grdpts_id_blocks_1_to_5_NS


         subroutine add_grdpts_id_block_6_NS(
     $     bf_grdpts_id,
     $     j_min, j_max)

           implicit none

           integer      , dimension(:,:)              , intent(out):: bf_grdpts_id
           integer(ikind)                             , intent(in) :: j_min
           integer(ikind)                             , intent(in) :: j_max

           
           integer(ikind) :: i,j


           do j=j_min, j_max
              do i=1, size(bf_grdpts_id,1)
                 bf_grdpts_id(i,j) = no_pt
              end do
           end do

         end subroutine add_grdpts_id_block_6_NS


         subroutine add_grdpts_id_blocks_1_to_5_E(
     $     bf_grdpts_id,
     $     border_S, border_N, interior_profile,
     $     outside_j_max1, outside_j_max2,
     $     interior_j_max1,
     $     j_min1, j_min3, j_min4)

           implicit none

           integer      , dimension(:,:)              , intent(out):: bf_grdpts_id
           integer      , dimension(2*bc_size,bc_size), intent(in) :: border_S
           integer      , dimension(2*bc_size,bc_size), intent(in) :: border_N
           integer      , dimension(2*bc_size)        , intent(in) :: interior_profile
           integer(ikind)                             , intent(in) :: outside_j_max1
           integer(ikind)                             , intent(in) :: outside_j_max2
           integer(ikind)                             , intent(in) :: interior_j_max1
           integer(ikind)                             , intent(in) :: j_min1
           integer(ikind)                             , intent(in) :: j_min3
           integer(ikind)                             , intent(in) :: j_min4

           integer(ikind) :: i,j,j_start

           if((interior_j_max1-bc_size).gt.0) then
              j_start = max(bc_size,interior_j_max1-bc_size)
           else
              j_start = interior_j_max1
           end if

           !block 1 + partially 6
           do j=1, outside_j_max1
              do i=1, size(bf_grdpts_id,1)
                 bf_grdpts_id(i,j) = no_pt
              end do
           end do

           !block 2 + partially 6
           do j=1, min(bc_size, interior_j_max1)
              do i=1, 2*bc_size
                 bf_grdpts_id(i,j_min1+j) = border_S(i,j)
              end do
              do i=2*bc_size+1, size(bf_grdpts_id,1)
                 bf_grdpts_id(i,j_min1+j) = no_pt
              end do
           end do

           !block 3 + partially 6
           do j=min(bc_size, interior_j_max1)+1, interior_j_max1-bc_size
              do i=1, 2*bc_size
                 bf_grdpts_id(i,j_min1+j) = interior_profile(i)
              end do
              do i=2*bc_size+1, size(bf_grdpts_id,1)
                 bf_grdpts_id(i,j_min1+j) = no_pt
              end do
           end do


           !block 4 + partially 6
           do j=1, interior_j_max1-j_start
              do i=1, 2*bc_size
                 bf_grdpts_id(i,j_min3+j) = border_N(i,j)
              end do
              do i=2*bc_size+1, size(bf_grdpts_id,1)
                 bf_grdpts_id(i,j_min3+j) = no_pt
              end do
           end do

           !block 5 + partially 6
           do j=1, outside_j_max2
              do i=1, size(bf_grdpts_id,1)
                 bf_grdpts_id(i,j_min4+j) = no_pt
              end do
           end do

         end subroutine add_grdpts_id_blocks_1_to_5_E


        subroutine add_grdpts_id_blocks_1_to_5_W(
     $     bf_grdpts_id,
     $     border_S, border_N, interior_profile,
     $     outside_j_max1, outside_j_max2,
     $     interior_j_max1,
     $     j_min1, j_min3, j_min4)

           implicit none

           integer      , dimension(:,:)              , intent(out):: bf_grdpts_id
           integer      , dimension(2*bc_size,bc_size), intent(in) :: border_S
           integer      , dimension(2*bc_size,bc_size), intent(in) :: border_N
           integer      , dimension(2*bc_size)        , intent(in) :: interior_profile
           integer(ikind)                             , intent(in) :: outside_j_max1
           integer(ikind)                             , intent(in) :: outside_j_max2
           integer(ikind)                             , intent(in) :: interior_j_max1
           integer(ikind)                             , intent(in) :: j_min1
           integer(ikind)                             , intent(in) :: j_min3
           integer(ikind)                             , intent(in) :: j_min4

           integer(ikind) :: i,j,i_match, j_start

           i_match = size(bf_grdpts_id,1) - 2*bc_size

           if((interior_j_max1-bc_size).gt.0) then
              j_start = max(bc_size,interior_j_max1-bc_size)
           else
              j_start = interior_j_max1
           end if


           !block 1 + partially 6
           do j=1, outside_j_max1
              do i=1, size(bf_grdpts_id,1)
                 bf_grdpts_id(i,j) = no_pt
              end do
           end do

           !block 2 + partially 6
           do j=1, min(bc_size, interior_j_max1)
              do i=1, i_match
                  bf_grdpts_id(i,j_min1+j) = no_pt
              end do
              do i=1, 2*bc_size
                 bf_grdpts_id(i_match+i,j_min1+j) = border_S(i,j)
              end do             
           end do

           !block 3 + partially 6
           do j=min(bc_size, interior_j_max1)+1, interior_j_max1-bc_size
              do i=1, i_match
                  bf_grdpts_id(i,j_min1+j) = no_pt
              end do
              do i=1, 2*bc_size
                 bf_grdpts_id(i_match+i,j_min1+j) = interior_profile(i)
              end do              
           end do

           !block 4 + partially 6
           do j=1, interior_j_max1-j_start
              do i=1, i_match
                 bf_grdpts_id(i,j_min3+j) = no_pt
              end do
              do i=1, 2*bc_size
                 bf_grdpts_id(i_match+i,j_min3+j) = border_N(i,j)
              end do
           end do

           !block 5 + partially 6
           do j=1, outside_j_max2
              do i=1, size(bf_grdpts_id,1)
                 bf_grdpts_id(i,j_min4+j) = no_pt
              end do
           end do

         end subroutine add_grdpts_id_blocks_1_to_5_W


         subroutine get_additional_blocks_N(
     $     i_border_W, i_border_E,
     $     border_W, border_E, interior_profile)

          implicit none

          integer(ikind)                              , intent(in)  :: i_border_W
          integer(ikind)                              , intent(in)  :: i_border_E
          integer       , dimension(bc_size,2*bc_size), intent(out) :: border_W
          integer       , dimension(bc_size,2*bc_size), intent(out) :: border_E
          integer       , dimension(2*bc_size)        , intent(out) :: interior_profile

          integer :: i,j
          integer(ikind), dimension(2*bc_size,2*bc_size) :: corner_E
          integer(ikind), dimension(2*bc_size,2*bc_size) :: corner_W

          !corner_W
          do j=1, bc_size
             corner_W(1,j)=no_pt
             corner_W(2,j)=bc_pt
             corner_W(3,j)=bc_interior_pt
             corner_W(4,j)=interior_pt
          end do

          j=bc_size+1
          corner_W(1,j)=no_pt
          corner_W(2,j)=bc_pt
          corner_W(3,j)=bc_interior_pt
          corner_W(4,j)=bc_interior_pt

          j=bc_size+2
          corner_W(1,j)=no_pt
          corner_W(2,j)=bc_pt
          corner_W(3,j)=bc_pt
          corner_W(4,j)=bc_pt


          !corner_E
          do j=1, bc_size
             corner_E(1,j)=interior_pt
             corner_E(2,j)=bc_interior_pt
             corner_E(3,j)=bc_pt
             corner_E(4,j)=no_pt
          end do

          j=bc_size+1
          corner_E(1,j)=bc_interior_pt
          corner_E(2,j)=bc_interior_pt
          corner_E(3,j)=bc_pt
          corner_E(4,j)=no_pt

          j=bc_size+2
          corner_E(1,j)=bc_pt
          corner_E(2,j)=bc_pt
          corner_E(3,j)=bc_pt
          corner_E(4,j)=no_pt


          !interior profile
          do i=1, bc_size
             interior_profile(i) = interior_pt
          end do
          interior_profile(bc_size+1) = bc_interior_pt
          do i=bc_size+2, 2*bc_size
             interior_profile(i) = bc_pt
          end do

          !border_W
          call fill_border_NS(i_border_W, border_W, corner_W, corner_E, interior_profile)

          !border_E
          call fill_border_NS(i_border_E, border_E, corner_W, corner_E, interior_profile)

        end subroutine get_additional_blocks_N                


        subroutine get_additional_blocks_S(
     $     i_border_W, i_border_E,
     $     border_W, border_E, interior_profile)

          implicit none

          integer(ikind)                              , intent(in)  :: i_border_W
          integer(ikind)                              , intent(in)  :: i_border_E
          integer       , dimension(bc_size,2*bc_size), intent(out) :: border_W
          integer       , dimension(bc_size,2*bc_size), intent(out) :: border_E
          integer       , dimension(2*bc_size)        , intent(out) :: interior_profile

          integer :: i,j
          integer(ikind), dimension(2*bc_size,2*bc_size) :: corner_E
          integer(ikind), dimension(2*bc_size,2*bc_size) :: corner_W


          !corner_W
          j=1
          corner_W(1,j)=no_pt
          corner_W(2,j)=bc_pt
          corner_W(3,j)=bc_pt
          corner_W(4,j)=bc_pt

          j=bc_size
          corner_W(1,j)=no_pt
          corner_W(2,j)=bc_pt
          corner_W(3,j)=bc_interior_pt
          corner_W(4,j)=bc_interior_pt

          do j=bc_size+1, 2*bc_size
             corner_W(1,j)=no_pt
             corner_W(2,j)=bc_pt
             corner_W(3,j)=bc_interior_pt
             corner_W(4,j)=interior_pt
          end do


          !corner_E
          j=1
          corner_E(1,j)=bc_pt
          corner_E(2,j)=bc_pt
          corner_E(3,j)=bc_pt
          corner_E(4,j)=no_pt

          j=bc_size
          corner_E(1,j)=bc_interior_pt
          corner_E(2,j)=bc_interior_pt
          corner_E(3,j)=bc_pt
          corner_E(4,j)=no_pt

          do j=bc_size+1, 2*bc_size
             corner_E(1,j)=interior_pt
             corner_E(2,j)=bc_interior_pt
             corner_E(3,j)=bc_pt
             corner_E(4,j)=no_pt
          end do


          !interior profile
          do i=1, bc_size-1
             interior_profile(i) = bc_pt
          end do
          interior_profile(bc_size) = bc_interior_pt
          do i=bc_size+1, 2*bc_size
             interior_profile(i) = interior_pt
          end do          
     

          !border_W
          call fill_border_NS(
     $         i_border_W, border_W,
     $         corner_W, corner_E, interior_profile)

          !border_E
          call fill_border_NS(
     $         i_border_E, border_E,
     $         corner_W, corner_E, interior_profile)

        end subroutine get_additional_blocks_S


        subroutine get_additional_blocks_E(
     $     j_border_S, j_border_N,
     $     border_S, border_N, interior_profile)

          implicit none

          integer(ikind)                              , intent(in)  :: j_border_S
          integer(ikind)                              , intent(in)  :: j_border_N          
          integer       , dimension(2*bc_size,bc_size), intent(out) :: border_S
          integer       , dimension(2*bc_size,bc_size), intent(out) :: border_N
          integer       , dimension(2*bc_size)        , intent(out) :: interior_profile

          integer :: i,j
          integer(ikind), dimension(2*bc_size,2*bc_size) :: corner_S
          integer(ikind), dimension(2*bc_size,2*bc_size) :: corner_N


          !corner_S
          j=1
          do i=1, 2*bc_size
             corner_S(i,j)=no_pt
          end do

          j=bc_size
          do i=1, 2*bc_size
             corner_S(i,j) = bc_pt
          end do

          j=bc_size+1
          corner_S(1,j)=bc_interior_pt
          corner_S(2,j)=bc_interior_pt
          corner_S(3,j)=bc_interior_pt
          corner_S(4,j)=bc_pt

          j=bc_size+2
          corner_S(1,j)=interior_pt
          corner_S(2,j)=interior_pt
          corner_S(3,j)=bc_interior_pt
          corner_S(4,j)=bc_pt          


          !corner_N
          j=1
          corner_N(1,j)=interior_pt
          corner_N(2,j)=interior_pt
          corner_N(3,j)=bc_interior_pt
          corner_N(4,j)=bc_pt

          j=bc_size
          corner_N(1,j)=bc_interior_pt
          corner_N(2,j)=bc_interior_pt
          corner_N(3,j)=bc_interior_pt
          corner_N(4,j)=bc_pt

          j=bc_size+1
          do i=1, 2*bc_size
             corner_N(i,j)=bc_pt
          end do

          j=bc_size+2
          do i=1, 2*bc_size
             corner_N(i,j)=no_pt
          end do


          !interior profile
          do i=1, bc_size
             interior_profile(i) = interior_pt
          end do
          interior_profile(bc_size+1) = bc_interior_pt
          do i=bc_size+2, 2*bc_size
             interior_profile(i) = bc_pt
          end do


          !border_S
          call fill_border_EW(
     $         j_border_S, border_S,
     $         corner_S, corner_N, interior_profile)


          !border_N
          call fill_border_EW(
     $         j_border_N, border_N,
     $         corner_S, corner_N, interior_profile)

        end subroutine get_additional_blocks_E


        !< compute the additional block W for Fig.4
        subroutine get_additional_blocks_W(
     $     j_border_S, j_border_N,
     $     border_S, border_N, interior_profile)

          implicit none

          integer(ikind)                              , intent(in)  :: j_border_S
          integer(ikind)                              , intent(in)  :: j_border_N
          integer       , dimension(2*bc_size,bc_size), intent(out) :: border_S
          integer       , dimension(2*bc_size,bc_size), intent(out) :: border_N
          integer       , dimension(2*bc_size)        , intent(out) :: interior_profile

          integer :: i,j
          integer(ikind), dimension(2*bc_size,2*bc_size) :: corner_S
          integer(ikind), dimension(2*bc_size,2*bc_size) :: corner_N

          
          !corner_S
          j=1
          do i=1, 2*bc_size
             corner_S(i,j)=no_pt
          end do

          j=bc_size
          do i=1, 2*bc_size
             corner_S(i,j) = bc_pt
          end do

          j=bc_size+1
          corner_S(1,j)=bc_pt
          corner_S(2,j)=bc_interior_pt
          corner_S(3,j)=bc_interior_pt
          corner_S(4,j)=bc_interior_pt

          j=bc_size+2
          corner_S(1,j)=bc_pt
          corner_S(2,j)=bc_interior_pt
          corner_S(3,j)=interior_pt
          corner_S(4,j)=interior_pt


          !corner_N
          j=1
          corner_N(1,j)=bc_pt
          corner_N(2,j)=bc_interior_pt
          corner_N(3,j)=interior_pt
          corner_N(4,j)=interior_pt

          j=bc_size
          corner_N(1,j)=bc_pt
          corner_N(2,j)=bc_interior_pt
          corner_N(3,j)=bc_interior_pt
          corner_N(4,j)=bc_interior_pt

          j=bc_size+1
          do i=1, 2*bc_size
             corner_N(i,j)=bc_pt
          end do

          j=bc_size+2
          do i=1, 2*bc_size
             corner_N(i,j)=no_pt
          end do



          !interior profile
          do i=1, bc_size-1
             interior_profile(i) = bc_pt
          end do
          i=bc_size
          interior_profile(bc_size) = bc_interior_pt
          do i=bc_size+1, 2*bc_size
             interior_profile(i) = interior_pt
          end do

          !border_S
          call fill_border_EW(
     $         j_border_S, border_S,
     $         corner_S, corner_N, interior_profile)


          !border_N
          call fill_border_EW(
     $         j_border_N, border_N,
     $         corner_S, corner_N, interior_profile)

        end subroutine get_additional_blocks_W


        subroutine fill_border_NS(i_border, border, corner_W, corner_E, interior_profile)

          implicit none

          integer(ikind)                                , intent(in) :: i_border
          integer(ikind), dimension(bc_size,2*bc_size)  , intent(out):: border
          integer(ikind), dimension(2*bc_size,2*bc_size), intent(in) :: corner_W
          integer(ikind), dimension(2*bc_size,2*bc_size), intent(in) :: corner_E
          integer(ikind), dimension(2*bc_size)          , intent(in) :: interior_profile


          integer(ikind) :: i,j


          if((i_border.ge.0).and.(i_border.le.bc_size)) then
             do j=1, 2*bc_size
                do i=1, bc_size
                   border(i,j) = corner_W(i+i_border,j)
                end do
             end do

          else
             if((i_border.ge.(nx-bc_size)).and.(i_border.le.nx)) then
                do j=1, 2*bc_size
                   do i=1, bc_size
                      border(i,j) = corner_E(i_border-(nx-bc_size)+i,j)
                   end do
                end do
             else
                do j=1, 2*bc_size
                   do i=1, bc_size
                      border(i,j) = interior_profile(j)
                   end do
                end do
             end if
          end if

        end subroutine fill_border_NS


        subroutine fill_border_EW(j_border, border, corner_S, corner_N, interior_profile)

          implicit none

          integer(ikind)                                , intent(in) :: j_border
          integer(ikind), dimension(2*bc_size,bc_size)  , intent(out):: border
          integer(ikind), dimension(2*bc_size,2*bc_size), intent(in) :: corner_S
          integer(ikind), dimension(2*bc_size,2*bc_size), intent(in) :: corner_N
          integer(ikind), dimension(2*bc_size)          , intent(in) :: interior_profile


          integer(ikind) :: i,j


          if((j_border.ge.0).and.(j_border.le.bc_size)) then
             do j=1, bc_size
                do i=1, 2*bc_size
                   border(i,j) = corner_S(i,j+j_border)
                end do
             end do

          else
             if((j_border.ge.(ny-bc_size)).and.(j_border.le.ny)) then
                do j=1, bc_size
                   do i=1, 2*bc_size
                      border(i,j) = corner_N(i,j_border-(ny-bc_size)+j)
                   end do
                end do
             else
                do j=1, bc_size
                   do i=1, 2*bc_size
                      border(i,j) = interior_profile(i)
                   end do
                end do
             end if
          end if

        end subroutine fill_border_EW        

      end module bf_layer_allocate_module
