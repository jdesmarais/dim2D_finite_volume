      !> @file
      !> module encapsulating the subroutines needed to allocate 
      !> the attributes of a bf_layer object
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating the subroutines needed to allocate 
      !> the attributes of a bf_layer object
      !> \image html  bf_layer_allocate_module.png
      !> \image latex bf_layer_allocate_module.eps
      !
      !> @warning
      !> the x_map and y_map are not correcly initialized if the
      !> buffer layer does not have any grid point in common with
      !> the interior domain
      !
      !> @date
      ! 27_06_2014 - documentation update        - J.L. Desmarais
      ! 23_10_2014 - addition of x_map and y_map - J.L. Desmarais
      !-----------------------------------------------------------------
      module bf_layer_allocate_module

        use parameters_bf_layer, only :
     $     no_pt, interior_pt,
     $     bc_interior_pt, bc_pt,
     $     align_N, align_S, align_E, align_W

        use parameters_constant, only :
     $       N,S,E,W,
     $       x_direction,
     $       y_direction

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none


        private
        public :: allocate_bf_layer_N,
     $            allocate_bf_layer_S,
     $            allocate_bf_layer_E,
     $            allocate_bf_layer_W,
     $            get_additional_blocks_N,
     $            get_additional_blocks_S,
     $            get_additional_blocks_E,
     $            get_additional_blocks_W,
     $            get_match_interior


        contains


        subroutine allocate_bf_layer_N(
     $       bf_x_map, interior_x_map,
     $       bf_y_map, interior_y_map,
     $       bf_nodes, interior_nodes,
     $       bf_grdpts_id,
     $       bf_alignment, final_alignment)

          implicit none

          real(rkind), dimension(:)    , allocatable, intent(inout) :: bf_x_map
          real(rkind), dimension(nx)                , intent(in)    :: interior_x_map
          real(rkind), dimension(:)    , allocatable, intent(inout) :: bf_y_map
          real(rkind), dimension(ny)                , intent(in)    :: interior_y_map
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
     $         bf_x_map, interior_x_map,
     $         bf_y_map, interior_y_map,
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
     $       bf_x_map, interior_x_map,
     $       bf_y_map, interior_y_map,
     $       bf_nodes, interior_nodes,
     $       bf_grdpts_id,
     $       bf_alignment, final_alignment)

          implicit none

          real(rkind), dimension(:)    , allocatable, intent(inout) :: bf_x_map
          real(rkind), dimension(nx)                , intent(in)    :: interior_x_map
          real(rkind), dimension(:)    , allocatable, intent(inout) :: bf_y_map
          real(rkind), dimension(ny)                , intent(in)    :: interior_y_map
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
     $         bf_x_map, interior_x_map,
     $         bf_y_map, interior_y_map,
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
     $     bf_x_map, interior_x_map,
     $     bf_y_map, interior_y_map,
     $     bf_nodes, interior_nodes,
     $     bf_grdpts_id,
     $     bf_alignment, final_alignment)

          implicit none

          real(rkind), dimension(:)    , allocatable, intent(inout) :: bf_x_map
          real(rkind), dimension(nx)                , intent(in)    :: interior_x_map
          real(rkind), dimension(:)    , allocatable, intent(inout) :: bf_y_map
          real(rkind), dimension(ny)                , intent(in)    :: interior_y_map
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
     $         bf_x_map, interior_x_map,
     $         bf_y_map, interior_y_map,
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
     $     bf_x_map, interior_x_map,
     $     bf_y_map, interior_y_map,
     $     bf_nodes, interior_nodes,
     $     bf_grdpts_id,
     $     bf_alignment, final_alignment)

          implicit none

          real(rkind), dimension(:)    , allocatable, intent(inout) :: bf_x_map
          real(rkind), dimension(nx)                , intent(in)    :: interior_x_map
          real(rkind), dimension(:)    , allocatable, intent(inout) :: bf_y_map
          real(rkind), dimension(ny)                , intent(in)    :: interior_y_map
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
     $         bf_x_map, interior_x_map,
     $         bf_y_map, interior_y_map,
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
     $     bf_x_map, interior_x_map,
     $     bf_y_map, interior_y_map,
     $     bf_nodes, interior_nodes,
     $     i_min1,
     $     interior_i_max1,
     $     bf_alignment, new_sizes)

          implicit none

          real(rkind), dimension(:)    , allocatable, intent(out) :: bf_x_map
          real(rkind), dimension(nx)                , intent(in)  :: interior_x_map
          real(rkind), dimension(:)    , allocatable, intent(out) :: bf_y_map
          real(rkind), dimension(ny)                , intent(in)  :: interior_y_map
          real(rkind), dimension(:,:,:), allocatable, intent(out) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in)  :: interior_nodes
          integer(ikind)                            , intent(in)  :: i_min1
          integer(ikind)                            , intent(in)  :: interior_i_max1
          integer(ikind), dimension(2,2)            , intent(in)  :: bf_alignment
          integer(ikind), dimension(2)              , intent(in)  :: new_sizes


          integer(ikind) :: i_match, j_match
          integer(ikind) :: i,j
          integer        :: k

          allocate(bf_x_map(new_sizes(1)))
          allocate(bf_y_map(new_sizes(2)))
          allocate(bf_nodes(new_sizes(1), new_sizes(2), ne))

          !match indices
          i_match = bf_alignment(1,1)-(bc_size+1)+i_min1
          j_match = ny-2*bc_size

          !x_map copy
          call create_map_from_interior(
     $         bf_x_map, interior_x_map,
     $         i_min1, interior_i_max1, i_match)
          
          !y_map copy
          call create_map_right(
     $         bf_y_map, interior_y_map,
     $         j_match)         

          !copy of grid points from the interior
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
     $     bf_x_map, interior_x_map,
     $     bf_y_map, interior_y_map,
     $     bf_nodes, interior_nodes,
     $     i_min1,
     $     interior_i_max1,
     $     bf_alignment, new_sizes)

          implicit none

          real(rkind), dimension(:)    , allocatable, intent(out) :: bf_x_map
          real(rkind), dimension(nx)                , intent(in)  :: interior_x_map
          real(rkind), dimension(:)    , allocatable, intent(out) :: bf_y_map
          real(rkind), dimension(ny)                , intent(in)  :: interior_y_map
          real(rkind), dimension(:,:,:), allocatable, intent(out) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in)  :: interior_nodes
          integer(ikind)                            , intent(in)  :: i_min1
          integer(ikind)                            , intent(in)  :: interior_i_max1
          integer(ikind), dimension(2,2)            , intent(in)  :: bf_alignment
          integer(ikind), dimension(2)              , intent(in)  :: new_sizes

          integer(ikind) :: i_match, j_match
          integer(ikind) :: i,j
          integer        :: k

          allocate(bf_x_map(new_sizes(1)))
          allocate(bf_y_map(new_sizes(2)))
          allocate(bf_nodes(new_sizes(1), new_sizes(2), ne))

          !match indices
          i_match = bf_alignment(1,1)-(bc_size+1)+i_min1
          j_match = 0

          !x_map copy
          call create_map_from_interior(
     $         bf_x_map, interior_x_map,
     $         i_min1, interior_i_max1, i_match)

          !y_map copy
          call create_map_left(
     $         bf_y_map, interior_y_map,
     $         j_match)

          !copy of grid points from the interior
          do k=1, ne
             do j=1, 2*bc_size
                do i=1, interior_i_max1
                   bf_nodes(i_min1+i,j+size(bf_nodes,2)-2*bc_size,k) =
     $                  interior_nodes(i_match+i,j_match+j,k)
                end do
             end do
          end do

        end subroutine allocate_nodes_S


        subroutine allocate_nodes_E(
     $     bf_x_map, interior_x_map,
     $     bf_y_map, interior_y_map,
     $     bf_nodes, interior_nodes,
     $     j_min1,
     $     interior_j_max1,
     $     bf_alignment, new_sizes)

          implicit none

          real(rkind), dimension(:)    , allocatable, intent(out) :: bf_x_map
          real(rkind), dimension(nx)                , intent(in)  :: interior_x_map
          real(rkind), dimension(:)    , allocatable, intent(out) :: bf_y_map
          real(rkind), dimension(ny)                , intent(in)  :: interior_y_map
          real(rkind), dimension(:,:,:), allocatable, intent(out) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in)  :: interior_nodes
          integer(ikind)                            , intent(in)  :: j_min1
          integer(ikind)                            , intent(in)  :: interior_j_max1
          integer(ikind), dimension(2,2)            , intent(in)  :: bf_alignment
          integer(ikind), dimension(2)              , intent(in)  :: new_sizes

          integer(ikind) :: i_match, j_match
          integer(ikind) :: i,j
          integer        :: k

          allocate(bf_x_map(new_sizes(1)))
          allocate(bf_y_map(new_sizes(2)))
          allocate(bf_nodes(new_sizes(1), new_sizes(2), ne))

          !match indices
          i_match = nx-2*bc_size
          j_match = bf_alignment(2,1)-(bc_size+1)+j_min1

          !x_map copy
          call create_map_right(
     $         bf_x_map, interior_x_map,
     $         i_match)

          !y_map copy
          call create_map_from_interior(
     $         bf_y_map, interior_y_map,
     $         j_min1, interior_j_max1, j_match)

          !copy of grid points from the interior
          do k=1, ne
             do j=1, interior_j_max1
                do i=1, 2*bc_size
                   bf_nodes(i,j_min1+j,k) = interior_nodes(i_match+i,j_match+j,k)
                end do
             end do
          end do

        end subroutine allocate_nodes_E


        subroutine allocate_nodes_W(
     $     bf_x_map, interior_x_map,
     $     bf_y_map, interior_y_map,
     $     bf_nodes, interior_nodes,
     $     j_min1,
     $     interior_j_max1,
     $     bf_alignment, new_sizes)

          implicit none

          real(rkind), dimension(:)    , allocatable, intent(out) :: bf_x_map
          real(rkind), dimension(nx)                , intent(in)  :: interior_x_map
          real(rkind), dimension(:)    , allocatable, intent(out) :: bf_y_map
          real(rkind), dimension(ny)                , intent(in)  :: interior_y_map
          real(rkind), dimension(:,:,:), allocatable, intent(out) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in)  :: interior_nodes
          integer(ikind)                            , intent(in)  :: j_min1
          integer(ikind)                            , intent(in)  :: interior_j_max1
          integer(ikind), dimension(2,2)            , intent(in)  :: bf_alignment
          integer(ikind), dimension(2)              , intent(in)  :: new_sizes

          integer(ikind) :: i_match, j_match
          integer(ikind) :: i,j
          integer        :: k

          allocate(bf_x_map(new_sizes(1)))
          allocate(bf_y_map(new_sizes(2)))
          allocate(bf_nodes(new_sizes(1), new_sizes(2), ne))

          !match indices
          i_match = 0
          j_match = bf_alignment(2,1)-(bc_size+1)+j_min1

          !x_map copy
          call create_map_left(
     $         bf_x_map, interior_x_map,
     $         i_match)

          !y_map copy
          call create_map_from_interior(
     $         bf_y_map, interior_y_map,
     $         j_min1, interior_j_max1, j_match)

          !copy of grid points from the interior
          do k=1, ne
             do j=1, interior_j_max1
                do i=1, 2*bc_size
                   bf_nodes(i+size(bf_nodes,1)-2*bc_size,j_min1+j,k) =
     $                  interior_nodes(i_match+i,j_match+j,k)
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
     $         size(bf_grdpts_id,2)-2*bc_size)       

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
          new_sizes(2) = final_alignment(2,2) - align_N +
     $                   2*bc_size + 1

          !compute the new alignment
          bf_alignment(1,1) = final_alignment(1,1)
          bf_alignment(2,1) = align_N
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
          new_sizes(2) = align_S - final_alignment(2,1) +
     $                   2*bc_size + 1

          !compute the new alignment
          bf_alignment(1,1) = final_alignment(1,1)
          bf_alignment(2,1) = final_alignment(2,1)
          bf_alignment(1,2) = final_alignment(1,2)
          bf_alignment(2,2) = align_S

        end subroutine get_new_sizes_and_alignment_S


        subroutine get_new_sizes_and_alignment_E(
     $     final_alignment,
     $     new_sizes, bf_alignment)

          implicit none

          integer(ikind), dimension(2,2), intent(in)  :: final_alignment
          integer(ikind), dimension(2)  , intent(out) :: new_sizes
          integer(ikind), dimension(2,2), intent(out) :: bf_alignment
          
          !compute the new_sizes
          new_sizes(1) = final_alignment(1,2) - align_E +
     $                   2*bc_size + 1
          new_sizes(2) = final_alignment(2,2) - final_alignment(2,1) +
     $                   2*bc_size + 1

          !compute the new alignment
          bf_alignment(1,1) = align_E
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
          new_sizes(1) = align_W - final_alignment(1,1) +
     $                   2*bc_size + 1
          new_sizes(2) = final_alignment(2,2) - final_alignment(2,1) +
     $                   2*bc_size + 1

          !compute the new alignment
          bf_alignment(1,1) = final_alignment(1,1)
          bf_alignment(2,1) = final_alignment(2,1)
          bf_alignment(1,2) = align_W
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


         !< decide the size of the sub layers inside
         !> an interior layer to know whether a special
         !> treatment for the corner grid points id is needed
         !> for example when dividing the first interior layer
         !> into 2+3+4 for the reallocation layers
         subroutine get_match_interior(
     $     dir, ndir,
     $     bf_alignment, i_start_interior, size_interior,
     $     size_layer1, size_layer2)

           implicit none
           
           integer                       , intent(in) :: dir
           integer(ikind)                , intent(in) :: ndir
           integer(ikind), dimension(2,2), intent(in) :: bf_alignment
           integer(ikind)                , intent(in) :: i_start_interior
           integer(ikind)                , intent(in) :: size_interior
           integer(ikind)                , intent(out):: size_layer1
           integer(ikind)                , intent(out):: size_layer2
           
           
           integer(ikind) :: layer_coord
           
           
           layer_coord = bf_alignment(dir,1)-bc_size+
     $          i_start_interior
           if((layer_coord.ge.1).and.(layer_coord.le.bc_size)) then
              size_layer1 = min(bc_size,size_interior)
           else
              size_layer1 = 0
           end if
           
           layer_coord = bf_alignment(dir,1)-(bc_size+1)+
     $          i_start_interior+
     $          size_interior
           if((layer_coord.ge.(ndir-1)).and.(layer_coord.le.ndir)) then
              size_layer2 = min(bc_size,size_interior)
           else
              size_layer2 = 0
           end if

         end subroutine get_match_interior


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

        subroutine create_map_from_interior(
     $     bf_map, interior_map,
     $     i_min1, interior_i_max1, i_match)

          implicit none

          real(rkind), dimension(:), intent(out) :: bf_map
          real(rkind), dimension(:), intent(in)  :: interior_map
          integer(ikind)           , intent(in)  :: i_min1
          integer(ikind)           , intent(in)  :: interior_i_max1
          integer(ikind)           , intent(in)  :: i_match

          integer(ikind) :: i
          integer(ikind) :: i_left
          integer(ikind) :: i_right
          real(rkind)    :: dx
          

          !x_map outside domain: left
          i_left = i_match+1
          dx = interior_map(i_left+1) - interior_map(i_left)
          do i=1, i_min1
             bf_map(i) = interior_map(i_left) - (i_min1+1-i)*dx
          end do

          !copy of x_map from the interior
          do i=1, interior_i_max1
             bf_map(i_min1+i) = interior_map(i_match+i)
          end do

          !x_map outside domain: right
          i_right = i_match+interior_i_max1
          dx = interior_map(i_right)-interior_map(i_right-1)
          do i=i_min1+interior_i_max1+1,size(bf_map,1)
             bf_map(i) = interior_map(i_right) + (i-(i_min1+interior_i_max1))*dx
          end do

        end subroutine create_map_from_interior


        subroutine create_map_right(
     $     bf_y_map, interior_y_map,
     $     j_match)

          implicit none
          
          real(rkind), dimension(:), intent(out) :: bf_y_map
          real(rkind), dimension(:), intent(in)  :: interior_y_map
          integer(ikind)           , intent(in)  :: j_match

          integer(ikind) :: j
          real(rkind)    :: dy

          !copy of y_map from the interior
          do j=1, 2*bc_size
             bf_y_map(j) = interior_y_map(j_match+j)
          end do

          !creation of y_map outside
          dy = interior_y_map(size(interior_y_map,1))-
     $         interior_y_map(size(interior_y_map,1)-1)
          do j=2*bc_size+1, size(bf_y_map,1)
             bf_y_map(j) = bf_y_map(2*bc_size) + (j-2*bc_size)*dy
          end do

        end subroutine create_map_right


        subroutine create_map_left(
     $     bf_y_map, interior_y_map,
     $     j_match)

          implicit none
          
          real(rkind), dimension(:), intent(out) :: bf_y_map
          real(rkind), dimension(:), intent(in)  :: interior_y_map
          integer(ikind)           , intent(in)  :: j_match

          integer(ikind) :: j
          real(rkind)    :: dy

          !create map outside
          dy = interior_y_map(2) - interior_y_map(1)
          do j=1, size(bf_y_map,1)-2*bc_size
             bf_y_map(j) = -(size(bf_y_map,1)-2*bc_size+1-j)*dy + interior_y_map(j_match+1)
          end do

          !copy from inside
          do j=1, 2*bc_size
             bf_y_map(j+size(bf_y_map,1)-2*bc_size) = interior_y_map(j_match+j)
          end do

        end subroutine create_map_left

      end module bf_layer_allocate_module
