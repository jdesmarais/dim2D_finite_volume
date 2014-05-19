      module bf_layer_allocate_module

        use parameters_bf_layer, only : no_pt, interior_pt,
     $                                  bc_interior_pt, bc_pt
        use parameters_constant, only : N,S,E,W
        use parameters_input   , only : nx,ny,ne,bc_size
        use parameters_kind    , only : ikind, rkind     

        implicit none


        private
        public :: allocate_bf_layer_N,
     $            allocate_bf_layer_S,
     $            allocate_bf_layer_E,
     $            allocate_bf_layer_W

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


          !compute the new sizes
          call get_new_sizes_and_alignment_N(
     $         final_alignment, new_sizes, bf_alignment)

          !allocate the nodes
          call allocate_nodes_N(
     $         bf_alignment, new_sizes, bf_nodes, interior_nodes)

          !allocate the grdpts_id
          call allocate_grdpts_id_N(
     $         bf_alignment, new_sizes, bf_grdpts_id)

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


          !compute the new sizes
          call get_new_sizes_and_alignment_S(
     $         final_alignment, new_sizes, bf_alignment)

          !allocate the nodes
          call allocate_nodes_S(
     $         bf_alignment, new_sizes, bf_nodes, interior_nodes)

          !allocate the grdpts_id
          call allocate_grdpts_id_S(
     $         bf_alignment, new_sizes, bf_grdpts_id)

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


          !compute the new sizes
          call get_new_sizes_and_alignment_E(
     $         final_alignment, new_sizes, bf_alignment)

          !allocate the nodes
          call allocate_nodes_E(
     $         bf_alignment, new_sizes, bf_nodes, interior_nodes)

          !allocate the grdpts_id
          call allocate_grdpts_id_E(
     $         bf_alignment, new_sizes, bf_grdpts_id)

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


          !compute the new sizes
          call get_new_sizes_and_alignment_W(
     $         final_alignment, new_sizes, bf_alignment)

          !allocate the nodes
          call allocate_nodes_W(
     $         bf_alignment, new_sizes, bf_nodes, interior_nodes)

          !allocate the grdpts_id
          call allocate_grdpts_id_W(
     $         bf_alignment, new_sizes, bf_grdpts_id)

        end subroutine allocate_bf_layer_W


        subroutine allocate_nodes_N(
     $     bf_alignment, new_size, bf_nodes, interior_nodes)

          implicit none

          integer(ikind), dimension(2,2)            , intent(in) :: bf_alignment
          integer(ikind), dimension(2)              , intent(in) :: new_size
          real(rkind), dimension(:,:,:), allocatable, intent(out) :: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in) :: interior_nodes

          integer(ikind) :: i_match, j_match
          integer(ikind) :: i,j,k

          allocate(bf_nodes(new_size(1), new_size(2), ne))

          !copy of grid points from the interior
          i_match = bf_alignment(1,1)-(bc_size+1)
          j_match = ny-2*bc_size
          do k=1, ne
             do j=1, 2*bc_size
                do i=1, size(bf_nodes,1)
                   bf_nodes(i,j,k) = interior_nodes(
     $                  i_match+i,j_match+j,k)
                end do
             end do
          end do

        end subroutine allocate_nodes_N


        subroutine allocate_nodes_S(
     $     bf_alignment, new_size, bf_nodes, interior_nodes)

          implicit none

          integer(ikind), dimension(2,2)            , intent(in) :: bf_alignment
          integer(ikind), dimension(2)              , intent(in) :: new_size
          real(rkind), dimension(:,:,:), allocatable, intent(out):: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in) :: interior_nodes

          integer(ikind) :: i_match, j_match
          integer(ikind) :: i,j,k

          allocate(bf_nodes(new_size(1), new_size(2), ne))

          !copy of grid points from the interior
          i_match = bf_alignment(1,1)-(bc_size+1)
          j_match = 0
          do k=1, ne
             do j=1, 2*bc_size
                do i=1, size(bf_nodes,1)
                   bf_nodes(i,j+1,k) = interior_nodes(i_match+i,j_match+j,k)
                end do
             end do
          end do

        end subroutine allocate_nodes_S


        subroutine allocate_nodes_E(
     $     bf_alignment, new_size, bf_nodes, interior_nodes)

          implicit none

          integer(ikind), dimension(2,2)            , intent(in) :: bf_alignment
          integer(ikind), dimension(2)              , intent(in) :: new_size
          real(rkind), dimension(:,:,:), allocatable, intent(out):: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in) :: interior_nodes

          integer(ikind) :: i_match, j_match
          integer(ikind) :: i,j,k

          allocate(bf_nodes(new_size(1), new_size(2), ne))

          !copy of grid points from the interior
          i_match = nx-2*bc_size
          j_match = bf_alignment(2,1)-(bc_size+1)
          do k=1, ne
             do j=1, size(bf_nodes,2)
                do i=1, 2*bc_size
                   bf_nodes(i,j,k) = interior_nodes(i_match+i,j_match+j,k)
                end do
             end do
          end do

        end subroutine allocate_nodes_E


        subroutine allocate_nodes_W(
     $     bf_alignment, new_size, bf_nodes, interior_nodes)

          implicit none

          integer(ikind), dimension(2,2)            , intent(in) :: bf_alignment
          integer(ikind), dimension(2)              , intent(in) :: new_size
          real(rkind), dimension(:,:,:), allocatable, intent(out):: bf_nodes
          real(rkind), dimension(nx,ny,ne)          , intent(in) :: interior_nodes

          integer(ikind) :: i_match, j_match
          integer(ikind) :: i,j,k

          allocate(bf_nodes(new_size(1), new_size(2), ne))

          !copy of grid points from the interior
          i_match = 0
          j_match = bf_alignment(2,1)-(bc_size+1)
          do k=1, ne
             do j=1, size(bf_nodes,2)
                do i=1, 2*bc_size
                   bf_nodes(i+1,j,k) = interior_nodes(i_match+i,j_match+j,k)
                end do
             end do
          end do

        end subroutine allocate_nodes_W

      
        subroutine allocate_grdpts_id_N(
     $     bf_alignment, new_size, bf_grdpts_id)

          implicit none

          integer(ikind), dimension(2,2)             , intent(in) :: bf_alignment
          integer(ikind), dimension(2)               , intent(in) :: new_size
          integer       , dimension(:,:), allocatable, intent(out):: bf_grdpts_id

          allocate(bf_grdpts_id(new_size(1), new_size(2)))

          call add_interior_layer_NS(
     $         bf_grdpts_id, 1, bc_size, bf_alignment)

          call add_bc_interior_layer_NS(
     $         bf_grdpts_id, bc_size+1, bf_alignment)

          call add_bc_layer_NS(
     $         bf_grdpts_id, bc_size+2)

          call add_no_pt_layer_NS(
     $         bf_grdpts_id, bc_size+3, size(bf_grdpts_id,2))

        end subroutine allocate_grdpts_id_N


        subroutine allocate_grdpts_id_S(
     $     bf_alignment, new_size, bf_grdpts_id)

          implicit none

          integer(ikind), dimension(2,2)             , intent(in) :: bf_alignment
          integer(ikind), dimension(2)               , intent(in) :: new_size
          integer       , dimension(:,:), allocatable, intent(out):: bf_grdpts_id

          allocate(bf_grdpts_id(new_size(1), new_size(2)))

          call add_no_pt_layer_NS(
     $         bf_grdpts_id, 1, size(bf_grdpts_id,2)-(bc_size+2))

          call add_bc_layer_NS(
     $         bf_grdpts_id, size(bf_grdpts_id,2)-(bc_size+2)+1)

          call add_bc_interior_layer_NS(
     $         bf_grdpts_id, size(bf_grdpts_id,2)-(bc_size+2)+2, bf_alignment)

          call add_interior_layer_NS(
     $         bf_grdpts_id,
     $         size(bf_grdpts_id,2)-(bc_size+2)+3,
     $         size(bf_grdpts_id,2), bf_alignment)

        end subroutine allocate_grdpts_id_S


        subroutine allocate_grdpts_id_E(
     $     bf_alignment, new_size, bf_grdpts_id)

          implicit none

          integer(ikind), dimension(2,2)             , intent(in) :: bf_alignment
          integer(ikind), dimension(2)               , intent(in) :: new_size
          integer       , dimension(:,:), allocatable, intent(out):: bf_grdpts_id

          allocate(bf_grdpts_id(new_size(1), new_size(2)))

          call add_edge_layer_E(
     $         bf_grdpts_id,
     $         1, bc_size,
     $         bf_alignment, .false.)

          call add_interior_layer_E(
     $         bf_grdpts_id, bc_size+1, size(bf_grdpts_id,2)-bc_size)

          call add_edge_layer_E(
     $         bf_grdpts_id,
     $         size(bf_grdpts_id,2)-bc_size+1, size(bf_grdpts_id,2),
     $         bf_alignment, .true.)

        end subroutine allocate_grdpts_id_E


        subroutine allocate_grdpts_id_W(
     $     bf_alignment, new_size, bf_grdpts_id)

          implicit none

          integer(ikind), dimension(2,2)             , intent(in) :: bf_alignment
          integer(ikind), dimension(2)               , intent(in) :: new_size
          integer       , dimension(:,:), allocatable, intent(out):: bf_grdpts_id

          allocate(bf_grdpts_id(new_size(1), new_size(2)))

          call add_edge_layer_W(
     $         bf_grdpts_id,
     $         1, bc_size,
     $         bf_alignment, .false.)

          call add_interior_layer_W(
     $         bf_grdpts_id, bc_size+1, size(bf_grdpts_id,2)-bc_size)

          call add_edge_layer_W(
     $         bf_grdpts_id,
     $         size(bf_grdpts_id,2)-bc_size+1, size(bf_grdpts_id,2),
     $         bf_alignment, .true.)

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

      
         subroutine add_interior_layer_NS(
     $     grdpts_id, j_min, j_max, alignment)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)                , intent(in)    :: j_min
          integer(ikind)                , intent(in)    :: j_max
          integer(ikind), dimension(2,2), intent(in)    :: alignment


          integer(ikind) :: i,j
          integer(ikind) :: i_min1, i_min2
          integer(ikind) :: i_max1, i_max2, i_max3

          
          if(alignment(1,1).gt.(bc_size+2)) then
             i_min1=0
             i_min2=i_min1+0
          else
             if(alignment(1,1).eq.(bc_size+2)) then
                i_min1=0
                i_min2=1
             else
                i_min1=bc_size-1
                i_min2=i_min1+1
             end if
          end if

          if(alignment(1,2).lt.(nx-bc_size-1)) then
             i_max1 = size(grdpts_id,1)
             i_max2 = 0
             i_max3 = 0
          else
             if(alignment(1,2).eq.(nx-bc_size-1)) then
                i_max1 = size(grdpts_id,1)-1
                i_max2 = size(grdpts_id,1)
                i_max3 = 0
             else
                i_max1 = size(grdpts_id,1)-bc_size
                i_max2 = size(grdpts_id,1)-bc_size+1
                i_max3 = size(grdpts_id,1)
             end if
          end if


          do j=j_min,j_max
             do i=1,i_min1
                grdpts_id(i,j) = bc_pt
             end do

             do i=i_min1+1,i_min2
                grdpts_id(i,j) = bc_interior_pt
             end do

             do i=i_min2+1, i_max1
                grdpts_id(i,j) = interior_pt
             end do

             do i=i_max1+1, i_max2
                grdpts_id(i,j) = bc_interior_pt
             end do

             do i=i_max2+1, i_max3
                grdpts_id(i,j) = bc_pt
             end do
          end do

        end subroutine add_interior_layer_NS


        subroutine add_bc_interior_layer_NS(
     $     grdpts_id, j, alignment)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)                , intent(in)    :: j
          integer(ikind), dimension(2,2), intent(in)    :: alignment


          integer(ikind) :: i
          integer(ikind) :: i_min1
          integer(ikind) :: i_max1, i_max2

          
          if(alignment(1,1).gt.(bc_size+1)) then
             i_min1=0
          else
             if(alignment(1,1).eq.(bc_size+1)) then
                i_min1=1
             end if
          end if

          if(alignment(1,2).lt.(nx-bc_size)) then
             i_max1 = size(grdpts_id,1)
             i_max2 = 0
          else
             if(alignment(1,2).eq.(nx-bc_size)) then
                i_max1 = size(grdpts_id,1)-1
                i_max2 = size(grdpts_id,1)
             end if
          end if

          do i=1,i_min1
             grdpts_id(i,j) = bc_pt
          end do

          do i=i_min1+1,i_max1
             grdpts_id(i,j) = bc_interior_pt
          end do

          do i=i_max1+1, i_max2
             grdpts_id(i,j) = bc_pt
          end do
        end subroutine add_bc_interior_layer_NS


        subroutine add_bc_layer_NS(
     $     grdpts_id, j)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)                , intent(in)    :: j

          integer(ikind) :: i

          do i=1,size(grdpts_id,1)
             grdpts_id(i,j) = bc_pt
          end do

        end subroutine add_bc_layer_NS


        subroutine add_no_pt_layer_NS(
     $     grdpts_id, j_min, j_max)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)                , intent(in)    :: j_min
          integer(ikind)                , intent(in)    :: j_max

          integer(ikind) :: i,j

          do j=j_min, j_max
             do i=1,size(grdpts_id,1)
                grdpts_id(i,j) = no_pt
             end do
          end do

        end subroutine add_no_pt_layer_NS


        !< add the interior layer for the initialization
        !> of gridpts_id for eastern buffer layers
        subroutine add_interior_layer_E(
     $     grdpts_id, j_min, j_max)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: j_min
          integer(ikind)         , intent(in)    :: j_max

          integer(ikind) :: i,j

          do j=j_min, j_max
             do i=1, bc_size
                grdpts_id(i,j) = interior_pt
             end do

             i=bc_size+1
             grdpts_id(i,j) = bc_interior_pt

             do i=bc_size+2, 2*bc_size
                grdpts_id(i,j) = bc_pt
             end do

             do i=2*bc_size+1, size(grdpts_id,1)
                grdpts_id(i,j) = no_pt
             end do
          end do          

        end subroutine add_interior_layer_E


        !< add the interior layer for the initialization
        !> of gridpts_id for western buffer layers
        subroutine add_interior_layer_W(
     $     grdpts_id, j_min, j_max)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)         , intent(in)    :: j_min
          integer(ikind)         , intent(in)    :: j_max

          integer(ikind) :: i,j

          do j=j_min, j_max
             do i=1, size(grdpts_id,1)-(2*bc_size)
                grdpts_id(i,j) = no_pt
             end do

             do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                grdpts_id(i,j) = bc_pt
             end do

             i=size(grdpts_id,1)-bc_size
             grdpts_id(i,j) = bc_interior_pt
             
             do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                grdpts_id(i,j) = interior_pt
             end do
             
          end do          

        end subroutine add_interior_layer_W


        !< add the bc_pt layer for the initialization
        !> of gridpts_id
        !> type=.true.  = upper layer
        !> type=.false. = lower layer
        subroutine add_edge_layer_E(
     $     grdpts_id, j_min, j_max,
     $     alignment, type)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)                , intent(in)    :: j_min
          integer(ikind)                , intent(in)    :: j_max
          integer(ikind), dimension(2,2), intent(in)    :: alignment
          logical                       , intent(in)    :: type

          integer(ikind) :: i,j

          if(type) then
             if(alignment(2,2).lt.(ny-bc_size-1)) then
                do j=j_min, j_max
                   do i=1, bc_size
                      grdpts_id(i,j) = interior_pt
                   end do
                   i=bc_size+1
                   grdpts_id(i,j) = bc_interior_pt

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do
                end do
             else

                if(alignment(2,2).eq.(ny-bc_size-1)) then
                   j=j_min
                   do i=1, bc_size
                      grdpts_id(i,j) = interior_pt
                   end do
                   i=bc_size+1
                   grdpts_id(i,j) = bc_interior_pt

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   j=j_min+1
                   do i=1,bc_size+1
                      grdpts_id(i,j)=bc_interior_pt
                   end do

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                else

                   j=j_min
                   do i=1, bc_size+1
                      grdpts_id(i,j) = bc_interior_pt
                   end do

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do

                   j=j_min+1
                   do i=1,2*bc_size
                      grdpts_id(i,j)=bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do                   
                   
                end if
             end if
          else
             if(alignment(2,1).gt.(bc_size+2)) then
                do j=j_min, j_max
                   do i=1, bc_size
                      grdpts_id(i,j) = interior_pt
                   end do
                   i=bc_size+1
                   grdpts_id(i,j) = bc_interior_pt

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do
                end do
             else

                if(alignment(2,1).eq.(bc_size+2)) then
                   j=j_min
                   do i=1,bc_size+1
                      grdpts_id(i,j)=bc_interior_pt
                   end do

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   j=j_min+1
                    do i=1, bc_size
                      grdpts_id(i,j) = interior_pt
                   end do
                   i=bc_size+1
                   grdpts_id(i,j) = bc_interior_pt

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                else

                   j=j_min
                   do i=1,2*bc_size
                      grdpts_id(i,j)=bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do

                   j=j_min+1
                   do i=1, bc_size+1
                      grdpts_id(i,j) = bc_interior_pt
                   end do

                   do i=bc_size+2, 2*bc_size
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=2*bc_size+1,size(grdpts_id,1)
                      grdpts_id(i,j) = no_pt
                   end do                   
                   
                end if
             end if
          end if

        end subroutine add_edge_layer_E


        !< add the bc_pt layer for the initialization
        !> of gridpts_id
        !> type=.true.  = upper layer
        !> type=.false. = lower layer
        subroutine add_edge_layer_W(
     $     grdpts_id, j_min, j_max, alignment, type)

          implicit none

          integer       , dimension(:,:), intent(inout) :: grdpts_id
          integer(ikind)                , intent(in)    :: j_min
          integer(ikind)                , intent(in)    :: j_max
          integer(ikind), dimension(2,2), intent(in)    :: alignment
          logical                       , intent(in)    :: type

          integer(ikind) :: i,j

          if(type) then

             if(alignment(2,2).lt.(ny-bc_size-1)) then

                do j=j_min, j_max
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   i=size(grdpts_id,1)-bc_size
                   grdpts_id(i,j) = bc_interior_pt
                   
                   do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = interior_pt
                   end do
                end do

             else

                if(alignment(2,2).eq.(ny-bc_size-1)) then
                   j=j_min
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   i=size(grdpts_id,1)-bc_size
                   grdpts_id(i,j) = bc_interior_pt
                   
                   do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = interior_pt
                   end do

                   j=j_min+1
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=size(grdpts_id,1)-bc_size, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_interior_pt
                   end do
                   
                   
                else

                   j=j_min
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=size(grdpts_id,1)-bc_size, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_interior_pt
                   end do

                   j=j_min+1
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_pt
                   end do

                end if
             end if

          else

             if(alignment(2,1).gt.(bc_size+2)) then
                do j=j_min, j_max
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   i=size(grdpts_id,1)-bc_size
                   grdpts_id(i,j) = bc_interior_pt
                   
                   do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = interior_pt
                   end do
                end do
             else

                if(alignment(2,1).eq.(bc_size+2)) then
                   j=j_min
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=size(grdpts_id,1)-bc_size, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_interior_pt
                   end do

                   j=j_min+1
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   i=size(grdpts_id,1)-bc_size
                   grdpts_id(i,j) = bc_interior_pt
                   
                   do i=size(grdpts_id,1)-bc_size+1, size(grdpts_id,1)
                      grdpts_id(i,j) = interior_pt
                   end do
                   
                else

                   j=j_min
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_pt
                   end do                   

                   j=j_min+1
                   do i=1, size(grdpts_id,1)-(2*bc_size)
                      grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=size(grdpts_id,1)-(2*bc_size)+1, size(grdpts_id,1)-bc_size-1
                      grdpts_id(i,j) = bc_pt
                   end do

                   do i=size(grdpts_id,1)-bc_size, size(grdpts_id,1)
                      grdpts_id(i,j) = bc_interior_pt
                   end do
                end if
             end if
          end if

        end subroutine add_edge_layer_W

      end module bf_layer_allocate_module
