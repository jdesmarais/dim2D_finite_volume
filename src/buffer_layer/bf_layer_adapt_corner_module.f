      module bf_layer_adapt_corner_module

        use bf_layer_merge_module, only : get_new_size

        use parameters_bf_layer  , only : exchange_pt, interior_pt,
     $                                    bc_interior_pt, bc_pt, no_pt
        use parameters_constant  , only : N,S,E,W,
     $                                    x_direction, y_direction,
     $                                    min_border, max_border
        use parameters_input     , only : nx,ny,ne,bc_size
        use parameters_kind      , only : ikind, rkind


        implicit none


        private
        public :: adapt_bf_layer_N_to_corner,
     $            adapt_bf_layer_S_to_corner,
     $            adapt_bf_layer_E_to_corner,
     $            adapt_bf_layer_W_to_corner


        contains


        !< adapt nodes and grdpts to corner for northern buffer layers
        subroutine adapt_bf_layer_N_to_corner(
     $     nodes, grdpts_id, alignment,
     $     interior_nodes, new_alignment, new_neighbors)


          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: nodes
          integer       , dimension(:,:)     , allocatable, intent(inout) :: grdpts_id
          integer(ikind), dimension(2,2)                  , intent(inout) :: alignment
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2)                  , intent(in)    :: new_alignment
          logical       , dimension(4)                    , intent(in)    :: new_neighbors

          logical :: reallocation_needed
          integer(ikind), dimension(2)                  :: new_size
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
          
          integer(ikind) :: i_min1
          integer(ikind) :: i_min2
          integer(ikind) :: j_min
          integer(ikind) :: interior_i_max1
          integer(ikind) :: interior_i_max2

          
          !check whether there is any change to be made
          reallocation_needed = check_reallocation(alignment, new_alignment)

          !if no reallocation is required, no nodes transfer is needed
          !so only the grdpts_id are modified locally for the neighbors
          if(.not.reallocation_needed) then
             
             call modify_grdpts_id_NS_for_corner(
     $            grdpts_id, new_neighbors,
     $            bc_size+1, 2*bc_size+1)

          else

             !determine the new size
             new_size = get_new_size(alignment, new_alignment)

             !determine the match indices
             call get_match_indices(
     $            y_direction,
     $            interior_i_max1, interior_i_max2,
     $            i_min1, i_min2,
     $            j_min,
     $            alignment, new_alignment)

             !adapt the nodes + move allocation
             allocate(new_nodes(new_size(1), new_size(2), ne))
             call adapt_nodes_N_to_corner(
     $            new_nodes,
     $            nodes, interior_nodes,
     $            new_alignment,
     $            interior_i_max1, interior_i_max2,
     $            i_min1, i_min2,
     $            j_min)
             call MOVE_ALLOC(new_nodes, nodes)

             !adapt the grdpts_id + move allocation
             allocate(new_grdpts_id(new_size(1), new_size(2)))
             call adapt_grdpts_id_N_to_corner(
     $            new_grdpts_id,
     $            grdpts_id,
     $            interior_i_max1, interior_i_max2,
     $            i_min1, i_min2,
     $            j_min,
     $            new_neighbors(E),
     $            new_neighbors(W))
             call MOVE_ALLOC(new_grdpts_id, grdpts_id)
     $            
             !change the alignment
             alignment = new_alignment

          end if

        end subroutine adapt_bf_layer_N_to_corner


        !< adapt nodes and grdpts to corner for southern buffer layers
        subroutine adapt_bf_layer_S_to_corner(
     $     nodes, grdpts_id, alignment,
     $     interior_nodes, new_alignment, new_neighbors)


          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: nodes
          integer       , dimension(:,:)     , allocatable, intent(inout) :: grdpts_id
          integer(ikind), dimension(2,2)                  , intent(inout) :: alignment
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2)                  , intent(in)    :: new_alignment
          logical       , dimension(4)                    , intent(in)    :: new_neighbors

          logical :: reallocation_needed
          integer(ikind), dimension(2)                  :: new_size
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
          
          integer(ikind) :: i_min1
          integer(ikind) :: i_min2
          integer(ikind) :: j_min
          integer(ikind) :: interior_i_max1
          integer(ikind) :: interior_i_max2

          
          !check whether there is any change to be made
          reallocation_needed = check_reallocation(alignment, new_alignment)

          !if no reallocation is required, no nodes transfer is needed
          !so only the grdpts_id are modified locally for the neighbors
          if(.not.reallocation_needed) then
             
             call modify_grdpts_id_NS_for_corner(
     $            grdpts_id, new_neighbors,
     $            size(grdpts_id,2)-2*bc_size, 
     $            size(grdpts_id,2)-bc_size)

          else

             !determine the new size
             new_size = get_new_size(alignment, new_alignment)

             !determine the match indices
             call get_match_indices(
     $            y_direction,
     $            interior_i_max1, interior_i_max2,
     $            i_min1, i_min2,
     $            j_min,
     $            alignment, new_alignment)

             !adapt the nodes + move allocation
             allocate(new_nodes(new_size(1), new_size(2), ne))
             call adapt_nodes_S_to_corner(
     $            new_nodes,
     $            nodes, interior_nodes,
     $            new_alignment,
     $            interior_i_max1, interior_i_max2,
     $            i_min1, i_min2,
     $            j_min, new_size(2))
             call MOVE_ALLOC(new_nodes, nodes)

             !adapt the grdpts_id + move allocation
             allocate(new_grdpts_id(new_size(1), new_size(2)))
             call adapt_grdpts_id_S_to_corner(
     $            new_grdpts_id,
     $            grdpts_id,
     $            interior_i_max1, interior_i_max2,
     $            i_min1, i_min2,
     $            j_min, new_size(2),
     $            new_neighbors(E),
     $            new_neighbors(W))
             call MOVE_ALLOC(new_grdpts_id, grdpts_id)
     $            
             !change the alignment
             alignment = new_alignment

          end if

        end subroutine adapt_bf_layer_S_to_corner


        !< adapt nodes and grdpts to corner for southern buffer layers
        subroutine adapt_bf_layer_E_to_corner(
     $     nodes, grdpts_id, alignment,
     $     interior_nodes, new_alignment, new_neighbors)


          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: nodes
          integer       , dimension(:,:)     , allocatable, intent(inout) :: grdpts_id
          integer(ikind), dimension(2,2)                  , intent(inout) :: alignment
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2)                  , intent(in)    :: new_alignment
          logical       , dimension(4)                    , intent(in)    :: new_neighbors

          logical :: reallocation_needed
          integer(ikind), dimension(2)                  :: new_size
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
          
          integer(ikind) :: j_min1
          integer(ikind) :: j_min2
          integer(ikind) :: i_min
          integer(ikind) :: interior_j_max1
          integer(ikind) :: interior_j_max2

          
          !check whether there is any change to be made
          reallocation_needed = check_reallocation(alignment, new_alignment)

          !if no reallocation is required, no nodes transfer is needed
          !so only the grdpts_id are modified locally for the neighbors
          if(.not.reallocation_needed) then
             
             call modify_grdpts_id_EW_for_corner(
     $            grdpts_id, new_neighbors,
     $            bc_size+1, 2*bc_size+1)

          else

             !determine the new size
             new_size = get_new_size(alignment, new_alignment)

             !determine the match indices
             call get_match_indices(
     $            x_direction,
     $            interior_j_max1, interior_j_max2,
     $            j_min1, j_min2,
     $            i_min,
     $            alignment, new_alignment)

             !adapt the nodes + move allocation
             allocate(new_nodes(new_size(1), new_size(2), ne))
             call adapt_nodes_E_to_corner(
     $            new_nodes,
     $            nodes, interior_nodes,
     $            new_alignment,
     $            interior_j_max1, interior_j_max2,
     $            j_min1, j_min2)
             call MOVE_ALLOC(new_nodes, nodes)

             !adapt the grdpts_id + move allocation
             allocate(new_grdpts_id(new_size(1), new_size(2)))
             call adapt_grdpts_id_E_to_corner(
     $            new_grdpts_id,
     $            grdpts_id,
     $            interior_j_max1, interior_j_max2,
     $            j_min1, j_min2,
     $            new_size(1),
     $            new_neighbors(N),
     $            new_neighbors(S))
             call MOVE_ALLOC(new_grdpts_id, grdpts_id)
     $            
             !change the alignment
             alignment = new_alignment

          end if

        end subroutine adapt_bf_layer_E_to_corner


        !< adapt nodes and grdpts to corner for southern buffer layers
        subroutine adapt_bf_layer_W_to_corner(
     $     nodes, grdpts_id, alignment,
     $     interior_nodes, new_alignment, new_neighbors)


          real(rkind)   , dimension(:,:,:)   , allocatable, intent(inout) :: nodes
          integer       , dimension(:,:)     , allocatable, intent(inout) :: grdpts_id
          integer(ikind), dimension(2,2)                  , intent(inout) :: alignment
          real(rkind)   , dimension(nx,ny,ne)             , intent(in)    :: interior_nodes
          integer(ikind), dimension(2,2)                  , intent(in)    :: new_alignment
          logical       , dimension(4)                    , intent(in)    :: new_neighbors

          logical :: reallocation_needed
          integer(ikind), dimension(2)                  :: new_size
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          integer       , dimension(:,:)  , allocatable :: new_grdpts_id
          
          integer(ikind) :: j_min1
          integer(ikind) :: j_min2
          integer(ikind) :: i_min
          integer(ikind) :: interior_j_max1
          integer(ikind) :: interior_j_max2

          
          !check whether there is any change to be made
          reallocation_needed = check_reallocation(alignment, new_alignment)

          !if no reallocation is required, no nodes transfer is needed
          !so only the grdpts_id are modified locally for the neighbors
          if(.not.reallocation_needed) then
             
             call modify_grdpts_id_EW_for_corner(
     $            grdpts_id, new_neighbors,
     $            size(grdpts_id,1)-(2*bc_size),
     $            size(grdpts_id,1)-bc_size)

          else

             !determine the new size
             new_size = get_new_size(alignment, new_alignment)

             !determine the match indices
             call get_match_indices(
     $            x_direction,
     $            interior_j_max1, interior_j_max2,
     $            j_min1, j_min2,
     $            i_min,
     $            alignment, new_alignment)

             !adapt the nodes + move allocation
             allocate(new_nodes(new_size(1), new_size(2), ne))
             call adapt_nodes_W_to_corner(
     $            new_nodes,
     $            nodes, interior_nodes,
     $            new_alignment,
     $            interior_j_max1, interior_j_max2,
     $            j_min1, j_min2)
             call MOVE_ALLOC(new_nodes, nodes)

             !adapt the grdpts_id + move allocation
             allocate(new_grdpts_id(new_size(1), new_size(2)))
             call adapt_grdpts_id_W_to_corner(
     $            new_grdpts_id,
     $            grdpts_id,
     $            interior_j_max1, interior_j_max2,
     $            j_min1, j_min2,
     $            new_size(1),
     $            new_neighbors(N),
     $            new_neighbors(S))
             call MOVE_ALLOC(new_grdpts_id, grdpts_id)
     $            
             !change the alignment
             alignment = new_alignment

          end if

        end subroutine adapt_bf_layer_W_to_corner


        !< check if a reallcoation is needed
        function check_reallocation(alignment, new_alignment)
     $     result(reallocation_needed)

          implicit none

          integer(ikind), dimension(2,2), intent(in) :: alignment
          integer(ikind), dimension(2,2), intent(in) :: new_alignment
          logical                                    :: reallocation_needed

          integer :: i,j

          reallocation_needed = .false.
          do j=1,2
             do i=1,2
                if(alignment(i,j).ne.new_alignment(i,j)) then
                   reallocation_needed = .true.
                   exit
                end if
             end do
             
             if(reallocation_needed) then
                exit
             end if
          end do

        end function check_reallocation


        !< modify the gridpts_id for the corner for northern
        !> or southern buffer layers
        subroutine modify_grdpts_id_NS_for_corner(
     $     gridpts_id, new_neighbors, j_min, j_max)

          implicit none

          integer, dimension(:,:), intent(inout) :: gridpts_id
          logical, dimension(4)  , intent(in)    :: new_neighbors
          integer                , intent(in)    :: j_min
          integer                , intent(in)    :: j_max


          integer(ikind) :: i,j, i_min, i_max


          if(new_neighbors(W)) then
             i_max = bc_size
          else
             i_max = 0
          end if

          if(new_neighbors(E)) then
             i_min = size(gridpts_id,1)-bc_size+1
          else
             i_min = size(gridpts_id,1)+1
          end if

          do j=j_min, j_max
             do i=1, i_max
                gridpts_id(i,j) = exchange_pt
             end do

             do i=i_min, size(gridpts_id,1)
                gridpts_id(i,j) = exchange_pt
             end do
          end do

        end subroutine modify_grdpts_id_NS_for_corner


        !< modify the gridpts_id for the corner for northern
        !> or southern buffer layers
        subroutine modify_grdpts_id_EW_for_corner(
     $     gridpts_id, new_neighbors, i_min, i_max)

          implicit none

          integer, dimension(:,:), intent(inout) :: gridpts_id
          logical, dimension(4)  , intent(in)    :: new_neighbors
          integer                , intent(in)    :: i_min
          integer                , intent(in)    :: i_max


          integer(ikind) :: i,j, j_min, j_max


          if(new_neighbors(S)) then
             j_max = bc_size
          else
             j_max = 0
          end if

          if(new_neighbors(N)) then
             j_min = size(gridpts_id,2)-bc_size+1
          else
             j_min = size(gridpts_id,2)+1
          end if

          do j=1, j_max
             do i=i_min, i_max
                gridpts_id(i,j) = exchange_pt
             end do
          end do

          do j=j_min, size(gridpts_id,2)
             do i=i_min, i_max
                gridpts_id(i,j) = exchange_pt
             end do
          end do

        end subroutine modify_grdpts_id_EW_for_corner


        !< get the indices to easily correspond the gridpoints
        !> from the newly allocated tables to the previously
        !> allocated tables
        subroutine get_match_indices(
     $     direction,
     $     interior_i_max1, interior_i_max2,
     $     i_min1, i_min2,
     $     j_min,
     $     alignment, new_alignment)

          implicit none

          integer                                 , intent(in)  :: direction
          integer                                 , intent(out) :: interior_i_max1
          integer                                 , intent(out) :: interior_i_max2
          integer                                 , intent(out) :: i_min1
          integer                                 , intent(out) :: i_min2
          integer                                 , intent(out) :: j_min
          integer(ikind), dimension(2,2)          , intent(in)  :: alignment
          integer(ikind), dimension(2,2)          , intent(in)  :: new_alignment

          integer :: dir1, dir2

          select case(direction)
            case(y_direction)
               dir1 = x_direction
               dir2 = y_direction
            case(x_direction)
               dir1 = y_direction
               dir2 = x_direction
            case default
               print '(''bf_layer_adapt_corner_module'')'
               print '(''get_match_indices'')'
               print '(''direction not recognized'')'
               print '(''direction: '',I2)', direction
               stop 'either use x_direction or y_direction'
          end select

          interior_i_max1 = max(alignment(dir1,min_border) - new_alignment(dir1,min_border),0)
          interior_i_max2 = max(new_alignment(dir1,max_border) - alignment(dir1,max_border),0)

          i_min1 = interior_i_max1
          i_min2 = i_min1 + alignment(dir1,max_border) - alignment(dir1,min_border) + 2*bc_size+1

          j_min  = alignment(dir2,max_border)-alignment(dir2,min_border)+2*bc_size+1

        end subroutine get_match_indices
        

        !< adapt the nodes of the northern buffer layers
        !> to the corner
        subroutine adapt_nodes_N_to_corner(
     $     new_nodes,
     $     nodes, interior_nodes,
     $     new_alignment,
     $     interior_i_max1, interior_i_max2,
     $     i_min1, i_min2,
     $     j_min)
        
          implicit none

          real(rkind), dimension(:,:,:)   , intent(out) :: new_nodes
          real(rkind), dimension(:,:,:)   , intent(in)  :: nodes
          real(rkind), dimension(nx,ny,ne), intent(in)  :: interior_nodes
          integer(ikind), dimension(2,2)  , intent(in)  :: new_alignment
          integer(ikind)                  , intent(in)  :: interior_i_max1
          integer(ikind)                  , intent(in)  :: interior_i_max2
          integer(ikind)                  , intent(in)  :: i_min1
          integer(ikind)                  , intent(in)  :: i_min2
          integer(ikind)                  , intent(in)  :: j_min
          
          integer(ikind) :: i,j
          integer        :: k

          
          do k=1, ne
             do j=1, 2*bc_size
                do i=1, interior_i_max1
                   new_nodes(i,j,k) = interior_nodes(
     $                  new_alignment(1,1)-(bc_size+1)+i,
     $                  ny-2*bc_size+j,
     $                  k)
                end do

                do i=1, size(nodes,1)
                   new_nodes(i_min1+i,j,k) = nodes(i,j,k)
                end do
                
                do i=1, interior_i_max2
                   new_nodes(i_min2+i,j,k) = interior_nodes(
     $                  new_alignment(1,2)-interior_i_max2+bc_size+i,
     $                  ny-2*bc_size+j,
     $                  k)
                end do
                
             end do

             do j=2*bc_size+1, j_min
                do i=1, size(nodes,1)
                   new_nodes(i_min1+i,j,k) = nodes(i,j,k)
                end do
             end do
          end do

        end subroutine adapt_nodes_N_to_corner


        !< adapt the nodes of the southern buffer layers
        !> to the corner
        subroutine adapt_nodes_S_to_corner(
     $     new_nodes,
     $     nodes, interior_nodes,
     $     new_alignment,
     $     interior_i_max1, interior_i_max2,
     $     i_min1, i_min2,
     $     j_min, j_max)
        
          implicit none

          real(rkind), dimension(:,:,:)   , intent(out) :: new_nodes
          real(rkind), dimension(:,:,:)   , intent(in)  :: nodes
          real(rkind), dimension(nx,ny,ne), intent(in)  :: interior_nodes
          integer(ikind), dimension(2,2)  , intent(in)  :: new_alignment
          integer(ikind)                  , intent(in)  :: interior_i_max1
          integer(ikind)                  , intent(in)  :: interior_i_max2
          integer(ikind)                  , intent(in)  :: i_min1
          integer(ikind)                  , intent(in)  :: i_min2
          integer(ikind)                  , intent(in)  :: j_min
          integer(ikind)                  , intent(in)  :: j_max
          
          integer(ikind) :: i,j
          integer        :: k

          
          do k=1, ne
             do j=j_max-j_min+1, j_max-(2*bc_size)
                do i=1, size(nodes,1)
                   new_nodes(i_min1+i,j,k) = nodes(i,j-(j_max-j_min),k)
                end do
             end do

             do j=j_max-(2*bc_size)+1, j_max
                do i=1, interior_i_max1
                   new_nodes(i,j,k) = interior_nodes(
     $                  new_alignment(1,1)-(bc_size+1)+i,
     $                  j-(j_max-2*bc_size),
     $                  k)
                end do

                do i=1, size(nodes,1)
                   new_nodes(i_min1+i,j,k) = nodes(i,j-(j_max-j_min),k)
                end do
                
                do i=1, interior_i_max2
                   new_nodes(i_min2+i,j,k) = interior_nodes(
     $                  new_alignment(1,2)-interior_i_max2+bc_size+i,
     $                  j-(j_max-2*bc_size),
     $                  k)
                end do
                
             end do

             
          end do

        end subroutine adapt_nodes_S_to_corner


        !< adapt the nodes of the eastern buffer layers
        !> to the corner
        subroutine adapt_nodes_E_to_corner(
     $     new_nodes,
     $     nodes, interior_nodes,
     $     new_alignment,
     $     interior_j_max1, interior_j_max2,
     $     j_min1, j_min2)
        
          implicit none

          real(rkind), dimension(:,:,:)   , intent(out) :: new_nodes
          real(rkind), dimension(:,:,:)   , intent(in)  :: nodes
          real(rkind), dimension(nx,ny,ne), intent(in)  :: interior_nodes
          integer(ikind), dimension(2,2)  , intent(in)  :: new_alignment
          integer(ikind)                  , intent(in)  :: interior_j_max1
          integer(ikind)                  , intent(in)  :: interior_j_max2
          integer(ikind)                  , intent(in)  :: j_min1
          integer(ikind)                  , intent(in)  :: j_min2
          
          integer(ikind) :: i,j
          integer        :: k

          
          do k=1, ne
             do j=1, interior_j_max1
                do i=1, 2*bc_size
                   new_nodes(i,j,k) = interior_nodes(
     $                  nx-2*bc_size+i,
     $                  new_alignment(2,1)-(bc_size+1)+j,
     $                  k)
                end do
             end do

             do j=1, size(nodes,2)
                do i=1, size(nodes,1)
                   new_nodes(i,j_min1+j,k) = nodes(i,j,k)
                end do
             end do

             do j=1, interior_j_max2
                do i=1, 2*bc_size
                   new_nodes(i,j_min2+j,k) = interior_nodes(
     $                  nx-2*bc_size+i,
     $                  new_alignment(2,2)-interior_j_max2+bc_size+j,
     $                  k)
                end do
             end do
          end do

        end subroutine adapt_nodes_E_to_corner


        !< adapt the nodes of the western buffer layers
        !> to the corner
        subroutine adapt_nodes_W_to_corner(
     $     new_nodes,
     $     nodes, interior_nodes,
     $     new_alignment,
     $     interior_j_max1, interior_j_max2,
     $     j_min1, j_min2)
        
          implicit none

          real(rkind), dimension(:,:,:)   , intent(out) :: new_nodes
          real(rkind), dimension(:,:,:)   , intent(in)  :: nodes
          real(rkind), dimension(nx,ny,ne), intent(in)  :: interior_nodes
          integer(ikind), dimension(2,2)  , intent(in)  :: new_alignment
          integer(ikind)                  , intent(in)  :: interior_j_max1
          integer(ikind)                  , intent(in)  :: interior_j_max2
          integer(ikind)                  , intent(in)  :: j_min1
          integer(ikind)                  , intent(in)  :: j_min2
          
          integer(ikind) :: i,j,i_min,i_max
          integer        :: k

          i_min = size(new_nodes,1)-2*bc_size
          i_max = size(new_nodes,1)-size(nodes,1)
          
          do k=1, ne
             do j=1, interior_j_max1
                do i=i_min+1, size(new_nodes,1)
                   new_nodes(i,j,k) = interior_nodes(
     $                  i-i_min,
     $                  new_alignment(2,1)-(bc_size+1)+j,
     $                  k)
                end do
             end do

             do j=1, size(nodes,2)
                do i=1, size(nodes,1)
                   new_nodes(i_max+i,j_min1+j,k) = nodes(i,j,k)
                end do
             end do

             do j=1, interior_j_max2
                do i=i_min+1, size(new_nodes,1)
                   new_nodes(i,j_min2+j,k) = interior_nodes(
     $                  i-i_min,
     $                  new_alignment(2,2)-interior_j_max2+bc_size+j,
     $                  k)
                end do
             end do
          end do

        end subroutine adapt_nodes_W_to_corner


        !< adapt the grdpts_id of the northern buffer layers
        !> to the corner
        subroutine adapt_grdpts_id_N_to_corner(
     $     new_grdpts_id,
     $     grdpts_id,
     $     interior_i_max1, interior_i_max2,
     $     i_min1, i_min2,
     $     j_min,
     $     neighbor_E,
     $     neighbor_W)
        
          implicit none

          integer       , dimension(:,:), intent(out) :: new_grdpts_id
          integer       , dimension(:,:), intent(in)  :: grdpts_id
          integer(ikind)                , intent(in)  :: interior_i_max1
          integer(ikind)                , intent(in)  :: interior_i_max2
          integer(ikind)                , intent(in)  :: i_min1
          integer(ikind)                , intent(in)  :: i_min2
          integer(ikind)                , intent(in)  :: j_min
          logical                       , intent(in)  :: neighbor_E
          logical                       , intent(in)  :: neighbor_W
          
          integer(ikind) :: i,j


          !exchange_layer
          call add_exchange_pt_layer_NS(
     $         new_grdpts_id,
     $         grdpts_id,
     $         1, bc_size,
     $         interior_i_max1,
     $         interior_i_max2,
     $         i_min1, 0,
     $         i_min2)

          !interior_layer
          call add_interior_pt_layer_NS(
     $         new_grdpts_id,
     $         grdpts_id,
     $         bc_size+1,
     $         interior_i_max1,
     $         interior_i_max2,
     $         i_min1, 0,
     $         i_min2,
     $         neighbor_E,
     $         neighbor_W)

          !bc_interior
          call add_bc_interior_pt_layer_NS(
     $         new_grdpts_id,
     $         grdpts_id,
     $         bc_size+2,
     $         interior_i_max1,
     $         interior_i_max2,
     $         i_min1, 0,
     $         i_min2,
     $         neighbor_E,
     $         neighbor_W)

          !bc_layer
          call add_bc_pt_layer_NS(
     $         new_grdpts_id,
     $         grdpts_id,
     $         bc_size+3,
     $         interior_i_max1,
     $         interior_i_max2,
     $         i_min1, 0,
     $         i_min2,
     $         neighbor_E,
     $         neighbor_W)

          !no_pt
          call add_no_pt_layer_NS(
     $         new_grdpts_id,
     $         grdpts_id,
     $         bc_size+4, j_min,
     $         interior_i_max1,
     $         interior_i_max2,
     $         i_min1, 0,
     $         i_min2)

          !last top layer
          do j=size(grdpts_id,2)+1, size(new_grdpts_id,2)
             do i=1, size(new_grdpts_id,1)
                new_grdpts_id(i,j) = no_pt
             end do
          end do

        end subroutine adapt_grdpts_id_N_to_corner


        !< adapt the grdpts_id of the northern buffer layers
        !> to the corner
        subroutine adapt_grdpts_id_S_to_corner(
     $     new_grdpts_id,
     $     grdpts_id,
     $     interior_i_max1, interior_i_max2,
     $     i_min1, i_min2,
     $     j_min, j_max,
     $     neighbor_E,
     $     neighbor_W)
        
          implicit none

          integer       , dimension(:,:), intent(out) :: new_grdpts_id
          integer       , dimension(:,:), intent(in)  :: grdpts_id
          integer(ikind)                , intent(in)  :: interior_i_max1
          integer(ikind)                , intent(in)  :: interior_i_max2
          integer(ikind)                , intent(in)  :: i_min1
          integer(ikind)                , intent(in)  :: i_min2
          integer(ikind)                , intent(in)  :: j_min
          integer(ikind)                , intent(in)  :: j_max
          logical                       , intent(in)  :: neighbor_E
          logical                       , intent(in)  :: neighbor_W
          
          integer(ikind) :: i,j


          !last top layer
          do j=1, j_max-j_min
             do i=1, size(new_grdpts_id,1)
                new_grdpts_id(i,j) = no_pt
             end do
          end do

          !no_pt
          call add_no_pt_layer_NS(
     $         new_grdpts_id,
     $         grdpts_id,
     $         j_max-j_min+1, j_max-bc_size-3,
     $         interior_i_max1,
     $         interior_i_max2,
     $         i_min1, -(j_max-j_min),
     $         i_min2)

          !bc_layer
          call add_bc_pt_layer_NS(
     $         new_grdpts_id,
     $         grdpts_id,
     $         j_max-bc_size-2,
     $         interior_i_max1,
     $         interior_i_max2,
     $         i_min1, -(j_max-j_min),
     $         i_min2,
     $         neighbor_E,
     $         neighbor_W)

          !bc_interior
          call add_bc_interior_pt_layer_NS(
     $         new_grdpts_id,
     $         grdpts_id,
     $         j_max-bc_size-1,
     $         interior_i_max1,
     $         interior_i_max2,
     $         i_min1, -(j_max-j_min),
     $         i_min2,
     $         neighbor_E,
     $         neighbor_W)

          !interior_layer
          call add_interior_pt_layer_NS(
     $         new_grdpts_id,
     $         grdpts_id,
     $         j_max-bc_size,
     $         interior_i_max1,
     $         interior_i_max2,
     $         i_min1, -(j_max-j_min),
     $         i_min2,
     $         neighbor_E,
     $         neighbor_W)

          !exchange_layer
          call add_exchange_pt_layer_NS(
     $         new_grdpts_id,
     $         grdpts_id,
     $         j_max-bc_size+1, j_max,
     $         interior_i_max1,
     $         interior_i_max2,
     $         i_min1, -(j_max-j_min),
     $         i_min2)                    

        end subroutine adapt_grdpts_id_S_to_corner


        !< adapt the grdpts_id of the easthern buffer layers
        !> to the corner
        subroutine adapt_grdpts_id_E_to_corner(
     $     new_grdpts_id,
     $     grdpts_id,
     $     interior_j_max1, interior_j_max2,
     $     j_min1, j_min2,
     $     i_max,
     $     neighbor_N,
     $     neighbor_S)
        
          implicit none

          integer       , dimension(:,:), intent(out) :: new_grdpts_id
          integer       , dimension(:,:), intent(in)  :: grdpts_id
          integer(ikind)                , intent(in)  :: interior_j_max1
          integer(ikind)                , intent(in)  :: interior_j_max2
          integer(ikind)                , intent(in)  :: j_min1
          integer(ikind)                , intent(in)  :: j_min2
          integer(ikind)                , intent(in)  :: i_max
          logical                       , intent(in)  :: neighbor_N
          logical                       , intent(in)  :: neighbor_S
          
          integer(ikind) :: i,j, j_min_block, j_max_block
          integer(ikind) :: j_max

          j_max = size(grdpts_id,2)
          

          !block1: nodes copied from the interior domain
          !        frontier between block1 and block2
          if(interior_j_max1.ge.2) then
             
             !> grdpt_id corresponding to nodes copied
             !> from the interior domain: 1st block
             
             !1st block : j=1, bc_size
             if(neighbor_S) then
                do j=1, bc_size
                   do i=1, 2*bc_size+1
                      new_grdpts_id(i,j) = exchange_pt
                   end do
                   do i=2*bc_size+2, i_max
                      new_grdpts_id(i,j) = no_pt
                   end do
                end do
             else
                call add_layer1_out_E(
     $               1, new_grdpts_id, i_max)
                
                call add_layer2_out_E(
     $               2, new_grdpts_id, i_max)
             end if
             
             !1st block : j>bc_size
             call add_layer_interior_E(
     $            3, interior_j_max1, new_grdpts_id, i_max)

             !> frontier between the 1st and 2nd blocks
             call add_layer1_in_E(
     $            interior_j_max1+1, new_grdpts_id,
     $            1, grdpts_id,
     $            i_max)

             call add_layer1_in_E(
     $            interior_j_max1+2, new_grdpts_id,
     $            2, grdpts_id,
     $            i_max)

             !> update the lower border for the next block copied
             j_min_block = bc_size+1

          else

             if(interior_j_max1.eq.1) then
                
                !1st block : j=1, bc_size
                if(neighbor_S) then
                   j=1
                   do i=1, 2*bc_size+1
                      new_grdpts_id(i,j) = exchange_pt
                   end do
                   do i=2*bc_size+2, i_max
                      new_grdpts_id(i,j) = no_pt
                   end do

                   j=2
                   do i=1, 2*bc_size+1
                      new_grdpts_id(i,j) = exchange_pt
                   end do

                   do i=2*bc_size+2, size(grdpts_id,1)
                      new_grdpts_id(i,j) = grdpts_id(i,j-1)
                   end do

                   do i=size(grdpts_id,1)+1, i_max
                      new_grdpts_id(i,j) = no_pt
                   end do
                   
                else

                   !> grdpt_id corresponding to nodes copied
                   !> from the interior domain: 1st block
                   call add_layer1_out_E(
     $                  1, new_grdpts_id, i_max)
                   
                   !> frontier between the 1st and 2nd blocks
                   call add_layer2_in_E(
     $                  2, new_grdpts_id,
     $                  1, grdpts_id,
     $                  i_max)
                end if                
                
                call add_layer1_in_E(
     $               3, new_grdpts_id,
     $               2, grdpts_id,
     $               i_max)

                !> update the lower border for the next block copied
                j_min_block = bc_size+1
                
             else
                
                if(neighbor_S) then
                   
                   do j=1, bc_size
                      do i=1, 2*bc_size+1
                         new_grdpts_id(i,j) = exchange_pt
                      end do

                      do i=2*bc_size+2, size(grdpts_id,1)
                         new_grdpts_id(i,j) = grdpts_id(i,j)
                      end do
                      
                      do i=size(grdpts_id,1)+1, i_max
                         new_grdpts_id(i,j) = no_pt
                      end do                      
                      
                   end do

                   !> update the lower border for the next block copied
                   j_min_block = bc_size+1
                else

                   !> update the lower border for the next block copied
                   j_min_block = 1
                end if

             end if
          end if          

          
          !block2: grdpt copied from the sublayer
          !block3: grdpt corresponding to nodes copied
          !        from interior to the buffer layer
          if(interior_j_max2.ge.2) then

             !> update the upper border for block2 copied
             j_max_block = j_max-bc_size

             !> copy the block corresponding to sublayer_
             !> first
             do j=j_min_block, j_max_block
                do i=1, size(grdpts_id,1)
                   new_grdpts_id(i,j_min1+j) = grdpts_id(i,j)
                end do
                
                do i=size(grdpts_id,1)+1, i_max
                   new_grdpts_id(i,j_min1+j) = no_pt
                end do
             end do

             !> frontier between block4 and block5
             call add_layer1_in_E(
     $            j_min1+j_max-1, new_grdpts_id,
     $            j_max_block+1, grdpts_id,
     $            i_max)

             call add_layer1_in_E(
     $            j_min1+j_max, new_grdpts_id,
     $            j_max_block+2, grdpts_id,
     $            i_max)

             call add_layer_interior_E(
     $            j_min2+1,
     $            j_min2+interior_j_max2-bc_size,
     $            new_grdpts_id,
     $            i_max)

             !> frontier at the edge
             if(neighbor_N) then
                do j=j_min2+interior_j_max2-1, j_min2+interior_j_max2
                   do i=1, 2*bc_size+1
                      new_grdpts_id(i,j) = exchange_pt
                   end do
                   do i=2*bc_size+2, i_max
                      new_grdpts_id(i,j) = no_pt
                   end do
                end do
             else
                call add_layer2_out_E(
     $               j_min2+interior_j_max2-1, new_grdpts_id, i_max)
                
                call add_layer1_out_E(
     $               j_min2+interior_j_max2  , new_grdpts_id, i_max)
             end if
          else

             if(interior_j_max2.eq.1) then
                
                !> update the upper border for block2 copied
                j_max_block = j_max-bc_size

                !> copy the block corresponding to sublayer_
                !> first
                do j=j_min_block, j_max_block
                   do i=1, size(grdpts_id,1)
                      new_grdpts_id(i,j_min1+j) = grdpts_id(i,j)
                   end do
                   
                   do i=size(grdpts_id,1)+1, i_max
                      new_grdpts_id(i,j_min1+j) = no_pt
                   end do
                end do

                !> frontier between block4 and the edge
                call add_layer1_in_E(
     $               j_min1+j_max_block+1, new_grdpts_id,
     $               j_max_block+1, grdpts_id,
     $               i_max)

                !block5: last two layers
                if(neighbor_N) then

                   j=j_min1+j_max_block+2

                   do i=1, 2*bc_size+1
                      new_grdpts_id(i,j) = exchange_pt
                   end do
                   
                   do i=2*bc_size+2, size(grdpts_id,1)
                      new_grdpts_id(i,j) = 
     $                     grdpts_id(i,j_max_block+2)
                   end do
                   
                   do i=size(grdpts_id,1)+1, i_max
                      new_grdpts_id(i,j) = no_pt
                   end do

                   j=j_min1+j_max_block+3

                   do i=1, 2*bc_size+1
                      new_grdpts_id(i,j) = exchange_pt
                   end do

                   do i=2*bc_size+2, i_max
                      new_grdpts_id(i,j) = no_pt
                   end do

                else

                   call add_layer2_in_E(
     $                  j_min1+j_max_block+2, new_grdpts_id,
     $                  j_max_block+2, grdpts_id,
     $                  i_max)

                   call add_layer1_out_E(
     $                  j_min1+j_max_block+3, new_grdpts_id, i_max)
                end if

             else

                !> update the upper border for block2 copied
                j_max_block = j_max

                !> copy the block corresponding to sublayer_
                !> first
                if(neighbor_N) then
                   
                   do j=j_min_block, j_max_block-bc_size
                      do i=1, size(grdpts_id,1)
                         new_grdpts_id(i,j_min1+j) = grdpts_id(i,j)
                      end do
                      
                      do i=size(grdpts_id,1)+1, i_max
                         new_grdpts_id(i,j_min1+j) = no_pt
                      end do
                   end do

                   do j=j_max_block-bc_size+1, j_max_block
                      do i=1, 2*bc_size+1
                         new_grdpts_id(i,j_min1+j) = exchange_pt
                      end do
                      
                      do i=2*bc_size+2, size(grdpts_id,1)
                         new_grdpts_id(i,j_min1+j) = grdpts_id(i,j)
                      end do
                      
                      do i=size(grdpts_id,1)+1, i_max
                         new_grdpts_id(i,j_min1+j) = no_pt
                      end do
                   end do

                else

                   do j=j_min_block, j_max_block
                      do i=1, size(grdpts_id,1)
                         new_grdpts_id(i,j_min1+j) = grdpts_id(i,j)
                      end do
                      
                      do i=size(grdpts_id,1)+1, i_max
                         new_grdpts_id(i,j_min1+j) = no_pt
                      end do
                   end do

                end if
                
             end if

          end if

        end subroutine adapt_grdpts_id_E_to_corner


        !< adapt the grdpts_id of the westhern buffer layers
        !> to the corner
        subroutine adapt_grdpts_id_W_to_corner(
     $     new_grdpts_id,
     $     grdpts_id,
     $     interior_j_max1, interior_j_max2,
     $     j_min1, j_min2,
     $     i_max,
     $     neighbor_N,
     $     neighbor_S)
        
          implicit none

          integer       , dimension(:,:), intent(out) :: new_grdpts_id
          integer       , dimension(:,:), intent(in)  :: grdpts_id
          integer(ikind)                , intent(in)  :: interior_j_max1
          integer(ikind)                , intent(in)  :: interior_j_max2
          integer(ikind)                , intent(in)  :: j_min1
          integer(ikind)                , intent(in)  :: j_min2
          integer(ikind)                , intent(in)  :: i_max
          logical                       , intent(in)  :: neighbor_N
          logical                       , intent(in)  :: neighbor_S
          
          integer(ikind) :: i,j, j_min_block, j_max_block
          integer(ikind) :: i_min, j_max

          i_min = size(new_grdpts_id,1) - size(grdpts_id,1)
          j_max = size(grdpts_id,2)
          

          !block1: nodes copied from the interior domain
          !        frontier between block1 and block2
          if(interior_j_max1.ge.2) then
             
             !> grdpt_id corresponding to nodes copied
             !> from the interior domain: 1st block
             !1st block : j=1, bc_size
             if(neighbor_S) then
                do j=1, bc_size
                   do i=1, i_max-(2*bc_size+1)
                      new_grdpts_id(i,j) = no_pt
                   end do
                   do i=i_max-(2*bc_size+1)+1, i_max
                      new_grdpts_id(i,j) = exchange_pt
                   end do
                end do
             else
                call add_layer1_out_W(
     $               1, new_grdpts_id, i_max)

                call add_layer2_out_W(
     $               2, new_grdpts_id, i_max)
             end if             

             call add_layer_interior_W(
     $            3, interior_j_max1, new_grdpts_id, i_max)

             !> frontier between the 1st and 2nd blocks
             call add_layer1_in_W(
     $            interior_j_max1+1, new_grdpts_id,
     $            1, grdpts_id,
     $            i_max)

             call add_layer1_in_W(
     $            interior_j_max1+2, new_grdpts_id,
     $            2, grdpts_id,
     $            i_max)

             !> update the lower border for the next block copied
             j_min_block = bc_size+1

          else

             if(interior_j_max1.eq.1) then
                
                !1st block : j=1, bc_size
                if(neighbor_S) then
                   j=1
                   do i=1, i_max-(2*bc_size+1)
                      new_grdpts_id(i,j) = no_pt
                   end do
                   do i=i_max-(2*bc_size+1)+1, i_max
                      new_grdpts_id(i,j) = exchange_pt
                   end do

                   j=2
                   !i_min = i_max-size(grdpts_id,1)

                   do i=1, i_min
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=i_min+1, i_max-(2*bc_size+1)
                      new_grdpts_id(i,j) = grdpts_id(i-i_min,j-1)
                   end do

                   do i=i_max-(2*bc_size+1)+1, i_max
                      new_grdpts_id(i,j) = exchange_pt
                   end do
                   
                else

                   !> grdpt_id corresponding to nodes copied
                   !> from the interior domain: 1st block
                   call add_layer1_out_W(
     $                  1, new_grdpts_id, i_max)
                   
                   !> frontier between the 1st and 2nd blocks
                   call add_layer2_in_W(
     $                  2, new_grdpts_id,
     $                  1, grdpts_id,
     $                  i_max)

                end if
                
                call add_layer1_in_W(
     $               3, new_grdpts_id,
     $               2, grdpts_id,
     $               i_max)

                !> update the lower border for the next block copied
                j_min_block = bc_size+1
                
             else

                if(neighbor_S) then
                   
                   !i_min = i_max-size(grdpts_id,1)

                   do j=1, bc_size
                      do i=1, i_min
                         new_grdpts_id(i,j) = no_pt
                      end do

                      do i=i_min+1, i_max-(2*bc_size+1)
                         new_grdpts_id(i,j) = grdpts_id(i-i_min,j)
                      end do
                      
                      do i=i_max-(2*bc_size+1)+1, i_max
                         new_grdpts_id(i,j) = exchange_pt
                      end do
                   end do

                   !> update the lower border for the next block copied
                   j_min_block = bc_size+1
                else
                
                   !> update the lower border for the next block copied
                   j_min_block = 1
                end if

             end if
          end if          

          
          !block4: grdpt copied from the sublayer_third
          !        frontier between block4 and block5
          !block5: grdpt corresponding to nodes copied
          !        from interior to the buffer layer
          !        frontier between block5 and the edge
          if(interior_j_max2.ge.2) then

             !> update the upper border for block2 copied
             j_max_block = j_max-bc_size

             !> copy the block corresponding to sublayer_
             !> first
             !i_min = i_max-size(grdpts_id,1)
             do j=j_min_block, j_max_block
                do i=1, i_min
                   new_grdpts_id(i,j_min1+j) = no_pt
                end do
                
                do i=i_min+1,i_max
                   new_grdpts_id(i,j_min1+j) = grdpts_id(i-i_min,j)
                end do
             end do

             !> frontier between block4 and block5
             call add_layer1_in_W(
     $            j_min1+j_max-1, new_grdpts_id,
     $            j_max_block+1, grdpts_id,
     $            i_max)

             call add_layer1_in_W(
     $            j_min1+j_max, new_grdpts_id,
     $            j_max_block+2, grdpts_id,
     $            i_max)

             call add_layer_interior_W(
     $            j_min2+1, j_min2+interior_j_max2-bc_size,
     $            new_grdpts_id,
     $            i_max)

             !> frontier at the edge
             if(neighbor_N) then
                do j=j_min2+interior_j_max2-bc_size+1, j_min2+interior_j_max2
                   do i=1, i_max-(2*bc_size+1)
                      new_grdpts_id(i,j) = no_pt
                   end do
                   do i=i_max-(2*bc_size+1)+1, i_max
                      new_grdpts_id(i,j) = exchange_pt
                   end do
                end do
             else
                call add_layer2_out_W(
     $               j_min2+interior_j_max2-1, new_grdpts_id, i_max)
                
                call add_layer1_out_W(
     $               j_min2+interior_j_max2  , new_grdpts_id, i_max)
             end if

          else

             if(interior_j_max2.eq.1) then
                
                !> update the upper border for block2 copied
                j_max_block = j_max-bc_size

                !> copy the block corresponding to sublayer_
                !> first
                !i_min = i_max-size(grdpts_id,1)
                do j=j_min_block, j_max_block
                   do i=1, i_min
                      new_grdpts_id(i,j_min1+j) = no_pt
                   end do

                   do i=i_min+1, i_max
                      new_grdpts_id(i,j_min1+j) = grdpts_id(i-i_min,j)
                   end do
                end do

                !> frontier between block4 and the edge
                call add_layer1_in_W(
     $               j_min1+j_max_block+1, new_grdpts_id,
     $               j_max_block+1, grdpts_id,
     $               i_max)

                !block5: last two layers
                if(neighbor_N) then

                   j=j_min1+j_max_block+2
                   !i_min = i_max-size(grdpts_id,1)

                   do i=1, i_min
                      new_grdpts_id(i,j) = no_pt
                   end do
                   
                   do i=i_min+1, i_max-(2*bc_size+1)
                      new_grdpts_id(i,j) = grdpts_id(i,j_max_block+2)
                   end do
                   
                   do i=i_max-(2*bc_size+1)+1, i_max
                      new_grdpts_id(i,j) = exchange_pt
                   end do

                   j=j_min1+j_max_block+3

                   do i=1, i_max-(2*bc_size+1)
                      new_grdpts_id(i,j) = no_pt
                   end do

                   do i=i_max-(2*bc_size+1)+1, i_max
                      new_grdpts_id(i,j) = exchange_pt
                   end do

                else

                   call add_layer2_in_W(
     $                  j_min1+j_max_block+2, new_grdpts_id,
     $                  j_max_block+2, grdpts_id,
     $                  i_max)

                   call add_layer1_out_W(
     $                  j_min1+j_max_block+3, new_grdpts_id, i_max)
                end if
                
             else

                !> update the upper border for block2 copied
                j_max_block = j_max

                !i_min = i_max-size(grdpts_id,1)

                if(neighbor_N) then

                   !> copy the block corresponding to sublayer_
                   !> third                   
                   do j=j_min_block, j_max_block-bc_size
                      do i=1, i_min
                         new_grdpts_id(i,j_min1+j) = no_pt
                      end do
                      
                      do i=i_min+1,i_max
                         new_grdpts_id(i,j_min1+j) = grdpts_id(i-i_min,j)
                      end do
                   end do

                   do j=j_max_block-bc_size+1, j_max_block
                      do i=1, i_min
                         new_grdpts_id(i,j_min1+j) = no_pt
                      end do
                      
                      do i=i_min+1,i_max-(2*bc_size+1)
                         new_grdpts_id(i,j_min1+j) = grdpts_id(i-i_min,j)
                      end do

                      do i=i_max-(2*bc_size+1)+1, i_max
                         new_grdpts_id(i,j_min1+j) = exchange_pt
                      end do
                   end do

                else                   
                   
                   !> copy the block corresponding to sublayer_
                   !> third                   
                   do j=j_min_block, j_max_block
                      do i=1, i_min
                         new_grdpts_id(i,j_min1+j) = no_pt
                      end do
                      
                      do i=i_min+1,i_max
                         new_grdpts_id(i,j_min1+j) = grdpts_id(i-i_min,j)
                      end do
                   end do

                end if
                
             end if

          end if

        end subroutine adapt_grdpts_id_W_to_corner


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param new_grdpts_id
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer1_out_E(
     $     j, temp_grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(out) :: temp_grdptid
          integer(ikind)         , intent(in)  :: i_max

          integer(ikind) :: i

          ! exchange
          !    |  bc_pt   no_pt
          !    |   _|_   ___|____
          !   / \ /   \ /        \
          !  ---------------------
          !  |0 0 2 2 2          |
          !  ---------------------

          do i=1, bc_size
             temp_grdptid(i,j) = exchange_pt
          end do

          do i=bc_size+1, 2*bc_size+1
             temp_grdptid(i,j) = bc_pt
          end do

          do i=2*bc_size+2, i_max
             temp_grdptid(i,j) = no_pt
          end do             

        end subroutine add_layer1_out_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer2_out_E(
     $     j, temp_grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(out) :: temp_grdptid
          integer(ikind)         , intent(in)  :: i_max

          integer(ikind) :: i

          !
          ! exchange
          !    |  bc_interior_pt
          !    |   |        
          !    |   | bc_pt no_pt
          !    |   |  |  ___|____
          !   / \ / \ | /        \
          !  ---------------------
          !  |0 0 1 1 2          |
          !  ---------------------

          do i=1, bc_size
             temp_grdptid(i,j) = exchange_pt
          end do

          do i=bc_size+1, 2*bc_size
             temp_grdptid(i,j) = bc_interior_pt
          end do

          i=2*bc_size+1
          temp_grdptid(i,j) = bc_pt

          do i=2*bc_size+2, i_max
             temp_grdptid(i,j) = no_pt
          end do             

        end subroutine add_layer2_out_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer:
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j_min
        !> y-index identifying the lower border of the layers
        !> initialized
        !
        !>@param j_max
        !> y-index identifying the upper border of the layers
        !> initialized
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer_interior_E(
     $     j_min, j_max, temp_grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in)  :: j_min
          integer(ikind)         , intent(in)  :: j_max
          integer, dimension(:,:), intent(out) :: temp_grdptid
          integer(ikind)         , intent(in)  :: i_max


          integer(ikind) :: i,j


          ! exchange
          !    | interior_pt
          !    |  | bc_interior_pt
          !    |  | | bc_pt
          !    |  | | |   no_pt
          !    |  | | |  ___|____
          !   / \ | | | /        \
          !  ---------------------
          !  |0 0 1 2 3          |
          !  ---------------------

          do j=j_min, j_max
             do i=1, bc_size
                temp_grdptid(i,j) = exchange_pt
             end do

             i=bc_size+1
             temp_grdptid(i,j) = interior_pt

             i=bc_size+2
             temp_grdptid(i,j) = bc_interior_pt

             i=2*bc_size+1
             temp_grdptid(i,j) = bc_pt

             do i=2*bc_size+2, i_max
                temp_grdptid(i,j) = no_pt
             end do
          end do

        end subroutine add_layer_interior_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer:
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param j_match
        !> y-index identifying the layer corresponding in the
        !> sublayer copied for merge
        !
        !>@param grdptid
        !> table for the gridpoint ID of the sublayer copied for
        !> merge
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer1_in_E(
     $     j, temp_grdptid, j_match, grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in) :: j
          integer, dimension(:,:), intent(out):: temp_grdptid
          integer(ikind)         , intent(in) :: j_match
          integer, dimension(:,:), intent(in) :: grdptid
          integer(ikind)         , intent(in) :: i_max

          integer(ikind) :: i


          ! exchange
          !    | interior
          !    |  | bc_interior
          !    |  | | 
          !    |  | |     grdptid   no_pt
          !    |  | |  _____|____  ___|___
          !   / \ | | /          \/       \
          !  ------------------------------
          !  |0 0 1 2                     |
          !  ------------------------------


          do i=1, bc_size
             temp_grdptid(i,j) = exchange_pt
          end do

          i=bc_size+1
          temp_grdptid(i,j) = interior_pt
          
          i=bc_size+2
          temp_grdptid(i,j) = bc_interior_pt
          
          do i=bc_size+3, size(grdptid,1)
             temp_grdptid(i,j) = grdptid(i,j_match)
          end do
          
          do i=size(grdptid,1)+1, i_max
             temp_grdptid(i,j) = no_pt
          end do
          
        end subroutine add_layer1_in_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param j_match
        !> y-index identifying the layer corresponding in the
        !> sublayer copied for merge
        !
        !>@param grdptid
        !> table for the gridpoint ID of the sublayer copied for
        !> merge
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer2_in_E(
     $     j, temp_grdptid, j_match, grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in) :: j
          integer, dimension(:,:), intent(out):: temp_grdptid
          integer(ikind)         , intent(in) :: j_match
          integer, dimension(:,:), intent(in) :: grdptid
          integer(ikind)         , intent(in) :: i_max

          integer(ikind) :: i


          ! exchange
          !    | 
          !    | bc_interior
          !    |   | 
          !    |   |      grdptid   no_pt
          !    |   |   _____|____  ___|___
          !   / \ / \ /          \/       \
          !  ------------------------------
          !  |0 0 1 1                     |
          !  ------------------------------

          do i=1, bc_size
             temp_grdptid(i,j) = exchange_pt
          end do

          do i=bc_size+1, bc_size+2
             temp_grdptid(i,j) = bc_interior_pt
          end do
          
          do i=bc_size+3, size(grdptid,1)
             temp_grdptid(i,j) = grdptid(i,j_match)
          end do
          
          do i=size(grdptid,1)+1, i_max
             temp_grdptid(i,j) = no_pt
          end do
          
        end subroutine add_layer2_in_E


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Western sublayer
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer1_out_W(
     $     j, temp_grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(out) :: temp_grdptid
          integer(ikind)         , intent(in)  :: i_max

          integer(ikind) :: i

          !            
          !                 exchange
          !    no_pt    bc_pt  |  
          !   ___|____   _|_   |  
          !  /        \ /   \ / \ 
          ! ----------------------
          ! |         | 2 2 2 0 0|
          ! ----------------------

          do i=1, i_max-(2*bc_size+1)
             temp_grdptid(i,j) = no_pt
          end do

          do i=i_max-(2*bc_size+1)+1, i_max-bc_size
             temp_grdptid(i,j) = bc_pt
          end do

          do i=i_max-bc_size+1, i_max
             temp_grdptid(i,j) = exchange_pt
          end do

        end subroutine add_layer1_out_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer2_out_W(
     $     j, temp_grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in)  :: j
          integer, dimension(:,:), intent(out) :: temp_grdptid
          integer(ikind)         , intent(in)  :: i_max

          integer(ikind) :: i


          !
          !        bc_interior_pt
          !              |
          !  no_pt  bc_pt|  exchange
          !  ___|____ |  |    |        
          ! /        \| / \  / \       
          ! ---------------------      
          ! |        |2 1 1 |0 0|     
          ! ---------------------

          do i=1, i_max-(2*bc_size+1)
             temp_grdptid(i,j) = no_pt
          end do

          i=i_max-(2*bc_size+1)+1
          temp_grdptid(i,j) = bc_pt

          do i=i_max-(2*bc_size+1)+2, i_max-bc_size
             temp_grdptid(i,j) = bc_interior_pt
          end do

          do i=i_max-bc_size+1, i_max
             temp_grdptid(i,j) = exchange_pt
          end do

        end subroutine add_layer2_out_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Eastern sublayer        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j_min
        !> y-index identifying the lower border of the layers
        !> initialized
        !
        !>@param j_max
        !> y-index identifying the upper border of the layers
        !> initialized
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer_interior_W(
     $     j_min, j_max, temp_grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in)  :: j_min
          integer(ikind)         , intent(in)  :: j_max
          integer, dimension(:,:), intent(out) :: temp_grdptid
          integer(ikind)         , intent(in)  :: i_max


          integer(ikind) :: i,j
          
          !                exchange_pt  
          !         interior_pt|        
          !                 |  |        
          !   bc_interior_pt|  |        
          !    no_pt      | |  |        
          !   ___|__ bc_pt  |  |        
          !  /        \ | | | / \       
          ! -----------------------     
          !           | 3 2 1 0 0 |     
          ! -----------------------     

          do j=j_min, j_max
             do i=1, i_max-bc_size-3
                temp_grdptid(i,j) = no_pt
             end do
             
             i=i_max-bc_size-2
             temp_grdptid(i,j) = bc_pt

             i=i_max-bc_size-1
             temp_grdptid(i,j) = bc_interior_pt

             i=i_max-bc_size
             temp_grdptid(i,j) = interior_pt

             do i=i_max-bc_size+1, i_max
                temp_grdptid(i,j) = exchange_pt
             end do
             
          end do

        end subroutine add_layer_interior_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Western sublayer        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param j_match
        !> y-index identifying the layer corresponding in the
        !> sublayer copied for merge
        !
        !>@param grdptid
        !> table for the gridpoint ID of the sublayer copied for
        !> merge
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer1_in_W(
     $     j, temp_grdptid, j_match, grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in) :: j
          integer, dimension(:,:), intent(out):: temp_grdptid
          integer(ikind)         , intent(in) :: j_match
          integer, dimension(:,:), intent(in) :: grdptid
          integer(ikind)         , intent(in) :: i_max

          integer(ikind) :: i, i_min


          !
          !                      exchange
          !                  interior | 
          !            bc_interior |  | 
          !                      | |  | 
          !   no_pt      grdptid | |  | 
          !  ___|___  _____|____ | |  | 
          ! /       \/          \| | / \
          ! ------------------------------
          ! |        |           |2 1 0 0|
          ! ------------------------------

          i_min = i_max-size(grdptid,1)

          do i=1, i_min
             temp_grdptid(i,j) = no_pt
          end do

          do i=i_min+1, i_max-bc_size-2
             temp_grdptid(i,j) = grdptid(i-i_min,j_match)
          end do

          i=i_max-bc_size-1
          temp_grdptid(i,j) = bc_interior_pt

          i=i_max-bc_size
          temp_grdptid(i,j) = interior_pt

          do i=i_max-bc_size+1, i_max
             temp_grdptid(i,j) = exchange_pt
          end do
          
        end subroutine add_layer1_in_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine initializing a layer for the gridpoint
        !> ID table of an Western sublayer
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param j
        !> y-index identifying the layer
        !
        !>@param temp_grdptid
        !> temporary table for the gridpoint ID of the merged
        !> sublayer
        !
        !>@param j_match
        !> y-index identifying the layer corresponding in the
        !> sublayer copied for merge
        !
        !>@param grdptid
        !> table for the gridpoint ID of the sublayer copied for
        !> merge
        !
        !>@param i_max
        !> x-size of the largest sublayer
        !--------------------------------------------------------------
        subroutine add_layer2_in_W(
     $     j, temp_grdptid, j_match, grdptid, i_max)

          implicit none

          integer(ikind)         , intent(in) :: j
          integer, dimension(:,:), intent(out):: temp_grdptid
          integer(ikind)         , intent(in) :: j_match
          integer, dimension(:,:), intent(in) :: grdptid
          integer(ikind)         , intent(in) :: i_max

          integer(ikind) :: i, i_min


          !                       exchange
          !                           |
          !              bc_interior  |       
          !                       |   |  
          !   no_pt      grdptid  |   |  
          !  ___|___  _____|____  |   |  
          ! /       \/          \/ \ / \ 
          ! -----------------------------
          ! |                   |1 1 0 0|
          ! -----------------------------


          i_min = i_max-size(grdptid,1)

          do i=1, i_min
             temp_grdptid(i,j) = no_pt
          end do

          do i=i_min+1, i_max-bc_size-2
             temp_grdptid(i,j) = grdptid(i-i_min,j_match)
          end do          

          do i=i_max-bc_size-1, i_max-bc_size
             temp_grdptid(i,j) = bc_interior_pt
          end do

          do i=i_max-bc_size+1, i_max
             temp_grdptid(i,j) = exchange_pt
          end do          
          
        end subroutine add_layer2_in_W


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine adding the exchange layer
        !> for Northern and Southern sublayers        
        !
        !> @date
        !> 15_05_2014 - initial version - J.L. Desmarais
        !
        !>@param new_grdpts_id
        !> table for the reallocated sublayer
        !
        !>@param grdpts_id
        !> table copied from the original sublayer
        !
        !>@param j_min_exchange
        !> y-lower border for the layer
        !
        !>@param j_max_exchange
        !> y-upper border for the layer
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the original
        !> sublayer (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the original
        !> sublayer (see Fig. 1)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param j_min1
        !> y-index identifying the correspondance between the merged
        !> sublayer and the first sublayer copied
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !--------------------------------------------------------------
        subroutine add_exchange_pt_layer_NS(
     $     new_grdpts_id,
     $     grdpts_id,
     $     j_min_exchange, j_max_exchange,
     $     interior_i_max1,
     $     interior_i_max2,
     $     i_min1, j_min1,
     $     i_min2)

          implicit none
          
          integer, dimension(:,:), intent(out) :: new_grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer(ikind)         , intent(in)  :: j_min_exchange
          integer(ikind)         , intent(in)  :: j_max_exchange
          integer(ikind)         , intent(in)  :: interior_i_max1
          integer(ikind)         , intent(in)  :: interior_i_max2
          integer(ikind)         , intent(in)  :: i_min1
          integer(ikind)         , intent(in)  :: j_min1
          integer(ikind)         , intent(in)  :: i_min2

          integer(ikind) :: i,j


          !interior_i_max1    interior_i_max2
          !  ___|____            ___|____  
          ! /        \          /        \ 
          ! ------------------------------ j_max_exchange
          ! |exchange| grdptid1 |exchange| 
          ! ------------------------------ j_min_exchange
          ! |        |          |
          ! 1       i_min1    i_min2

          do j=j_min_exchange, j_max_exchange
             do i=1, interior_i_max1
                new_grdpts_id(i,j) =
     $               exchange_pt
             end do
             
             do i=1, size(grdpts_id,1)
                new_grdpts_id(i_min1+i,j) =
     $               grdpts_id(i,j_min1+j)
             end do

             do i=1, interior_i_max2
                new_grdpts_id(i_min2+i,j) =
     $               exchange_pt
             end do
          end do

        end subroutine add_exchange_pt_layer_NS


      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine adding the interior layer
        !> for Northern and Southern sublayers
        !  
        !> @date
        !> 30_04_2014 - initial version    - J.L. Desmarais
        !> 15_05_2014 - adapted for corner - J.L. Desmarais
        !
        !>@param new_grdpts_id
        !> table for the sublayer adapted to the corner
        !
        !>@param grdpts_id
        !> table copied from the original sublayer
        !
        !>@param j
        !> y-index identifying the layer modified in new_grdpts_id
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the original
        !> sublayer (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the original
        !> sublayer (see Fig. 1)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param j_min1
        !> y-index identifying the correspondance between the merged
        !> sublayer and the first sublayer copied
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !--------------------------------------------------------------
        subroutine add_interior_pt_layer_NS(
     $     new_grdpts_id,
     $     grdpts_id,
     $     j,
     $     interior_i_max1,
     $     interior_i_max2,
     $     i_min1, j_min1,
     $     i_min2,
     $     neighbor_E,
     $     neighbor_W)

          implicit none

          integer, dimension(:,:), intent(out) :: new_grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer(ikind)         , intent(in)  :: j
          integer(ikind)         , intent(in)  :: interior_i_max1
          integer(ikind)         , intent(in)  :: interior_i_max2
          integer(ikind)         , intent(in)  :: i_min1, j_min1
          integer(ikind)         , intent(in)  :: i_min2
          logical                , intent(in)  :: neighbor_E
          logical                , intent(in)  :: neighbor_W


          integer(ikind) :: i


          !interior_i_max1  interior_i_max2
          !  ___|____          ___|____  
          ! /        \ i_max1 /        \
          ! ----------------------------
          ! |interior|grdptid1|interior|
          ! ----------------------------
          ! |        | |    | |
          ! 1   i_min1 |    i_min2
          !            |    |            
          !    i_max_first_layer       
          !                 |
          !         i_max_second_layer
          !--------------------------------
          !            Fig. 5
          !--------------------------------


          if(neighbor_W) then
             new_grdpts_id(1,j) = exchange_pt
             new_grdpts_id(2,j) = exchange_pt
          else
             new_grdpts_id(1,j) = bc_pt
             new_grdpts_id(2,j) = bc_interior_pt
          end if

          do i=3, interior_i_max1+2
             new_grdpts_id(i,j) = interior_pt
          end do

          do i=3, size(grdpts_id,1)-2
             new_grdpts_id(i_min1+i,j) =
     $            grdpts_id(i,j_min1+j)
          end do

          do i=-1, interior_i_max2-2
             new_grdpts_id(i_min2+i,j) =
     $            interior_pt
          end do

          if(neighbor_E) then
             new_grdpts_id(i_min2+interior_i_max2-1,j) = exchange_pt
             new_grdpts_id(i_min2+interior_i_max2  ,j) = exchange_pt
          else
             new_grdpts_id(i_min2+interior_i_max2-1,j) = bc_interior_pt
             new_grdpts_id(i_min2+interior_i_max2  ,j) = bc_pt
          end if

        end subroutine add_interior_pt_layer_NS


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine adding the bc_interior layer
        !> for Northern and Southern sublayers        
        !  
        !> @date
        !> 30_04_2014 - initial version  - J.L. Desmarais
        !> 15_05_2014 - adapt for corner - J.L. Desmarais
        !
        !>@param new_grdpts_id
        !> table for the sublayer adapted to the corner
        !
        !>@param grdpts_id
        !> table copied from the original sublayer
        !
        !>@param j
        !> y-index identifying the layer modified in new_grdpts_id
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the original
        !> sublayer (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the original
        !> sublayer (see Fig. 1)
        !
        !>@param i_min_second_layer
        !> index indicating the lower border for the second layer
        !> above the exchange layer
        !>  (see Fig. 1)
        !
        !>@param i_max_second_layer
        !> index indicating the upper border for the second layer
        !> above the exchange layer
        !>  (see Fig. 1)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param j_min1
        !> y-index identifying the correspondance between the adapted
        !> sublayer and the original sublayer copied
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !--------------------------------------------------------------
        subroutine add_bc_interior_pt_layer_NS(
     $     new_grdpts_id,
     $     grdpts_id,
     $     j,
     $     interior_i_max1,
     $     interior_i_max2,
     $     i_min1, j_min1,
     $     i_min2,
     $     neighbor_E,
     $     neighbor_W)

          implicit none

          integer, dimension(:,:), intent(out) :: new_grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer(ikind)         , intent(in)  :: j
          integer(ikind)         , intent(in)  :: interior_i_max1
          integer(ikind)         , intent(in)  :: interior_i_max2
          integer(ikind)         , intent(in)  :: i_min1, j_min1
          integer(ikind)         , intent(in)  :: i_min2
          logical                , intent(in)  :: neighbor_E
          logical                , intent(in)  :: neighbor_W

          integer(ikind) :: i


          ! interior_i_max1      interior_i_max2
          !  _____|_____          _____|_____ 
          ! /           \ i_max1 /           \
          ! ----------------------------------
          ! |bc_interior|grdptid1|bc_interior|
          ! ----------------------------------
          ! |           |        |         
          ! 1      i_min1        i_min2      
          !------------------------------------
          !               Fig. 6
          !------------------------------------

          !block1 + block2
          if(neighbor_W) then

             do i=1, bc_size
                new_grdpts_id(i,j) = exchange_pt
             end do

             do i=bc_size+1, interior_i_max1+bc_size-1
                new_grdpts_id(i,j) = bc_interior_pt
             end do

             if(interior_i_max1.ge.1) then
                do i=bc_size, size(grdpts_id,1)-bc_size+1
                   new_grdpts_id(i_min1+i,j)=
     $                  grdpts_id(i,j_min1+j)
                end do
             else
                do i=bc_size+1, size(grdpts_id,1)-bc_size+1
                   new_grdpts_id(i_min1+i,j)=
     $                  grdpts_id(i,j_min1+j)
                end do
             end if
                
          else

             do i=1, bc_size-1
                new_grdpts_id(i,j) = bc_pt
             end do
             
             do i=bc_size, interior_i_max1+bc_size-1
                new_grdpts_id(i,j) =
     $               bc_interior_pt
             end do

             do i=bc_size, size(grdpts_id,1)-bc_size+1
                new_grdpts_id(i_min1+i,j)=
     $               grdpts_id(i,j_min1+j)
             end do
          end if

          !block4 + block5
          if(neighbor_E) then

             do i=2-bc_size, interior_i_max2-bc_size
                new_grdpts_id(i_min2+i,j) =
     $               bc_interior_pt
             end do

             do i=1, bc_size
                new_grdpts_id(i_min2+interior_i_max2-bc_size+i,j) = exchange_pt
             end do

          else
             do i=2-bc_size, interior_i_max2-bc_size+1
                new_grdpts_id(i_min2+i,j) =
     $               bc_interior_pt
             end do

             do i=2, bc_size
                new_grdpts_id(i_min2+interior_i_max2-bc_size+i,j) = bc_pt
             end do
          end if

        end subroutine add_bc_interior_pt_layer_NS


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine adding the bc_interior layer
        !> for Northern and Southern sublayers        
        !
        !> @date
        !> 30_04_2014 - initial version - J.L. Desmarais
        !
        !>@param new_grdpts_id
        !> table for the sublayer adapted to the corner
        !
        !>@param grdpts_id
        !> table copied from the original sublayer
        !
        !>@param j
        !> y-index identifying the layer modified in grdpts_id
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the
        !> original sublayer (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the
        !> original sublayer (see Fig. 1)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param j_min1
        !> y-index identifying the correspondance between the merged
        !> sublayer and the first sublayer copied
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !--------------------------------------------------------------
        subroutine add_bc_pt_layer_NS(
     $     new_grdpts_id,
     $     grdpts_id,
     $     j,
     $     interior_i_max1,
     $     interior_i_max2,
     $     i_min1, j_min1,
     $     i_min2,
     $     neighbor_E,
     $     neighbor_W)

          implicit none

          integer, dimension(:,:), intent(out) :: new_grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer(ikind)         , intent(in)  :: j
          integer(ikind)         , intent(in)  :: interior_i_max1
          integer(ikind)         , intent(in)  :: interior_i_max2
          integer(ikind)         , intent(in)  :: i_min1, j_min1
          integer(ikind)         , intent(in)  :: i_min2
          logical                , intent(in)  :: neighbor_E
          logical                , intent(in)  :: neighbor_W

          integer(ikind) :: i


          ! interior_i_max1      interior_i_max2
          !  _____|_____          _____|_____   
          ! /           \ i_max1 /           \
          ! ----------------------------------
          ! |   bc_pt   |grdptid1|    bc_pt  |
          ! ----------------------------------
          ! |           |        |           
          ! 1      i_min1      i_min2      

          if(neighbor_W) then
             do i=1, bc_size
                new_grdpts_id(i,j) = exchange_pt
             end do

             do i=bc_size+1, interior_i_max1
                new_grdpts_id(i,j) = bc_pt
             end do

             if(interior_i_max1.ge.1) then
                if(interior_i_max1.eq.1) then
                   do i=2, size(grdpts_id,1)
                      new_grdpts_id(i_min1+i,j)=
     $                     grdpts_id(i,j_min1+j)
                   end do
                else
                   do i=1, size(grdpts_id,1)
                      new_grdpts_id(i_min1+i,j)=
     $                     grdpts_id(i,j_min1+j)
                   end do
                end if
             else
                do i=bc_size+1, size (grdpts_id,1)
                   new_grdpts_id(i_min1+i,j)=
     $                  grdpts_id(i,j_min1+j)
                end do
             end if
          else
             do i=1, interior_i_max1
                new_grdpts_id(i,j) = bc_pt
             end do

             do i=1, size(grdpts_id,1)
                new_grdpts_id(i_min1+i,j)=
     $               grdpts_id(i,j_min1+j)
             end do
          end if          

          if(neighbor_E) then
             do i=1, interior_i_max2-bc_size
                new_grdpts_id(i_min2+i,j) = bc_pt
             end do

             do i=1, bc_size
                new_grdpts_id(i_min2+interior_i_max2-bc_size+i,j) = exchange_pt
             end do
          else
             do i=1, interior_i_max2
                new_grdpts_id(i_min2+i,j) = bc_pt
             end do
          end if

        end subroutine add_bc_pt_layer_NS


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine adding the no_pt
        !> for Northern and Southern sublayers
        !
        !> @date
        !> 30_04_2014 - initial version    - J.L. Desmarais
        !> 15_05_2014 - adapted for corner - J.L. Desmarais
        !
        !>@param new_grdpts_id
        !> table for the adapted sublayer
        !
        !>@param grdpts_id
        !> table copied from the original sublayer
        !
        !>@param j
        !> y-index identifying the layer modified in grdpts_id
        !
        !>@param interior_i_max1
        !> size of the first block copied from the inside if the
        !> alignment of the new sublayer is larger than the
        !> adapted sublayer (see Fig. 1)
        !
        !>@param interior_i_max2
        !> size of the second block copied from the inside if the
        !> alignment of the new sublayer is larger than the
        !> adapted sublayer (see Fig. 1)
        !
        !>@param i_min1
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block2 (see Fig. 1)
        !
        !>@param j_min1
        !> y-index identifying the correspondance between the merged
        !> sublayer and the first sublayer copied
        !
        !>@param i_min2
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block3 (see Fig. 1)
        !
        !>@param i_min3
        !> index indicating the gridpoint in the new sublayer nodes
        !> table corresponding to the start of block4 (see Fig. 1)
        !--------------------------------------------------------------
        subroutine add_no_pt_layer_NS(
     $     new_grdpts_id,
     $     grdpts_id,
     $     j_min_no_pt, j_max_no_pt,
     $     interior_i_max1,
     $     interior_i_max2,
     $     i_min1, j_min1,
     $     i_min2)

          implicit none

          integer, dimension(:,:), intent(out) :: new_grdpts_id
          integer, dimension(:,:), intent(in)  :: grdpts_id
          integer(ikind)         , intent(in)  :: j_min_no_pt, j_max_no_pt
          integer(ikind)         , intent(in)  :: interior_i_max1
          integer(ikind)         , intent(in)  :: interior_i_max2
          integer(ikind)         , intent(in)  :: i_min1, j_min1
          integer(ikind)         , intent(in)  :: i_min2

          integer(ikind) :: i,j


          !
          ! interior_i_max1      interior_i_max2
          !  _____|_____          _____|_____ 
          ! /           \        /           \
          ! ----------------------------------
          ! |   no_pt   |grdptid1|   no_pt   |
          ! ----------------------------------
          ! |           |        |
          ! 1         i_min1   i_min2

          do j=j_min_no_pt, j_max_no_pt
             do i=1, interior_i_max1
                new_grdpts_id(i,j)=
     $               no_pt
             end do
             
             do i=1, size(grdpts_id,1)
                new_grdpts_id(i_min1+i,j)=
     $               grdpts_id(i,j_min1+j)
             end do
             
             do i=1, interior_i_max2
                new_grdpts_id(i_min2+i,j) =
     $               no_pt
             end do
          end do

        end subroutine add_no_pt_layer_NS


      end module bf_layer_adapt_corner_module
