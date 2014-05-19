c$$$          integer(ikind), dimension(2,2) :: border_changes
c$$$          integer(ikind) :: i,j
c$$$          integer        :: k
c$$$          integer(ikind) :: i_min, i_max, j_min, j_max
c$$$          integer(ikind) :: new_size_x, new_size_y
c$$$          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
c$$$          integer       , dimension(:,:)  , allocatable :: new_grdptid
c$$$          integer(ikind), dimension(2) :: match_table
c$$$
c$$$
c$$$          !< determine the border changes corresponding to the
c$$$          !> alignment asked
c$$$          border_changes(1,1) = alignment(1,1) - this%alignment(1,1)
c$$$          border_changes(2,1) = alignment(2,1) - this%alignment(2,1)
c$$$          border_changes(1,2) = alignment(1,2) - this%alignment(1,2)
c$$$          border_changes(2,2) = alignment(2,2) - this%alignment(2,2)
c$$$
c$$$
c$$$          !< determine the new alignment between the interior and 
c$$$          !> buffer layer tables
c$$$          select case(this%localization)
c$$$            case(N)
c$$$               border_changes(1,1) = alignment(1,1) - this%alignment(1,1)
c$$$               border_changes(2,1) = 0 !alignment(2,1) - this%alignment(2,1)
c$$$               border_changes(1,2) = alignment(1,2) - this%alignment(1,2)
c$$$               border_changes(2,2) = alignment(2,2) - this%alignment(2,2)
c$$$
c$$$               this%alignment(1,1) = alignment(1,1)
c$$$               !this%alignment(2,1) = alignment(2,1)
c$$$               this%alignment(1,2) = alignment(1,2)
c$$$               this%alignment(2,2) = alignment(2,2)
c$$$           
c$$$            case(S)
c$$$               border_changes(1,1) = alignment(1,1) - this%alignment(1,1)
c$$$               border_changes(2,1) = alignment(2,1) - this%alignment(2,1)
c$$$               border_changes(1,2) = alignment(1,2) - this%alignment(1,2)
c$$$               border_changes(2,2) = 0 !alignment(2,2) - this%alignment(2,2)
c$$$
c$$$               this%alignment(1,1) = alignment(1,1)
c$$$               this%alignment(2,1) = alignment(2,1)
c$$$               this%alignment(1,2) = alignment(1,2)
c$$$               !this%alignment(2,2) = alignment(2,2)
c$$$
c$$$            case(E)
c$$$               border_changes(1,1) = 0 !alignment(1,1) - this%alignment(1,1)
c$$$               border_changes(2,1) = alignment(2,1) - this%alignment(2,1)
c$$$               border_changes(1,2) = alignment(1,2) - this%alignment(1,2)
c$$$               border_changes(2,2) = alignment(2,2) - this%alignment(2,2)
c$$$
c$$$               !this%alignment(1,1) = alignment(1,1)
c$$$               this%alignment(2,1) = alignment(2,1)
c$$$               this%alignment(1,2) = alignment(1,2)
c$$$               this%alignment(2,2) = alignment(2,2)
c$$$
c$$$            case(W)
c$$$               border_changes(1,1) = alignment(1,1) - this%alignment(1,1)
c$$$               border_changes(2,1) = alignment(2,1) - this%alignment(2,1)
c$$$               border_changes(1,2) = 0 !alignment(1,2) - this%alignment(1,2)
c$$$               border_changes(2,2) = alignment(2,2) - this%alignment(2,2)
c$$$
c$$$               this%alignment(1,1) = alignment(1,1)
c$$$               this%alignment(2,1) = alignment(2,1)
c$$$               !this%alignment(1,2) = alignment(1,2)
c$$$               this%alignment(2,2) = alignment(2,2)
c$$$
c$$$            case default
c$$$               print '(''bf_layer_class'')'
c$$$               print '(''reallocate_bf_layer'')'
c$$$               print '(''localization not recognized: '')', this%localization
c$$$               stop 'wrong initialization'
c$$$
c$$$          end select
c$$$
c$$$
c$$$          !< determine the borders when filling the new table with the
c$$$          !> old data
c$$$          i_min = 1 + max(0,border_changes(1,1)) - border_changes(1,1)
c$$$          i_max = size(this%nodes,1) + min(0,border_changes(1,2)) - border_changes(1,1)
c$$$          j_min = 1 + max(0,border_changes(2,1)) - border_changes(2,1)
c$$$          j_max = size(this%nodes,2) + min(0,border_changes(2,2)) - border_changes(2,1)
c$$$
c$$$          match_table(1) = border_changes(1,1)
c$$$          match_table(2) = border_changes(2,1)
c$$$
c$$$
c$$$          !< determine the new size of the nodes and grdptid tables
c$$$          new_size_x = size(this%nodes,1) - border_changes(1,1) + border_changes(1,2)
c$$$          new_size_y = size(this%nodes,2) - border_changes(2,1) + border_changes(2,2)
c$$$          
c$$$
c$$$          !< allocate the new tables
c$$$          allocate(new_nodes(new_size_x,new_size_y,ne))
c$$$          allocate(new_grdptid(new_size_x,new_size_y))
c$$$
c$$$
c$$$          !fill the new nodes table with the previous data
c$$$          !and transfer the new table in the buffer layer
c$$$          do k=1, ne
c$$$             do j=j_min, j_max
c$$$                do i=i_min, i_max
c$$$                   new_nodes(i,j,k) = this%nodes(
c$$$     $                  match_table(1)+i,
c$$$     $                  match_table(2)+j,
c$$$     $                  k)
c$$$                end do
c$$$             end do
c$$$          end do
c$$$          call MOVE_ALLOC(new_nodes,this%nodes)
c$$$
c$$$
c$$$          !fill the new grdptid with the previous data
c$$$          !and transfer the new table in the buffer layer
c$$$          do j=1, j_min-1
c$$$             do i=1, size(this%nodes,1)
c$$$                new_grdptid(i,j) = no_pt
c$$$             end do
c$$$          end do
c$$$          
c$$$          do j=j_min, j_max
c$$$
c$$$             do i=1,i_min-1
c$$$                new_grdptid(i,j) = no_pt
c$$$             end do
c$$$
c$$$             do i=i_min, i_max
c$$$                new_grdptid(i,j) = this%grdpts_id(
c$$$     $               match_table(1)+i,
c$$$     $               match_table(2)+j)
c$$$             end do
c$$$
c$$$             do i=i_max+1, size(this%nodes,1)
c$$$                new_grdptid(i,j) = no_pt
c$$$             end do
c$$$
c$$$          end do          
c$$$
c$$$          do j=j_max+1,size(this%nodes,2)
c$$$             do i=1, size(this%nodes,1)
c$$$                new_grdptid(i,j) = no_pt
c$$$             end do
c$$$          end do
c$$$             
c$$$          call MOVE_ALLOC(new_grdptid,this%grdpts_id)
