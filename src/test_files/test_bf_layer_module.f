      module test_bf_layer_module

        use bf_layer_class               , only : bf_layer
        use bf_layer_update_grdpts_module, only : update_grdpts
        use parameters_constant          , only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_kind              , only : rkind, ikind
        use parameters_input             , only : nx,ny,ne,bc_size

        private
        public ::
     $       print_nodes,
     $       print_sizes,
     $       bf_layer_test_allocation,
     $       bf_layer_test_reallocation,
     $       bf_layer_test_update_grdpts,
     $       test_bf_layer_local_coord,
     $       ini_nodes,
     $       ini_alignment_table,
     $       ini_neighbors_table

        contains


        subroutine print_nodes(nodes, filename)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          character(*)                 , intent(in) :: filename

          integer :: ios
          
          open(unit=1,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) nodes
              close(unit=1)
           else
              stop 'file opening pb'
           end if

        end subroutine print_nodes


        subroutine print_sizes(nodes, filename)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          character(*)                 , intent(in) :: filename

          integer :: ios
          
          open(unit=1,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) size(nodes,1),
     $             size(nodes,2),
     $             size(nodes,3)
              close(unit=1)
           else
              stop 'file opening pb'
           end if

        end subroutine print_sizes


        subroutine bf_layer_test_allocation(
     $     bf_layer_tested,
     $     bf_layer_loc,
     $     alignment,
     $     nodes,
     $     neighbors,
     $     sizes_filename,
     $     nodes_filename,
     $     grdid_filename)

          implicit none

          class(bf_layer)                  , intent(inout) :: bf_layer_tested
          integer                          , intent(in)    :: bf_layer_loc
          integer        , dimension(2,2)  , intent(in)    :: alignment
          logical        , dimension(4)    , intent(in)    :: neighbors
          real(rkind)    , dimension(:,:,:), intent(in)    :: nodes
          character(*)                     , intent(in)    :: sizes_filename
          character(*)                     , intent(in)    :: nodes_filename
          character(*)                     , intent(in)    :: grdid_filename

          call bf_layer_tested%ini(bf_layer_loc)
          call bf_layer_tested%allocate_bf_layer(
     $         alignment, nodes, neighbors)
          call bf_layer_tested%print_sizes(sizes_filename)
          call bf_layer_tested%print_nodes(nodes_filename)
          call bf_layer_tested%print_grdpts_id(grdid_filename)

        end subroutine bf_layer_test_allocation


        subroutine bf_layer_test_reallocation(
     $     bf_layer_tested,
     $     border_changes,
     $     match_table,
     $     sizes_filename,
     $     nodes_filename,
     $     grdid_filename)

          implicit none

          class(bf_layer)                  , intent(inout) :: bf_layer_tested
          integer     , dimension(2,2)     , intent(in)    :: border_changes
          integer     , dimension(2)       , intent(inout) :: match_table
          character(*)                     , intent(in)    :: sizes_filename
          character(*)                     , intent(in)    :: nodes_filename
          character(*)                     , intent(in)    :: grdid_filename

          call bf_layer_tested%reallocate_bf_layer(
     $         border_changes,
     $         match_table)
          call bf_layer_tested%print_sizes(sizes_filename)
          call bf_layer_tested%print_nodes(nodes_filename)
          call bf_layer_tested%print_grdpts_id(grdid_filename)
          
        end subroutine bf_layer_test_reallocation


        subroutine bf_layer_test_update_grdpts(
     $     bf_layer_tested,
     $     selected_grdpts,
     $     match_table,
     $     sizes_filename,
     $     nodes_filename,
     $     grdid_filename)

          implicit none

          class(bf_layer)                  , intent(inout) :: bf_layer_tested
          integer     , dimension(:,:)     , intent(in)    :: selected_grdpts
          integer     , dimension(2)       , intent(in)    :: match_table
          character(*)                     , intent(in)    :: sizes_filename
          character(*)                     , intent(in)    :: nodes_filename
          character(*)                     , intent(in)    :: grdid_filename

          call update_grdpts(bf_layer_tested,
     $         selected_grdpts,
     $         match_table)
          call bf_layer_tested%print_sizes(sizes_filename)
          call bf_layer_tested%print_nodes(nodes_filename)
          call bf_layer_tested%print_grdpts_id(grdid_filename)
          
        end subroutine bf_layer_test_update_grdpts


        subroutine ini_nodes(nodes)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes

          integer :: i,j,k


          call srand(10)

          do k=1, ne
             do j=1, ny
                do i=1, nx
                   nodes(i,j,k) = RAND() !(k-1)*1000 + (j-1)*100 + (i-1)
                end do
             end do
          end do

        end subroutine ini_nodes


        subroutine ini_alignment_table(test_alignment)
        
          implicit none

          integer, dimension(:,:,:), intent(inout) :: test_alignment

          test_alignment(1,1,1) = bc_size+9
          test_alignment(1,1,2) = bc_size+10
          test_alignment(1,2,1) = bc_size+1
          test_alignment(1,2,2) = bc_size+2

          test_alignment(2,1,1) = bc_size+1
          test_alignment(2,1,2) = bc_size+2
          test_alignment(2,2,1) = bc_size+16
          test_alignment(2,2,2) = bc_size+16
          
          test_alignment(3,1,1) = bc_size+16
          test_alignment(3,1,2) = bc_size+16
          test_alignment(3,2,1) = bc_size+9
          test_alignment(3,2,2) = bc_size+10
          
        end subroutine ini_alignment_table


        subroutine ini_neighbors_table(test_neighbors)
        
          implicit none

          logical, dimension(:,:), intent(inout) :: test_neighbors

          test_neighbors(1,:) = [.true.,.true.,.true.,.true.]
          test_neighbors(2,:) = [.true.,.true.,.true.,.true.]
          test_neighbors(3,:) = [.true.,.true.,.true.,.true.]
          
        end subroutine ini_neighbors_table


        subroutine test_bf_layer_local_coord(
     $     bf_layer_tested,
     $     nodes,
     $     general_coord,
     $     sizes_filename,
     $     nodes_filename,
     $     grdid_filename)

          implicit none

          type(bf_layer)                  , intent(inout) :: bf_layer_tested
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer(ikind), dimension(2)    , intent(in)    :: general_coord
          character(*)                    , intent(in)    :: sizes_filename
          character(*)                    , intent(in)    :: nodes_filename
          character(*)                    , intent(in)    :: grdid_filename

          integer(ikind) :: i,j,k
          real(rkind)    :: id
          integer(ikind), dimension(2) :: local_coord

          
          !re-initialize the nodes of the buffer layer tested
          do k=1, ne
             do j=1, size(bf_layer_tested%nodes,2)
                do i=1, size(bf_layer_tested%nodes,1)
                   bf_layer_tested%nodes(i,j,k) = 1.0
                end do
             end do
          end do

          !print the general coordinate tested in the nodes
          !and the buffer layer
          id = bf_layer_tested%localization/10.0
          nodes(general_coord(1), general_coord(2), 1) = id
          local_coord = bf_layer_tested%get_local_coord(general_coord)
          bf_layer_tested%nodes(local_coord(1),local_coord(2),1) = id

          !print tables
          call bf_layer_tested%print_sizes(sizes_filename)
          call bf_layer_tested%print_nodes(nodes_filename)
          call bf_layer_tested%print_grdpts_id(grdid_filename)

        end subroutine test_bf_layer_local_coord


      end module test_bf_layer_module
