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
     $       ini_nodes

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

          call bf_layer_tested%ini([bf_layer_loc,1])
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

      end module test_bf_layer_module
