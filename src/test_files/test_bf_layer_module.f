      module test_bf_layer_module

        !use ifport

        use bf_layer_class               , only : bf_layer
c$$$        use bf_layer_update_grdpts_module, only : update_grdpts
        use parameters_bf_layer          , only : interior_pt, bc_interior_pt, bc_pt
        use parameters_constant          , only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_kind              , only : rkind, ikind
        use parameters_input             , only : nx,ny,ne,bc_size

        private
        public ::
     $       print_interior_data,
     $       bf_layer_test_allocation,
     $       bf_layer_test_reallocation,
     $       bf_layer_test_merge,
     $       bf_layer_test_copy_neighbors,
     $       test_bf_layer_local_coord,
     $       ini_nodes,
     $       ini_grdpts_id,
     $       ini_alignment_table,
     $       ini_neighbors_table,
     $       ini_general_coord

        contains


        subroutine print_interior_data(
     $       nodes, grdpts_id,
     $       filename_nodes, filename_grdpts_id, filename_sizes)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer    , dimension(:,:)  , intent(in) :: grdpts_id 
          character(*)                 , intent(in) :: filename_nodes
          character(*)                 , intent(in) :: filename_grdpts_id
          character(*)                 , intent(in) :: filename_sizes

          integer :: ios
          
          !nodes
          open(unit=1,
     $          file=filename_nodes,
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

           !gridpts_id
           open(unit=1,
     $          file=filename_grdpts_id,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) grdpts_id
              close(unit=1)
           else
              stop 'file opening pb'
           end if

           !sizes
           open(unit=1,
     $          file=filename_sizes,
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

        end subroutine print_interior_data


        subroutine bf_layer_test_allocation(
     $     bf_layer_tested,
     $     bf_layer_loc,
     $     alignment,
     $     nodes,
     $     sizes_filename,
     $     nodes_filename,
     $     grdid_filename)

          implicit none

          class(bf_layer)                  , intent(inout) :: bf_layer_tested
          integer                          , intent(in)    :: bf_layer_loc
          integer        , dimension(2,2)  , intent(in)    :: alignment
          real(rkind)    , dimension(:,:,:), intent(in)    :: nodes
          character(*)                     , intent(in)    :: sizes_filename
          character(*)                     , intent(in)    :: nodes_filename
          character(*)                     , intent(in)    :: grdid_filename

          call bf_layer_tested%ini(bf_layer_loc)
          call bf_layer_tested%allocate_bf_layer(nodes,
     $                                           alignment)
          call bf_layer_tested%print_binary(nodes_filename,
     $                                      grdid_filename,
     $                                      sizes_filename)

        end subroutine bf_layer_test_allocation


        subroutine bf_layer_test_reallocation(
     $     bf_layer_tested,
     $     interior_nodes,
     $     new_alignment,
     $     sizes_filename,
     $     nodes_filename,
     $     grdid_filename)

          implicit none

          class(bf_layer)                  , intent(inout) :: bf_layer_tested
          real(rkind) , dimension(nx,ny,ne), intent(in)    :: interior_nodes
          integer     , dimension(2,2)     , intent(in)    :: new_alignment
          character(*)                     , intent(in)    :: sizes_filename
          character(*)                     , intent(in)    :: nodes_filename
          character(*)                     , intent(in)    :: grdid_filename

          call bf_layer_tested%reallocate_bf_layer(interior_nodes, new_alignment)
          call bf_layer_tested%print_binary(nodes_filename,
     $                                      grdid_filename,
     $                                      sizes_filename)
          
        end subroutine bf_layer_test_reallocation


        subroutine bf_layer_test_merge(
     $     bf_layer_tested,
     $     bf_layer_tested2,
     $     nodes,
     $     new_alignment,
     $     sizes_filename,
     $     nodes_filename,
     $     grdid_filename)

          implicit none

          class(bf_layer)                 , intent(inout) :: bf_layer_tested
          class(bf_layer)                 , intent(inout) :: bf_layer_tested2
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer    , dimension(2,2)     , intent(in)    :: new_alignment
          character(*)                    , intent(in)    :: sizes_filename
          character(*)                    , intent(in)    :: nodes_filename
          character(*)                    , intent(in)    :: grdid_filename

          call bf_layer_tested%merge_bf_layer(bf_layer_tested2,
     $                                        nodes,
     $                                        new_alignment)

          call bf_layer_tested%print_binary(nodes_filename,
     $                                      grdid_filename,
     $                                      sizes_filename)
          
        end subroutine bf_layer_test_merge



        subroutine bf_layer_test_copy_neighbors(
     $     bf_layer_tested, neighbor1, neighbor2,
     $     bf_layer_char, neighbor1_char, neighbor2_char,
     $     file_index)

          implicit none

          type(bf_layer), intent(inout) :: bf_layer_tested
          type(bf_layer), intent(inout) :: neighbor1
          type(bf_layer), intent(inout) :: neighbor2
          character(2)  , intent(in)    :: bf_layer_char
          character(2)  , intent(in)    :: neighbor1_char
          character(2)  , intent(in)    :: neighbor2_char
          integer       , intent(in)    :: file_index

          real(rkind) :: mainlayer_id
          real(rkind) :: neighbor1_id
          real(rkind) :: neighbor2_id

          character(len=21) :: sizes_filename
          character(len=21) :: nodes_filename
          character(len=21) :: grdid_filename


          !id for the buffer layers
          mainlayer_id = 0.2
          neighbor1_id = 0.5
          neighbor2_id = 0.9

          !initialize the nodes of the main layer+neighbors1+neighbors2
          call ini_cst_nodes(bf_layer_tested, mainlayer_id)
          call ini_cst_nodes(neighbor1, neighbor1_id)
          call ini_cst_nodes(neighbor2, neighbor2_id)

          !exchange between the main buffer layer and the neighbors
          call bf_layer_tested%copy_from_neighbor1(neighbor1)
          call bf_layer_tested%copy_to_neighbor1(neighbor1)
          call bf_layer_tested%copy_from_neighbor2(neighbor2)
          call bf_layer_tested%copy_to_neighbor2(neighbor2)


          !print the content of the buffer layers
          !print main buffer layer
          write(sizes_filename,'(A2,''1_sizes'',I1,''.dat'')') bf_layer_char, file_index
          write(nodes_filename,'(A2,''1_nodes'',I1,''.dat'')') bf_layer_char, file_index
          write(grdid_filename,'(A2,''1_grdpt_id'',I1,''.dat'')') bf_layer_char, file_index

          call bf_layer_tested%print_binary(
     $         nodes_filename,
     $         grdid_filename,
     $         sizes_filename)


          !print neighbor1
          write(sizes_filename,'(A2,''1_sizes'',I1,''.dat'')') neighbor1_char, file_index
          write(nodes_filename,'(A2,''1_nodes'',I1,''.dat'')') neighbor1_char, file_index
          write(grdid_filename,'(A2,''1_grdpt_id'',I1,''.dat'')') neighbor1_char, file_index

          call neighbor1%print_binary(
     $         nodes_filename,
     $         grdid_filename,
     $         sizes_filename)


          !print neighbor2
          write(sizes_filename,'(A2,''1_sizes'',I1,''.dat'')') neighbor2_char, file_index
          write(nodes_filename,'(A2,''1_nodes'',I1,''.dat'')') neighbor2_char, file_index
          write(grdid_filename,'(A2,''1_grdpt_id'',I1,''.dat'')') neighbor2_char, file_index

          call neighbor2%print_binary(
     $         nodes_filename,
     $         grdid_filename,
     $         sizes_filename)          

        end subroutine bf_layer_test_copy_neighbors



c$$$        subroutine bf_layer_test_update_grdpts(
c$$$     $     bf_layer_tested,
c$$$     $     selected_grdpts,
c$$$     $     match_table,
c$$$     $     sizes_filename,
c$$$     $     nodes_filename,
c$$$     $     grdid_filename)
c$$$
c$$$          implicit none
c$$$
c$$$          class(bf_layer)                  , intent(inout) :: bf_layer_tested
c$$$          integer     , dimension(:,:)     , intent(in)    :: selected_grdpts
c$$$          integer     , dimension(2)       , intent(in)    :: match_table
c$$$          character(*)                     , intent(in)    :: sizes_filename
c$$$          character(*)                     , intent(in)    :: nodes_filename
c$$$          character(*)                     , intent(in)    :: grdid_filename
c$$$
c$$$          call update_grdpts(bf_layer_tested,
c$$$     $         selected_grdpts,
c$$$     $         match_table)
c$$$          call bf_layer_tested%print_sizes(sizes_filename)
c$$$          call bf_layer_tested%print_nodes(nodes_filename)
c$$$          call bf_layer_tested%print_grdpts_id(grdid_filename)
c$$$          
c$$$        end subroutine bf_layer_test_update_grdpts


        subroutine ini_cst_nodes(bf_layer_initialized, cst)

          implicit none

          type(bf_layer), intent(inout) :: bf_layer_initialized
          real(rkind)   , intent(in)    :: cst
          
          integer(ikind), dimension(2) :: sizes
          integer(ikind)               :: i,j

          sizes = bf_layer_initialized%get_sizes()
          
          do j=1, sizes(2)
             do i=1, sizes(1)
                call bf_layer_initialized%set_nodes_pt(
     $               i,j,1,
     $               cst+0.1*real(i-1)/real(sizes(1)-1))
             end do
          end do          

        end subroutine ini_cst_nodes


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


        subroutine ini_grdpts_id(grdpts_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpts_id

          integer :: i,j

          j=1
          do i=1, nx
             grdpts_id(i,j) = bc_pt
          end do

          j=2
          grdpts_id(1,j) = bc_pt
          do i=2, nx-1
             grdpts_id(i,j) = bc_interior_pt
          end do
          grdpts_id(nx,j) = bc_pt


          do j=3, ny-bc_size
             i=1
             grdpts_id(i,j) = bc_pt

             i=2
             grdpts_id(i,j) = bc_interior_pt

             do i=3, nx-bc_size
                grdpts_id(i,j) = interior_pt
             end do

             i=nx-1
             grdpts_id(i,j) = bc_interior_pt

             i=nx
             grdpts_id(i,j) = bc_pt

          end do

          j=ny-1
          grdpts_id(1,j) = bc_pt
          do i=2, nx-1
             grdpts_id(i,j) = bc_interior_pt
          end do
          grdpts_id(nx,j) = bc_pt

          j=ny
          do i=1, nx
             grdpts_id(i,j) = bc_pt
          end do

        end subroutine ini_grdpts_id


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

          
          integer(ikind), dimension(2)                  :: new_sizes
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          
          
          !get the sizes of the nodes of the buffer layer
          new_sizes = bf_layer_tested%get_sizes()

          !allocate the new nodes to be created
          allocate(new_nodes(new_sizes(1), new_sizes(2), ne))

          !re-initialize the nodes of the buffer layer tested
          do k=1, ne
             do j=1, new_sizes(2)
                do i=1, new_sizes(1)
                   new_nodes(i,j,k) = 1.0
                end do
             end do
          end do

          !print the general coordinate tested in the nodes
          !and the buffer layer
          id = bf_layer_tested%get_localization()/10.0
          nodes(general_coord(1), general_coord(2), 1) = id
          local_coord = bf_layer_tested%get_local_coord(general_coord)
          new_nodes(local_coord(1),local_coord(2),1) = id
          call bf_layer_tested%set_nodes(new_nodes)

          !print tables
          call bf_layer_tested%print_binary(
     $         nodes_filename, grdid_filename, sizes_filename)

        end subroutine test_bf_layer_local_coord


        subroutine ini_general_coord(mainlayer_id, alignment, general_coord)

          implicit none

          integer                       , intent(in)  :: mainlayer_id
          integer       , dimension(2,2), intent(in)  :: alignment
          integer(ikind), dimension(2)  , intent(out) :: general_coord

          select case(mainlayer_id)
            case(N)
               general_coord(1) = alignment(1,1)
               general_coord(2) = ny-1
            case(S)
               general_coord(1) = alignment(1,1)
               general_coord(2) = 2
            case(E)
               general_coord(1) = nx-1
               general_coord(2) = alignment(1,2)
            case(W)
               general_coord(1) = 2
               general_coord(2) = alignment(1,2)
            case(N_E)
               general_coord(1) = nx-1
               general_coord(2) = ny-1
            case(N_W)
               general_coord(1) = 2
               general_coord(2) = ny-1
            case(S_E)
               general_coord(1) = nx-1
               general_coord(2) = 2
            case(S_W)
               general_coord(1) = 2
               general_coord(2) = 2
            case default
               print '(''test_bf_layer_prog'')'
               print '(''ini_general_coord'')'
               print '(''localization not recognized:'',I2)',
     $              mainlayer_id
               stop 'was the buffer layer initialized ?'
          end select

        end subroutine ini_general_coord

      end module test_bf_layer_module
