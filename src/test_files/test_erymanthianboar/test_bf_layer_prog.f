      program test_bf_layer_prog

        use ifport

        use bf_layer_class, only :
     $     bf_layer

c$$$        use bf_layer_update_grdpts_module, only : update_grdpts

        use parameters_bf_layer, only :
     $       align_N,
     $       align_S,
     $       align_E,
     $       align_W

        use parameters_constant, only :
     $       N,
     $       S,
     $       E,
     $       W

        use parameters_kind, only :
     $       rkind,
     $       ikind

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       bc_size

        use test_bf_layer_module, only :
     $       print_interior_data,
     $       bf_layer_test_allocation,
     $       bf_layer_test_reallocation,
     $       bf_layer_test_merge,
     $       bf_layer_test_copy_neighbors,
     $       test_bf_layer_local_coord,
     $       ini_x_map,
     $       ini_y_map,
     $       ini_nodes,
     $       ini_grdpts_id,
     $       ini_general_coord

        implicit none

        !integer, parameter :: neighbor_case = 2
        integer, parameter :: size_case = 1
        integer, parameter :: distance_case = 1
        integer, parameter :: random_seed = 86456
        integer, parameter :: over_alignment_case_x = 1
        integer, parameter :: over_alignment_case_y = 1
        integer, parameter :: inverse_case = 1
        integer, parameter :: inverse_size_case = 1
        integer, parameter :: test_first_bf_layer_align_case = 0
        integer, parameter :: test_second_bf_layer_align_case = 1
        logical, parameter :: test_reallocation = .true.
        logical, parameter :: test_merge = .true.
        logical, parameter :: test_copy_with_neighbors=.true.


        type(bf_layer), dimension(8)   :: table_bf_layer_tested
        type(bf_layer), dimension(4)   :: table_1st_bf_layer_tested
        type(bf_layer), dimension(4)   :: table_2nd_bf_layer_tested
        type(bf_layer), dimension(4,3) :: table_bf_layer_copy_tested

        real(rkind)   , dimension(nx)       :: x_map
        real(rkind)   , dimension(ny)       :: y_map
        real(rkind)   , dimension(nx,ny,ne) :: nodes
        integer       , dimension(nx,ny)    :: grdpts_id
        integer(ikind), dimension(2,2)      :: alignment
c$$$        integer(ikind), dimension(2,2)      :: alignment_2nd_bf_layer
c$$$        integer(ikind), dimension(2,2)      :: alignment_after_merge
        integer(ikind)                      :: i !,j,k
        integer       , dimension(4)        :: bf_layer_loc
        character(2)  , dimension(4)        :: bf_layer_char
        character(len=21)                   :: sizes_filename
        character(len=21)                   :: x_map_filename
        character(len=21)                   :: y_map_filename
        character(len=21)                   :: nodes_filename
        character(len=21)                   :: grdid_filename
c$$$        integer(ikind), dimension(2,2)      :: border_changes
c$$$        logical       , dimension(4)        :: neighbors
c$$$        integer(ikind), dimension(2,2)      :: selected_grdpts
c$$$        integer       , dimension(2)        :: match_table
c$$$        integer(ikind), dimension(2)        :: general_coord

        integer(ikind), dimension(4,2,2)   :: test_alignment_reallocation
        integer(ikind), dimension(2,2)     :: alignment_reallocation
        integer(ikind), dimension(4,2,2)   :: test_alignment_merge
        integer(ikind), dimension(4,2,2)   :: test_alignment_2nd
        integer(ikind), dimension(4,3,2,2) :: test_alignment_copy_neighbors
        integer(ikind), dimension(2,2)     :: alignment_merge
        integer(ikind), dimension(2,2)     :: alignment_2nd

        integer :: relative_distance
        integer :: relative_size
        integer :: over_alignment_x
        integer :: over_alignment_y

        integer, dimension(4,2) :: neighbors

        neighbors(N,1) = W
        neighbors(N,2) = E
        neighbors(S,1) = W
        neighbors(S,2) = E
        neighbors(E,1) = S
        neighbors(E,2) = N
        neighbors(W,1) = S
        neighbors(W,2) = N

c$$$        integer, dimension(8,2,2) :: test_selected_grdpts

        call ini_x_map(x_map)
        call ini_y_map(y_map)
        call ini_nodes(nodes)
        call ini_grdpts_id(grdpts_id)

        !print the nodes
        call print_interior_data(
     $       x_map, y_map,
     $       nodes, grdpts_id,
     $       'interior_x_map.dat',
     $       'interior_y_map.dat',
     $       'interior_nodes.dat',
     $       'interior_grdpts_id.dat',
     $       'interior_sizes.dat')

        !buffer layers tested
        bf_layer_loc  = [N,S,E,W]
        bf_layer_char = ['N_','S_','E_','W_']


        !over alignments
        over_alignment_x = over_alignment_case_x
        over_alignment_y = over_alignment_case_y

        
        !alignment after reallocation
        call ini_alignment_after_reallocation(
     $       test_first_bf_layer_align_case,
     $       over_alignment_x,
     $       over_alignment_y,
     $       test_alignment_reallocation)

        !alignment after merge
        if(distance_case.ge.0) then
           relative_distance = distance_case
        else
           print '(''distance_case < 0'')'
           stop 'not possible'
        end if
        !call ini_using_case(distance_case, random_seed, relative_distance)
        call ini_using_case(size_case, random_seed, relative_size)
        
        !call ini_using_case(over_alignment_case_x, random_seed, over_alignment_x)
        !call ini_using_case(over_alignment_case_y, random_seed, over_alignment_y)

        call ini_alignment_2nd_bf_layer(
     $       test_first_bf_layer_align_case,
     $       test_second_bf_layer_align_case,
     $       relative_distance, relative_size, 
     $       over_alignment_x, over_alignment_y,
     $       test_alignment_2nd, test_alignment_merge)


        !neighbors
        !call ini_neighbors(neighbor_case, neighbors)
        
        !selected grid points
        !call ini_selected_grdpts(test_selected_grdpts)


        !tests on all the buffer layers
        do i=1, size(bf_layer_loc,1)

           !test allocation
           write(sizes_filename,'(A2,''1_sizes1.dat'')') bf_layer_char(i)
           write(x_map_filename,'(A2,''1_x_map1.dat'')') bf_layer_char(i)
           write(y_map_filename,'(A2,''1_y_map1.dat'')') bf_layer_char(i)
           write(nodes_filename,'(A2,''1_nodes1.dat'')') bf_layer_char(i)
           write(grdid_filename,'(A2,''1_grdpt_id1.dat'')') bf_layer_char(i)
        
           call ini_alignment(
     $          test_first_bf_layer_align_case,
     $          i, alignment)

           call bf_layer_test_allocation(
     $          table_bf_layer_tested(i),
     $          bf_layer_loc(i),
     $          alignment,
     $          x_map,
     $          y_map,
     $          nodes,
     $          sizes_filename,
     $          x_map_filename,
     $          y_map_filename,
     $          nodes_filename,
     $          grdid_filename)
           
           !test reallocation
           if(test_reallocation) then
              write(sizes_filename,'(A2,''1_sizes2.dat'')') bf_layer_char(i)
              write(x_map_filename,'(A2,''1_x_map2.dat'')') bf_layer_char(i)
              write(y_map_filename,'(A2,''1_y_map2.dat'')') bf_layer_char(i)
              write(nodes_filename,'(A2,''1_nodes2.dat'')') bf_layer_char(i)
              write(grdid_filename,'(A2,''1_grdpt_id2.dat'')') bf_layer_char(i)

              alignment_reallocation(1,1) = test_alignment_reallocation(bf_layer_loc(i),1,1)
              alignment_reallocation(1,2) = test_alignment_reallocation(bf_layer_loc(i),1,2)
              alignment_reallocation(2,1) = test_alignment_reallocation(bf_layer_loc(i),2,1)
              alignment_reallocation(2,2) = test_alignment_reallocation(bf_layer_loc(i),2,2)

              call bf_layer_test_reallocation(
     $             table_bf_layer_tested(i),
     $             x_map,
     $             y_map,
     $             nodes,
     $             alignment_reallocation,
     $             sizes_filename,
     $             x_map_filename,
     $             y_map_filename,
     $             nodes_filename,
     $             grdid_filename)
           end if

        end do
       


        !test merge
        if(test_merge) then
           do i=1,4

              !print the data of the first buffer layer
              write(sizes_filename,'(A2,''1_sizes3.dat'')') bf_layer_char(i)
              write(x_map_filename,'(A2,''1_x_map3.dat'')') bf_layer_char(i)
              write(y_map_filename,'(A2,''1_y_map3.dat'')') bf_layer_char(i)
              write(nodes_filename,'(A2,''1_nodes3.dat'')') bf_layer_char(i)
              write(grdid_filename,'(A2,''1_grdpt_id3.dat'')') bf_layer_char(i)

              call ini_alignment(
     $             test_first_bf_layer_align_case, i, alignment)

              call bf_layer_test_allocation(
     $             table_1st_bf_layer_tested(i),
     $             bf_layer_loc(i),
     $             alignment,
     $             x_map,
     $             y_map,
     $             nodes,
     $             sizes_filename,
     $             x_map_filename,
     $             y_map_filename,
     $             nodes_filename,
     $             grdid_filename)

              if(inverse_size_case.eq.1) then
                 if(relative_size.ne.0) then
                    call update_alignment_to_increase_size(
     $                   i,
     $                   relative_size,
     $                   alignment)

                    call bf_layer_test_reallocation(
     $                   table_1st_bf_layer_tested(i),
     $                   x_map,
     $                   y_map,
     $                   nodes,
     $                   alignment,
     $                   sizes_filename,
     $                   x_map_filename,
     $                   y_map_filename,
     $                   nodes_filename,
     $                   grdid_filename)
                 end if
              end if

              !create additional bf layer and print its data
              write(sizes_filename,'(A2,''2_sizes3.dat'')') bf_layer_char(i)
              write(x_map_filename,'(A2,''2_x_map3.dat'')') bf_layer_char(i)
              write(y_map_filename,'(A2,''2_y_map3.dat'')') bf_layer_char(i)
              write(nodes_filename,'(A2,''2_nodes3.dat'')') bf_layer_char(i)
              write(grdid_filename,'(A2,''2_grdpt_id3.dat'')') bf_layer_char(i)

              alignment_2nd(1,1) = test_alignment_2nd(bf_layer_loc(i),1,1)
              alignment_2nd(1,2) = test_alignment_2nd(bf_layer_loc(i),1,2)
              alignment_2nd(2,1) = test_alignment_2nd(bf_layer_loc(i),2,1)
              alignment_2nd(2,2) = test_alignment_2nd(bf_layer_loc(i),2,2)

              call bf_layer_test_allocation(
     $             table_2nd_bf_layer_tested(i),
     $             bf_layer_loc(i),
     $             alignment_2nd,
     $             x_map,
     $             y_map,
     $             nodes,
     $             sizes_filename,
     $             x_map_filename,
     $             y_map_filename,
     $             nodes_filename,
     $             grdid_filename)

              if(inverse_size_case.ne.1) then
                 if(relative_size.ne.0) then
                    call update_alignment_to_increase_size(
     $                   i,
     $                   relative_size,
     $                   alignment_2nd)

                    call bf_layer_test_reallocation(
     $                   table_2nd_bf_layer_tested(i),
     $                   x_map,
     $                   y_map,
     $                   nodes,
     $                   alignment_2nd,
     $                   sizes_filename,
     $                   x_map_filename,
     $                   y_map_filename,
     $                   nodes_filename,
     $                   grdid_filename)
                 end if
              end if

              !test the merging procedure
              write(sizes_filename,'(A2,''1_sizes4.dat'')') bf_layer_char(i)
              write(x_map_filename,'(A2,''1_x_map4.dat'')') bf_layer_char(i)
              write(y_map_filename,'(A2,''1_y_map4.dat'')') bf_layer_char(i)
              write(nodes_filename,'(A2,''1_nodes4.dat'')') bf_layer_char(i)
              write(grdid_filename,'(A2,''1_grdpt_id4.dat'')') bf_layer_char(i)

              alignment_merge(1,1) = test_alignment_merge(bf_layer_loc(i),1,1)
              alignment_merge(1,2) = test_alignment_merge(bf_layer_loc(i),1,2)
              alignment_merge(2,1) = test_alignment_merge(bf_layer_loc(i),2,1)
              alignment_merge(2,2) = test_alignment_merge(bf_layer_loc(i),2,2)

              if(inverse_case.eq.1) then
                 call bf_layer_test_merge(
     $                table_1st_bf_layer_tested(i),
     $                table_2nd_bf_layer_tested(i),
     $                x_map,
     $                y_map,
     $                nodes,
     $                alignment_merge,
     $                sizes_filename,
     $                x_map_filename,
     $                y_map_filename,
     $                nodes_filename,
     $                grdid_filename)

              else
                 call bf_layer_test_merge(
     $                table_2nd_bf_layer_tested(i),
     $                table_1st_bf_layer_tested(i),
     $                x_map,
     $                y_map,
     $                nodes,
     $                alignment_merge,
     $                sizes_filename,
     $                x_map_filename,
     $                y_map_filename,
     $                nodes_filename,
     $                grdid_filename)
              end if

           end do
        end if


        !test copy from the neighbors
        if(test_copy_with_neighbors) then

           call ini_alignment_for_neighbor_copy(
     $          test_alignment_copy_neighbors)

           do i=1,4


              !buffer layer
              alignment(1,1) = test_alignment_copy_neighbors(bf_layer_loc(i),1,1,1)
              alignment(1,2) = test_alignment_copy_neighbors(bf_layer_loc(i),1,1,2)
              alignment(2,1) = test_alignment_copy_neighbors(bf_layer_loc(i),1,2,1)
              alignment(2,2) = test_alignment_copy_neighbors(bf_layer_loc(i),1,2,2)

              call table_bf_layer_copy_tested(i,1)%ini(bf_layer_loc(i))
              call table_bf_layer_copy_tested(i,1)%allocate_bf_layer(x_map,y_map,nodes,alignment)


              !neighbor1
              alignment(1,1) = test_alignment_copy_neighbors(bf_layer_loc(i),2,1,1)
              alignment(1,2) = test_alignment_copy_neighbors(bf_layer_loc(i),2,1,2)
              alignment(2,1) = test_alignment_copy_neighbors(bf_layer_loc(i),2,2,1)
              alignment(2,2) = test_alignment_copy_neighbors(bf_layer_loc(i),2,2,2)

              call table_bf_layer_copy_tested(i,2)%ini(neighbors(i,1))
              call table_bf_layer_copy_tested(i,2)%allocate_bf_layer(x_map,y_map,nodes,alignment)


              !neighbor2
              alignment(1,1) = test_alignment_copy_neighbors(bf_layer_loc(i),3,1,1)
              alignment(1,2) = test_alignment_copy_neighbors(bf_layer_loc(i),3,1,2)
              alignment(2,1) = test_alignment_copy_neighbors(bf_layer_loc(i),3,2,1)
              alignment(2,2) = test_alignment_copy_neighbors(bf_layer_loc(i),3,2,2)

              call table_bf_layer_copy_tested(i,3)%ini(neighbors(i,2))
              call table_bf_layer_copy_tested(i,3)%allocate_bf_layer(x_map,y_map,nodes,alignment)


              !test copy neighbors
              call bf_layer_test_copy_neighbors(
     $             table_bf_layer_copy_tested(i,1),
     $             table_bf_layer_copy_tested(i,2),
     $             table_bf_layer_copy_tested(i,3),
     $             bf_layer_char(i),
     $             bf_layer_char(neighbors(i,1)),
     $             bf_layer_char(neighbors(i,2)),
     $             4+i)

           end do

        end if


c$$$        !test the local coordinates
c$$$        !--------------------------------
c$$$        !re-initialize the interior nodes
c$$$        do k=1, ne
c$$$           do j=1, ny
c$$$              do i=1, nx
c$$$                 nodes(i,j,k) = 1.0
c$$$              end do
c$$$           end do
c$$$        end do
c$$$
c$$$        !test on all the buffer layers
c$$$        do i=1, size(bf_layer_loc,1)
c$$$
c$$$           write(sizes_filename,'(A2,''1_sizes4.dat'')') bf_layer_char(i)
c$$$           write(nodes_filename,'(A2,''1_nodes4.dat'')') bf_layer_char(i)
c$$$           write(grdid_filename,'(A2,''1_grdpt_id4.dat'')') bf_layer_char(i)
c$$$
c$$$           call ini_general_coord(
c$$$     $          table_bf_layer_tested(i)%localization,
c$$$     $          table_bf_layer_tested(i)%alignment,
c$$$     $          general_coord)
c$$$
c$$$           call test_bf_layer_local_coord(
c$$$     $          table_bf_layer_tested(i),
c$$$     $          nodes,
c$$$     $          general_coord,
c$$$     $          sizes_filename,
c$$$     $          nodes_filename,
c$$$     $          grdid_filename)
c$$$           
c$$$        end do
c$$$
c$$$        !for the test local coordinates: the nodes should be
c$$$        !printed once all the buffer layers were tested
c$$$        call print_interior_data(
c$$$     $       nodes, grdpts_id,
c$$$     $       'interior_nodes4.dat',
c$$$     $       'interior_grdpts_id4.dat',
c$$$     $       'interior_sizes4.dat')

        contains

        subroutine ini_alignment(
     $       test_first_bf_layer_align_case, mainlayer_id, alignment)

          implicit none

          integer                       , intent(in)  :: test_first_bf_layer_align_case
          integer                       , intent(in)  :: mainlayer_id
          integer(ikind), dimension(2,2), intent(out) :: alignment


          integer(ikind) :: border_min,border_max

          
          !basic alignment depending on cardinal point
          select case(mainlayer_id)
            case(N)
               alignment(2,1) = align_N
               alignment(2,2) = align_N
            case(S)
               alignment(2,1) = align_S
               alignment(2,2) = align_S
            case(E)
               alignment(1,1) = align_E
               alignment(1,2) = align_E
            case(W)
               alignment(1,1) = align_W
               alignment(1,2) = align_W
          end select

          border_min = bc_size+1+test_first_bf_layer_align_case
          border_max = bc_size+7+test_first_bf_layer_align_case

          !distance from corner
c$$$          if(test_first_bf_layer_align_case.ge.-4) then
c$$$             border_min = bc_size+1+test_first_bf_layer_align_case
c$$$             border_max = bc_size+5+test_first_bf_layer_align_case
c$$$          else
c$$$             print '(''test_bf_layer_prog'')'
c$$$             print '(''ini_alignment'')'
c$$$             print '(''test_first_bf_layer_align_case not correct'')'
c$$$          end if
             
          !borders
          select case(mainlayer_id)
            case(N,S)
               alignment(1,1) = border_min
               alignment(1,2) = border_max
            case(E,W)
               alignment(2,1) = border_min
               alignment(2,2) = border_max
          end select
            
        end subroutine ini_alignment


        subroutine ini_alignment_after_reallocation(
     $     test_first_bf_layer_align_case,
     $     over_alignment_x,
     $     over_alignment_y,
     $     alignment)
        
          implicit none

          integer                  , intent(in)    :: test_first_bf_layer_align_case
          integer                  , intent(in)    :: over_alignment_x
          integer                  , intent(in)    :: over_alignment_y
          integer, dimension(:,:,:), intent(inout) :: alignment

          alignment(N,1,1) = bc_size+1+test_first_bf_layer_align_case-over_alignment_x
          alignment(N,1,2) = bc_size+7+test_first_bf_layer_align_case+over_alignment_x
          alignment(N,2,1) = align_N
          alignment(N,2,2) = align_N+over_alignment_y

          alignment(S,1,1) = bc_size+1+test_first_bf_layer_align_case-over_alignment_x
          alignment(S,1,2) = bc_size+7+test_first_bf_layer_align_case+over_alignment_x
          alignment(S,2,1) = align_S-over_alignment_y
          alignment(S,2,2) = align_S

          alignment(E,1,1) = align_E
          alignment(E,1,2) = align_E+over_alignment_y
          alignment(E,2,1) = bc_size+1+test_first_bf_layer_align_case-over_alignment_x
          alignment(E,2,2) = bc_size+7+test_first_bf_layer_align_case+over_alignment_x

          alignment(W,1,1) = align_W-over_alignment_y
          alignment(W,1,2) = align_W
          alignment(W,2,1) = bc_size+1+test_first_bf_layer_align_case-over_alignment_x
          alignment(W,2,2) = bc_size+7+test_first_bf_layer_align_case+over_alignment_x

        end subroutine ini_alignment_after_reallocation


        subroutine ini_alignment_2nd_bf_layer(
     $     test_first_bf_layer_alignment,
     $     test_second_bf_layer_alignment,
     $     relative_distance, relative_size, 
     $     over_alignment_x, over_alignment_y,
     $     alignment, final_alignment)

          implicit none

          integer                  , intent(in)    :: test_first_bf_layer_alignment
          integer                  , intent(in)    :: test_second_bf_layer_alignment
          integer                  , intent(in)    :: relative_distance
          integer                  , intent(in)    :: relative_size
          integer                  , intent(in)    :: over_alignment_x
          integer                  , intent(in)    :: over_alignment_y
          integer, dimension(:,:,:), intent(inout) :: alignment
          integer, dimension(:,:,:), intent(inout) :: final_alignment


          !basic alignment requirements needed by the 2nd buffer layer
          !depending on the cardinal point
          alignment(N,2,1) = align_N
          alignment(N,2,2) = align_N
          
          alignment(S,2,1) = align_S
          alignment(S,2,2) = align_S

          alignment(E,1,1) = align_E
          alignment(E,1,2) = align_E

          alignment(W,1,1) = align_W
          alignment(W,1,2) = align_W

          select case(test_second_bf_layer_alignment)

            !a bit further from the first buffer layer
            case(1)
               alignment(N,1,1) = bc_size+1+7+
     $              test_first_bf_layer_alignment+
     $              2*bc_size+relative_distance
               alignment(N,1,2) = alignment(N,1,1) + relative_size
          
               alignment(S,1,1) = bc_size+1+7+
     $              test_first_bf_layer_alignment+
     $              2*bc_size+relative_distance
               alignment(S,1,2) = alignment(S,1,1) + relative_size    

               alignment(E,2,1) = bc_size+7+1+
     $              test_first_bf_layer_alignment+
     $              2*bc_size+relative_distance
               alignment(E,2,2) = alignment(E,2,1) + relative_size

               alignment(W,2,1) = bc_size+7+1+
     $              test_first_bf_layer_alignment+
     $              2*bc_size+relative_distance
               alignment(W,2,2) = alignment(W,2,1) + relative_size

            !on the other corner - 2
            case(2)
               alignment(N,1,2) = nx-bc_size-2
               alignment(N,1,1) = alignment(N,1,2) - relative_size
               
               alignment(S,1,2) = nx-bc_size-2
               alignment(S,1,1) = alignment(S,1,2) - relative_size

               alignment(E,2,2) = ny-bc_size-2
               alignment(E,2,1) = alignment(E,2,2)-relative_size

               alignment(W,2,2) = ny-bc_size-2
               alignment(W,2,1) = alignment(E,2,2)-relative_size

            !on the other corner - 1
            case(3)
               alignment(N,1,2) = nx-bc_size-1
               alignment(N,1,1) = alignment(N,1,2) - relative_size
               
               alignment(S,1,2) = nx-bc_size-1
               alignment(S,1,1) = alignment(S,1,2) - relative_size

               alignment(E,2,2) = ny-bc_size-1
               alignment(E,2,1) = alignment(E,2,2)-relative_size

               alignment(W,2,2) = ny-bc_size-1
               alignment(W,2,1) = alignment(E,2,2)-relative_size

            !on the other corner
            case(4)
                alignment(N,1,2) = nx-bc_size
                alignment(N,1,1) = alignment(N,1,2) - relative_size
               
                alignment(S,1,2) = nx-bc_size
                alignment(S,1,1) = alignment(S,1,2) - relative_size
                
                alignment(E,2,2) = ny-bc_size
                alignment(E,2,1) = alignment(E,2,2)-relative_size
                
                alignment(W,2,2) = ny-bc_size
                alignment(W,2,1) = alignment(E,2,2)-relative_size

            !on the other corner + 1
            case(5)
                alignment(N,1,2) = nx-bc_size+1
                alignment(N,1,1) = alignment(N,1,2) - relative_size
               
                alignment(S,1,2) = nx-bc_size+1
                alignment(S,1,1) = alignment(S,1,2) - relative_size
                
                alignment(E,2,2) = ny-bc_size+1
                alignment(E,2,1) = alignment(E,2,2)-relative_size
                
                alignment(W,2,2) = ny-bc_size+1
                alignment(W,2,1) = alignment(E,2,2)-relative_size

            !on the other corner + 2
            case(6)
                alignment(N,1,2) = nx-bc_size+2
                alignment(N,1,1) = alignment(N,1,2) - relative_size
               
                alignment(S,1,2) = nx-bc_size+2
                alignment(S,1,1) = alignment(S,1,2) - relative_size
                
                alignment(E,2,2) = ny-bc_size+2
                alignment(E,2,1) = alignment(E,2,2)-relative_size
                
                alignment(W,2,2) = ny-bc_size+2
                alignment(W,2,1) = alignment(E,2,2)-relative_size

            case default
               print '(''test_bf_layer_prog'')'
               print '(''ini_alignment_2nd_bf_layer'')'
               print '(''second_bf_layer_alignment not recognized'')'

          end select

          final_alignment(N,1,1) = bc_size+1+
     $                             test_first_bf_layer_alignment -
     $                             over_alignment_x
          final_alignment(N,1,2) = alignment(N,1,2) + over_alignment_x
          final_alignment(N,2,1) = align_N
          final_alignment(N,2,2) = align_N + over_alignment_y

          final_alignment(S,1,1) = bc_size+1+
     $                             test_first_bf_layer_alignment -
     $                             over_alignment_x
          final_alignment(S,1,2) = alignment(S,1,2) + over_alignment_x
          final_alignment(S,2,1) = align_S - over_alignment_y
          final_alignment(S,2,2) = align_S
          
          final_alignment(E,1,1) = align_E
          final_alignment(E,1,2) = align_E + over_alignment_y
          final_alignment(E,2,1) = bc_size+1+
     $                             test_first_bf_layer_alignment -
     $                             over_alignment_x
          final_alignment(E,2,2) = alignment(E,2,2) + over_alignment_x
          
          final_alignment(W,1,1) = align_W - over_alignment_y
          final_alignment(W,1,2) = align_W
          final_alignment(W,2,1) = bc_size+1+
     $                             test_first_bf_layer_alignment -
     $                             over_alignment_x
          final_alignment(W,2,2) = alignment(W,2,2) + over_alignment_x

        end subroutine ini_alignment_2nd_bf_layer


        subroutine ini_alignment_for_neighbor_copy(
     $     test_alignment_copy_neighbors)

          implicit none

          integer(ikind), dimension(4,3,2,2), intent(out) :: test_alignment_copy_neighbors
          
          integer(ikind) :: large_layer
          integer(ikind) :: small_layer

          large_layer = 7
          small_layer = 3

          !north
          test_alignment_copy_neighbors(N,1,1,1) = -large_layer
          test_alignment_copy_neighbors(N,1,1,2) = nx+large_layer
          test_alignment_copy_neighbors(N,1,2,1) = align_N
          test_alignment_copy_neighbors(N,1,2,2) = align_N

          test_alignment_copy_neighbors(N,2,1,1) = align_W-small_layer
          test_alignment_copy_neighbors(N,2,1,2) = align_W
          test_alignment_copy_neighbors(N,2,2,1) = ny-bc_size-small_layer
          test_alignment_copy_neighbors(N,2,2,2) = ny-bc_size

          test_alignment_copy_neighbors(N,3,1,1) = align_E
          test_alignment_copy_neighbors(N,3,1,2) = align_E+small_layer
          test_alignment_copy_neighbors(N,3,2,1) = ny-bc_size-small_layer
          test_alignment_copy_neighbors(N,3,2,2) = ny-bc_size

          !south
          test_alignment_copy_neighbors(S,1,1,1) = -large_layer
          test_alignment_copy_neighbors(S,1,1,2) = nx+large_layer
          test_alignment_copy_neighbors(S,1,2,1) = align_S
          test_alignment_copy_neighbors(S,1,2,2) = align_S

          test_alignment_copy_neighbors(S,2,1,1) = align_W-small_layer
          test_alignment_copy_neighbors(S,2,1,2) = align_W
          test_alignment_copy_neighbors(S,2,2,1) = 1+bc_size
          test_alignment_copy_neighbors(S,2,2,2) = 1+bc_size+small_layer

          test_alignment_copy_neighbors(S,3,1,1) = align_E
          test_alignment_copy_neighbors(S,3,1,2) = align_E+small_layer
          test_alignment_copy_neighbors(S,3,2,1) = 1+bc_size
          test_alignment_copy_neighbors(S,3,2,2) = 1+bc_size+small_layer
          
          !east
          test_alignment_copy_neighbors(E,1,1,1) = align_E
          test_alignment_copy_neighbors(E,1,1,2) = align_E+large_layer
          test_alignment_copy_neighbors(E,1,2,1) = bc_size+1
          test_alignment_copy_neighbors(E,1,2,2) = ny-bc_size

          test_alignment_copy_neighbors(E,2,1,1) = nx-small_layer
          test_alignment_copy_neighbors(E,2,1,2) = nx+small_layer
          test_alignment_copy_neighbors(E,2,2,1) = align_S
          test_alignment_copy_neighbors(E,2,2,2) = align_S

          test_alignment_copy_neighbors(E,3,1,1) = nx-small_layer
          test_alignment_copy_neighbors(E,3,1,2) = nx+small_layer
          test_alignment_copy_neighbors(E,3,2,1) = align_N
          test_alignment_copy_neighbors(E,3,2,2) = align_N

          !west
          test_alignment_copy_neighbors(W,1,1,1) = align_W-large_layer
          test_alignment_copy_neighbors(W,1,1,2) = align_W
          test_alignment_copy_neighbors(W,1,2,1) = bc_size+1
          test_alignment_copy_neighbors(W,1,2,2) = ny-bc_size

          test_alignment_copy_neighbors(W,2,1,1) = -small_layer
          test_alignment_copy_neighbors(W,2,1,2) = small_layer
          test_alignment_copy_neighbors(W,2,2,1) = align_S
          test_alignment_copy_neighbors(W,2,2,2) = align_S

          test_alignment_copy_neighbors(W,3,1,1) = -small_layer
          test_alignment_copy_neighbors(W,3,1,2) = small_layer
          test_alignment_copy_neighbors(W,3,2,1) = align_N
          test_alignment_copy_neighbors(W,3,2,2) = align_N

        end subroutine ini_alignment_for_neighbor_copy

        
        subroutine update_alignment_to_increase_size(
     $     mainlayer_id,
     $     relative_size,
     $     alignment)

          implicit none

          integer                       , intent(in)   :: mainlayer_id
          integer                       , intent(in)   :: relative_size
          integer(ikind), dimension(2,2), intent(inout):: alignment

          select case(mainlayer_id)
            case(N)
               alignment(2,2) = alignment(2,2) + relative_size
            case(S)
               alignment(2,1) = alignment(2,1) - relative_size
            case(E)
               alignment(1,2) = alignment(1,2) + relative_size
            case(W)
               alignment(1,1) = alignment(1,1) - relative_size
            case default
               print '(''test_bf_layer_prog'')'
               print '(''update_alignemnt_to_increase_size'')'
               print '(''mainlayer_id not recognized'')'
               print '(''mainlayer : '', I2)', mainlayer_id
          end select

        end subroutine update_alignment_to_increase_size


        subroutine ini_selected_grdpts(selected_grdpts)

          implicit none

          integer, dimension(:,:,:), intent(inout) :: selected_grdpts

          selected_grdpts(N,1,1) = 4
          selected_grdpts(N,1,2) = 3
          selected_grdpts(N,2,1) = 4
          selected_grdpts(N,2,2) = 4

          selected_grdpts(S,1,1) = 4
          selected_grdpts(S,1,2) = 3
          selected_grdpts(S,2,1) = 4
          selected_grdpts(S,2,2) = 2

          selected_grdpts(E,1,1) = 3
          selected_grdpts(E,1,2) = 3
          selected_grdpts(E,2,1) = 4
          selected_grdpts(E,2,2) = 3

          selected_grdpts(W,1,1) = 3
          selected_grdpts(W,1,2) = 3
          selected_grdpts(W,2,1) = 2
          selected_grdpts(W,2,2) = 4

        end subroutine ini_selected_grdpts

        subroutine ini_neighbors(
     $     neighbor_case, neighbors)

          implicit none

          integer              , intent(in)  :: neighbor_case
          logical, dimension(4), intent(out) :: neighbors

          select case(neighbor_case)
            case(1)
               neighbors = [.false.,.false.,.false.,.false.]
            case(2)
               neighbors = [.true.,.false.,.true.,.false.]
            case(3)
               neighbors = [.false.,.true.,.false.,.true.]
            case(4)
               neighbors = [.true.,.true.,.true.,.true.]
          end select

        end subroutine ini_neighbors

        subroutine ini_using_case(
     $     distance_case, seed, relative_distance)

          implicit none

          integer, intent(in) :: distance_case
          integer, intent(in) :: seed
          integer, intent(out):: relative_distance          

          select case(distance_case)
            case(1)
               relative_distance=0
            case(2)
               relative_distance=1
            case(3)
               call srand(seed)
               relative_distance = nint(5.0*RAND())
            case default
               print '(''test_bf_layer_prog'')'
               print '(''ini_using_case'')'
               print '(''distance_case not recognized'')'
               print '(''distance_case'', I2)', distance_case
               stop 'change distance_case'               
          end select               

        end subroutine ini_using_case

      end program test_bf_layer_prog
