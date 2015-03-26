      program test_bf_layer_exchange

        use bf_layer_exchange_module, only :
     $     do_grdpts_overlap_along_x_dir,
     $     get_match_indices_for_exchange_with_neighbor1,
     $     get_match_indices_for_exchange_with_neighbor2,
     $     copy_from_bf1_to_bf2,
     $     get_sync_indices_with_interior,
     $     sync_nodes,
     $     get_sync_indices_with_neighbor1,
     $     get_sync_indices_with_neighbor2

        use check_data_module, only :
     $       is_int_vector_validated,
     $       is_int_matrix_validated,
     $       is_real_vector_validated,
     $       is_real_matrix3D_validated

        use parameters_bf_layer, only :
     $       align_N,
     $       align_S,
     $       align_E,
     $       align_W

        use parameters_constant, only :
     $       N,S,E,W,
     $       x_direction,y_direction

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        interface
           subroutine match_interface(
     $          bf_alignment,
     $          nbf_alignment,
     $          bf_i_min,
     $          bf_j_min,
     $          nbf_i_min,
     $          nbf_j_min,
     $          bf_copy_size_x,
     $          bf_copy_size_y)

           import ikind

           integer(ikind), dimension(2,2), intent(in)  :: bf_alignment
           integer(ikind), dimension(2,2), intent(in)  :: nbf_alignment
           integer(ikind)                , intent(out) :: bf_i_min
           integer(ikind)                , intent(out) :: bf_j_min
           integer(ikind)                , intent(out) :: nbf_i_min
           integer(ikind)                , intent(out) :: nbf_j_min
           integer(ikind)                , intent(out) :: bf_copy_size_x
           integer(ikind)                , intent(out) :: bf_copy_size_y
           
          end subroutine match_interface

        end interface

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled      = .false.
        test_validated = .true.


        test_loc = test_do_grdpts_overlap_along_x_dir(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_do_grdpts_overlap_along_x_dir: '',L1)', test_loc
        print '()'

        test_loc = test_get_match_indices_for_exchange_with_neighbor1(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_match_indices_for_exchange_with_neighbor1: '',L1)', test_loc
        print '()'

        test_loc = test_get_match_indices_for_exchange_with_neighbor2(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_match_indices_for_exchange_with_neighbor2: '',L1)', test_loc
        print '()'

        test_loc = test_copy_from_bf1_to_bf2(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_copy_from_bf1_to_bf2: '',L1)', test_loc
        print '()'

        test_loc = test_get_sync_indices_with_interior(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_sync_indices_with_interior: '',L1)', test_loc
        print '()'

        test_loc = test_sync_nodes(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_sync_nodes: '',L1)', test_loc
        print '()'

        test_loc = test_get_sync_indices_with_neighbor1(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_sync_indices_with_neighbor1: '',L1)', test_loc
        print '()'

        detailled = .true.

        test_loc = test_get_sync_indices_with_neighbor2(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_sync_indices_with_neighbor2: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_do_grdpts_overlap_along_x_dir(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer(ikind), dimension(2,2,4) :: bf_alignment
          integer(ikind), dimension(2,2,4) :: nbf_alignment
          logical       , dimension(4)     :: overlap_test
          logical                          :: overlap
          logical                          :: test_loc

          integer :: k

          !input
          bf_alignment(1,1,1) = 2
          bf_alignment(1,2,1) = 8

          nbf_alignment(1,1,1) = 1
          nbf_alignment(1,2,1) = 6

          overlap_test(1) = .true.


          bf_alignment(1,1,2) = 2
          bf_alignment(1,2,2) = 8

          nbf_alignment(1,1,2) = -3
          nbf_alignment(1,2,2) = -1

          overlap_test(2) = .true.
          

          bf_alignment(1,1,3) = 2
          bf_alignment(1,2,3) = 8

          nbf_alignment(1,1,3) = 12
          nbf_alignment(1,2,3) = 15

          overlap_test(3) = .true.


          bf_alignment(1,1,4) = 2
          bf_alignment(1,2,4) = 8

          nbf_alignment(1,1,4) = 13
          nbf_alignment(1,2,4) = 16

          overlap_test(4) = .false.

          
          test_validated = .true.

          do k=1,4
             
             !output
             overlap = do_grdpts_overlap_along_x_dir(
     $            bf_alignment(:,:,k),
     $            nbf_alignment(:,:,k))

             !validation
             test_loc = (overlap.eqv.overlap_test(k))
             test_validated = test_validated.and.test_loc

             !detailled
             if((.not.test_loc).and.detailled) then

                print '(''overlap('',I1,''): '',L1,'' -> '',L1)',
     $               k,test_loc, overlap_test(k)

             end if

          end do

        end function test_do_grdpts_overlap_along_x_dir

      
        function test_get_match_indices_for_exchange_with_neighbor(
     $     detailled,
     $     get_match_indices_for_exchange_with_neighbor)
     $     result(test_validated)

          implicit none

          logical, intent(in)        :: detailled
          procedure(match_interface) :: get_match_indices_for_exchange_with_neighbor
          logical                    :: test_validated

          integer(ikind), dimension(2,2,1) :: bf_alignment
          integer(ikind), dimension(2,2,1) :: nbf_alignment
          integer                          :: k
          logical                          :: test_loc
          integer, dimension(6)            :: output
          integer, dimension(6,1)          :: output_test


          bf_alignment(1,1,1) = -4
          bf_alignment(1,2,1) =  5
          bf_alignment(2,1,1) =  align_S-4
          bf_alignment(2,2,1) =  align_S

          nbf_alignment(1,1,1) =  align_W-3
          nbf_alignment(1,2,1) =  align_W
          nbf_alignment(2,1,1) =  align_S+1
          nbf_alignment(2,2,1) =  align_S+5

          output_test(1,1) = 4
          output_test(2,1) = 6
          output_test(3,1) = 1
          output_test(4,1) = 1
          output_test(5,1) = 8
          output_test(6,1) = 4

          test_validated = .true.

          do k=1,1

             !output
             call get_match_indices_for_exchange_with_neighbor(
     $            bf_alignment(:,:,k),
     $            nbf_alignment(:,:,k),
     $            output(1),
     $            output(2),
     $            output(3),
     $            output(4),
     $            output(5),
     $            output(6))

             !validation
             test_loc = is_int_vector_validated(output,output_test(:,k))
             test_validated = test_validated.and.test_loc

             !detailled
             if((.not.test_loc).and.detailled) then

                print '(''test '',I1,'' : '',L1)', k, test_loc
                print '(''  - bf_i_min:       '',I3,'' -> '',I3)', output(1), output_test(1,k)
                print '(''  - bf_j_min:       '',I3,'' -> '',I3)', output(2), output_test(2,k)
                print '(''  - nbf_i_min:      '',I3,'' -> '',I3)', output(3), output_test(3,k)
                print '(''  - nbf_j_min:      '',I3,'' -> '',I3)', output(4), output_test(4,k)
                print '(''  - bf_copy_size_x: '',I3,'' -> '',I3)', output(5), output_test(5,k)
                print '(''  - bf_copy_size_y: '',I3,'' -> '',I3)', output(6), output_test(6,k)
                
             end if

          end do

        end function test_get_match_indices_for_exchange_with_neighbor


        function test_get_match_indices_for_exchange_with_neighbor1(
     $     detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = test_get_match_indices_for_exchange_with_neighbor(
     $         detailled,
     $         get_match_indices_for_exchange_with_neighbor1)

        end function test_get_match_indices_for_exchange_with_neighbor1


        function test_get_match_indices_for_exchange_with_neighbor2(
     $     detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = test_get_match_indices_for_exchange_with_neighbor(
     $         detailled,
     $         get_match_indices_for_exchange_with_neighbor2)

        end function test_get_match_indices_for_exchange_with_neighbor2


        function test_copy_from_bf1_to_bf2(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer :: bf1_i_min
          integer :: bf1_j_min
          integer :: bf2_i_min
          integer :: bf2_j_min
          integer :: bf_copy_size_x
          integer :: bf_copy_size_y

          real(rkind), dimension(6,8,ne) :: bf1_nodes
          real(rkind), dimension(7,3,ne) :: bf2_nodes
          integer    , dimension(6,8)    :: bf1_grdpts_id
          integer    , dimension(7,3)    :: bf2_grdpts_id
          real(rkind), dimension(6,8,ne) :: bf1_nodes_test
          real(rkind), dimension(7,3,ne) :: bf2_nodes_test
          integer    , dimension(6,8)    :: bf1_grdpts_id_test
          integer    , dimension(7,3)    :: bf2_grdpts_id_test

          integer(ikind) :: i,j
          integer        :: k


          !input
          bf1_i_min = 3
          bf1_j_min = 4
          bf2_i_min = 5
          bf2_j_min = 2
          bf_copy_size_x = 3
          bf_copy_size_y = 2

          do k=1,ne
             do j=1,size(bf1_nodes,2)
                do i=1,size(bf1_nodes,1)
                   bf1_nodes(i,j,k) = (i-1) + 10*(j-1) + 100*k
                end do
             end do
          end do

          do j=1,size(bf1_nodes,2)
             do i=1,size(bf1_nodes,1)
                bf1_grdpts_id(i,j) = (i-1) + 10*(j-1)
             end do
          end do

          do k=1,ne
             do j=1,size(bf2_nodes,2)
                do i=1,size(bf2_nodes,1)
                   bf2_nodes(i,j,k) = 5*(i-1) + 50*(j-1) + 500*k
                end do
             end do
          end do

          do j=1,size(bf2_nodes,2)
             do i=1,size(bf2_nodes,1)
                bf2_grdpts_id(i,j) = 5*(i-1) + 50*(j-1)
             end do
          end do

          bf1_nodes_test     = bf1_nodes
          bf2_nodes_test     = bf2_nodes
          bf1_grdpts_id_test = bf1_grdpts_id
          bf2_grdpts_id_test = bf2_grdpts_id

          bf2_nodes_test(
     $         bf2_i_min:bf2_i_min+bf_copy_size_x-1,
     $         bf2_j_min:bf2_j_min+bf_copy_size_y-1,
     $         :) =
     $         bf1_nodes_test(
     $         bf1_i_min:bf1_i_min+bf_copy_size_x-1,
     $         bf1_j_min:bf1_j_min+bf_copy_size_y-1,
     $         :)
          
          bf2_grdpts_id_test(
     $         bf2_i_min:bf2_i_min+bf_copy_size_x-1,
     $         bf2_j_min:bf2_j_min+bf_copy_size_y-1) =
     $         bf1_grdpts_id_test(
     $         bf1_i_min:bf1_i_min+bf_copy_size_x-1,
     $         bf1_j_min:bf1_j_min+bf_copy_size_y-1)

          !output
          call copy_from_bf1_to_bf2(
     $         bf1_i_min, bf1_j_min, bf2_i_min, bf2_j_min,
     $         bf_copy_size_x, bf_copy_size_y,
     $         bf1_nodes, bf1_grdpts_id,
     $         bf2_nodes, bf2_grdpts_id)

          !validation
          test_validated = .true.

          test_loc = is_real_matrix3D_validated(bf1_nodes,bf1_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_real_matrix3D_validated(bf2_nodes,bf2_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_int_matrix_validated(bf1_grdpts_id,bf1_grdpts_id_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_int_matrix_validated(bf2_grdpts_id,bf2_grdpts_id_test,detailled)
          test_validated = test_validated.and.test_loc

        end function test_copy_from_bf1_to_bf2


        function test_get_sync_indices_with_interior(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          integer       , dimension(4)     :: localization
          integer(ikind), dimension(2,2,4) :: bf_alignment
          integer(ikind), dimension(2,4)   :: bf_size
          integer(ikind), dimension(10)    :: output
          integer(ikind), dimension(10,4)  :: output_test

          logical :: test_loc
          integer :: k


          !input
          localization(1)     = N
          bf_alignment(1,1,1) = 4
          bf_alignment(1,2,1) = 5
          bf_alignment(2,1,1) = align_N
          bf_alignment(2,2,1) = align_N+1
          bf_size(:,1)        = [6,6]

          output_test(1:2,1)  = [2,ny-2*bc_size+1]
          output_test(3:4,1)  = [2,ny-bc_size+1]
          output_test(5:6,1)  = [1,3]
          output_test(7:8,1)  = [1,1]
          output_test(9:10,1) = [6,bc_size]

          localization(2)     = S
          bf_alignment(1,1,2) = 4
          bf_alignment(1,2,2) = 5
          bf_alignment(2,1,2) = align_S-1
          bf_alignment(2,2,2) = align_S
          bf_size(:,2)        = [6,6]

          output_test(1:2,2)  = [2,3]
          output_test(3:4,2)  = [2,1]
          output_test(5:6,2)  = [1,3]
          output_test(7:8,2)  = [1,5]
          output_test(9:10,2) = [6,bc_size]

          localization(3)     = E
          bf_alignment(1,1,3) = align_E
          bf_alignment(1,2,3) = align_E+1
          bf_alignment(2,1,3) = 4
          bf_alignment(2,2,3) = 5
          bf_size(:,3)        = [6,6]

          output_test(1:2,3)  = [nx-2*bc_size+1,2]
          output_test(3:4,3)  = [nx-bc_size+1  ,2]
          output_test(5:6,3)  = [3,1]
          output_test(7:8,3)  = [1,1]
          output_test(9:10,3) = [bc_size,6]

          localization(4)     = W
          bf_alignment(1,1,4) = align_W-1
          bf_alignment(1,2,4) = align_W
          bf_alignment(2,1,4) = 4
          bf_alignment(2,2,4) = 5
          bf_size(:,4)        = [6,6]

          output_test(1:2,4)  = [bc_size+1,2]
          output_test(3:4,4)  = [1,2]
          output_test(5:6,4)  = [3,1]
          output_test(7:8,4)  = [5,1]
          output_test(9:10,4) = [bc_size,6]


          test_validated = .true.

          do k=1,4

             !output
             call get_sync_indices_with_interior(
     $            localization(k),
     $            bf_alignment(:,:,k),
     $            bf_size(:,k),
     $            output(1:2),
     $            output(3:4),
     $            output(5:6),
     $            output(7:8),
     $            output(9:10))

             !validation
             test_loc = is_int_vector_validated(output,output_test(:,k),detailled)
             test_validated = test_validated.and.test_loc
             
             !detailled
             if((.not.test_loc).and.detailled) then

                print '(''test '',I2,'': '',L1)', k, test_loc
                print '(''   - in_send: '',2I3,'' -> '',2I3)', output(1:2) , output_test(1:2,k)
                print '(''   - in_recv: '',2I3,'' -> '',2I3)', output(3:4) , output_test(3:4,k)
                print '(''   - bf_send: '',2I3,'' -> '',2I3)', output(5:6) , output_test(5:6,k)
                print '(''   - bf_recv: '',2I3,'' -> '',2I3)', output(7:8) , output_test(7:8,k)
                print '(''   - ex_size: '',2I3,'' -> '',2I3)', output(9:10), output_test(9:10,k)
                print '()'

             end if

          end do

        end function test_get_sync_indices_with_interior


        function test_sync_nodes(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer(ikind), dimension(2) :: bf1_send
          integer(ikind), dimension(2) :: bf1_recv
          integer(ikind), dimension(2) :: bf2_send
          integer(ikind), dimension(2) :: bf2_recv
          integer(ikind), dimension(2) :: ex_size

          real(rkind), dimension(10,9,ne) :: bf1_nodes
          real(rkind), dimension(7,6,ne)  :: bf2_nodes
          real(rkind), dimension(10,9,ne) :: bf1_nodes_test
          real(rkind), dimension(7,6,ne)  :: bf2_nodes_test

          integer(ikind) :: i,j
          integer        :: k


          !input
          bf1_send = [3,4]
          bf1_recv = [7,7]

          bf2_send = [5,2]
          bf2_recv = [1,5]

          ex_size  = [3,2]

          do k=1,ne
             do j=1,size(bf1_nodes,2)
                do i=1,size(bf1_nodes,1)
                   bf1_nodes(i,j,k) = (i-1) + 10*(j-1) + 100*k
                end do
             end do
          end do

          do k=1,ne
             do j=1,size(bf2_nodes,2)
                do i=1,size(bf2_nodes,1)
                   bf2_nodes(i,j,k) = 5*(i-1) + 50*(j-1) + 500*k
                end do
             end do
          end do

          bf1_nodes_test = bf1_nodes
          bf2_nodes_test = bf2_nodes

          bf2_nodes_test(
     $         bf2_recv(1):bf2_recv(1)+ex_size(1)-1,
     $         bf2_recv(2):bf2_recv(2)+ex_size(2)-1,
     $         :) =
     $         bf1_nodes_test(
     $         bf1_send(1):bf1_send(1)+ex_size(1)-1,
     $         bf1_send(2):bf1_send(2)+ex_size(2)-1,
     $         :)

          bf1_nodes_test(
     $         bf1_recv(1):bf1_recv(1)+ex_size(1)-1,
     $         bf1_recv(2):bf1_recv(2)+ex_size(2)-1,
     $         :) =
     $         bf2_nodes_test(
     $         bf2_send(1):bf2_send(1)+ex_size(1)-1,
     $         bf2_send(2):bf2_send(2)+ex_size(2)-1,
     $         :)
          

          !output
          call sync_nodes(
     $         bf1_nodes,
     $         bf1_send,
     $         bf1_recv,
     $         bf2_nodes,
     $         bf2_send,
     $         bf2_recv,
     $         ex_size)


          !validation
          test_validated = .true.

          test_loc = is_real_matrix3D_validated(bf1_nodes,bf1_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_real_matrix3D_validated(bf2_nodes,bf2_nodes_test,detailled)
          test_validated = test_validated.and.test_loc


        end function test_sync_nodes


        function test_get_sync_indices_with_neighbor1(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer       , dimension(4)     :: bf1_localization
          integer(ikind), dimension(2,2,4) :: bf1_alignment
          integer(ikind), dimension(2,4)   :: bf1_size

          integer(ikind), dimension(2,2,4) :: bf2_alignment
          integer(ikind), dimension(2,4)   :: bf2_size

          integer(ikind), dimension(10,4)  :: output_test
          integer(ikind), dimension(10)    :: output

          integer :: k
          logical :: test_loc
          

          test_validated = .true.


          !input
          !North
          bf1_localization(1) = N

          bf1_alignment(1,1,1) = align_W-10
          bf1_alignment(1,2,1) = align_W+4
          bf1_alignment(2,1,1) = align_N
          bf1_alignment(2,2,1) = align_N+1

          bf1_size(:,1) = [19,6]

          output_test(1:2,1)  = [6,3]
          output_test(3:4,1)  = [6,1]

          bf2_alignment(1,1,1) = align_W-5
          bf2_alignment(1,2,1) = align_W
          bf2_alignment(2,1,1) = align_N-4
          bf2_alignment(2,2,1) = align_N-1
          
          bf2_size(:,1) = [10,8]

          output_test(5:6,1)  = [1,5]
          output_test(7:8,1)  = [1,7]
          output_test(9:10,1) = [10,2]


          !South
          bf1_localization(2) = S

          bf1_alignment(1,1,2) = align_W-10
          bf1_alignment(1,2,2) = align_W+4
          bf1_alignment(2,1,2) = align_S-1
          bf1_alignment(2,2,2) = align_S

          bf1_size(:,2) = [19,6]

          output_test(1:2,2)  = [6,3]
          output_test(3:4,2)  = [6,5]

          bf2_alignment(1,1,2) = align_W-5
          bf2_alignment(1,2,2) = align_W
          bf2_alignment(2,1,2) = align_S+1
          bf2_alignment(2,2,2) = align_S+5
          
          bf2_size(:,2) = [10,8]

          output_test(5:6,2)  = [1,3]
          output_test(7:8,2)  = [1,1]
          output_test(9:10,2) = [10,2]


          !East
          bf1_localization(3) = E

          bf1_alignment(1,1,3) = align_E
          bf1_alignment(1,2,3) = align_E+1
          bf1_alignment(2,1,3) = align_S+1
          bf1_alignment(2,2,3) = align_S+15

          bf1_size(:,3) = [6,19]

          output_test(1:2,3)  = [1,3]
          output_test(3:4,3)  = [1,1]

          bf2_alignment(1,1,3) = align_E-5
          bf2_alignment(1,2,3) = align_E+10
          bf2_alignment(2,1,3) = align_S-4
          bf2_alignment(2,2,3) = align_S
          
          bf2_size(:,3) = [20,9]

          output_test(5:6,3)  = [6,6]
          output_test(7:8,3)  = [6,8]
          output_test(9:10,3) = [6,2]


          !West
          bf1_localization(4) = W

          bf1_alignment(1,1,4) = align_W-1
          bf1_alignment(1,2,4) = align_W
          bf1_alignment(2,1,4) = align_S+1
          bf1_alignment(2,2,4) = align_S+15

          bf1_size(:,4) = [6,19]

          output_test(1:2,4)  = [1,3]
          output_test(3:4,4)  = [1,1]

          bf2_alignment(1,1,4) = align_W-10
          bf2_alignment(1,2,4) = align_W+4
          bf2_alignment(2,1,4) = align_S-4
          bf2_alignment(2,2,4) = align_S
          
          bf2_size(:,4) = [19,9]

          output_test(5:6,4)  = [10,6]
          output_test(7:8,4)  = [10,8]
          output_test(9:10,4) = [6,2]


          do k=1,4

             !output
             call get_sync_indices_with_neighbor1(
     $            bf1_localization(k),
     $            bf1_alignment(:,:,k),
     $            bf1_size(:,k),
     $            output(1:2),
     $            output(3:4),
     $            bf2_alignment(:,:,k),
     $            bf2_size(:,k),
     $            output(5:6),
     $            output(7:8),
     $            output(9:10))

             !validation
             test_loc = is_int_vector_validated(output,output_test(:,k),detailled)
             test_validated = test_validated.and.test_loc

             !detailled
             if((.not.test_loc).and.detailled) then

                print '(''test '',I2,'': '',L1)', k,test_loc
                print '(''   - bf1_send: '',2I3,'' -> '',2I3)', output(1:2) , output_test(1:2,k)
                print '(''   - bf1_recv: '',2I3,'' -> '',2I3)', output(3:4) , output_test(3:4,k)
                print '(''   - bf2_send: '',2I3,'' -> '',2I3)', output(5:6) , output_test(5:6,k)
                print '(''   - bf2_recv: '',2I3,'' -> '',2I3)', output(7:8) , output_test(7:8,k)
                print '(''   - ex_size : '',2I3,'' -> '',2I3)', output(9:10), output_test(9:10,k)
                print '()'

             end if

          end do

        end function test_get_sync_indices_with_neighbor1


        function test_get_sync_indices_with_neighbor2(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer       , dimension(4)     :: bf1_localization
          integer(ikind), dimension(2,2,4) :: bf1_alignment
          integer(ikind), dimension(2,4)   :: bf1_size

          integer(ikind), dimension(2,2,4) :: bf2_alignment
          integer(ikind), dimension(2,4)   :: bf2_size

          integer(ikind), dimension(10,4)  :: output_test
          integer(ikind), dimension(10)    :: output

          integer :: k
          logical :: test_loc
          

          test_validated = .true.


          !input
          !North
          bf1_localization(1) = N

          bf1_alignment(1,1,1) = align_E-10
          bf1_alignment(1,2,1) = align_E+4
          bf1_alignment(2,1,1) = align_N
          bf1_alignment(2,2,1) = align_N+1

          bf1_size(:,1) = [19,6]

          output_test(1:2,1)  = [11,3]
          output_test(3:4,1)  = [11,1]

          bf2_alignment(1,1,1) = align_E
          bf2_alignment(1,2,1) = align_E+5
          bf2_alignment(2,1,1) = align_N-4
          bf2_alignment(2,2,1) = align_N-1
          
          bf2_size(:,1) = [10,8]

          output_test(5:6,1)  = [1,5]
          output_test(7:8,1)  = [1,7]
          output_test(9:10,1) = [9,2]


          !South
          bf1_localization(2) = S

          bf1_alignment(1,1,2) = align_E-10
          bf1_alignment(1,2,2) = align_E+4
          bf1_alignment(2,1,2) = align_S-1
          bf1_alignment(2,2,2) = align_S

          bf1_size(:,2) = [19,6]

          output_test(1:2,2)  = [11,3]
          output_test(3:4,2)  = [11,5]

          bf2_alignment(1,1,2) = align_E
          bf2_alignment(1,2,2) = align_E+5
          bf2_alignment(2,1,2) = align_S+1
          bf2_alignment(2,2,2) = align_S+5
          
          bf2_size(:,2) = [10,8]

          output_test(5:6,2)  = [1,3]
          output_test(7:8,2)  = [1,1]
          output_test(9:10,2) = [9,2]


          !East
          bf1_localization(3) = E

          bf1_alignment(1,1,3) = align_E
          bf1_alignment(1,2,3) = align_E+1
          bf1_alignment(2,1,3) = align_N-15
          bf1_alignment(2,2,3) = align_N-1

          bf1_size(:,3) = [6,19]

          output_test(1:2,3)  = [1,16]
          output_test(3:4,3)  = [1,18]

          bf2_alignment(1,1,3) = align_E-5
          bf2_alignment(1,2,3) = align_E+10
          bf2_alignment(2,1,3) = align_N
          bf2_alignment(2,2,3) = align_N+4
          
          bf2_size(:,3) = [20,9]

          output_test(5:6,3)  = [6,3]
          output_test(7:8,3)  = [6,1]
          output_test(9:10,3) = [6,2]


          !West
          bf1_localization(4) = W

          bf1_alignment(1,1,4) = align_W-1
          bf1_alignment(1,2,4) = align_W
          bf1_alignment(2,1,4) = align_N-15
          bf1_alignment(2,2,4) = align_N-1

          bf1_size(:,4) = [6,19]

          output_test(1:2,4)  = [1,16]
          output_test(3:4,4)  = [1,18]

          bf2_alignment(1,1,4) = align_W-10
          bf2_alignment(1,2,4) = align_W+4
          bf2_alignment(2,1,4) = align_N
          bf2_alignment(2,2,4) = align_N+4
          
          bf2_size(:,4) = [19,9]

          output_test(5:6,4)  = [10,3]
          output_test(7:8,4)  = [10,1]
          output_test(9:10,4) = [6,2]


          do k=1,4

             !output
             call get_sync_indices_with_neighbor2(
     $            bf1_localization(k),
     $            bf1_alignment(:,:,k),
     $            bf1_size(:,k),
     $            output(1:2),
     $            output(3:4),
     $            bf2_alignment(:,:,k),
     $            bf2_size(:,k),
     $            output(5:6),
     $            output(7:8),
     $            output(9:10))

             !validation
             test_loc = is_int_vector_validated(output,output_test(:,k),detailled)
             test_validated = test_validated.and.test_loc

             !detailled
             if((.not.test_loc).and.detailled) then

                print '(''test '',I2,'': '',L1)', k,test_loc
                print '(''   - bf1_send: '',2I3,'' -> '',2I3)', output(1:2) , output_test(1:2,k)
                print '(''   - bf1_recv: '',2I3,'' -> '',2I3)', output(3:4) , output_test(3:4,k)
                print '(''   - bf2_send: '',2I3,'' -> '',2I3)', output(5:6) , output_test(5:6,k)
                print '(''   - bf2_recv: '',2I3,'' -> '',2I3)', output(7:8) , output_test(7:8,k)
                print '(''   - ex_size : '',2I3,'' -> '',2I3)', output(9:10), output_test(9:10,k)
                print '()'

             end if

          end do

        end function test_get_sync_indices_with_neighbor2


      end program test_bf_layer_exchange
