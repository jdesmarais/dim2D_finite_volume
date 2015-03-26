      program test_bf_layer_sync

        use bf_layer_sync_class, only :
     $       bf_layer_sync

        use check_data_module, only :
     $       is_int_matrix_validated,
     $       is_real_matrix3D_validated

        use parameters_bf_layer, only :
     $       N,S,E,W,
     $       align_N,
     $       align_S,
     $       align_E,
     $       align_W

        use parameters_input, only :
     $       bc_size,
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        integer, parameter :: test_id_copy_to_neighbor1   = 1
        integer, parameter :: test_id_copy_from_neighbor1 = 2
        integer, parameter :: test_id_copy_to_neighbor2   = 3
        integer, parameter :: test_id_copy_from_neighbor2 = 4
        integer, parameter :: test_id_sync_with_neighbor  = 5

        
        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.

        test_loc = test_set_neighbor1_share(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_set_neighbor1_share: '',L1)', test_loc
        print '()'


        test_loc = test_set_neighbor2_share(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_set_neighbor2_share: '',L1)', test_loc
        print '()'


        test_loc = test_get_neighbor1_id(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_neighbor1_id: '',L1)', test_loc
        print '()'


        test_loc = test_get_neighbor2_id(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_neighbor2_id: '',L1)', test_loc
        print '()'


        test_loc = test_copy_from_neighbor1(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_copy_from_neighbor1: '',L1)', test_loc
        print '()'


        test_loc = test_copy_to_neighbor1(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_copy_to_neighbor1: '',L1)', test_loc
        print '()'


        test_loc = test_copy_from_neighbor2(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_copy_from_neighbor2: '',L1)', test_loc
        print '()'


        test_loc = test_copy_to_neighbor2(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_copy_to_neighbor2: '',L1)', test_loc
        print '()'


        test_loc = test_sync_nodes_with_interior(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_sync_nodes_with_interior: '',L1)', test_loc
        print '()'


        test_loc = test_sync_nodes_with_neighbor1(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_sync_nodes_with_neighbor1: '',L1)', test_loc
        print '()'

        
        test_loc = test_sync_nodes_with_neighbor2(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_sync_nodes_with_neighbor2: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_set_neighbor1_share(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_sync)              :: bf_layer_used
          integer       , dimension(8)     :: bf_localization
          integer(ikind), dimension(2,2,8) :: bf_alignment
          logical       , dimension(8)     :: neighbor_share_test
          logical                          :: neighbor_share

          integer :: k
          logical :: test_loc

          test_validated = .true.


          !input
          !N/S
          bf_localization(1)     = N
          bf_alignment(1,1,1)    = align_W-3
          bf_alignment(1,2,1)    = align_W+2
          bf_alignment(2,1,1)    = align_N
          bf_alignment(2,2,1)    = align_N+1
          neighbor_share_test(1) = .true.

          bf_localization(2)     = N
          bf_alignment(1,1,2)    = align_W+2*bc_size
          bf_alignment(1,2,2)    = align_W+2*bc_size+1
          bf_alignment(2,1,2)    = align_N
          bf_alignment(2,2,2)    = align_N+1
          neighbor_share_test(2) = .true.

          bf_localization(3)     = N
          bf_alignment(1,1,3)    = align_W+2*bc_size+1
          bf_alignment(1,2,3)    = align_W+2*bc_size+2
          bf_alignment(2,1,3)    = align_N
          bf_alignment(2,2,3)    = align_N+1
          neighbor_share_test(3) = .false.

          bf_localization(4)     = N
          bf_alignment(1,1,4)    = align_W+2*bc_size+5
          bf_alignment(1,2,4)    = align_W+2*bc_size+6
          bf_alignment(2,1,4)    = align_N
          bf_alignment(2,2,4)    = align_N+1
          neighbor_share_test(4) = .false.

          !E/W
          bf_localization(5)     = E
          bf_alignment(1,1,5)    = align_E
          bf_alignment(1,2,5)    = align_E+2
          bf_alignment(2,1,5)    = align_S+bc_size+1
          bf_alignment(2,2,5)    = align_S+bc_size+2
          neighbor_share_test(5) = .true.

          bf_localization(6)     = E
          bf_alignment(1,1,6)    = align_E
          bf_alignment(1,2,6)    = align_E+2
          bf_alignment(2,1,6)    = align_S+bc_size+2
          bf_alignment(2,2,6)    = align_S+bc_size+3
          neighbor_share_test(6) = .true.

          bf_localization(7)     = E
          bf_alignment(1,1,7)    = align_E
          bf_alignment(1,2,7)    = align_E+2
          bf_alignment(2,1,7)    = align_S+bc_size+3
          bf_alignment(2,2,7)    = align_S+bc_size+4
          neighbor_share_test(7) = .false.

          bf_localization(8)     = E
          bf_alignment(1,1,8)    = align_E
          bf_alignment(1,2,8)    = align_E+2
          bf_alignment(2,1,8)    = align_S+bc_size+4
          bf_alignment(2,2,8)    = align_S+bc_size+5
          neighbor_share_test(8) = .false.

          do k=1,size(bf_localization,1)

             !output
             call bf_layer_used%ini(bf_localization(k))
             call bf_layer_used%set_alignment_tab(bf_alignment(:,:,k))

             call bf_layer_used%set_neighbor1_share()
             
             neighbor_share = bf_layer_used%can_exchange_with_neighbor1()


             !validation
             test_loc = neighbor_share.eqv.neighbor_share_test(k)
             test_validated = test_validated.and.test_loc


             !detailled
             if((.not.test_loc).and.detailled) then
                print '(''test '',I2, '' : '',L1,'' -> '',L1)', 
     $               k,
     $               neighbor_share,
     $               neighbor_share_test(k)
             end if

          end do

        end function test_set_neighbor1_share


        function test_set_neighbor2_share(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_sync)              :: bf_layer_used
          integer       , dimension(8)     :: bf_localization
          integer(ikind), dimension(2,2,8) :: bf_alignment
          logical       , dimension(8)     :: neighbor_share_test
          logical                          :: neighbor_share

          integer :: k
          logical :: test_loc

          test_validated = .true.


          !input
          !N/S
          bf_localization(1)     = N
          bf_alignment(1,1,1)    = align_E-2
          bf_alignment(1,2,1)    = align_E+3
          bf_alignment(2,1,1)    = align_N
          bf_alignment(2,2,1)    = align_N+1
          neighbor_share_test(1) = .true.

          bf_localization(2)     = N
          bf_alignment(1,1,2)    = align_E-(2*bc_size+1)
          bf_alignment(1,2,2)    = align_E-2*bc_size
          bf_alignment(2,1,2)    = align_N
          bf_alignment(2,2,2)    = align_N+1
          neighbor_share_test(2) = .true.

          bf_localization(3)     = N
          bf_alignment(1,1,3)    = align_E-(2*bc_size+2)
          bf_alignment(1,2,3)    = align_E-(2*bc_size+1)
          bf_alignment(2,1,3)    = align_N
          bf_alignment(2,2,3)    = align_N+1
          neighbor_share_test(3) = .false.

          bf_localization(4)     = N
          bf_alignment(1,1,4)    = align_E-(2*bc_size+3)
          bf_alignment(1,2,4)    = align_E-(2*bc_size+2)
          bf_alignment(2,1,4)    = align_N
          bf_alignment(2,2,4)    = align_N+1
          neighbor_share_test(4) = .false.

          !E/W
          bf_localization(5)     = E
          bf_alignment(1,1,5)    = align_E
          bf_alignment(1,2,5)    = align_E+2
          bf_alignment(2,1,5)    = align_N-(bc_size+2)
          bf_alignment(2,2,5)    = align_N-(bc_size+1)
          neighbor_share_test(5) = .true.

          bf_localization(6)     = E
          bf_alignment(1,1,6)    = align_E
          bf_alignment(1,2,6)    = align_E+2
          bf_alignment(2,1,6)    = align_N-(bc_size+3)
          bf_alignment(2,2,6)    = align_N-(bc_size+2)
          neighbor_share_test(6) = .true.

          bf_localization(7)     = E
          bf_alignment(1,1,7)    = align_E
          bf_alignment(1,2,7)    = align_E+2
          bf_alignment(2,1,7)    = align_N-(bc_size+4)
          bf_alignment(2,2,7)    = align_N-(bc_size+3)
          neighbor_share_test(7) = .false.

          bf_localization(8)     = E
          bf_alignment(1,1,8)    = align_E
          bf_alignment(1,2,8)    = align_E+2
          bf_alignment(2,1,8)    = align_N-(bc_size+5)
          bf_alignment(2,2,8)    = align_N-(bc_size+4)
          neighbor_share_test(8) = .false.

          do k=1,size(bf_localization,1)

             !output
             call bf_layer_used%ini(bf_localization(k))
             call bf_layer_used%set_alignment_tab(bf_alignment(:,:,k))

             call bf_layer_used%set_neighbor2_share()
             
             neighbor_share = bf_layer_used%can_exchange_with_neighbor2()


             !validation
             test_loc = neighbor_share.eqv.neighbor_share_test(k)
             test_validated = test_validated.and.test_loc


             !detailled
             if((.not.test_loc).and.detailled) then
                print '(''test '',I2, '' : '',L1,'' -> '',L1)', 
     $               k,
     $               neighbor_share,
     $               neighbor_share_test(k)
             end if

          end do

        end function test_set_neighbor2_share        


        function test_get_neighbor1_id(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_sync)          :: bf_layer_used
          integer       , dimension(4) :: bf_localization
          integer       , dimension(4) :: neighbor1_id_test
          integer       , dimension(4) :: neighbor_index_test
          integer                      :: neighbor1_id
          integer                      :: neighbor_index

          logical :: test_loc
          integer :: k

          test_validated = .true.


          bf_localization(1)     = N
          neighbor1_id_test(1)   = W
          neighbor_index_test(1) = 2

          bf_localization(2)     = S
          neighbor1_id_test(2)   = W
          neighbor_index_test(2) = 1

          bf_localization(3)     = E
          neighbor1_id_test(3)   = S
          neighbor_index_test(3) = 2

          bf_localization(4)     = W
          neighbor1_id_test(4)   = S
          neighbor_index_test(4) = 1


          do k=1,4

             !output
             call bf_layer_used%ini(bf_localization(k))
             call bf_layer_used%get_neighbor1_id(
     $            neighbor1_id,
     $            neighbor_index)

             !validation
             test_loc =
     $            (neighbor1_id.eq.neighbor1_id_test(k)).and.
     $            (neighbor_index.eq.neighbor_index_test(k))

             !detailled
             if((.not.test_loc).and.detailled) then
                print '(''test '',I2,'' : '',L1)', k, test_loc
                print '(''  - neighbor1_id  : '',I2,'' -> '',I2)', neighbor1_id, neighbor1_id_test(k)
                print '(''  - neighbor_index: '',I2,'' -> '',I2)', neighbor_index, neighbor_index_test(k)
                print '()'
             end if

          end do

        end function test_get_neighbor1_id


        function test_get_neighbor2_id(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_sync)          :: bf_layer_used
          integer       , dimension(4) :: bf_localization
          integer       , dimension(4) :: neighbor2_id_test
          integer       , dimension(4) :: neighbor_index_test
          integer                      :: neighbor2_id
          integer                      :: neighbor_index

          logical :: test_loc
          integer :: k

          test_validated = .true.


          bf_localization(1)     = N
          neighbor2_id_test(1)   = E
          neighbor_index_test(1) = 2

          bf_localization(2)     = S
          neighbor2_id_test(2)   = E
          neighbor_index_test(2) = 1

          bf_localization(3)     = E
          neighbor2_id_test(3)   = N
          neighbor_index_test(3) = 2

          bf_localization(4)     = W
          neighbor2_id_test(4)   = N
          neighbor_index_test(4) = 1


          do k=1,4

             !output
             call bf_layer_used%ini(bf_localization(k))
             call bf_layer_used%get_neighbor2_id(
     $            neighbor2_id,
     $            neighbor_index)

             !validation
             test_loc =
     $            (neighbor2_id.eq.neighbor2_id_test(k)).and.
     $            (neighbor_index.eq.neighbor_index_test(k))

             !detailled
             if((.not.test_loc).and.detailled) then
                print '(''test '',I2,'' : '',L1)', k, test_loc
                print '(''  - neighbor2_id  : '',I2,'' -> '',I2)', neighbor2_id, neighbor2_id_test(k)
                print '(''  - neighbor_index: '',I2,'' -> '',I2)', neighbor_index, neighbor_index_test(k)
                print '()'
             end if

          end do

        end function test_get_neighbor2_id


        function test_copy_from_neighbor1(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_sync) :: bf_layer_used
          type(bf_layer_sync) :: nbf_layer_used

          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes
          real(rkind)   , dimension(:,:,:), allocatable :: nbf_nodes
          integer       , dimension(:,:)  , allocatable :: bf_grdpts_id
          integer       , dimension(:,:)  , allocatable :: nbf_grdpts_id
          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes_test
          real(rkind)   , dimension(:,:,:), allocatable :: nbf_nodes_test
          integer       , dimension(:,:)  , allocatable :: bf_grdpts_id_test
          integer       , dimension(:,:)  , allocatable :: nbf_grdpts_id_test
          
          logical        :: test_loc

          test_validated = .true.

          !input
          call ini_for_neighbor_test(
     $         bf_layer_used,
     $         nbf_layer_used,
     $         bf_nodes,
     $         bf_nodes_test,
     $         bf_grdpts_id,
     $         bf_grdpts_id_test,
     $         nbf_nodes,
     $         nbf_nodes_test,
     $         nbf_grdpts_id,
     $         nbf_grdpts_id_test,
     $         test_id_copy_from_neighbor1)

          !output
          call bf_layer_used%copy_from_neighbor1(nbf_layer_used)
          call bf_layer_used%get_nodes_array(bf_nodes)
          call bf_layer_used%get_grdpts_id(bf_grdpts_id)
          call nbf_layer_used%get_nodes_array(nbf_nodes)
          call nbf_layer_used%get_grdpts_id(nbf_grdpts_id)

          !validation
          test_loc = is_real_matrix3D_validated(bf_nodes,bf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_real_matrix3D_validated(nbf_nodes,nbf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_int_matrix_validated(bf_grdpts_id,bf_grdpts_id_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_int_matrix_validated(nbf_grdpts_id,nbf_grdpts_id_test,detailled)
          test_validated = test_validated.and.test_loc

       end function test_copy_from_neighbor1


       function test_copy_to_neighbor1(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_sync)          :: bf_layer_used
          type(bf_layer_sync)          :: nbf_layer_used

          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes
          real(rkind)   , dimension(:,:,:), allocatable :: nbf_nodes
          integer       , dimension(:,:)  , allocatable :: bf_grdpts_id
          integer       , dimension(:,:)  , allocatable :: nbf_grdpts_id
          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes_test
          real(rkind)   , dimension(:,:,:), allocatable :: nbf_nodes_test
          integer       , dimension(:,:)  , allocatable :: bf_grdpts_id_test
          integer       , dimension(:,:)  , allocatable :: nbf_grdpts_id_test
          
          logical        :: test_loc          

          test_validated = .true.

          !input
          call ini_for_neighbor_test(
     $         bf_layer_used,
     $         nbf_layer_used,
     $         bf_nodes,
     $         bf_nodes_test,
     $         bf_grdpts_id,
     $         bf_grdpts_id_test,
     $         nbf_nodes,
     $         nbf_nodes_test,
     $         nbf_grdpts_id,
     $         nbf_grdpts_id_test,
     $         test_id_copy_to_neighbor1)

          !output
          call bf_layer_used%copy_to_neighbor1(nbf_layer_used)
          call bf_layer_used%get_nodes_array(bf_nodes)
          call bf_layer_used%get_grdpts_id(bf_grdpts_id)
          call nbf_layer_used%get_nodes_array(nbf_nodes)
          call nbf_layer_used%get_grdpts_id(nbf_grdpts_id)

          !validation
          test_loc = is_real_matrix3D_validated(bf_nodes,bf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_real_matrix3D_validated(nbf_nodes,nbf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_int_matrix_validated(bf_grdpts_id,bf_grdpts_id_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_int_matrix_validated(nbf_grdpts_id,nbf_grdpts_id_test,detailled)
          test_validated = test_validated.and.test_loc

       end function test_copy_to_neighbor1


       subroutine ini_for_neighbor_test(
     $     bf_layer_used,
     $     nbf_layer_used,
     $     bf_nodes,
     $     bf_nodes_test,
     $     bf_grdpts_id,
     $     bf_grdpts_id_test,
     $     nbf_nodes,
     $     nbf_nodes_test,
     $     nbf_grdpts_id,
     $     nbf_grdpts_id_test,
     $     test_id)

         implicit none

         type(bf_layer_sync)                          , intent(inout) :: bf_layer_used
         type(bf_layer_sync)                          , intent(inout) :: nbf_layer_used
         real(rkind)   , dimension(:,:,:), allocatable, intent(inout) :: bf_nodes
         real(rkind)   , dimension(:,:,:), allocatable, intent(out)   :: bf_nodes_test
         integer       , dimension(:,:)  , allocatable, intent(out)   :: bf_grdpts_id
         integer       , dimension(:,:)  , allocatable, intent(out)   :: bf_grdpts_id_test
         real(rkind)   , dimension(:,:,:), allocatable, intent(out)   :: nbf_nodes
         real(rkind)   , dimension(:,:,:), allocatable, intent(out)   :: nbf_nodes_test
         integer       , dimension(:,:)  , allocatable, intent(out)   :: nbf_grdpts_id
         integer       , dimension(:,:)  , allocatable, intent(out)   :: nbf_grdpts_id_test
         integer                                      , intent(in)    :: test_id

         integer                                       :: bf_localization
         integer                                       :: nbf_localization
         integer(ikind), dimension(2,2)                :: bf_alignment
         integer(ikind), dimension(2,2)                :: nbf_alignment

         integer(ikind) :: i,j
         integer        :: k
         
         integer(ikind), dimension(2) :: bf_recv
         integer(ikind), dimension(2) :: bf_send
         integer(ikind), dimension(2) :: nbf_send
         integer(ikind), dimension(2) :: nbf_recv
         integer(ikind), dimension(2) :: ex_size

         bf_localization  = N
         nbf_localization = W

         bf_alignment(1,1) = align_W-3
         bf_alignment(1,2) = align_W+4
         bf_alignment(2,1) = align_N
         bf_alignment(2,2) = align_N+1

         nbf_alignment(1,1) = align_W-1
         nbf_alignment(1,2) = align_W
         nbf_alignment(2,1) = align_N-2
         nbf_alignment(2,2) = align_N-1


         allocate(bf_nodes(
     $        bf_alignment(1,2)-bf_alignment(1,1)+1+2*bc_size,
     $        bf_alignment(2,2)-bf_alignment(2,1)+1+2*bc_size,
     $        ne))

         allocate(nbf_nodes(
     $        nbf_alignment(1,2)-nbf_alignment(1,1)+1+2*bc_size,
     $        nbf_alignment(2,2)-nbf_alignment(2,1)+1+2*bc_size,
     $        ne))

         allocate(bf_grdpts_id(
     $        bf_alignment(1,2)-bf_alignment(1,1)+1+2*bc_size,
     $        bf_alignment(2,2)-bf_alignment(2,1)+1+2*bc_size))

         allocate(nbf_grdpts_id(
     $        nbf_alignment(1,2)-nbf_alignment(1,1)+1+2*bc_size,
     $        nbf_alignment(2,2)-nbf_alignment(2,1)+1+2*bc_size))
         

         do k=1,ne
            do j=1, size(bf_nodes,2)
               do i=1, size(bf_nodes,1)
                  bf_nodes(i,j,k) = (i-1) + 10*(j-1) + 100*k
               end do
            end do
         end do

         do j=1, size(bf_nodes,2)
            do i=1, size(bf_nodes,1)
               bf_grdpts_id(i,j) = -(i-1) - 10*(j-1)
            end do
         end do

         do k=1,ne
            do j=1, size(nbf_nodes,2)
               do i=1, size(nbf_nodes,1)
                  nbf_nodes(i,j,k) = (i-1) + 50*(j-1) + 500*k
               end do
            end do
         end do

         do j=1, size(nbf_nodes,2)
            do i=1, size(nbf_nodes,1)
               nbf_grdpts_id(i,j) = -(i-1) - 50*(j-1)
            end do
         end do

         allocate(bf_nodes_test,source=bf_nodes)
         allocate(nbf_nodes_test,source=nbf_nodes)
         allocate(bf_grdpts_id_test,source=bf_grdpts_id)
         allocate(nbf_grdpts_id_test,source=nbf_grdpts_id)

         select case(test_id)
           case(test_id_copy_from_neighbor1)
              nbf_send = [1,3]
              ex_size  = [6,4]
              bf_recv  = [3,1]

              bf_nodes_test(
     $             bf_recv(1):bf_recv(1)+ex_size(1)-1,
     $             bf_recv(2):bf_recv(2)+ex_size(2)-1,
     $             :)
     $             =
     $             nbf_nodes_test(
     $             nbf_send(1):nbf_send(1)+ex_size(1)-1,
     $             nbf_send(2):nbf_send(2)+ex_size(2)-1,
     $             :)
              
              bf_grdpts_id_test(
     $             bf_recv(1):bf_recv(1)+ex_size(1)-1,
     $             bf_recv(2):bf_recv(2)+ex_size(2)-1)
     $             =
     $             nbf_grdpts_id_test(
     $             nbf_send(1):nbf_send(1)+ex_size(1)-1,
     $             nbf_send(2):nbf_send(2)+ex_size(2)-1)

           case(test_id_copy_to_neighbor1)
              nbf_recv = [1,3]
              ex_size  = [6,4]
              bf_send  = [3,1]

              nbf_nodes_test(
     $             nbf_recv(1):nbf_recv(1)+ex_size(1)-1,
     $             nbf_recv(2):nbf_recv(2)+ex_size(2)-1,
     $             :)
     $             =
     $             bf_nodes_test(
     $             bf_send(1):bf_send(1)+ex_size(1)-1,
     $             bf_send(2):bf_send(2)+ex_size(2)-1,
     $             :)
              
              nbf_grdpts_id_test(
     $             nbf_recv(1):nbf_recv(1)+ex_size(1)-1,
     $             nbf_recv(2):nbf_recv(2)+ex_size(2)-1)
     $             =
     $             bf_grdpts_id_test(
     $             bf_send(1):bf_send(1)+ex_size(1)-1,
     $             bf_send(2):bf_send(2)+ex_size(2)-1)

           case(test_id_copy_from_neighbor2)
              nbf_recv = [1,3]
              ex_size  = [6,4]
              bf_send  = [3,1]

              nbf_nodes_test(
     $             nbf_recv(1):nbf_recv(1)+ex_size(1)-1,
     $             nbf_recv(2):nbf_recv(2)+ex_size(2)-1,
     $             :)
     $             =
     $             bf_nodes_test(
     $             bf_send(1):bf_send(1)+ex_size(1)-1,
     $             bf_send(2):bf_send(2)+ex_size(2)-1,
     $             :)
              
              nbf_grdpts_id_test(
     $             nbf_recv(1):nbf_recv(1)+ex_size(1)-1,
     $             nbf_recv(2):nbf_recv(2)+ex_size(2)-1)
     $             =
     $             bf_grdpts_id_test(
     $             bf_send(1):bf_send(1)+ex_size(1)-1,
     $             bf_send(2):bf_send(2)+ex_size(2)-1)

           case(test_id_copy_to_neighbor2)
              nbf_send = [1,3]
              ex_size  = [6,4]
              bf_recv  = [3,1]

              bf_nodes_test(
     $             bf_recv(1):bf_recv(1)+ex_size(1)-1,
     $             bf_recv(2):bf_recv(2)+ex_size(2)-1,
     $             :)
     $             =
     $             nbf_nodes_test(
     $             nbf_send(1):nbf_send(1)+ex_size(1)-1,
     $             nbf_send(2):nbf_send(2)+ex_size(2)-1,
     $             :)
              
              bf_grdpts_id_test(
     $             bf_recv(1):bf_recv(1)+ex_size(1)-1,
     $             bf_recv(2):bf_recv(2)+ex_size(2)-1)
     $             =
     $             nbf_grdpts_id_test(
     $             nbf_send(1):nbf_send(1)+ex_size(1)-1,
     $             nbf_send(2):nbf_send(2)+ex_size(2)-1)

           case(test_id_sync_with_neighbor)
              bf_send  = [3,3]
              bf_recv  = [3,1]
              nbf_send = [1,3]
              nbf_recv = [1,5]
              ex_size  = [6,2]

              bf_nodes_test(
     $             bf_recv(1):bf_recv(1)+ex_size(1)-1,
     $             bf_recv(2):bf_recv(2)+ex_size(2)-1,
     $             :)
     $             =
     $             nbf_nodes_test(
     $             nbf_send(1):nbf_send(1)+ex_size(1)-1,
     $             nbf_send(2):nbf_send(2)+ex_size(2)-1,
     $             :)

              nbf_nodes_test(
     $             nbf_recv(1):nbf_recv(1)+ex_size(1)-1,
     $             nbf_recv(2):nbf_recv(2)+ex_size(2)-1,
     $             :)
     $             =
     $             bf_nodes_test(
     $             bf_send(1):bf_send(1)+ex_size(1)-1,
     $             bf_send(2):bf_send(2)+ex_size(2)-1,
     $             :)
              
              bf_grdpts_id_test(
     $             bf_recv(1):bf_recv(1)+ex_size(1)-1,
     $             bf_recv(2):bf_recv(2)+ex_size(2)-1)
     $             =
     $             nbf_grdpts_id_test(
     $             nbf_send(1):nbf_send(1)+ex_size(1)-1,
     $             nbf_send(2):nbf_send(2)+ex_size(2)-1)

              nbf_grdpts_id_test(
     $             nbf_recv(1):nbf_recv(1)+ex_size(1)-1,
     $             nbf_recv(2):nbf_recv(2)+ex_size(2)-1)
     $             =
     $             bf_grdpts_id_test(
     $             bf_send(1):bf_send(1)+ex_size(1)-1,
     $             bf_send(2):bf_send(2)+ex_size(2)-1)

         end select
              
         call bf_layer_used%ini(bf_localization)
         call nbf_layer_used%ini(nbf_localization)

         call bf_layer_used%set_alignment_tab(bf_alignment)
         call nbf_layer_used%set_alignment_tab(nbf_alignment)

         call bf_layer_used%set_nodes(bf_nodes)
         call nbf_layer_used%set_nodes(nbf_nodes)
         call bf_layer_used%set_grdpts_id(bf_grdpts_id)
         call nbf_layer_used%set_grdpts_id(nbf_grdpts_id)

       end subroutine ini_for_neighbor_test


       function test_copy_from_neighbor2(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_sync) :: bf_layer_used
          type(bf_layer_sync) :: nbf_layer_used

          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes
          real(rkind)   , dimension(:,:,:), allocatable :: nbf_nodes
          integer       , dimension(:,:)  , allocatable :: bf_grdpts_id
          integer       , dimension(:,:)  , allocatable :: nbf_grdpts_id
          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes_test
          real(rkind)   , dimension(:,:,:), allocatable :: nbf_nodes_test
          integer       , dimension(:,:)  , allocatable :: bf_grdpts_id_test
          integer       , dimension(:,:)  , allocatable :: nbf_grdpts_id_test
          
          logical        :: test_loc

          test_validated = .true.

          !input
          call ini_for_neighbor_test(
     $         bf_layer_used,
     $         nbf_layer_used,
     $         bf_nodes,
     $         bf_nodes_test,
     $         bf_grdpts_id,
     $         bf_grdpts_id_test,
     $         nbf_nodes,
     $         nbf_nodes_test,
     $         nbf_grdpts_id,
     $         nbf_grdpts_id_test,
     $         test_id_copy_from_neighbor2)

          !output
          call nbf_layer_used%copy_from_neighbor2(bf_layer_used)
          call bf_layer_used%get_nodes_array(bf_nodes)
          call bf_layer_used%get_grdpts_id(bf_grdpts_id)
          call nbf_layer_used%get_nodes_array(nbf_nodes)
          call nbf_layer_used%get_grdpts_id(nbf_grdpts_id)

          !validation
          test_loc = is_real_matrix3D_validated(bf_nodes,bf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_real_matrix3D_validated(nbf_nodes,nbf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_int_matrix_validated(bf_grdpts_id,bf_grdpts_id_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_int_matrix_validated(nbf_grdpts_id,nbf_grdpts_id_test,detailled)
          test_validated = test_validated.and.test_loc

       end function test_copy_from_neighbor2


       function test_copy_to_neighbor2(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_sync)          :: bf_layer_used
          type(bf_layer_sync)          :: nbf_layer_used

          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes
          real(rkind)   , dimension(:,:,:), allocatable :: nbf_nodes
          integer       , dimension(:,:)  , allocatable :: bf_grdpts_id
          integer       , dimension(:,:)  , allocatable :: nbf_grdpts_id
          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes_test
          real(rkind)   , dimension(:,:,:), allocatable :: nbf_nodes_test
          integer       , dimension(:,:)  , allocatable :: bf_grdpts_id_test
          integer       , dimension(:,:)  , allocatable :: nbf_grdpts_id_test
          
          logical        :: test_loc          

          test_validated = .true.

          !input
          call ini_for_neighbor_test(
     $         bf_layer_used,
     $         nbf_layer_used,
     $         bf_nodes,
     $         bf_nodes_test,
     $         bf_grdpts_id,
     $         bf_grdpts_id_test,
     $         nbf_nodes,
     $         nbf_nodes_test,
     $         nbf_grdpts_id,
     $         nbf_grdpts_id_test,
     $         test_id_copy_to_neighbor2)

          !output
          call nbf_layer_used%copy_to_neighbor2(bf_layer_used)
          call bf_layer_used%get_nodes_array(bf_nodes)
          call bf_layer_used%get_grdpts_id(bf_grdpts_id)
          call nbf_layer_used%get_nodes_array(nbf_nodes)
          call nbf_layer_used%get_grdpts_id(nbf_grdpts_id)

          !validation
          test_loc = is_real_matrix3D_validated(bf_nodes,bf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_real_matrix3D_validated(nbf_nodes,nbf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_int_matrix_validated(bf_grdpts_id,bf_grdpts_id_test,detailled)
          test_validated = test_validated.and.test_loc
          
          test_loc = is_int_matrix_validated(nbf_grdpts_id,nbf_grdpts_id_test,detailled)
          test_validated = test_validated.and.test_loc

       end function test_copy_to_neighbor2


       function test_sync_nodes_with_interior(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_sync)                           :: bf_layer_used
          integer                                       :: bf_localization
          integer(ikind), dimension(2,2)                :: bf_alignment
          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne)           :: interior_nodes
          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes_test
          real(rkind)   , dimension(nx,ny,ne)           :: interior_nodes_test

          integer(ikind) :: i,j
          integer        :: k
         
          integer(ikind), dimension(2) :: bf_recv
          integer(ikind), dimension(2) :: bf_send
          integer(ikind), dimension(2) :: int_recv
          integer(ikind), dimension(2) :: int_send
          integer(ikind), dimension(2) :: ex_size

          logical :: test_loc

          test_validated = .true.


          !input
          bf_localization  = W

          bf_alignment(1,1) = align_W-1
          bf_alignment(1,2) = align_W
          bf_alignment(2,1) = align_S+5
          bf_alignment(2,2) = align_S+8

          allocate(bf_nodes(
     $         bf_alignment(1,2)-bf_alignment(1,1)+1+2*bc_size,
     $         bf_alignment(2,2)-bf_alignment(2,1)+1+2*bc_size,
     $         ne))

          do k=1,ne
             do j=1, size(bf_nodes,2)
                do i=1, size(bf_nodes,1)
                   bf_nodes(i,j,k) = (i-1) + 10*(j-1) + 100*k
                end do
             end do
          end do

          do k=1,ne
             do j=1, size(interior_nodes,2)
                do i=1, size(interior_nodes,1)
                   interior_nodes(i,j,k) = -(i-1) - 10*(j-1) - 100*k
                end do
             end do
          end do

          allocate(bf_nodes_test,source=bf_nodes)
          interior_nodes_test = interior_nodes
          
          if((align_S+10)>ny) then
             print '(''in test_sync_nodes_with_interior, align_S+10>ny'')'
             stop 'change inputs to have align_S+10.le.ny'
          end if

          bf_send  = [3,1]
          bf_recv  = [5,1]
          int_send = [3,align_S+3]
          int_recv = [1,align_S+3]
          ex_size  = [bc_size,size(bf_nodes,2)]

          bf_nodes_test(
     $         bf_recv(1):bf_recv(1)+ex_size(1)-1,
     $         bf_recv(2):bf_recv(2)+ex_size(2)-1,
     $         :)
     $         =
     $         interior_nodes_test(
     $         int_send(1):int_send(1)+ex_size(1)-1,
     $         int_send(2):int_send(2)+ex_size(2)-1,
     $         :)

          interior_nodes_test(
     $         int_recv(1):int_recv(1)+ex_size(1)-1,
     $         int_recv(2):int_recv(2)+ex_size(2)-1,
     $         :)
     $         =
     $         bf_nodes_test(
     $         bf_send(1):bf_send(1)+ex_size(1)-1,
     $         bf_send(2):bf_send(2)+ex_size(2)-1,
     $         :)

          call bf_layer_used%ini(bf_localization)
          call bf_layer_used%set_alignment_tab(bf_alignment)
          call bf_layer_used%set_nodes(bf_nodes)
          

          !output
          call bf_layer_used%sync_nodes_with_interior(interior_nodes)
          call bf_layer_used%get_nodes_array(bf_nodes)

          !validation
          test_loc = is_real_matrix3D_validated(interior_nodes,interior_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_real_matrix3D_validated(bf_nodes,bf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

        end function test_sync_nodes_with_interior


        function test_sync_nodes_with_neighbor1(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_sync) :: bf_layer_used
          type(bf_layer_sync) :: nbf_layer_used

          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes
          real(rkind)   , dimension(:,:,:), allocatable :: nbf_nodes
          integer       , dimension(:,:)  , allocatable :: bf_grdpts_id
          integer       , dimension(:,:)  , allocatable :: nbf_grdpts_id
          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes_test
          real(rkind)   , dimension(:,:,:), allocatable :: nbf_nodes_test
          integer       , dimension(:,:)  , allocatable :: bf_grdpts_id_test
          integer       , dimension(:,:)  , allocatable :: nbf_grdpts_id_test
          
          logical        :: test_loc

          test_validated = .true.

          !input
          call ini_for_neighbor_test(
     $         bf_layer_used,
     $         nbf_layer_used,
     $         bf_nodes,
     $         bf_nodes_test,
     $         bf_grdpts_id,
     $         bf_grdpts_id_test,
     $         nbf_nodes,
     $         nbf_nodes_test,
     $         nbf_grdpts_id,
     $         nbf_grdpts_id_test,
     $         test_id_sync_with_neighbor)

          !output
          call bf_layer_used%sync_nodes_with_neighbor1(nbf_layer_used)
          call bf_layer_used%get_nodes_array(bf_nodes)
          call bf_layer_used%get_grdpts_id(bf_grdpts_id)
          call nbf_layer_used%get_nodes_array(nbf_nodes)
          call nbf_layer_used%get_grdpts_id(nbf_grdpts_id)

          !validation
          test_loc = is_real_matrix3D_validated(bf_nodes,bf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_real_matrix3D_validated(nbf_nodes,nbf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

       end function test_sync_nodes_with_neighbor1


       function test_sync_nodes_with_neighbor2(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_sync) :: bf_layer_used
          type(bf_layer_sync) :: nbf_layer_used

          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes
          real(rkind)   , dimension(:,:,:), allocatable :: nbf_nodes
          integer       , dimension(:,:)  , allocatable :: bf_grdpts_id
          integer       , dimension(:,:)  , allocatable :: nbf_grdpts_id
          real(rkind)   , dimension(:,:,:), allocatable :: bf_nodes_test
          real(rkind)   , dimension(:,:,:), allocatable :: nbf_nodes_test
          integer       , dimension(:,:)  , allocatable :: bf_grdpts_id_test
          integer       , dimension(:,:)  , allocatable :: nbf_grdpts_id_test
          
          logical        :: test_loc

          test_validated = .true.

          !input
          call ini_for_neighbor_test(
     $         bf_layer_used,
     $         nbf_layer_used,
     $         bf_nodes,
     $         bf_nodes_test,
     $         bf_grdpts_id,
     $         bf_grdpts_id_test,
     $         nbf_nodes,
     $         nbf_nodes_test,
     $         nbf_grdpts_id,
     $         nbf_grdpts_id_test,
     $         test_id_sync_with_neighbor)

          !output
          call nbf_layer_used%sync_nodes_with_neighbor2(bf_layer_used)
          call bf_layer_used%get_nodes_array(bf_nodes)
          call bf_layer_used%get_grdpts_id(bf_grdpts_id)
          call nbf_layer_used%get_nodes_array(nbf_nodes)
          call nbf_layer_used%get_grdpts_id(nbf_grdpts_id)

          !validation
          test_loc = is_real_matrix3D_validated(bf_nodes,bf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_real_matrix3D_validated(nbf_nodes,nbf_nodes_test,detailled)
          test_validated = test_validated.and.test_loc

       end function test_sync_nodes_with_neighbor2

      end program test_bf_layer_sync
