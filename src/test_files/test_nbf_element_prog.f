      program test_nbf_element_prog

        use bf_mainlayer_class  , only : bf_mainlayer
        use bf_sublayer_class   , only : bf_sublayer
        use nbf_element_class   , only : nbf_element
        use parameters_bf_layer , only : align_N, align_S,
     $                                   align_E, align_W
        use parameters_constant , only : N,S,E,W
        use parameters_input    , only : nx,ny,ne,bc_size
        use parameters_kind     , only : ikind, rkind
        use test_bf_layer_module, only : print_interior_data,
     $                                   ini_nodes,
     $                                   ini_grdpts_id,
     $                                   ini_cst_nodes


        implicit none


        real(rkind)       , dimension(nx,ny,ne) :: nodes
        integer           , dimension(nx,ny)    :: grdpts_id
        integer           , dimension(4)        :: bf_layer_loc
        character(2)      , dimension(4)        :: bf_layer_char

        type(bf_mainlayer), dimension(4)        :: bf_mainlayers
        type(bf_sublayer) , pointer             :: added_sublayer
        type(nbf_element) , dimension(4,3)      :: nbf_elements
        integer(ikind)    , dimension(4,3,2,2)  :: test_alignment
        integer(ikind)    , dimension(2,2)      :: alignment
        logical                                 :: test
        integer           , dimension(4,2)      :: neighbors
        real(rkind)                             :: scale

        integer(ikind) :: g,h,i,j
        

        !initialize the nodes and the grdpts_id
        call ini_nodes(nodes)
        call ini_grdpts_id(grdpts_id)

        !print the nodes
        call print_interior_data(
     $       nodes, grdpts_id,
     $       'interior_nodes.dat',
     $       'interior_grdpts_id.dat',
     $       'interior_sizes.dat')

        !buffer layers tested
        bf_layer_loc  = [N,S,E,W]
        bf_layer_char = ['N_','S_','E_','W_']
        
        neighbors(N,1) = W
        neighbors(N,2) = E
        neighbors(S,1) = W
        neighbors(S,2) = E
        neighbors(W,1) = S
        neighbors(W,2) = N
        neighbors(E,1) = S
        neighbors(E,2) = N


        !initialize the alignment for all the buffer layers
        call ini_alignment(test_alignment)


        !initialize the main layers
        do j=1,4
           call bf_mainlayers(j)%ini(j)
        end do


        !add the buffer layers to the mainlayers
        !loop over the cardinal coordinates
        print '()'
        print '(''buffer layer alignments'')'
        print '(''-----------------------'')'
        do j=1,4

           !loop over the number of buffer layers
           !per cardinal coordinate
           do i=1,3

              !alignment for the buffer layer
              do g=1,2
                 do h=1,2
                    alignment(h,g) = test_alignment(j,i,h,g)
                 end do
              end do

              !print the alignment
              print '(A2,I1,'' ('',2I3,'')  ('',2I3,'')'')',
     $             bf_layer_char(j), i,
     $             alignment(1,1), alignment(1,2),
     $             alignment(2,1), alignment(2,2)

              !allocate the buffer layer
              added_sublayer => bf_mainlayers(j)%add_sublayer(
     $             nodes,alignment)

              !initialize the corresponding nbf_element
              !with the reference to the buffer sublayer
              call nbf_elements(j,i)%ini(added_sublayer)

           end do
        end do
        print '(''-----------------------'')'
        print '()'
              

        !print the main layers before the tests
        do j=1,4
           call bf_mainlayers(j)%print_binary(
     $          'nodes1.dat',
     $          'grdpt_id1.dat',
     $          'sizes1.dat')

        end do


        !test the fct is_before()
        print '()'
        print '(''test is_before() '')'
        print '(''--------------------------'')'
        do j=1,4
           do i=2,3
              test = nbf_elements(j,1)%is_before(nbf_elements(j,i))

              if(test) then
                 
                 print '(A1,'': ('',2I3,'') .le. ('',2I3,'')'')',
     $                bf_layer_char(j),
     $                test_alignment(j,1,1,1), test_alignment(j,1,1,2),
     $                test_alignment(j,i,1,1), test_alignment(j,i,1,2)

              else
                 print '(A1,'': ('',2I3,'') .gt. ('',2I3,'')'')',
     $                bf_layer_char(j),
     $                test_alignment(j,1,1,1), test_alignment(j,1,1,2),
     $                test_alignment(j,i,1,1), test_alignment(j,i,1,2)

              end if

           end do
        end do
        print '(''--------------------------'')'
        print '()'        

        
        !test the fct refers_to()
        print '()'
        print '(''test refers_to() '')'
        print '(''----------------------'')'
        do j=1,4
           do i=1,3
              
              added_sublayer => bf_mainlayers(j)%get_head_sublayer()
              do h=1,3
                 test = nbf_elements(j,i)%refers_to(added_sublayer)
                 added_sublayer => added_sublayer%get_next()

                 if(test) then
                    print '(''nbf('',A1,I2,'') ---> bf('',A1,I2,'')'')',
     $                   bf_layer_char(j),i,
     $                   bf_layer_char(j),h
                 else
                    print '(''nbf('',A1,I2,'') -X-> bf('',A1,I2,'')'')',
     $                   bf_layer_char(j),i,
     $                   bf_layer_char(j),h
                 end if

              end do
           end do
        end do
        print '(''----------------------'')'
        print '()'


        !test the fct can_exchange_with()
        print '()'
        print '(''can_exchange_with() '')'
        print '(''----------------------'')'
        do j=1,4
           do i=1,3,2
              
              !loop over the neighboring buffer layers
              do g=1,2
                 added_sublayer => bf_mainlayers(neighbors(j,g))%get_head_sublayer()

                 !loop over the number of buffer layers in the main layer
                 do h=1,3

                    test = nbf_elements(j,i)%can_exchange_with(added_sublayer)
                    added_sublayer => added_sublayer%get_next()
                    
                    if(test) then
                       call print_exchange(
     $                      bf_layer_char(j),i,
     $                      bf_layer_char(neighbors(j,g)),h)
                             
                    end if
                 end do

              end do
           end do
        end do
        print '(''----------------------'')'
        print '()'


        !initialization of the nodes of the different layers
        scale = 0.07
        do j=1,4
           added_sublayer => bf_mainlayers(j)%get_head_sublayer()

           do i=1,3
              
              !initialize the nodes of the buffer layer with a constant value
              call ini_cst_nodes(added_sublayer, (3*(j-1)+i)*scale)
              added_sublayer => added_sublayer%get_next()

           end do
        end do

        !test the fct copy_from_neighbor1_to()
        !test the fct copy_from_neighbor2_to()
        print '()'
        print '(''copy_from_neighbor1_to() '')'
        print '(''copy_from_neighbor2_to() '')'
        print '(''------------------------'')'
        do j=1,4
           added_sublayer => bf_mainlayers(j)%get_head_sublayer()
           do i=1,3
              
              !potential neighbor1 buffer layers
              g=1

              !loop over the number of buffer layers in the main layer
              do h=1,3,2

                 call nbf_elements(neighbors(j,g),h)%copy_from_neighbor1_to(added_sublayer)
                                  
                 call print_exchange_from_to(
     $                bf_layer_char(neighbors(j,g)),h,
     $                bf_layer_char(j),i)
                 
              end do

              !potential neighbor2 buffer layers
              g=2

              !loop over the number of buffer layers in the main layer
              do h=1,3,2

                 call nbf_elements(neighbors(j,g),h)%copy_from_neighbor2_to(added_sublayer)
                 
                 call print_exchange_from_to(
     $                bf_layer_char(neighbors(j,g)),h,
     $                bf_layer_char(j),i)

              end do
              added_sublayer => added_sublayer%get_next()

           end do
        end do
        print '(''------------------------'')'
        print '()'

        
        !print the main layers after copy from the neighbors
        do j=1,4
           call bf_mainlayers(j)%print_binary(
     $          'nodes2.dat',
     $          'grdpt_id2.dat',
     $          'sizes2.dat')

        end do


        !test the fct copy_to_neighbor1_from()
        !test the fct copy_to_neighbor2_from()
        print '()'
        print '(''copy_to_neighbor1_from() '')'
        print '(''copy_to_neighbor2_from() '')'
        print '(''------------------------'')'
        do j=1,4
           added_sublayer => bf_mainlayers(j)%get_head_sublayer()
           do i=1,3
              
              !potential neighbor1 buffer layers
              g=1

              !loop over the number of buffer layers in the main layer
              do h=1,3,2

                 call nbf_elements(neighbors(j,g),h)%copy_to_neighbor1_from(added_sublayer)
                                  
                 call print_exchange_from_to(
     $                bf_layer_char(neighbors(j,g)),h,
     $                bf_layer_char(j),i)
                 
              end do

              !potential neighbor2 buffer layers
              g=2

              !loop over the number of buffer layers in the main layer
              do h=1,3,2

                 call nbf_elements(neighbors(j,g),h)%copy_to_neighbor2_from(added_sublayer)
                 
                 call print_exchange_from_to(
     $                bf_layer_char(neighbors(j,g)),h,
     $                bf_layer_char(j),i)

              end do
              added_sublayer => added_sublayer%get_next()

           end do
        end do
        print '(''------------------------'')'
        print '()'


        !print the main layers after copy to the neighbors
        do j=1,4
           call bf_mainlayers(j)%print_binary(
     $          'nodes3.dat',
     $          'grdpt_id3.dat',
     $          'sizes3.dat')

        end do


        contains

        subroutine ini_alignment(
     $     test_alignment)

          implicit none

          integer(ikind), dimension(4,3,2,2), intent(out) :: test_alignment
          
          integer(ikind) :: large_layer
          integer(ikind) :: small_layer

          large_layer = 5
          small_layer = 3

          !north
          test_alignment(N,1,1,1) = align_W+1
          test_alignment(N,1,1,2) = test_alignment(N,1,1,1)+large_layer
          test_alignment(N,1,2,1) = align_N
          test_alignment(N,1,2,2) = align_N

          test_alignment(N,2,1,1) = test_alignment(N,1,1,2)+2*bc_size+small_layer
          test_alignment(N,2,1,2) = test_alignment(N,2,1,1)+large_layer
          test_alignment(N,2,2,1) = align_N
          test_alignment(N,2,2,2) = align_N

          test_alignment(N,3,1,2) = align_E-1
          test_alignment(N,3,1,1) = test_alignment(N,3,1,2)-2*bc_size-large_layer
          test_alignment(N,3,2,1) = align_N
          test_alignment(N,3,2,2) = align_N

          !south
          test_alignment(S,1,1,1) = align_W+1
          test_alignment(S,1,1,2) = test_alignment(S,1,1,1)+large_layer
          test_alignment(S,1,2,1) = align_S
          test_alignment(S,1,2,2) = align_S

          test_alignment(S,2,1,1) = test_alignment(S,1,1,2)+2*bc_size+small_layer
          test_alignment(S,2,1,2) = test_alignment(S,2,1,1)+large_layer
          test_alignment(S,2,2,1) = align_S
          test_alignment(S,2,2,2) = align_S

          test_alignment(S,3,1,2) = align_E-1
          test_alignment(S,3,1,1) = test_alignment(S,3,1,2)-2*bc_size-large_layer
          test_alignment(S,3,2,1) = align_S
          test_alignment(S,3,2,2) = align_S
          
          !east
          test_alignment(E,1,1,1) = align_E
          test_alignment(E,1,1,2) = align_E
          test_alignment(E,1,2,1) = align_S+1
          test_alignment(E,1,2,2) = test_alignment(E,1,2,1)+large_layer

          test_alignment(E,2,1,1) = align_E
          test_alignment(E,2,1,2) = align_E
          test_alignment(E,2,2,1) = test_alignment(E,1,2,2)+2*bc_size+small_layer
          test_alignment(E,2,2,2) = test_alignment(E,2,2,1)+large_layer

          test_alignment(E,3,1,1) = align_E
          test_alignment(E,3,1,2) = align_E
          test_alignment(E,3,2,2) = align_N-1
          test_alignment(E,3,2,1) = test_alignment(E,3,2,2)-2*bc_size-large_layer


          !west
          test_alignment(W,1,1,1) = align_W
          test_alignment(W,1,1,2) = align_W
          test_alignment(W,1,2,1) = align_S+1
          test_alignment(W,1,2,2) = test_alignment(W,1,2,1)+large_layer

          test_alignment(W,2,1,1) = align_W
          test_alignment(W,2,1,2) = align_W
          test_alignment(W,2,2,1) = test_alignment(W,1,2,2)+2*bc_size+small_layer
          test_alignment(W,2,2,2) = test_alignment(W,2,2,1)+large_layer

          test_alignment(W,3,1,1) = align_W
          test_alignment(W,3,1,2) = align_W
          test_alignment(W,3,2,2) = align_N-1
          test_alignment(W,3,2,1) = test_alignment(W,3,2,2)-2*bc_size-large_layer

        end subroutine ini_alignment


        subroutine print_exchange(
     $     mainlayer_id1, sublayer_id1,
     $     mainlayer_id2, sublayer_id2)

          implicit none

          character(*), intent(in) :: mainlayer_id1
          integer     , intent(in) :: sublayer_id1
          character(*), intent(in) :: mainlayer_id2
          integer     , intent(in) :: sublayer_id2

          print '(''nbf('',A1,I2,'') <---> bf('',A1,I2,'')'')',
     $         mainlayer_id1, sublayer_id1,
     $         mainlayer_id2, sublayer_id2


        end subroutine print_exchange

      
        subroutine print_exchange_from_to(
     $     mainlayer_id1, sublayer_id1,
     $     mainlayer_id2, sublayer_id2)

          implicit none

          character(*), intent(in) :: mainlayer_id1
          integer     , intent(in) :: sublayer_id1
          character(*), intent(in) :: mainlayer_id2
          integer     , intent(in) :: sublayer_id2

          print '(''nbf('',A1,I2,'') ---> bf('',A1,I2,'')'')',
     $         mainlayer_id1, sublayer_id1,
     $         mainlayer_id2, sublayer_id2


        end subroutine print_exchange_from_to

      end program test_nbf_element_prog
