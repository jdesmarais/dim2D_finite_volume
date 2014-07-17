      program test_nbf_interface_prog

        use nbf_interface_class , only : nbf_interface
        use bf_mainlayer_class  , only : bf_mainlayer
        use bf_sublayer_class   , only : bf_sublayer
        use parameters_constant , only : N,S,E,W
        use parameters_input    , only : nx,ny,ne
        use parameters_kind     , only : ikind, rkind
        use test_bf_layer_module, only : print_interior_data,
     $                                   ini_nodes,
     $                                   ini_grdpts_id,
     $                                   ini_cst_nodes
        use test_nbf_list_module, only : ini_alignment
        
        implicit none

        
        real(rkind)        , dimension(nx,ny,ne) :: nodes
        integer            , dimension(nx,ny)    :: grdpts_id
        integer            , dimension(4)        :: bf_layer_loc
        character(2)       , dimension(4)        :: bf_layer_char
                           
        type(bf_mainlayer) , dimension(4)        :: bf_mainlayers
        type(nbf_interface)                      :: nbf_interface_used
        type(bf_sublayer)  , pointer             :: added_sublayer
        !type(bf_sublayer)  , pointer             :: added_sublayer2
                           
        integer(ikind)     , dimension(4,5,2,2)  :: test_alignment
        integer(ikind)     , dimension(2,2)      :: alignment
        integer            , dimension(4,2)      :: neighbors
        real(rkind)                              :: scale

        integer :: g,h,i,j,i_max

        real(rkind) :: dx,dy

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
        
        !initialize the nbf_interface
        call nbf_interface_used%ini()


        !add the buffer layers to the mainlayers
        !+ add the links to the nbf_lists
        !loop over the cardinal coordinates
        print '()'
        print '(''buffer layer alignments'')'
        print '(''-----------------------'')'
        do j=1,4

           select case(j)
             case(N,S)
                i_max = 5
             case(E,W)
                i_max = 3
           end select

           !loop over the number of buffer layers
           !per main layer
           do i=1,i_max

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
     $             nodes,alignment,dx,dy)

              !set whether the buffe rlayer can exchange
              !with neighboring layers
              call added_sublayer%set_neighbor1_share()
              call added_sublayer%set_neighbor2_share()

              !if the buffer layer is a potential neigboring buffer
              !layer, we add the link to the nbf_interface
              if(added_sublayer%can_exchange_with_neighbor1()) then
                 call nbf_interface_used%link_neighbor1_to_bf_sublayer(
     $                added_sublayer)
              end if

              if(added_sublayer%can_exchange_with_neighbor2()) then
                 call nbf_interface_used%link_neighbor2_to_bf_sublayer(
     $                added_sublayer)
              end if

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

        
        !print the content of the nbf_interface
        print '()'
        print '(''nbf_interface content'')'
        print '(''-----------------------'')'
        call nbf_interface_used%print_on_screen()
        print '(''-----------------------'')'
        print '()'

        
        !initialization of the nodes of the different layers
        scale = 0.07
        do j=1,4
           added_sublayer => bf_mainlayers(j)%get_head_sublayer()

           do i=1, bf_mainlayers(j)%get_nb_sublayers()
              
              !initialize the nodes of the buffer layer with a constant value
              call ini_cst_nodes(added_sublayer, (3*(j-1)+i)*scale)
              added_sublayer => added_sublayer%get_next()

           end do
        end do


        !test the fct update_grdpts_from_neighbors()
        do j=1,4
           
           added_sublayer => bf_mainlayers(j)%get_head_sublayer()

           do i=1, bf_mainlayers(j)%get_nb_sublayers()

              call nbf_interface_used%update_grdpts_from_neighbors(added_sublayer)

              added_sublayer => added_sublayer%get_next()
           end do

        end do


        !print the main layers after copy from neighbors
        do j=1,4
           call bf_mainlayers(j)%print_binary(
     $          'nodes2.dat',
     $          'grdpt_id2.dat',
     $          'sizes2.dat')

        end do


        !reinitialization of the nodes of the different layers
        scale = 0.07
        do j=1,4
           added_sublayer => bf_mainlayers(j)%get_head_sublayer()

           do i=1, bf_mainlayers(j)%get_nb_sublayers()
              
              !initialize the nodes of the buffer layer with a constant value
              call ini_cst_nodes(added_sublayer, (4*(i-1)+j)*scale)
              added_sublayer => added_sublayer%get_next()

           end do
        end do


        !test the fct update_neighbors_grdpts()
        do j=1,4
           
           added_sublayer => bf_mainlayers(j)%get_head_sublayer()

           do i=1, bf_mainlayers(j)%get_nb_sublayers()

              call nbf_interface_used%update_neighbor_grdpts(added_sublayer)

              added_sublayer => added_sublayer%get_next()
           end do

        end do

        !print the main layers after copy to neighbors
        do j=1,4
           call bf_mainlayers(j)%print_binary(
     $          'nodes3.dat',
     $          'grdpt_id3.dat',
     $          'sizes3.dat')

        end do

      end program test_nbf_interface_prog
