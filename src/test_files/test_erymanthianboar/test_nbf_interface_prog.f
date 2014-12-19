      program test_nbf_interface_prog

        use nbf_interface_class, only :
     $       nbf_interface

        use bf_mainlayer_class, only :
     $       bf_mainlayer

        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_constant, only :
     $       N,S,E,W,
     $       left,right

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size

        use parameters_kind, only :
     $       ikind, rkind

        use test_bf_layer_module, only :
     $       print_interior_data_wo_maps,
     $       ini_nodes,
     $       ini_grdpts_id,
     $       ini_cst_nodes

        use test_nbf_list_module, only :
     $       ini_alignment
        
        implicit none

c$$$        real(rkind)        , dimension(nx)       :: x_map
c$$$        real(rkind)        , dimension(ny)       :: y_map
c$$$        real(rkind)        , dimension(nx,ny,ne) :: nodes
c$$$        integer            , dimension(nx,ny)    :: grdpts_id
c$$$        integer            , dimension(4)        :: bf_layer_loc
c$$$        character(2)       , dimension(4)        :: bf_layer_char
c$$$                           
c$$$        type(bf_mainlayer) , dimension(4)        :: bf_mainlayers
c$$$        type(nbf_interface)                      :: nbf_interface_used
c$$$        type(bf_sublayer)  , pointer             :: added_sublayer
c$$$        !type(bf_sublayer)  , pointer             :: added_sublayer2
c$$$                           
c$$$        integer(ikind)     , dimension(4,5,2,2)  :: test_alignment
c$$$        integer(ikind)     , dimension(2,2)      :: alignment
c$$$        integer            , dimension(4,2)      :: neighbors
c$$$        real(rkind)                              :: scale
c$$$
c$$$        integer :: g,h,i,j,i_max
c$$$
        logical :: detailled
        logical :: test_loc
        logical :: test_validated        
c$$$
c$$$
c$$$        !initialize the nodes and the grdpts_id
c$$$        call ini_nodes(nodes)
c$$$        call ini_grdpts_id(grdpts_id)
c$$$
c$$$        !print the nodes
c$$$        call print_interior_data_wo_maps(
c$$$     $       nodes, grdpts_id,
c$$$     $       'interior_nodes.dat',
c$$$     $       'interior_grdpts_id.dat',
c$$$     $       'interior_sizes.dat')
c$$$
c$$$        !buffer layers tested
c$$$        bf_layer_loc  = [N,S,E,W]
c$$$        bf_layer_char = ['N_','S_','E_','W_']
c$$$        
c$$$        neighbors(N,1) = W
c$$$        neighbors(N,2) = E
c$$$        neighbors(S,1) = W
c$$$        neighbors(S,2) = E
c$$$        neighbors(W,1) = S
c$$$        neighbors(W,2) = N
c$$$        neighbors(E,1) = S
c$$$        neighbors(E,2) = N
c$$$
c$$$
c$$$        !initialize the alignment for all the buffer layers
c$$$        call ini_alignment(test_alignment)
c$$$
c$$$        !initialize the main layers
c$$$        do j=1,4
c$$$           call bf_mainlayers(j)%ini(j)           
c$$$        end do
c$$$        
c$$$        !initialize the nbf_interface
c$$$        call nbf_interface_used%ini()
c$$$
c$$$
c$$$        !add the buffer layers to the mainlayers
c$$$        !+ add the links to the nbf_lists
c$$$        !loop over the cardinal coordinates
c$$$        print '()'
c$$$        print '(''buffer layer alignments'')'
c$$$        print '(''-----------------------'')'
c$$$        do j=1,4
c$$$
c$$$           select case(j)
c$$$             case(N,S)
c$$$                i_max = 5
c$$$             case(E,W)
c$$$                i_max = 3
c$$$           end select
c$$$
c$$$           !loop over the number of buffer layers
c$$$           !per main layer
c$$$           do i=1,i_max
c$$$
c$$$              !alignment for the buffer layer
c$$$              do g=1,2
c$$$                 do h=1,2
c$$$                    alignment(h,g) = test_alignment(j,i,h,g)
c$$$                 end do
c$$$              end do
c$$$
c$$$              !print the alignment
c$$$              print '(A2,I1,'' ('',2I3,'')  ('',2I3,'')'')',
c$$$     $             bf_layer_char(j), i,
c$$$     $             alignment(1,1), alignment(1,2),
c$$$     $             alignment(2,1), alignment(2,2)
c$$$
c$$$              !allocate the buffer layer
c$$$              added_sublayer => bf_mainlayers(j)%add_sublayer(
c$$$     $             x_map,y_map,nodes,alignment)
c$$$
c$$$              !set whether the buffe rlayer can exchange
c$$$              !with neighboring layers
c$$$              call added_sublayer%set_neighbor1_share()
c$$$              call added_sublayer%set_neighbor2_share()
c$$$
c$$$              !if the buffer layer is a potential neigboring buffer
c$$$              !layer, we add the link to the nbf_interface
c$$$              if(added_sublayer%can_exchange_with_neighbor1()) then
c$$$                 call nbf_interface_used%link_neighbor1_to_bf_sublayer(
c$$$     $                added_sublayer)
c$$$              end if
c$$$
c$$$              if(added_sublayer%can_exchange_with_neighbor2()) then
c$$$                 call nbf_interface_used%link_neighbor2_to_bf_sublayer(
c$$$     $                added_sublayer)
c$$$              end if
c$$$
c$$$           end do
c$$$        end do
c$$$        print '(''-----------------------'')'
c$$$        print '()'
c$$$
c$$$
c$$$        !print the main layers before the tests
c$$$        do j=1,4
c$$$           call bf_mainlayers(j)%print_binary(
c$$$     $          'xmap1.dat',
c$$$     $          'ymap1.dat',
c$$$     $          'nodes1.dat',
c$$$     $          'grdpt_id1.dat',
c$$$     $          'sizes1.dat')
c$$$
c$$$        end do
c$$$
c$$$        
c$$$        !print the content of the nbf_interface
c$$$        print '()'
c$$$        print '(''nbf_interface content'')'
c$$$        print '(''-----------------------'')'
c$$$        call nbf_interface_used%print_on_screen()
c$$$        print '(''-----------------------'')'
c$$$        print '()'
c$$$
c$$$        
c$$$        !initialization of the nodes of the different layers
c$$$        scale = 0.07
c$$$        do j=1,4
c$$$           added_sublayer => bf_mainlayers(j)%get_head_sublayer()
c$$$
c$$$           do i=1, bf_mainlayers(j)%get_nb_sublayers()
c$$$              
c$$$              !initialize the nodes of the buffer layer with a constant value
c$$$              call ini_cst_nodes(added_sublayer, (3*(j-1)+i)*scale)
c$$$              added_sublayer => added_sublayer%get_next()
c$$$
c$$$           end do
c$$$        end do
c$$$
c$$$
c$$$        !test the fct update_grdpts_from_neighbors()
c$$$        do j=1,4
c$$$           
c$$$           added_sublayer => bf_mainlayers(j)%get_head_sublayer()
c$$$
c$$$           do i=1, bf_mainlayers(j)%get_nb_sublayers()
c$$$
c$$$              call nbf_interface_used%update_grdpts_from_neighbors(added_sublayer)
c$$$
c$$$              added_sublayer => added_sublayer%get_next()
c$$$           end do
c$$$
c$$$        end do
c$$$
c$$$
c$$$        !print the main layers after copy from neighbors
c$$$        do j=1,4
c$$$           call bf_mainlayers(j)%print_binary(
c$$$     $          'xmap2.dat',
c$$$     $          'ymap2.dat',
c$$$     $          'nodes2.dat',
c$$$     $          'grdpt_id2.dat',
c$$$     $          'sizes2.dat')
c$$$
c$$$        end do
c$$$
c$$$
c$$$        !reinitialization of the nodes of the different layers
c$$$        scale = 0.07
c$$$        do j=1,4
c$$$           added_sublayer => bf_mainlayers(j)%get_head_sublayer()
c$$$
c$$$           do i=1, bf_mainlayers(j)%get_nb_sublayers()
c$$$              
c$$$              !initialize the nodes of the buffer layer with a constant value
c$$$              call ini_cst_nodes(added_sublayer, (4*(i-1)+j)*scale)
c$$$              added_sublayer => added_sublayer%get_next()
c$$$
c$$$           end do
c$$$        end do
c$$$
c$$$
c$$$        !test the fct update_neighbors_grdpts()
c$$$        do j=1,4
c$$$           
c$$$           added_sublayer => bf_mainlayers(j)%get_head_sublayer()
c$$$
c$$$           do i=1, bf_mainlayers(j)%get_nb_sublayers()
c$$$
c$$$              call nbf_interface_used%update_neighbor_grdpts(added_sublayer)
c$$$
c$$$              added_sublayer => added_sublayer%get_next()
c$$$           end do
c$$$
c$$$        end do
c$$$
c$$$        !print the main layers after copy to neighbors
c$$$        do j=1,4
c$$$           call bf_mainlayers(j)%print_binary(
c$$$     $          'xmap3.dat',
c$$$     $          'ymap3.dat',
c$$$     $          'nodes3.dat',
c$$$     $          'grdpt_id3.dat',
c$$$     $          'sizes3.dat')
c$$$
c$$$        end do

        
        !test: ask_neighbor_for_bc_overlap
        detailled      = .true.
        test_loc       = test_ask_neighbors_for_bc_overlap(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_ask_neighbors_for_bc_overlap: '',L1)', test_loc
        print '()'

        contains


        function test_ask_neighbors_for_bc_overlap(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(nbf_interface)        :: nbf_interface_used
          type(bf_sublayer), pointer :: bf_layer_N_ptr
          type(bf_sublayer), pointer :: bf_layer_E_ptr
          
          integer       , dimension(:,:), allocatable :: grdpts_id
          integer(ikind), dimension(2,2)              :: bf_alignment
          integer(ikind)                              :: x_border_data
          integer(ikind)                              :: x_border_comp
          logical                                     :: err

          !test case:
          !2 buffer layers: N and E at the interface NE
          !the North buffer layer should be increased to
          !be able to compute all the boundary grid points
          !        
          !      North buffer layer
          !
          !                  ___ nx-3
          !                 |  ___ nx-2
          !                 | |  ___ nx-1
          !                 | | |  __  nx
          !                 | | | |
          !        ________________ 
          !         3 3 3 3 3 3 3 3|___
          ! ny   -  2 2 2 2 2 2 2 3|3  |
          ! ny-1 - _ _ _ _ _ _ _2_2|3_3|
          ! ny-2 -         |   |  2|2 3|
          ! ny-3 - ________|___|__ |2 3| East buffer layer
          !                |   |   |2 3|
          !       interior |1 1|2 2|2 3|
          !                |1 1|2 3|3 3|
          !------------------------------------------------
          !initialize the nbf_interface
          call nbf_interface_used%ini()


          !initialize the N buffer layer
          allocate(grdpts_id(10,5))

          grdpts_id = reshape((/
     $         1,1,1,1,1,1,1,1,1,1,
     $         1,1,1,1,1,1,1,1,1,2,
     $         1,1,1,1,1,1,1,1,2,2,
     $         2,2,2,2,2,2,2,2,2,3,
     $         3,3,3,3,3,3,3,3,3,3/),
     $         (/10,5/))

          bf_alignment(1,2) = nx-2
          bf_alignment(1,1) = bf_alignment(1,2) - size(grdpts_id,1) + (2*bc_size+1)
          bf_alignment(2,1) = ny-1
          bf_alignment(2,2) = bf_alignment(2,1) + size(grdpts_id,2) - (2*bc_size+1)

          
          allocate(bf_layer_N_ptr)
          call bf_layer_N_ptr%ini(N)
          call bf_layer_N_ptr%set_grdpts_id(grdpts_id)
          call bf_layer_N_ptr%set_alignment_tab(bf_alignment)

          call nbf_interface_used%link_neighbor2_to_bf_sublayer(bf_layer_N_ptr)


          !initialize the E buffer layer
          allocate(grdpts_id(6,7))

          grdpts_id = reshape((/
     $         1,1,2,3,3,3,
     $         1,1,2,2,2,3,
     $         1,1,1,1,2,3,
     $         1,1,1,1,2,3,
     $         1,1,1,2,2,3,
     $         1,1,2,2,3,3,
     $         2,2,2,3,3,0/),
     $         (/6,7/))

          bf_alignment(1,1) = nx-1
          bf_alignment(1,2) = bf_alignment(1,1) + size(grdpts_id,1) - (2*bc_size+1)
          bf_alignment(2,2) = ny-2
          bf_alignment(2,1) = bf_alignment(2,2) - size(grdpts_id,2) + (2*bc_size+1)

          allocate(bf_layer_E_ptr)
          call bf_layer_E_ptr%ini(E)
          call bf_layer_E_ptr%set_grdpts_id(grdpts_id)
          call bf_layer_E_ptr%set_alignment_tab(bf_alignment)

          call nbf_interface_used%link_neighbor2_to_bf_sublayer(bf_layer_E_ptr)


          !test: ask_neighbors_for_bc_overlap
          x_border_data = nx+1
          x_border_comp = nbf_interface_used%ask_neighbors_for_bc_overlap(
     $         N,
     $         2,
     $         [nx,ny-1],
     $         right,
     $         err)

          test_validated = x_border_data.eq.x_border_comp
          if(detailled.and.(.not.test_validated)) then
             print '(''x_border: '',I3,'' -> '',I3)', x_border_comp, x_border_data
          end if

        end function test_ask_neighbors_for_bc_overlap

      end program test_nbf_interface_prog
