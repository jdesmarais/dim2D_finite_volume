      !the idea with this test file is to recreate the
      !same results as presented in test_bf_layer_prog
      !but using the interface_abstract as the object
      !encapsulating the buffer layer objects      
      program test_bf_interface_prog

        use bf_sublayer_class   , only : bf_sublayer
        use bf_interface_class  , only : bf_interface
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

        type(bf_interface)                  :: interface_tested
        real(rkind)   , dimension(nx,ny,ne) :: nodes
        integer       , dimension(nx,ny)    :: grdpts_id
        integer       , dimension(2,2)      :: alignment
        integer                             :: mainlayer_id
        type(bf_sublayer), pointer          :: added_sublayer
        type(bf_sublayer), pointer          :: bf_sublayer_merged1
        type(bf_sublayer), pointer          :: bf_sublayer_merged2
        type(bf_sublayer), pointer          :: bf_sublayer_reallocated
        real(rkind)                         :: scale


        !initialize the nodes and print them
        call ini_nodes(nodes)
        call ini_grdpts_id(grdpts_id)
        call print_interior_data(nodes,
     $                           grdpts_id, 
     $                           'interior_nodes.dat',
     $                           'interior_grdpts_id.dat',
     $                           'interior_sizes.dat')
        
        !initialize the interface
        call interface_tested%ini()

        
        !add the N_E buffer layer
        call get_alignment(1, alignment, mainlayer_id)
        added_sublayer => interface_tested%allocate_sublayer(
     $       mainlayer_id, nodes, alignment)
        scale = 0.1
        call ini_cst_nodes(added_sublayer, scale)
        
        !add the W buffer layer under the N_W corner
        call get_alignment(4, alignment, mainlayer_id)
        added_sublayer => interface_tested%allocate_sublayer(
     $       mainlayer_id, nodes, alignment)
        scale = 0.2
        call ini_cst_nodes(added_sublayer, scale)


        !print the interface
        call interface_tested%print_binary(
     $       'nodes1.dat',
     $       'grdpt_id1.dat',
     $       'sizes1.dat',
     $       '1.dat')

        
        !add the E buffer layer under the N_E corner
        call get_alignment(3, alignment, mainlayer_id)
        bf_sublayer_merged1 => interface_tested%allocate_sublayer(
     $       mainlayer_id, nodes, alignment)
        
        !add the N_W buffer layer
        call get_alignment(2, alignment, mainlayer_id)
        added_sublayer => interface_tested%allocate_sublayer(
     $       mainlayer_id, nodes, alignment)

        
        !print the interface
        call interface_tested%print_binary(
     $       'nodes2.dat',
     $       'grdpt_id2.dat',
     $       'sizes2.dat',
     $       '2.dat')


        !add the E buffer layer
        call get_alignment(5, alignment, mainlayer_id)
        bf_sublayer_merged2 => interface_tested%allocate_sublayer(
     $       mainlayer_id, nodes, alignment)

        !add the W buffer layer
        call get_alignment(6, alignment, mainlayer_id)
        bf_sublayer_reallocated => interface_tested%allocate_sublayer(
     $       mainlayer_id, nodes, alignment)

        
        !print the interface
        call interface_tested%print_binary(
     $       'nodes3.dat',
     $       'grdpt_id3.dat',
     $       'sizes3.dat',
     $       '3.dat')


        !add the S_E buffer layer
        call get_alignment(7, alignment, mainlayer_id)
        added_sublayer => interface_tested%allocate_sublayer(
     $       mainlayer_id, nodes, alignment)
        scale = 0.3
        call ini_cst_nodes(added_sublayer, scale)
        

        !add the S_W buffer layer
        call get_alignment(8, alignment, mainlayer_id)
        added_sublayer => interface_tested%allocate_sublayer(
     $       mainlayer_id, nodes, alignment)
        scale = 0.4
        call ini_cst_nodes(added_sublayer, scale)


        !print the interface
        call interface_tested%print_binary(
     $       'nodes4.dat',
     $       'grdpt_id4.dat',
     $       'sizes4.dat',
     $       '4.dat')


        !merge the E buffer layers
        call get_alignment(9, alignment, mainlayer_id)
        added_sublayer => interface_tested%merge_sublayers(
     $       bf_sublayer_merged1, bf_sublayer_merged2,
     $       nodes, alignment)
        

        !reallocate the W buffer layer
        call get_alignment(10, alignment, mainlayer_id)
        call interface_tested%reallocate_sublayer(
     $       bf_sublayer_reallocated,
     $       nodes, alignment)


        !print the interface
        call interface_tested%print_binary(
     $       'nodes5.dat',
     $       'grdpt_id5.dat',
     $       'sizes5.dat',
     $       '5.dat')

        contains

        subroutine get_alignment(
     $       index, alignment, mainlayer_id)

          implicit none

          integer                       , intent(in)  :: index
          integer(ikind), dimension(2,2), intent(out) :: alignment
          integer                       , intent(out) :: mainlayer_id


          integer :: small_layer
          integer :: large_layer

          small_layer = 3
          large_layer = 7

          select case(index)

            !N_E layer at the corner
            case(1)
               mainlayer_id   = N
               alignment(1,1) = align_E-small_layer
               alignment(1,2) = align_E+small_layer
               alignment(2,1) = align_N
               alignment(2,2) = align_N

            !N_W layer at the corner
            case(2)
               mainlayer_id   = N
               alignment(1,1) = align_W-small_layer
               alignment(1,2) = align_W+small_layer
               alignment(2,1) = align_N
               alignment(2,2) = align_N

            !E layer just below the N_E corner
            case(3)
               mainlayer_id   = E
               alignment(1,1) = align_E
               alignment(1,2) = align_E
               alignment(2,1) = align_N-1-small_layer
               alignment(2,2) = align_N-1

            !W layer just below the N_E corner
            case(4)
               mainlayer_id   = W
               alignment(1,1) = align_W
               alignment(1,2) = align_W
               alignment(2,1) = align_N-bc_size-small_layer
               alignment(2,2) = align_N-bc_size

            !S_E layer above the corner
            case(5)
               mainlayer_id   = E
               alignment(1,1) = align_E
               alignment(1,2) = align_E
               alignment(2,1) = align_S+large_layer
               alignment(2,2) = align_S+large_layer+small_layer

            !S_W layer above the corner
            case(6)
               mainlayer_id   = W
               alignment(1,1) = align_W
               alignment(1,2) = align_W
               alignment(2,1) = align_S+large_layer
               alignment(2,2) = align_S+large_layer+small_layer

            !S_W layer
            case(7)
               mainlayer_id   = S
               alignment(1,1) = align_E-small_layer
               alignment(1,2) = align_E+small_layer
               alignment(2,1) = align_S
               alignment(2,2) = align_S

            !S_E layer
            case(8)
               mainlayer_id   = S
               alignment(1,1) = align_W-small_layer
               alignment(1,2) = align_W+small_layer
               alignment(2,1) = align_S
               alignment(2,2) = align_S

            !east merge
            case(9)
               mainlayer_id   = E
               alignment(1,1) = align_E
               alignment(1,2) = align_E
               alignment(2,1) = align_S+bc_size
               alignment(2,2) = align_N-1
               
            !west reallocation
            case(10)
               mainlayer_id   = W
               alignment(1,1) = align_W
               alignment(1,2) = align_W
               alignment(2,1) = align_S+1
               alignment(2,2) = align_S+large_layer+small_layer

            case default
               print '(''case not implemented'')'
               stop ''
          end select               

        end subroutine get_alignment

      end program test_bf_interface_prog
