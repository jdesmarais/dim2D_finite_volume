      program test_bf_sublayers_merge_prog

        use bf_sublayer_class          , only : bf_sublayer

        use bf_sublayers_merge_module  , only : merge_sublayers_N,
     $                                          merge_sublayers_S,
     $                                          merge_sublayers_E,
     $                                          merge_sublayers_W
        use interface_abstract_class   , only : interface_abstract
        
        use parameters_constant        , only : N,S,E,W
        use parameters_input           , only : nx,ny,ne,bc_size
        use parameters_kind            , only : ikind, rkind
                                       
        use test_bf_layer_module       , only : print_nodes,
     $                                          print_grdpts_id,
     $                                          print_sizes,
     $                                          ini_grdpts_id,
     $                                          ini_nodes
        use test_cases_interface_module, only : ini_interface
        use interface_print_module     , only : print_interface
        use sublayers_ini_module       , only : ini_relative_sizes,
     $                                          ini_relative_distance,
     $                                          ini_alignment,
     $                                          ini_neighbors

        implicit none

        !parameters for the test
        integer, parameter :: no_alignment = 0
        integer, parameter :: test_case_id = 3
        integer, parameter :: mainlayer_id = 1
        integer, parameter :: nb_sublayers = 2
        integer, parameter :: size_case = 1
        integer, parameter :: distance_case = 1
        integer, parameter :: merge_id = 1
        integer, parameter :: merge_inverse = 0
        integer, parameter :: random_seed = 86456
        integer, parameter :: neighbor_case = 1

        integer, parameter :: interface_before = 0
        integer, parameter :: interface_after = 1


        !local variables
        real(rkind), dimension(nx,ny,ne)   :: nodes
        integer    , dimension(nx,ny)      :: grdpts_id
        type(interface_abstract)           :: interface_used
        integer, dimension(2,nb_sublayers) :: relative_sizes
        integer, dimension(nb_sublayers+1) :: relative_distance
        integer(ikind), dimension(2,2)     :: alignment
        logical                            :: neighbor_log1
        logical                            :: neighbor_log2

        integer :: i

        !variables for the test
        type(bf_sublayer), pointer :: sublayer1
        type(bf_sublayer), pointer :: sublayer2
        type(bf_sublayer), pointer :: merged_sublayer


        !verify the test parameters
        call verify_test_parameters()


        !initialize the nodes for the test
        call ini_nodes(nodes)
        call ini_grdpts_id(grdpts_id)
        call print_sizes(nodes,'interior_sizes.dat')
        call print_grdpts_id(grdpts_id,'interior_grdpts_id.dat')
        call print_nodes(nodes,'interior_nodes.dat')


        !initialize the relative_size and relative_distance
        call ini_relative_sizes(
     $       size_case,
     $       random_seed,
     $       relative_sizes)

        call ini_relative_distance(
     $       distance_case,
     $       random_seed,
     $       relative_distance)

        call ini_alignment(
     $       alignment,
     $       mainlayer_id,
     $       merge_id,
     $       relative_sizes,
     $       relative_distance)

        call ini_neighbors(
     $       neighbor_case,
     $       neighbor_log1, neighbor_log2)


        !initialize the interface with two sublayers
        !in the same mainlayer
        call ini_interface(
     $       interface_used,
     $       test_case_id,
     $       nodes,
     $       mainlayer_id=mainlayer_id,
     $       nb_sublayers=nb_sublayers,
     $       relative_distance=relative_distance,
     $       relative_sizes=relative_sizes)
        

        !select the two sublayers for the merge
        if(merge_inverse.eq.0) then
           sublayer1 => interface_used%mainlayer_pointers(mainlayer_id)%ptr%head_sublayer
           do i=1, merge_id-1
              sublayer1 => sublayer1%next
           end do
           sublayer2 => sublayer1%next
        else
           sublayer2 => interface_used%mainlayer_pointers(mainlayer_id)%ptr%head_sublayer
           do i=1, merge_id-1
              sublayer2 => sublayer2%next
           end do
           sublayer1 => sublayer2%next
        end if


        !print interface before
        call print_interface(interface_used, interface_before)


        !merge the sublayers
        select case(mainlayer_id)

          !north main layer
          case(N)
             if(no_alignment.eq.1) then
                merged_sublayer => merge_sublayers_N(
     $               interface_used%mainlayer_pointers(N)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               neighbor_E_i=neighbor_log2,
     $               neighbor_W_i=neighbor_log1)
             else
                merged_sublayer => merge_sublayers_N(
     $               interface_used%mainlayer_pointers(N)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               alignment=alignment,
     $               neighbor_E_i=neighbor_log2,
     $               neighbor_W_i=neighbor_log1)
             end if
             interface_used%mainlayer_pointers(N)%ptr%nb_sublayers = 
     $            interface_used%mainlayer_pointers(N)%ptr%nb_sublayers-1

          !south main layer
          case(S)
             if(no_alignment.eq.1) then
                merged_sublayer => merge_sublayers_S(
     $               interface_used%mainlayer_pointers(S)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               neighbor_E_i=neighbor_log2,
     $               neighbor_W_i=neighbor_log1)
             else
                merged_sublayer => merge_sublayers_S(
     $               interface_used%mainlayer_pointers(S)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               alignment=alignment,
     $               neighbor_E_i=neighbor_log2,
     $               neighbor_W_i=neighbor_log1)
             end if
             interface_used%mainlayer_pointers(S)%ptr%nb_sublayers = 
     $            interface_used%mainlayer_pointers(S)%ptr%nb_sublayers-1

          !east main layer
          case(E)
             if(no_alignment.eq.1) then
                merged_sublayer => merge_sublayers_E(
     $               interface_used%mainlayer_pointers(E)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               neighbor_N_i=neighbor_log2,
     $               neighbor_S_i=neighbor_log1)
             else
                merged_sublayer => merge_sublayers_E(
     $               interface_used%mainlayer_pointers(E)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               alignment=alignment,
     $               neighbor_N_i=neighbor_log2,
     $               neighbor_S_i=neighbor_log1)
             end if
             interface_used%mainlayer_pointers(E)%ptr%nb_sublayers = 
     $            interface_used%mainlayer_pointers(E)%ptr%nb_sublayers-1

          !west main layer
          case(W)
             if(no_alignment.eq.1) then
                merged_sublayer => merge_sublayers_W(
     $               interface_used%mainlayer_pointers(W)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               neighbor_N_i=neighbor_log2,
     $               neighbor_S_i=neighbor_log1)
             else
                merged_sublayer => merge_sublayers_W(
     $               interface_used%mainlayer_pointers(W)%ptr,
     $               sublayer1,
     $               sublayer2,
     $               nodes,
     $               alignment=alignment,
     $               neighbor_N_i=neighbor_log2,
     $               neighbor_S_i=neighbor_log1)
             end if
             interface_used%mainlayer_pointers(W)%ptr%nb_sublayers = 
     $            interface_used%mainlayer_pointers(W)%ptr%nb_sublayers-1

          case default
             print '(''mainlayer not recognized'')'
             print '(''mainlayer_id: '',I2)', mainlayer_id
             stop 'test not implemented yet'
        end select



        !print interface after
        call print_interface(interface_used, interface_after)


        contains

        subroutine verify_test_parameters()

          if((no_alignment.ne.0).and.(no_alignment.ne.1)) then
             print '(''no_alignment should equal 0 or 1'')'
             print '(''it determines whether the alignment should'')'
             print '(''be taken into account'')'
             stop 'change no_alignment'
          end if
          
          if(test_case_id.ne.3) then
             print '(''only one initialization is implemented'')'
             print '(''for this test'')'
             print '(''please select test_case_id = 3'')'
             stop 'change test_case_id'
          end if

          if(  (mainlayer_id.ne.N).and.
     $         (mainlayer_id.ne.S).and.
     $         (mainlayer_id.ne.E).and.
     $         (mainlayer_id.ne.W)) then
             print '(''mainlayer_id not recognized'')'
             print '(''please select N, S, E, or W'')'
             stop 'change mainlayer_id'
          end if

          if((size_case.le.0).or.(size_case.ge.5)) then
             print '(''size_case not recognized'')'
             print '(''please select between 1 and 4'')'
             stop 'change size_case'
          end if

          if((distance_case.le.0).or.(distance_case.ge.5)) then
             print '(''distance_case not recognized'')'
             print '(''please select between 1 and 4'')'
             stop 'change distance_case'
          end if

          if(merge_id.gt.(nb_sublayers-1)) then
             print '(''merge_id not recognized'')'
             print '(''please select merge_id<nb_sublayers'')'
             stop 'change merge_id'
          end if

          if((merge_inverse.ne.0).and.(merge_inverse.ne.1)) then
             print '(''merge_inverse not recognized'')'
             print '(''please select 0 or 1'')'
             stop 'change merge_inverse'
          end if

          if((neighbor_case.le.0).or.(neighbor_case.ge.5)) then
             print '(''neighbor_case not recognized'')'
             print '(''please select between 1 and 4'')'
             stop 'change neighbor_case'
          end if

        end subroutine verify_test_parameters

      end program test_bf_sublayers_merge_prog
