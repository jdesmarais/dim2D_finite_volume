      program test_bf_mainlayer_prog
      
        use bf_mainlayer_class  , only : bf_mainlayer
        use bf_sublayer_class   , only : bf_sublayer
        use parameters_constant , only : N, E
        use parameters_input    , only : nx,ny,ne
        use parameters_kind     , only : rkind
        use test_bf_layer_module, only : ini_nodes,
     $                                   ini_alignment_table,
     $                                   ini_neighbors_table
        use test_bf_mainlayer_module, only : print_mainlayer

        implicit none

        type(bf_mainlayer) :: bf_mainlayer_N_tested
        type(bf_mainlayer) :: bf_mainlayer_E_tested
        integer, parameter :: nb_sublayers = 10
        type(bf_sublayer), pointer :: new_sublayer

        integer, dimension(nb_sublayers,2,2) :: test_alignment
        logical, dimension(nb_sublayers,4)   :: test_neighbors
        integer, dimension(2,2)              :: alignment
        real(rkind), dimension(nx,ny,ne)     :: nodes
        logical, dimension(4)                :: neighbors
        logical, parameter                   :: random_test = .true.
        integer                              :: random_seed

        integer :: i


        if(random_test) then
           random_seed = 8945423
           call srand(random_seed)
        end if


        !test initialize the mainlayer
        print '()'
        print '(''test ini() : initialization of bf_mainlayer'')'
        print '(''-------------------------------------------'')'
        call bf_mainlayer_N_tested%ini(N)
        print '(''nb_sublayers     : '',I1)', bf_mainlayer_N_tested%get_nb_sublayers()
        print '(''head associated? : '',L1)', associated(bf_mainlayer_N_tested%get_head_sublayer())
        print '(''tail associated? : '',L1)', associated(bf_mainlayer_N_tested%get_tail_sublayer())
        print '()'


        !initialize the alignment required for the buffer layers
        call ini_alignment_table(test_alignment)
        call ini_neighbors_table(test_neighbors)


        !create sublayers for the North main layer
        print '()'
        print '(''test add_sublayer():add sublayers to mainlayer'')'
        print '(''----------------------------------------------'')'
        do i=1, nb_sublayers

           if(.not.random_test) then
              alignment(1,1) = test_alignment(i,1,1)
              alignment(1,2) = test_alignment(i,1,2)
              alignment(2,1) = test_alignment(i,2,1)
              alignment(2,2) = test_alignment(i,2,2)
           else
              alignment(1,1) = nint(nx*RAND())
              alignment(1,2) = alignment(1,1)+(i-1)
              alignment(2,1) = ny+1
              alignment(2,2) = ny+1+(i-1)
           end if

           new_sublayer => bf_mainlayer_N_tested%add_sublayer(
     $          nodes, alignment, neighbors)

        end do

        call print_mainlayer(bf_mainlayer_N_tested)


        !create sublayers for the East main layer
        call bf_mainlayer_E_tested%ini(E)
        print '()'
        print '(''test add_sublayer():add sublayers to mainlayer'')'
        print '(''----------------------------------------------'')'
        do i=1, nb_sublayers

           if(.not.random_test) then
              alignment(1,1) = test_alignment(i,1,1)
              alignment(1,2) = test_alignment(i,1,2)
              alignment(2,1) = test_alignment(i,2,1)
              alignment(2,2) = test_alignment(i,2,2)
           else
              alignment(1,1) = nx+1
              alignment(1,2) = nx+1
              alignment(2,1) = nint(nx*RAND())
              alignment(2,2) = alignment(2,1)+(i-1)
           end if

           new_sublayer => bf_mainlayer_E_tested%add_sublayer(
     $          nodes, alignment, neighbors)

        end do

        call print_mainlayer(bf_mainlayer_E_tested)

      end program test_bf_mainlayer_prog
