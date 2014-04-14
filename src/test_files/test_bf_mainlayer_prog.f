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
        integer, parameter :: nb_sublayers = 3
        type(bf_sublayer), pointer :: new_sublayer

        integer, dimension(nb_sublayers,2,2) :: test_alignment
        logical, dimension(nb_sublayers,4)   :: test_neighbors
        integer, dimension(2,2)              :: alignment
        real(rkind), dimension(nx,ny,ne)     :: nodes
        logical, dimension(4)                :: neighbors
        integer :: i


        !test initialize the mainlayer
        print '()'
        print '(''test ini() : initialization of bf_mainlayer'')'
        print '(''-------------------------------------------'')'
        call bf_mainlayer_N_tested%ini(N)
        print '(''nb_sublayers     : '',I1)', bf_mainlayer_N_tested%nb_sublayers
        print '(''head associated? : '',L1)', associated(bf_mainlayer_N_tested%head_sublayer)
        print '(''tail associated? : '',L1)', associated(bf_mainlayer_N_tested%tail_sublayer)
        print '()'


        !initialize the alignment required for the buffer layers
        call ini_alignment_table(test_alignment)
        call ini_neighbors_table(test_neighbors)


        !create sublayers for the North main layer
        print '()'
        print '(''test add_sublayer():add sublayers to mainlayer'')'
        print '(''----------------------------------------------'')'
        do i=1, nb_sublayers

           alignment(1,1) = test_alignment(i,1,1)
           alignment(1,2) = test_alignment(i,1,2)
           alignment(2,1) = test_alignment(i,2,1)
           alignment(2,2) = test_alignment(i,2,2)

           new_sublayer => bf_mainlayer_N_tested%add_sublayer(
     $          alignment)

           call new_sublayer%element%ini(N)
           call new_sublayer%element%allocate_bf_layer(
     $          alignment,
     $          nodes,
     $          neighbors)

        end do

        call print_mainlayer(bf_mainlayer_N_tested)


        !create sublayers for the East main layer
        call bf_mainlayer_E_tested%ini(E)
        print '()'
        print '(''test add_sublayer():add sublayers to mainlayer'')'
        print '(''----------------------------------------------'')'
        do i=1, nb_sublayers

           alignment(1,1) = test_alignment(i,1,1)
           alignment(1,2) = test_alignment(i,1,2)
           alignment(2,1) = test_alignment(i,2,1)
           alignment(2,2) = test_alignment(i,2,2)

           new_sublayer => bf_mainlayer_E_tested%add_sublayer(
     $          alignment)

           call new_sublayer%element%ini(E)
           call new_sublayer%element%allocate_bf_layer(
     $          alignment,
     $          nodes,
     $          neighbors)

        end do

        call print_mainlayer(bf_mainlayer_E_tested)

      end program test_bf_mainlayer_prog
