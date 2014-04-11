      !the idea with this test file is to recreate the
      !same results as presented in test_bf_layer_prog
      !but using the interface_abstract as the object
      !encapsulating the buffer layer objects      
      program test_interface_abstract_prog

        use bf_sublayer_class       , only : bf_sublayer
        use interface_abstract_class, only : interface_abstract
        use parameters_constant     , only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input        , only : nx,ny,ne
        use parameters_kind         , only : rkind

        use test_bf_layer_module    , only : print_nodes,
     $                                       print_sizes,
     $                                       bf_layer_test_allocation,
     $                                       bf_layer_test_reallocation,
     $                                       bf_layer_test_update_grdpts,
     $                                       ini_nodes,
     $                                       ini_alignment_table,
     $                                       ini_neighbors_table

        implicit none

        integer, parameter :: nb_sublayers = 3
        
        type(interface_abstract)             :: interface_tested

        integer, dimension(nb_sublayers,2,2) :: test_alignment
        logical, dimension(nb_sublayers,4)   :: test_neighbors

        real(rkind), dimension(nx,ny,ne)     :: nodes
        integer, dimension(2,2)              :: alignment
        logical, dimension(4)                :: neighbors

        character(len=21)                    :: sizes_filename
        character(len=21)                    :: nodes_filename
        character(len=21)                    :: grdid_filename

        integer     , dimension(8) :: bf_layer_loc
        character(2), dimension(8) :: bf_layer_char
        integer                    :: i,j
        type(bf_sublayer), pointer :: added_sublayer


        !initialize the nodes and print them
        call ini_nodes(nodes)
        call print_sizes(nodes,'interior_sizes.dat')
        call print_nodes(nodes,'interior_nodes.dat')
        
        !buffer layer tested
        bf_layer_loc  = [N,S,E,W,N_E,N_W,S_E,S_W]
        bf_layer_char = ['N_','S_','E_','W_','NE','NW','SE','SW']

        !alignment and neighbors
        call ini_alignment_table(test_alignment)
        call ini_neighbors_table(test_neighbors)

        !create all the main layers
        !per main layer
        call print_nb_sublayers('sublayers_nb.dat')

        call interface_tested%ini()
        do i=1,4

           do j=1, nb_sublayers

              !initialize the alignment
              alignment(1,1) = test_alignment(j,1,1)
              alignment(1,2) = test_alignment(j,1,2)
              alignment(2,1) = test_alignment(j,2,1)
              alignment(2,2) = test_alignment(j,2,2)

              !initialize the neighbors
              neighbors(1) = test_neighbors(j,1)
              neighbors(2) = test_neighbors(j,2)
              neighbors(3) = test_neighbors(j,3)
              neighbors(4) = test_neighbors(j,4)
              
              !add the sublayer with this alignement and neighbors
              !and allocate the buffer layer using the nodes
              added_sublayer => interface_tested%add_sublayer(
     $             i, alignment, nodes, neighbors)
           	
              !print the allocated tables on files
              write(sizes_filename,'(A2,I1,''_sizes.dat'')') bf_layer_char(i),j
              write(nodes_filename,'(A2,I1,''_nodes.dat'')') bf_layer_char(i),j
              write(grdid_filename,'(A2,I1,''_grdpt_id.dat'')') bf_layer_char(i),j
              call added_sublayer%element%print_sizes(sizes_filename)
              call added_sublayer%element%print_nodes(nodes_filename)
              call added_sublayer%element%print_grdpts_id(grdid_filename)

           end do

        end do

        
        contains


        subroutine print_nb_sublayers(filename)

          implicit none

          character(*), intent(in) :: filename

          integer :: ios
          
          open(unit=1,
     $          file=filename,
     $          action="write", 
     $          status="unknown",
     $          form='unformatted',
     $          access='sequential',
     $          position='rewind',
     $          iostat=ios)

           if(ios.eq.0) then
              write(unit=1, iostat=ios) nb_sublayers
              close(unit=1)
           else
              stop 'file opening pb'
           end if

        end subroutine print_nb_sublayers
        
      end program test_interface_abstract_prog
