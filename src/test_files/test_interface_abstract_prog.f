      !the idea with this test file is to recreate the
      !same results as presented in test_bf_layer_prog
      !but using the interface_abstract as the object
      !encapsulating the buffer layer objects      
      program test_interface_abstract_prog

        use bf_layer_class          , only : bf_layer
        use interface_abstract_class, only : interface_abstract
        use parameters_constant     , only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input        , only : nx,ny,ne, bc_size
        use parameters_kind         , only : rkind

        use test_bf_layer_module    , only : print_nodes,
     $                                       print_sizes,
     $                                       bf_layer_test_allocation,
     $                                       bf_layer_test_reallocation,
     $                                       bf_layer_test_update_grdpts,
     $                                       ini_nodes
        implicit none

        integer, parameter :: nb_sublayers = 3
        
        type(interface_abstract)             :: interface_tested
        real(rkind), dimension(nx,ny,ne)     :: nodes
        integer, dimension(nb_sublayers,2,2) :: test_alignment
        integer, dimension(2,2)              :: alignment
        logical, dimension(4)                :: neighbors
        character(len=21)                    :: sizes_filename
        character(len=21)                    :: nodes_filename
        character(len=21)                    :: grdid_filename

        integer     , dimension(8) :: bf_layer_loc
        character(2), dimension(8) :: bf_layer_char
        integer                    :: i,j
        type(bf_layer), pointer    :: bf_layer_tested


        !initialize the nodes and print them
        call ini_nodes(nodes)
        call print_sizes(nodes,'interior_sizes.dat')
        call print_nodes(nodes,'interior_nodes.dat')
        
        !buffer layer tested
        bf_layer_loc  = [N,S,E,W,N_E,N_W,S_E,S_W]
        bf_layer_char = ['N_','S_','E_','W_','NE','NW','SE','SW']

        !alignment and neighbors
        test_alignment(1,1,1) = bc_size+1
        test_alignment(1,1,2) = bc_size+2
        test_alignment(1,2,1) = bc_size+1
        test_alignment(1,2,2) = bc_size+2

        test_alignment(2,1,1) = bc_size+9
        test_alignment(2,1,2) = bc_size+10
        test_alignment(2,2,1) = bc_size+9
        test_alignment(2,2,2) = bc_size+10

        test_alignment(3,1,1) = bc_size+16
        test_alignment(3,1,2) = bc_size+16
        test_alignment(3,2,1) = bc_size+16
        test_alignment(3,2,2) = bc_size+16

        neighbors = [.false.,.false.,.true.,.true.]


        !create all the main layers and one sublayer
        !per main layer
        call print_nb_sublayers('sublayers_nb.dat')
        do i=1, size(bf_layer_loc,1)

           !allocate space for the sublayer
           call interface_tested%allocate_bf_mainlayer(bf_layer_loc(i),nb_sublayers)

           do j=1, nb_sublayers

              !ini for the alignment
              alignment(1,1) = test_alignment(j,1,1)
              alignment(1,2) = test_alignment(j,1,2)
              alignment(2,1) = test_alignment(j,2,1)
              alignment(2,2) = test_alignment(j,2,2)

              !allocate the sublayer tables: nodes, gridpt_id...
              bf_layer_tested => interface_tested%get_bf_layer(bf_layer_loc(i),j)
              call bf_layer_tested%ini([bf_layer_loc(i),j])
              call bf_layer_tested%allocate_bf_layer(
     $             alignment, nodes, neighbors)
           	
              !print the allocated tables on files
              write(sizes_filename,'(A2,I1,''_sizes.dat'')') bf_layer_char(i),j
              write(nodes_filename,'(A2,I1,''_nodes.dat'')') bf_layer_char(i), j
              write(grdid_filename,'(A2,I1,''_grdpt_id.dat'')') bf_layer_char(i), j
              call bf_layer_tested%print_sizes(sizes_filename)
              call bf_layer_tested%print_nodes(nodes_filename)
              call bf_layer_tested%print_grdpts_id(grdid_filename)

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
