      !the idea with this test file is to recreate the
      !same results as presented in test_bf_layer_prog
      !but using the interface_abstract as the object
      !encapsulating the buffer layer objects      
      program test_interface_abstract_prog

        use bf_sublayer_class       , only : bf_sublayer
        use bf_mainlayer_class      , only : bf_mainlayer
        use interface_abstract_class, only : interface_abstract
        use parameters_constant     , only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input        , only : nx,ny,ne
        use parameters_kind         , only : ikind, rkind

        use test_bf_layer_module    , only : print_nodes,
     $                                       print_grdpts_id,
     $                                       print_sizes,
     $                                       bf_layer_test_allocation,
     $                                       bf_layer_test_reallocation,
     $                                       bf_layer_test_update_grdpts,
     $                                       ini_nodes,
     $                                       ini_grdpts_id,
     $                                       ini_alignment_table,
     $                                       ini_neighbors_table,
     $                                       ini_general_coord
        use test_bf_mainlayer_module, only : print_mainlayer

        implicit none

        integer, parameter :: nb_sublayers = 3
        
        type(interface_abstract)             :: interface_tested

        integer, dimension(nb_sublayers,2,2) :: test_alignment
        logical, dimension(nb_sublayers,4)   :: test_neighbors

        real(rkind), dimension(nx,ny,ne)     :: nodes
        integer    , dimension(nx,ny)        :: grdpts_id
        integer, dimension(2,2)              :: alignment
        logical, dimension(4)                :: neighbors
        integer, dimension(2)                :: general_coord
        integer, dimension(2)                :: local_coord

        character(len=22)                    :: sizes_filename
        character(len=22)                    :: nodes_filename
        character(len=22)                    :: grdid_filename

        integer     , dimension(8) :: bf_layer_loc
        character(2), dimension(8) :: bf_layer_char
        integer                    :: i,j,k
        type(bf_sublayer), pointer :: added_sublayer
        type(bf_sublayer), pointer :: sublayer_reached

        integer :: mainlayer_id
        integer :: sublayer_id

        real :: random


        !initialize the nodes and print them
        call ini_nodes(nodes)
        call ini_grdpts_id(grdpts_id)
        call print_sizes(nodes,'interior_sizes.dat')
        call print_grdpts_id(grdpts_id,'interior_grdpts_id.dat')
        call print_nodes(nodes,'interior_nodes.dat')
        
        !buffer layer tested
        bf_layer_loc  = [N,S,E,W,N_E,N_W,S_E,S_W]
        bf_layer_char = ['N_','S_','E_','W_','NE','NW','SE','SW']

        !alignment and neighbors
        call ini_alignment_table(test_alignment)
        call ini_neighbors_table(test_neighbors)

        !test add_sublayer
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

        !print the mainlayers of the interface_tested
        call print_mainlayers_interface(interface_tested)


        !test get_sublayer
        !reinitialize the nodes
        do k=1, ne
           do j=1,ny
              do i=1,nx
                 nodes(i,j,k) = 1.0
              end do
           end do
        end do


        !process the sublayers
        !initialize the random generation number
        call srand(100)

        do mainlayer_id=1,4

           do sublayer_id=1, nb_sublayers

              !get a random number
              random = RAND()

              !get the alignment for the sublayer
              alignment(1,1) = test_alignment(sublayer_id,1,1)
              alignment(1,2) = test_alignment(sublayer_id,1,2)
              alignment(2,1) = test_alignment(sublayer_id,2,1)
              alignment(2,2) = test_alignment(sublayer_id,2,2)

              !initialize the general coordinate asked for nodes
              call ini_general_coord(mainlayer_id,alignment,general_coord)

              !modify the corresponding element in the interior nodes
              nodes(general_coord(1), general_coord(2), 1) =
     $             random
!     $             real(mainlayer_id)/10.0

              !initialize the general coordinate asked
              call ini_general_coord_sublayer(
     $             mainlayer_id,alignment,general_coord)
              
              !get the corresponding sublayer
              sublayer_reached => interface_tested%get_sublayer(
     $             general_coord, local_coord)

              !reinitialize the nodes of the sublayer
              do k=1, ne
                 do j=1,size(sublayer_reached%element%nodes,2)
                    do i=1,size(sublayer_reached%element%nodes,1)
                       sublayer_reached%element%nodes(i,j,k) = 1.0
                    end do
                 end do
              end do

              !change the sublayer element corresponding to the general coord
              sublayer_reached%element%nodes(local_coord(1), local_coord(2),1)=
     $             random
!     $             sublayer_reached%element%localization/10.0
           	
              !print the allocated tables on files
              write(sizes_filename,'(A2,I1,''_sizes2.dat'')') bf_layer_char(mainlayer_id),sublayer_id
              write(nodes_filename,'(A2,I1,''_nodes2.dat'')') bf_layer_char(mainlayer_id),sublayer_id
              write(grdid_filename,'(A2,I1,''_grdpt_id2.dat'')') bf_layer_char(mainlayer_id),sublayer_id
              call sublayer_reached%element%print_sizes(sizes_filename)
              call sublayer_reached%element%print_nodes(nodes_filename)
              call sublayer_reached%element%print_grdpts_id(grdid_filename)

           end do

        end do
        
        call print_sizes(nodes,'interior_sizes2.dat')
        call print_grdpts_id(grdpts_id,'interior_grdpts_id2.dat')
        call print_nodes(nodes,'interior_nodes2.dat')
        
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


        subroutine ini_general_coord_sublayer(mainlayer_id, alignment, general_coord)

          implicit none

          integer                       , intent(in)  :: mainlayer_id
          integer       , dimension(2,2), intent(in)  :: alignment
          integer(ikind), dimension(2)  , intent(out) :: general_coord

          select case(mainlayer_id)
            case(N)
               general_coord(1) = alignment(1,1)
               general_coord(2) = ny+1
            case(S)
               general_coord(1) = alignment(1,1)
               general_coord(2) = 0
            case(E)
               general_coord(1) = nx+1
               general_coord(2) = alignment(1,2)
            case(W)
               general_coord(1) = 0
               general_coord(2) = alignment(1,2)
            case(N_E)
               general_coord(1) = nx+1
               general_coord(2) = ny+1
            case(N_W)
               general_coord(1) = 0
               general_coord(2) = ny+1
            case(S_E)
               general_coord(1) = nx+1
               general_coord(2) = 0
            case(S_W)
               general_coord(1) = 0
               general_coord(2) = 0
            case default
               print '(''test_bf_layer_prog'')'
               print '(''ini_general_coord'')'
               print '(''localization not recognized:'',I2)',
     $              mainlayer_id
               stop 'was the buffer layer initialized ?'
          end select

        end subroutine ini_general_coord_sublayer

        subroutine print_mainlayers_interface(
     $     interface_tested)

          implicit none

          type(interface_abstract), intent(in) :: interface_tested
          integer :: i
          type(bf_mainlayer), pointer :: mainlayer_printed

          do i=1, size(interface_tested%mainlayer_pointers,1)

             if(associated(interface_tested%mainlayer_pointers(i)%ptr)) then
                mainlayer_printed => interface_tested%mainlayer_pointers(i)%ptr
                call print_mainlayer(mainlayer_printed)
             end if

          end do

        end subroutine print_mainlayers_interface

      end program test_interface_abstract_prog
