      program test_bf_corner_prog

        use bf_sublayer_class       , only : bf_sublayer
        use bf_interface_class      , only : bf_interface
        use parameters_constant     , only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input        , only : nx,ny,ne
        use parameters_kind         , only : ikind, rkind

        use test_bf_layer_module    , only : print_interior_data,
     $                                       ini_nodes,
     $                                       ini_grdpts_id,
     $                                       ini_alignment_table,
     $                                       ini_neighbors_table

        implicit none

        integer, parameter :: nb_sublayers = 3
        
        type(bf_interface)                   :: interface_tested

        integer, dimension(nb_sublayers,2,2) :: test_alignment
        logical, dimension(nb_sublayers,4)   :: test_neighbors

        real(rkind), dimension(nx,ny,ne)     :: nodes
        integer    , dimension(nx,ny)        :: grdpts_id
        integer, dimension(2,2)              :: alignment
        logical, dimension(4)                :: neighbors
c$$$        integer, dimension(2)                :: general_coord
c$$$        integer, dimension(2)                :: local_coord

c$$$        character(len=22)                    :: sizes_filename
c$$$        character(len=22)                    :: nodes_filename
c$$$        character(len=22)                    :: grdid_filename

        integer     , dimension(8) :: bf_layer_loc
c$$$        character(2), dimension(8) :: bf_layer_char
        integer                    :: i,j
        type(bf_sublayer), pointer :: added_sublayer
c$$$        type(bf_sublayer), pointer :: sublayer_reached
c$$$ 
c$$$        integer :: mainlayer_id
c$$$        integer :: sublayer_id
c$$$
c$$$        real :: random

        !for the get_neighboring_sublayers test
        integer, dimension(4) :: corner_id_tab
        type(bf_sublayer), pointer :: sublayer1
        type(bf_sublayer), pointer :: sublayer2
        real(rkind)                :: color


        !initialize the nodes and print them
        call ini_nodes(nodes)
        call ini_grdpts_id(grdpts_id)
        call print_interior_data(nodes,
     $                           grdpts_id, 
     $                           'interior_nodes.dat',
     $                           'interior_grdpts_id.dat',
     $                           'interior_sizes.dat')
        
        !buffer layer tested
        bf_layer_loc  = [N,S,E,W,N_E,N_W,S_E,S_W]

        !alignment and neighbors
        call ini_alignment_table(test_alignment)
        call ini_neighbors_table(test_neighbors)

        !test add_sublayer
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
     $             i, nodes, alignment, neighbors)
           	
           end do

        end do

        !print the content of the interface on binary files
        call interface_tested%print_binary(
     $       'nodes.dat', 'grdpt_id.dat', 'sizes.dat',
     $       '.dat')


        !test get_neighboring_layers
        corner_id_tab = [N_E,S_E,S_W,N_W]

        do i=1, size(corner_id_tab,1)

           !get the neighboring sublayers corresponding
           !to the corner
           call interface_tested%get_neighboring_sublayers(
     $          corner_id_tab(i), sublayer1, sublayer2)

           !turn the nodes of the first neighboring sublayer
           !with one colour
           color = (i-1)*0.2 + 0.1
           call set_color(sublayer1, color)

           
           !turn the nodes of the second neighboring sublayer
           !with another color
           color = (i-1)*0.2 + 0.2
           call set_color(sublayer2, color)

        end do

        !print the content of the nodes and grdpts_id
        call print_interior_data(nodes,
     $                           grdpts_id, 
     $                           'interior_nodes2.dat',
     $                           'interior_grdpts_id2.dat',
     $                           'interior_sizes2.dat')

        !print the content of the interface on binary files
        call interface_tested%print_binary(
     $       'nodes2.dat', 'grdpt_id2.dat', 'sizes2.dat',
     $       '2.dat')




c$$$        !test get_sublayer
c$$$        !reinitialize the nodes
c$$$        do k=1, ne
c$$$           do j=1,ny
c$$$              do i=1,nx
c$$$                 nodes(i,j,k) = 1.0
c$$$              end do
c$$$           end do
c$$$        end do
c$$$
c$$$
c$$$        !process the sublayers
c$$$        !initialize the random generation number
c$$$        call srand(100)
c$$$
c$$$        do mainlayer_id=1,4
c$$$
c$$$           do sublayer_id=1, nb_sublayers
c$$$
c$$$              !get a random number
c$$$              random = RAND()
c$$$
c$$$              !get the alignment for the sublayer
c$$$              alignment(1,1) = test_alignment(sublayer_id,1,1)
c$$$              alignment(1,2) = test_alignment(sublayer_id,1,2)
c$$$              alignment(2,1) = test_alignment(sublayer_id,2,1)
c$$$              alignment(2,2) = test_alignment(sublayer_id,2,2)
c$$$
c$$$              !initialize the general coordinate asked for nodes
c$$$              call ini_general_coord(mainlayer_id,alignment,general_coord)
c$$$
c$$$              !modify the corresponding element in the interior nodes
c$$$              nodes(general_coord(1), general_coord(2), 1) =
c$$$     $             random
c$$$!     $             real(mainlayer_id)/10.0
c$$$
c$$$              !initialize the general coordinate asked
c$$$              call ini_general_coord_sublayer(
c$$$     $             mainlayer_id,alignment,general_coord)
c$$$              
c$$$              !get the corresponding sublayer
c$$$              sublayer_reached => interface_tested%get_sublayer(
c$$$     $             general_coord, local_coord)
c$$$
c$$$              !reinitialize the nodes of the sublayer
c$$$              do k=1, ne
c$$$                 do j=1,size(sublayer_reached%element%nodes,2)
c$$$                    do i=1,size(sublayer_reached%element%nodes,1)
c$$$                       sublayer_reached%element%nodes(i,j,k) = 1.0
c$$$                    end do
c$$$                 end do
c$$$              end do
c$$$
c$$$              !change the sublayer element corresponding to the general coord
c$$$              sublayer_reached%element%nodes(local_coord(1), local_coord(2),1)=
c$$$     $             random
c$$$!     $             sublayer_reached%element%localization/10.0
c$$$           	
c$$$              !print the allocated tables on files
c$$$              write(sizes_filename,'(A2,I1,''_sizes2.dat'')') bf_layer_char(mainlayer_id),sublayer_id
c$$$              write(nodes_filename,'(A2,I1,''_nodes2.dat'')') bf_layer_char(mainlayer_id),sublayer_id
c$$$              write(grdid_filename,'(A2,I1,''_grdpt_id2.dat'')') bf_layer_char(mainlayer_id),sublayer_id
c$$$              call sublayer_reached%element%print_sizes(sizes_filename)
c$$$              call sublayer_reached%element%print_nodes(nodes_filename)
c$$$              call sublayer_reached%element%print_grdpts_id(grdid_filename)
c$$$
c$$$           end do
c$$$
c$$$        end do
c$$$        
c$$$        call print_sizes(nodes,'interior_sizes2.dat')
c$$$        call print_grdpts_id(grdpts_id,'interior_grdpts_id2.dat')
c$$$        call print_nodes(nodes,'interior_nodes2.dat')
c$$$        
        contains

        subroutine set_color(sublayer, color)

          implicit none

          type(bf_sublayer), pointer, intent(inout) :: sublayer
          real(rkind)               , intent(in)    :: color


          integer(ikind), dimension(2)                  :: sizes_new_nodes
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes

          integer(ikind) :: i,j
          integer :: k


          !get the sizes of the nodes
          sizes_new_nodes = sublayer%get_sizes()


          !allocate the new_nodes to have the same shape
          allocate(new_nodes(
     $         sizes_new_nodes(1),
     $         sizes_new_nodes(2),
     $         ne))

          
          !initialize the new_nodes with the color
          do k=1, ne
             do j=1, sizes_new_nodes(2)
                do i=1, sizes_new_nodes(1)
                   new_nodes(i,j,k) = color
                end do
             end do
          end do


          !set the new_nodes as replacement of the
          !nodes of the sublayer
          call sublayer%set_nodes(new_nodes)

        end subroutine set_color


      end program test_bf_corner_prog
