      !the idea with this test file is to recreate the
      !same results as presented in test_bf_layer_prog
      !but using the interface_abstract as the object
      !encapsulating the buffer layer objects      
      program test_bf_interface_prog

        use bf_layer_bc_procedure_module, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       SW_corner_type,
     $       SE_corner_type,
     $       NW_corner_type,
     $       NE_corner_type

        use bf_mainlayer_class, only :
     $       bf_mainlayer

        use bf_sublayer_class, only :
     $       bf_sublayer

        use bf_interface_class, only :
     $       bf_interface

        use parameters_bf_layer, only :
     $       align_N,
     $       align_S,
     $       align_E,
     $       align_W

        use parameters_constant, only :
     $       N,
     $       S,
     $       E,
     $       W

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use test_bf_layer_module, only :
     $       print_interior_data,
     $       ini_x_map,
     $       ini_y_map,
     $       ini_nodes,
     $       ini_grdpts_id,
     $       ini_cst_nodes

        use sbf_list_class, only :
     $       sbf_list


        implicit none

c$$$        type(bf_interface)                  :: interface_tested
c$$$        real(rkind)   , dimension(nx)       :: x_map
c$$$        real(rkind)   , dimension(ny)       :: y_map
c$$$        real(rkind)   , dimension(nx,ny,ne) :: nodes
c$$$        integer       , dimension(nx,ny)    :: grdpts_id
c$$$        integer       , dimension(2,2)      :: alignment
c$$$        integer                             :: mainlayer_id
c$$$        type(bf_sublayer), pointer          :: added_sublayer
c$$$        type(bf_sublayer), pointer          :: bf_sublayer_merged1
c$$$        type(bf_sublayer), pointer          :: bf_sublayer_merged2
c$$$        type(bf_sublayer), pointer          :: bf_sublayer_reallocated
c$$$        type(bf_sublayer), pointer          :: merged_sublayer
c$$$        real(rkind)                         :: scale
c$$$
c$$$
c$$$        integer :: i, index
c$$$
        logical :: detailled
        logical :: test_loc
        logical :: test_validated
c$$$
c$$$
c$$$        !initialize the nodes and print them
c$$$        call ini_x_map(x_map)
c$$$        call ini_y_map(y_map)
c$$$        call ini_nodes(nodes)
c$$$        call ini_grdpts_id(grdpts_id)
c$$$        call print_interior_data(
c$$$     $       x_map,
c$$$     $       y_map,
c$$$     $       nodes,
c$$$     $       grdpts_id, 
c$$$     $       'interior_x_map.dat',
c$$$     $       'interior_y_map.dat',
c$$$     $       'interior_nodes.dat',
c$$$     $       'interior_grdpts_id.dat',
c$$$     $       'interior_sizes.dat')
c$$$        
c$$$        !initialize the interface
c$$$        call interface_tested%ini(x_map,y_map)
c$$$
c$$$        
c$$$        !add the N_E buffer layer
c$$$        call get_alignment(1, alignment, mainlayer_id)
c$$$        added_sublayer => interface_tested%allocate_sublayer(
c$$$     $       mainlayer_id, x_map, y_map, nodes, alignment)
c$$$        scale = 0.1
c$$$        call ini_cst_nodes(added_sublayer, scale)
c$$$        
c$$$        !add the W buffer layer under the N_W corner
c$$$        call get_alignment(4, alignment, mainlayer_id)
c$$$        added_sublayer => interface_tested%allocate_sublayer(
c$$$     $       mainlayer_id, x_map, y_map, nodes, alignment)
c$$$        scale = 0.2
c$$$        call ini_cst_nodes(added_sublayer, scale)
c$$$
c$$$
c$$$        !print the interface
c$$$        call interface_tested%print_binary(
c$$$     $       'x_map1.dat',
c$$$     $       'y_map1.dat',
c$$$     $       'nodes1.dat',
c$$$     $       'grdpt_id1.dat',
c$$$     $       'sizes1.dat',
c$$$     $       '1.dat')
c$$$
c$$$        
c$$$        !add the E buffer layer under the N_E corner
c$$$        call get_alignment(3, alignment, mainlayer_id)
c$$$        bf_sublayer_merged1 => interface_tested%allocate_sublayer(
c$$$     $       mainlayer_id, x_map, y_map, nodes, alignment)
c$$$        
c$$$        !add the N_W buffer layer
c$$$        call get_alignment(2, alignment, mainlayer_id)
c$$$        added_sublayer => interface_tested%allocate_sublayer(
c$$$     $       mainlayer_id, x_map, y_map, nodes, alignment)
c$$$
c$$$        
c$$$        !print the interface
c$$$        call interface_tested%print_binary(
c$$$     $       'x_map2.dat',
c$$$     $       'y_map2.dat',
c$$$     $       'nodes2.dat',
c$$$     $       'grdpt_id2.dat',
c$$$     $       'sizes2.dat',
c$$$     $       '2.dat')
c$$$
c$$$
c$$$        !add the E buffer layer
c$$$        call get_alignment(5, alignment, mainlayer_id)
c$$$        bf_sublayer_merged2 => interface_tested%allocate_sublayer(
c$$$     $       mainlayer_id, x_map, y_map, nodes, alignment)
c$$$
c$$$        !add the W buffer layer
c$$$        call get_alignment(6, alignment, mainlayer_id)
c$$$        bf_sublayer_reallocated => interface_tested%allocate_sublayer(
c$$$     $       mainlayer_id, x_map, y_map, nodes, alignment)
c$$$
c$$$        
c$$$        !print the interface
c$$$        call interface_tested%print_binary(
c$$$     $       'x_map3.dat',
c$$$     $       'y_map3.dat',
c$$$     $       'nodes3.dat',
c$$$     $       'grdpt_id3.dat',
c$$$     $       'sizes3.dat',
c$$$     $       '3.dat')
c$$$
c$$$
c$$$        !add the S_E buffer layer
c$$$        call get_alignment(7, alignment, mainlayer_id)
c$$$        added_sublayer => interface_tested%allocate_sublayer(
c$$$     $       mainlayer_id, x_map, y_map, nodes, alignment)
c$$$        scale = 0.3
c$$$        call ini_cst_nodes(added_sublayer, scale)
c$$$        
c$$$
c$$$        !add the S_W buffer layer
c$$$        call get_alignment(8, alignment, mainlayer_id)
c$$$        added_sublayer => interface_tested%allocate_sublayer(
c$$$     $       mainlayer_id, x_map, y_map, nodes, alignment)
c$$$        scale = 0.4
c$$$        call ini_cst_nodes(added_sublayer, scale)
c$$$
c$$$
c$$$        !print the interface
c$$$        call interface_tested%print_binary(
c$$$     $       'x_map4.dat',
c$$$     $       'y_map4.dat',
c$$$     $       'nodes4.dat',
c$$$     $       'grdpt_id4.dat',
c$$$     $       'sizes4.dat',
c$$$     $       '4.dat')
c$$$
c$$$
c$$$        !merge the E buffer layers
c$$$        call get_alignment(9, alignment, mainlayer_id)
c$$$        merged_sublayer => interface_tested%merge_sublayers(
c$$$     $       bf_sublayer_merged1, bf_sublayer_merged2,
c$$$     $       x_map, y_map, nodes, alignment)
c$$$        
c$$$
c$$$        !reallocate the W buffer layer
c$$$        call get_alignment(10, alignment, mainlayer_id)
c$$$        call interface_tested%reallocate_sublayer(
c$$$     $       bf_sublayer_reallocated,
c$$$     $       x_map, y_map, nodes, alignment)
c$$$
c$$$
c$$$        !print the interface
c$$$        call interface_tested%print_binary(
c$$$     $       'x_map5.dat',
c$$$     $       'y_map5.dat',
c$$$     $       'nodes5.dat',
c$$$     $       'grdpt_id5.dat',
c$$$     $       'sizes5.dat',
c$$$     $       '5.dat')
c$$$
c$$$        
c$$$        !add more sublayers on N,S
c$$$        do i=11,16
c$$$           call get_alignment(i, alignment, mainlayer_id)
c$$$           added_sublayer => interface_tested%allocate_sublayer(
c$$$     $          mainlayer_id, x_map, y_map, nodes, alignment)
c$$$        end do
c$$$
c$$$        call interface_tested%print_binary(
c$$$     $       'x_map6.dat',
c$$$     $       'y_map6.dat',
c$$$     $       'nodes6.dat',
c$$$     $       'grdpt_id6.dat',
c$$$     $       'sizes6.dat',
c$$$     $       '6.dat')
c$$$
c$$$        index = 7
c$$$
c$$$
c$$$        !test bf_depends_on_neighbors
c$$$        !-------------------------------------------------------
c$$$        !for each sublayer in each mainlayer, test if the buffer layer
c$$$        !depends on neighbors
c$$$        call test_bf_layer_depends_on_neighbors(interface_tested)
c$$$        
c$$$
c$$$        !test get_nbf_layers_sharing_grdpts_with
c$$$        !-------------------------------------------------------
c$$$        !for each sublayer in each mainlayer, test if the buffer layer
c$$$        !has neighbor dependencies and colorize the neighbors
c$$$        call test_get_nbf_layers_sharing_grdpts_with(interface_tested, index)
c$$$
c$$$        
c$$$        !test if the neighboring buffer layers will be removed
c$$$        !-------------------------------------------------------
c$$$        call test_does_a_neighbor_remains(interface_tested)
c$$$
c$$$
c$$$        !test the removal of sublayer
c$$$        !----------------------------
c$$$        call test_remove_sublayer(
c$$$     $       interface_tested,
c$$$     $       merged_sublayer,
c$$$     $       index,
c$$$     $       x_map,
c$$$     $       y_map)
c$$$
c$$$
c$$$        !re-test if the neighboring buffer layers will be removed
c$$$        !--------------------------------------------------------
c$$$        call test_does_a_neighbor_remains(interface_tested)
c$$$
c$$$        
c$$$        !test the netcdf writing
c$$$        !--------------------------------------------------------
c$$$        call test_print_netcdf(
c$$$     $       interface_tested,1)
c$$$
c$$$
c$$$        !test the determination of the bc_sections
c$$$        !--------------------------------------------------------
c$$$        call test_determine_interior_bc_layers(
c$$$     $       interface_tested,
c$$$     $       x_map,
c$$$     $       y_map,
c$$$     $       nodes,
c$$$     $       grdpts_id,
c$$$     $       index)
c$$$
c$$$
c$$$        !test the determination of the bc_procedures
c$$$        !--------------------------------------------------------
c$$$        call test_determine_interior_bc_procedures(
c$$$     $       interface_tested,
c$$$     $       x_map,
c$$$     $       y_map,
c$$$     $       nodes,
c$$$     $       grdpts_id,
c$$$     $       index)
c$$$
c$$$
c$$$        !test the node synchronization with the interior domain
c$$$        !--------------------------------------------------------
c$$$        call test_sync_nodes_with_interior(
c$$$     $       interface_tested,
c$$$     $       x_map,
c$$$     $       y_map,
c$$$     $       nodes,
c$$$     $       grdpts_id,
c$$$     $       index)
c$$$
c$$$        !test the node synchronization at the interface between
c$$$        !buffer main layers
c$$$        !--------------------------------------------------------
c$$$        call test_sync_nodes_at_mainlayer_interfaces(
c$$$     $       interface_tested,
c$$$     $       x_map,
c$$$     $       y_map,
c$$$     $       nodes,
c$$$     $       grdpts_id,
c$$$     $       index)
c$$$
c$$$        !test the integration borders for the buffer layers
c$$$        !--------------------------------------------------------
c$$$        call test_bf_integration_borders(
c$$$     $       interface_tested,
c$$$     $       x_map,
c$$$     $       y_map,
c$$$     $       nodes,
c$$$     $       grdpts_id,
c$$$     $       index)

        !test the resolving of bc overlap conflicts for buffer layers
        !------------------------------------------------------------
        detailled = .true.
        test_loc = test_resolve_bc_overlap_conflicts1(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_resolve_bc_overlap_conflicts1: '',L1)', test_loc
        print '()'

        detailled = .true.
        test_loc = test_resolve_bc_overlap_conflicts2(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_resolve_bc_overlap_conflicts2: '',L1)', test_loc
        print '()'


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

            !north outside allocation N_W
            case(11)
               mainlayer_id   = N
               alignment(1,1) = align_W-2*large_layer-small_layer
               alignment(1,2) = alignment(1,1)+small_layer
               alignment(2,1) = align_N
               alignment(2,2) = align_N               

            !north inside allocation
            case(12)
               mainlayer_id   = N
               alignment(1,1) = align_W+2*large_layer
               alignment(1,2) = alignment(1,1)+small_layer
               alignment(2,1) = align_N
               alignment(2,2) = align_N

            !north outside allocation N_E
            case(13)
               mainlayer_id   = N
               alignment(1,1) = align_E+2*large_layer
               alignment(1,2) = alignment(1,1)+small_layer
               alignment(2,1) = align_N
               alignment(2,2) = align_N

            !south outside allocation S_W
            case(14)
               mainlayer_id   = S
               alignment(1,1) = align_W-2*large_layer-small_layer
               alignment(1,2) = alignment(1,1)+small_layer
               alignment(2,1) = align_S
               alignment(2,2) = align_S               

            !south inside allocation
            case(15)
               mainlayer_id   = S
               alignment(1,1) = align_W+2*large_layer
               alignment(1,2) = alignment(1,1)+small_layer
               alignment(2,1) = align_S
               alignment(2,2) = align_S

            !south outside allocation S_E
            case(16)
               mainlayer_id   = S
               alignment(1,1) = align_E+2*large_layer
               alignment(1,2) = alignment(1,1)+small_layer
               alignment(2,1) = align_S
               alignment(2,2) = align_S

            case default
               print '(''case not implemented'')'
               stop ''
          end select               

        end subroutine get_alignment


        !< initialize the nodes using a constant variable
        subroutine ini_cst(nodes, color_i)

          implicit none

          real(rkind), dimension(:,:,:), intent(out) :: nodes
          real(rkind)        , optional, intent(in)  :: color_i

          integer(ikind) :: i,j
          real(rkind) :: color

          if(present(color_i)) then
             color = color_i
          else
             color = 1.0
          end if

          do j=1, size(nodes,2)
             do i=1, size(nodes,1)
                nodes(i,j,:)  = color
             end do
          end do

        end subroutine ini_cst


        !re-initialize the nodes of all sublayers
        subroutine reinitialize_nodes(interface_used)

          implicit none

          class(bf_interface), intent(inout) :: interface_used

          integer :: k,l
          type(bf_mainlayer), pointer :: mainlayer_ptr
          type(bf_sublayer) , pointer :: sublayer_ptr
          integer                     :: nb_sublayers
          integer(ikind), dimension(2):: sizes
          real(rkind), dimension(:,:,:), allocatable :: new_nodes
          
          
          do k=1,4

             mainlayer_ptr => interface_used%get_mainlayer(k)

             if(associated(mainlayer_ptr)) then
                nb_sublayers = mainlayer_ptr%get_nb_sublayers()

                sublayer_ptr => mainlayer_ptr%get_head_sublayer()

                do l=1, nb_sublayers
                   sizes = sublayer_ptr%get_sizes()
                   allocate(new_nodes(sizes(1), sizes(2), ne))
                   call ini_cst(new_nodes)
                   call sublayer_ptr%set_nodes(new_nodes)
                   sublayer_ptr => sublayer_ptr%get_next()
                end do

             end if

          end do          

        end subroutine reinitialize_nodes


        !re-initialize the nodes of all sublayers
        subroutine colorize_integration_borders(interface_used)

          implicit none

          class(bf_interface), intent(inout) :: interface_used

          integer :: k,l,m,i,j
          type(bf_mainlayer), pointer :: mainlayer_ptr
          type(bf_sublayer) , pointer :: sublayer_ptr
          integer                     :: nb_sublayers
          integer(ikind), dimension(2):: sizes
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes
          real(rkind)   , dimension(2)                  :: x_borders
          real(rkind)   , dimension(2)                  :: y_borders
          integer(ikind), dimension(:,:)  , allocatable :: N_bc_sections
          integer(ikind), dimension(:,:)  , allocatable :: S_bc_sections
          real(rkind) :: color_S_bc_sec
          real(rkind) :: color_integration
          real(rkind) :: color_N_bc_sec

          color_S_bc_sec    = 0.7
          color_integration = 0.4
          color_N_bc_sec    = 0.9

          
          do k=1,4

             mainlayer_ptr => interface_used%get_mainlayer(k)

             if(associated(mainlayer_ptr)) then
                nb_sublayers = mainlayer_ptr%get_nb_sublayers()

                sublayer_ptr => mainlayer_ptr%get_head_sublayer()

                do l=1, nb_sublayers
                   sizes = sublayer_ptr%get_sizes()

                   allocate(new_nodes(sizes(1), sizes(2), ne))

                   !get the x_borders and y_borders
                   x_borders = sublayer_ptr%get_x_borders()
                   y_borders = sublayer_ptr%get_y_borders()

                   !get the N and S bc_sections
                   call sublayer_ptr%get_S_bc_sections(S_bc_sections)
                   call sublayer_ptr%get_N_bc_sections(N_bc_sections)

                   
                   do j=1, sizes(2)
                      do i=1, sizes(1)
                         new_nodes(i,j,:) = 1.0
                      end do
                   end do

                   !colorize S_bc_sections
                   if(allocated(S_bc_sections)) then
                      do m=1, size(S_bc_sections,2)
                         do j=1,bc_size
                            do i=S_bc_sections(1,m), S_bc_sections(2,m)
                               new_nodes(i,j,:) = color_S_bc_sec
                            end do
                         end do
                      end do
                      deallocate(S_bc_sections)
                   end if

                   !colorize interior domain for integration
                   do j=y_borders(1), y_borders(2)
                      do i=x_borders(1), x_borders(2)
                         new_nodes(i,j,:) = color_integration
                      end do
                   end do

                   !colorize N_bc_sections
                   if(allocated(N_bc_sections)) then
                      do m=1, size(N_bc_sections,2)
                         do j=sizes(2)-bc_size+1,sizes(2)
                            do i=N_bc_sections(1,m), N_bc_sections(2,m)
                               new_nodes(i,j,:) = color_N_bc_sec
                            end do
                         end do
                      end do
                      deallocate(N_bc_sections)
                   end if

                   call sublayer_ptr%set_nodes(new_nodes)
                   sublayer_ptr => sublayer_ptr%get_next()
                end do

             end if

          end do          

        end subroutine colorize_integration_borders


        !re-initialize the nodes of all sublayers
        subroutine reinitialize_nodes_with_gradient(interface_used)

          implicit none

          class(bf_interface), intent(inout) :: interface_used

          integer :: k,l
          type(bf_mainlayer), pointer :: mainlayer_ptr
          type(bf_sublayer) , pointer :: sublayer_ptr
          integer                     :: nb_sublayers
          
          real(rkind), dimension(4) :: color
          real(rkind), dimension(4) :: scale

          color(N) = 0.1
          scale(N) = 0.3

          color(S) = 0.6
          scale(S) = 0.2

          color(E) = 0.4
          scale(E) = 0.2

          color(W) = 0.8
          scale(W) = 0.2

          do k=1,4

             mainlayer_ptr => interface_used%get_mainlayer(k)

             if(associated(mainlayer_ptr)) then
                nb_sublayers = mainlayer_ptr%get_nb_sublayers()

                sublayer_ptr => mainlayer_ptr%get_head_sublayer()

                do l=1, nb_sublayers
                   call ini_cst_nodes(sublayer_ptr,color(k),scale(k))
                   sublayer_ptr => sublayer_ptr%get_next()
                end do

             end if

          end do          

        end subroutine reinitialize_nodes_with_gradient


        !test bf_layer
        subroutine test_bf_layer_depends_on_neighbors(interface_used)
        
          implicit none

          class(bf_interface), intent(inout) :: interface_used

          integer :: k,l
          type(bf_mainlayer), pointer    :: mainlayer_ptr
          type(bf_sublayer) , pointer    :: sublayer_ptr
          integer                        :: nb_sublayers
          logical                        :: dependent
          integer(ikind), dimension(2,2) :: alignment


          do k=1,4
             
             print '(''mainlayer_id: '',I1)', k
             print '(''----------------'')'

             mainlayer_ptr => interface_used%get_mainlayer(k)
             
             if(associated(mainlayer_ptr)) then
                nb_sublayers = mainlayer_ptr%get_nb_sublayers()
                
                sublayer_ptr => mainlayer_ptr%get_head_sublayer()
                
                do l=1, nb_sublayers
                   
                   dependent = interface_used%bf_layer_depends_on_neighbors(sublayer_ptr)
                   alignment = sublayer_ptr%get_alignment_tab()
                   
                   print '(''sublayer '', ''('',4I3,'') dep: '',L1)',
     $                  alignment(1,1), alignment(1,2),
     $                  alignment(2,1), alignment(2,2),
     $                  dependent

                   sublayer_ptr => sublayer_ptr%get_next()
                end do
                
             end if
             
             print '()'
             
          end do

        end subroutine test_bf_layer_depends_on_neighbors


        !> test for the dependencies
        subroutine test_get_nbf_layers_sharing_grdpts_with(interface_used, index)

          implicit none

          class(bf_interface), intent(inout) :: interface_used
          integer            , intent(inout) :: index


          integer :: k,l
          type(bf_mainlayer), pointer    :: mainlayer_ptr
          type(bf_sublayer) , pointer    :: sublayer_ptr
          integer                        :: nb_sublayers
          type(sbf_list) :: nbf1_list
          type(sbf_list) :: nbf2_list

          
          do k=1,4           

             mainlayer_ptr => interface_used%get_mainlayer(k)
             
             if(associated(mainlayer_ptr)) then
                nb_sublayers = mainlayer_ptr%get_nb_sublayers()
                
                sublayer_ptr => mainlayer_ptr%get_head_sublayer()
                
                do l=1, nb_sublayers

                   call reinitialize_nodes(interface_used)

                   call nbf1_list%ini(8)
                   call nbf2_list%ini(8)

                   call interface_used%get_nbf_layers_sharing_grdpts_with(
     $                  sublayer_ptr, nbf1_list, nbf2_list)

                   !colorize the mainlayer in red
                   call colorize_main(sublayer_ptr)

                   !colorize the neighbors in green
                   call colorize_neighbors(nbf1_list)
                   call colorize_neighbors(nbf2_list)

                   call nbf1_list%destroy()
                   call nbf2_list%destroy()

                   !print the result on output files
                   call print_output(interface_used, index)

                   index = index+1

                   sublayer_ptr => sublayer_ptr%get_next()
                end do

             end if

          end do

        end subroutine test_get_nbf_layers_sharing_grdpts_with


        !< test whether the neighboring buffer layers remain
        subroutine test_does_a_neighbor_remains(interface_used)

          implicit none

          class(bf_interface), intent(in) :: interface_used

          integer :: k,l
          type(bf_mainlayer), pointer    :: mainlayer_ptr
          type(bf_sublayer) , pointer    :: sublayer_ptr
          integer                        :: nb_sublayers

          logical :: a_neighbor_remains          

          
          !set the remain status of the buffer layers to .true.
          do k=1,4           

             mainlayer_ptr => interface_used%get_mainlayer(k)
             
             if(associated(mainlayer_ptr)) then
                nb_sublayers = mainlayer_ptr%get_nb_sublayers()
                
                sublayer_ptr => mainlayer_ptr%get_head_sublayer()
                
                do l=1, nb_sublayers

                   call sublayer_ptr%set_remain_status(.true.)
                   sublayer_ptr => sublayer_ptr%get_next()

                end do

             end if

          end do

          
          !test if removal is approved
          do k=1,4           

             mainlayer_ptr => interface_used%get_mainlayer(k)
             
             print '(''mainlayer '',I2)', k
             print '(''---------------'')'

             if(associated(mainlayer_ptr)) then
                nb_sublayers = mainlayer_ptr%get_nb_sublayers()
                
                sublayer_ptr => mainlayer_ptr%get_head_sublayer()
                
                do l=1, nb_sublayers

                   a_neighbor_remains = interface_used%does_a_neighbor_remains(
     $                  sublayer_ptr,k)

                   print '(''sublayer '',I2,'': '',L1)', l, a_neighbor_remains

                   sublayer_ptr => sublayer_ptr%get_next()

                end do

             end if

             print '()'

          end do

        end subroutine test_does_a_neighbor_remains


        !test the removal of a sublayer
        subroutine test_remove_sublayer(
     $     interface_used,
     $     sublayer_ptr,
     $     index,
     $     interior_x_map,
     $     interior_y_map)
        
          implicit none

          class(bf_interface)       , intent(inout) :: interface_used
          type(bf_sublayer), pointer, intent(inout) :: sublayer_ptr
          integer                   , intent(inout) :: index
          real(rkind), dimension(nx), intent(in)    :: interior_x_map
          real(rkind), dimension(ny), intent(in)    :: interior_y_map

          print '(''test remove_sublayer'')'
          print '(''--------------------'')'

          call interface_used%remove_sublayer(
     $         sublayer_ptr,
     $         interior_x_map,
     $         interior_y_map)

          call print_output(interface_used, index)

          index = index+1

          print '()'

        end subroutine test_remove_sublayer


        !test the determination of interior boundary layers
        subroutine test_determine_interior_bc_layers(
     $     interface_used,
     $     x_map,
     $     y_map,
     $     nodes,
     $     grdpts_id,
     $     index)
        
          implicit none
          
          class(bf_interface)             , intent(inout) :: interface_used
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer    , dimension(nx,ny)   , intent(in)    :: grdpts_id
          integer                         , intent(inout) :: index

          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_N
          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_S
          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_E
          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_W

          integer(ikind) :: i,j
          integer        :: k

          !0) reinitialize the nodes of the boundary layers
          call interface_used%determine_interior_bc_layers(
     $         interior_bc_sections_N,
     $         interior_bc_sections_S,
     $         interior_bc_sections_E,
     $         interior_bc_sections_W)

          !display the boundary layers
          !1) reinitialize the nodes
          call reinitialize_nodes(interface_used)

          !2) for each piece of boundary layer, colorize the
          !   corresponding piece on the graph
          do k=1,ne
             do j=1,ny
                do i=1,nx
                   nodes(i,j,k) = 1.0
                end do
             end do
          end do

          if(allocated(interior_bc_sections_S)) then
             do k=1, size(interior_bc_sections_S,2)
                do j=1,bc_size
                   do i=interior_bc_sections_S(1,k),interior_bc_sections_S(2,k)
                      nodes(i,j,:) = 0.3
                   end do
                end do
             end do
          end if

          if(allocated(interior_bc_sections_E)) then
             do k=1, size(interior_bc_sections_E,2)
                do j=interior_bc_sections_E(1,k),interior_bc_sections_E(2,k)
                   do i=nx-bc_size+1,nx
                      nodes(i,j,:) = 0.5
                   end do
                end do
             end do
          end if

          if(allocated(interior_bc_sections_W)) then
             do k=1, size(interior_bc_sections_W,2)
                do j=interior_bc_sections_W(1,k),interior_bc_sections_W(2,k)
                   do i=1,bc_size
                      nodes(i,j,:) = 0.7
                   end do
                end do
             end do
          end if

          if(allocated(interior_bc_sections_N)) then
             do k=1, size(interior_bc_sections_N,2)
                do j=ny-bc_size+1,ny
                   do i=interior_bc_sections_N(1,k),interior_bc_sections_N(2,k)
                      nodes(i,j,:) = 0.1
                   end do
                end do
             end do
          end if  

          !print interface
          call print_output(interface_used, index)

          !print interior nodes
          call print_interior_data(
     $         x_map,
     $         y_map,
     $         nodes,
     $         grdpts_id,
     $         'interior_x_map_bc_sec.dat',
     $         'interior_y_map_bc_sec.dat',
     $         'interior_nodes_bc_sec.dat',
     $         'interior_grdpts_id_bc_sec.dat',
     $         'interior_sizes_bc_sec.dat')

          index = index+1

        end subroutine test_determine_interior_bc_layers


        !test the determination of interior boundary layers
        subroutine test_determine_interior_bc_procedures(
     $     interface_used,
     $     x_map,
     $     y_map,
     $     nodes,
     $     grdpts_id,
     $     index)
        
          implicit none
          
          class(bf_interface)             , intent(inout) :: interface_used
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer    , dimension(nx,ny)   , intent(in)    :: grdpts_id
          integer                         , intent(inout) :: index

          integer(ikind), dimension(:,:), allocatable :: bc_procedures

          integer(ikind) :: i,j
          integer        :: k

          real(rkind), dimension(8) :: color

          !0) reinitialize the nodes of the boundary layers
          call interface_used%determine_interior_bc_procedures(
     $         bc_procedures)

          !display the boundary layers
          !1) reinitialize the nodes of interface
          call reinitialize_nodes(interface_used)

          !2) reinitialize the interior nodes
          do k=1,ne
             do j=1,ny
                do i=1,nx
                   nodes(i,j,k) = 1.0
                end do
             end do
          end do

          !3) for each procedure, colorize the
          !   corresponding piece on the graph
          color(SW_corner_type) = 0.1
          color(SE_corner_type) = 0.2
          color(NW_corner_type) = 0.3
          color(NE_corner_type) = 0.4
          color(N_edge_type)    = 0.5
          color(S_edge_type)    = 0.6
          color(E_edge_type)    = 0.7
          color(W_edge_type)    = 0.8


          if(allocated(bc_procedures)) then
             do k=1, size(bc_procedures,2)

                select case(bc_procedures(1,k))

                  case(N_edge_type,S_edge_type)

                     do j = bc_procedures(3,k), bc_procedures(3,k)+1
                        do i=bc_procedures(2,k), bc_procedures(4,k)
                           nodes(i,j,:) = color(bc_procedures(1,k))
                        end do
                     end do

                  case(E_edge_type,W_edge_type)
                     
                     do j = bc_procedures(3,k), bc_procedures(4,k)
                        do i=bc_procedures(2,k), bc_procedures(2,k)+1
                           nodes(i,j,:) = color(bc_procedures(1,k))
                        end do
                     end do

                  case(SW_corner_type,SE_corner_type,
     $                 NW_corner_type,NE_corner_type)
                  
                     do j = bc_procedures(3,k), bc_procedures(3,k)+1
                        do i=bc_procedures(2,k), bc_procedures(2,k)+1
                           nodes(i,j,:) = color(bc_procedures(1,k))
                        end do
                     end do

                 case default
                    print '(''test_bf_interface_prog'')'
                    print '(''determine_interior_bc_procedures'')'
                    stop 'case not recognized'

               end select

             end do
          end if  

          !print interface
          call print_output(interface_used, index)

          !print interior nodes
          call print_interior_data(
     $         x_map,
     $         y_map,
     $         nodes,
     $         grdpts_id,
     $         'interior_x_map_bc_proc.dat',
     $         'interior_y_map_bc_proc.dat',
     $         'interior_nodes_bc_proc.dat',
     $         'interior_grdpts_id_bc_proc.dat',
     $         'interior_sizes_bc_proc.dat')

          index = index+1

        end subroutine test_determine_interior_bc_procedures


        !test the determination of interior boundary layers
        subroutine test_sync_nodes_with_interior(
     $     interface_used,
     $     x_map,
     $     y_map,
     $     nodes,
     $     grdpts_id,
     $     index)
        
          implicit none
          
          class(bf_interface)             , intent(inout) :: interface_used
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer    , dimension(nx,ny)   , intent(in)    :: grdpts_id
          integer                         , intent(inout) :: index

          integer(ikind) :: i,j
          integer        :: k

          !0) reinitialize the nodes of the boundary layers
          call reinitialize_nodes_with_gradient(interface_used)

          !2) reinitialize the interior nodes
          do k=1,ne
             do j=1,ny
                do i=1,nx
                   nodes(i,j,k) = real(i+j)/real(nx+ny)
                end do
             end do
          end do

          !3) exchange with interior
          call interface_used%sync_nodes_with_interior(nodes)

          !print interface
          call print_output(interface_used, index)

          !print interior nodes
          call print_interior_data(
     $         x_map,
     $         y_map,
     $         nodes,
     $         grdpts_id,
     $         'interior_x_map_sync_int.dat',
     $         'interior_y_map_sync_int.dat',
     $         'interior_nodes_sync_int.dat',
     $         'interior_grdpts_id_sync_int.dat',
     $         'interior_sizes_sync_int.dat')

          index = index+1

        end subroutine test_sync_nodes_with_interior


        !test the determination of interior boundary layers
        subroutine test_sync_nodes_at_mainlayer_interfaces(
     $     interface_used,
     $     x_map,
     $     y_map,
     $     nodes,
     $     grdpts_id,
     $     index)
        
          implicit none
          
          class(bf_interface)             , intent(inout) :: interface_used
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer    , dimension(nx,ny)   , intent(in)    :: grdpts_id
          integer                         , intent(inout) :: index

          integer(ikind) :: i,j
          integer        :: k

          !0) reinitialize the nodes of the boundary layers
          call reinitialize_nodes_with_gradient(interface_used)

          !2) reinitialize the interior nodes
          do k=1,ne
             do j=1,ny
                do i=1,nx
                   nodes(i,j,k) = real(i+j)/real(nx+ny)
                end do
             end do
          end do

          !3) exchange with interior
          call interface_used%sync_nodes_at_mainlayer_interfaces()

          !print interface
          call print_output(interface_used, index)

          !print interior nodes
          call print_interior_data(
     $         x_map,
     $         y_map,
     $         nodes,
     $         grdpts_id,
     $         'interior_x_map_sync_bml.dat',
     $         'interior_y_map_sync_bml.dat',
     $         'interior_nodes_sync_bml.dat',
     $         'interior_grdpts_id_sync_bml.dat',
     $         'interior_sizes_sync_bml.dat')

          index = index+1

        end subroutine test_sync_nodes_at_mainlayer_interfaces


        !test the determination of the integration borders
        subroutine test_bf_integration_borders(
     $     interface_used,
     $     x_map,
     $     y_map,
     $     nodes,
     $     grdpts_id,
     $     index)
        
          implicit none
          
          class(bf_interface)             , intent(inout) :: interface_used
          real(rkind), dimension(nx)      , intent(in)    :: x_map
          real(rkind), dimension(ny)      , intent(in)    :: y_map
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          integer    , dimension(nx,ny)   , intent(in)    :: grdpts_id
          integer                         , intent(inout) :: index

          integer(ikind) :: i,j
          integer        :: k

          !1) reinitialize the nodes of the boundary layers
          call reinitialize_nodes(interface_used)

          !2) reinitialize the interior nodes
          do k=1,ne
             do j=1,ny
                do i=1,nx
                   nodes(i,j,k) = 1.0
                end do
             end do
          end do

          !3) colorize the integration borders in the buffer layers
          call colorize_integration_borders(interface_used)

          !4) print interface
          call print_output(interface_used, index)

          !5) print interior nodes
          call print_interior_data(
     $         x_map,
     $         y_map,
     $         nodes,
     $         grdpts_id,
     $         'interior_x_map_int_borders.dat',
     $         'interior_y_map_int_borders.dat',
     $         'interior_nodes_int_borders.dat',
     $         'interior_grdpts_id_int_borders.dat',
     $         'interior_sizes_int_borders.dat')

          index = index+1

        end subroutine test_bf_integration_borders

      
        !< colorize the main layer whose dependencies are deterimed
        subroutine colorize_main(sublayer_used)

          implicit none

          type(bf_sublayer), intent(inout) :: sublayer_used

          integer(ikind), dimension(2) :: sizes
          real(rkind), dimension(:,:,:), allocatable :: new_nodes

          sizes = sublayer_used%get_sizes()
          allocate(new_nodes(sizes(1), sizes(2), ne))
          
          call ini_cst(new_nodes, 0.9d0)

          call sublayer_used%set_nodes(new_nodes)

        end subroutine colorize_main


        !< colorize the main layer whose dependencies are deterimed
        subroutine colorize_neighbors(sbf_list_used)

          implicit none

          type(sbf_list), intent(inout) :: sbf_list_used

          
          integer                      :: nb_ele
          integer                      :: k
          type(bf_sublayer), pointer   :: sublayer_ptr
          integer(ikind), dimension(2) :: sizes
          real(rkind)   , dimension(:,:,:), allocatable :: new_nodes

          
          nb_ele = sbf_list_used%get_nb_ele()

          do k=1, nb_ele

             sublayer_ptr => sbf_list_used%get_ele(k)

             sizes = sublayer_ptr%get_sizes()
             allocate(new_nodes(sizes(1), sizes(2), ne))
          
             call ini_cst(new_nodes, 0.7d0)
             
             call sublayer_ptr%set_nodes(new_nodes)

          end do

        end subroutine colorize_neighbors


        !print outputs
        subroutine print_output(interface_used, index)

          implicit none

          class(bf_interface), intent(in) :: interface_used
          integer            , intent(in) :: index


          integer :: format_index

          character(len=15) :: format_nodes
          character(len=18) :: format_grdpt
          character(len=10) :: format_nbsbl
          
          character(len=11) :: filename_x_map
          character(len=11) :: filename_y_map
          character(len=11) :: filename_nodes
          character(len=14) :: filename_grdpt
          character(len=11) :: filename_sizes
          character(len=6 ) :: filename_nbsbl


          !determine the number of integer needed to write the
          !file index
          if(index.le.9) then
             format_index = 1
          else
             if((index.ge.10).and.(index.le.99)) then
                format_index = 2
             else
                print '(''test_bf_interface_prog'')'
                print '(''print_output'')'
                stop 'file_index not supported'
             end if
          end if
          

          !determine the format for the name of the output files
          write(format_nodes, '(''(A5,I'',I1,'',A4)'')') format_index
          write(format_grdpt, '(''(A8,I'',I1,'',A4)'')') format_index
          write(format_nbsbl, '(''(I'',I1,'',A4)'')'  ) format_index
          

          !determine the name of the output files
          write(filename_x_map, format_nodes) 'x_map', index, '.dat'
          write(filename_y_map, format_nodes) 'y_map', index, '.dat'
          write(filename_nodes, format_nodes) 'nodes', index, '.dat'
          write(filename_grdpt, format_grdpt) 'grdpt_id', index, '.dat'
          write(filename_sizes, format_nodes) 'sizes', index, '.dat'
          write(filename_nbsbl, format_nbsbl) index, '.dat'


          !write the outputs
          call interface_used%print_binary(
     $         filename_x_map,
     $         filename_y_map,
     $         filename_nodes,
     $         filename_grdpt,
     $         filename_sizes,
     $         filename_nbsbl)
          
        end subroutine print_output


        subroutine test_print_netcdf(interface_used, index)

          implicit none

          class(bf_interface), intent(in) :: interface_used
          integer            , intent(in) :: index

          
          character*(*), dimension(4) :: name_var
          character*(*), dimension(4) :: longname_var
          character*(*), dimension(4) :: unit_var

          real(rkind) :: time

          parameter (name_var = [
     $         'mass      ',
     $         'momentum-x',
     $         'momentum-y',
     $         'energy    '])

          parameter (longname_var = [
     $         'non-dimensional mass      ',
     $         'non-dimensional momentum-x',
     $         'non-dimensional momentum-y',
     $         'non-dimensional energy    '])

          parameter (unit_var = [
     $         '(kg/m3)/(kg/m3)        ',
     $         '(kg/(m2.s))/(kg/(m2.s))',
     $         '(kg/(m2.s))/(kg/(m2.s))',
     $         '(J/m3)/(J/m3)          '])

          time = 2.0

          call interface_used%print_netcdf(
     $         index,
     $         name_var,
     $         longname_var,
     $         unit_var,
     $         time)

        end subroutine test_print_netcdf

      
        !test : resolve_bc_overlap_conflicts
        function test_resolve_bc_overlap_conflicts1(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_interface)                   :: bf_interface_used
          integer, dimension(:,:), allocatable :: grdpts_id
          integer(ikind), dimension(2,2)       :: bf_alignment
          type(bf_sublayer), pointer           :: bf_layer_N_ptr
          type(bf_sublayer), pointer           :: bf_layer_E_ptr
          real(rkind), dimension(nx)           :: interior_x_map
          real(rkind), dimension(ny)           :: interior_y_map
          real(rkind), dimension(nx,ny,ne)     :: interior_nodes
          type(bf_mainlayer), pointer          :: mainlayer_ptr


          !test case:
          !2 buffer layers: N and E at the interface NE
          !the North buffer layer should be increased to
          !be able to compute all the boundary grid points
          !        
          !      North buffer layer
          !                        
          !             nx-2___       ___ nx+1
          !                    |     |  ___ nx+2
          !             nx-3_  |     | |  ___ nx+3
          !                  | |     | | |  __  nx+4
          !                  | |     | | | |
          !        _________________________ 
          !       |3 3 3 3 3 3 3 3 3 3 3 3 3|___
          ! ny   -|3 2 2 2 2 2 2 2 2 2 2 2 3|3  |
          ! ny-1 -|2 2_ _ _ _ _ _ _ _ _ _2_2|3_3|
          ! ny-2 -          |   |          2|2 3|
          ! ny-3 - _________|___|__________ |2 3| East buffer layer
          !                 |   |           |2 3|
          !        interior |1 1|2 2 2 2 2 2|2 3|
          !                 |1 1|2 3 3 3 3 3|3 3|
          !                 -------------
          !------------------------------------------------
          call bf_interface_used%ini(interior_x_map,interior_y_map)


          !initialize the N buffer layer
          allocate(grdpts_id(13,5))

          grdpts_id = reshape((/
     $         1,1,1,1,1,1,1,1,1,1,1,1,1,
     $         1,1,1,1,1,1,1,1,1,1,1,1,2,
     $         1,1,1,1,1,1,1,1,1,1,1,2,2,
     $         2,2,2,2,2,2,2,2,2,2,2,2,3,
     $         3,3,3,3,3,3,3,3,3,3,3,3,3/),
     $         (/13,5/))

          bf_alignment(1,2) = nx+2
          bf_alignment(1,1) = bf_alignment(1,2) - size(grdpts_id,1) + (2*bc_size+1)
          bf_alignment(2,1) = ny-1
          bf_alignment(2,2) = bf_alignment(2,1) + size(grdpts_id,2) - (2*bc_size+1)
          
          if(detailled) then
             print '(''bf_alignment_N: '',4I4)', bf_alignment
          end if

          bf_layer_N_ptr => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          call bf_layer_N_ptr%set_grdpts_id(grdpts_id)


          !initialize the E buffer layer
          allocate(grdpts_id(10,7))

          grdpts_id = reshape((/
     $         1,1,2,3,3,3,3,3,3,3,
     $         1,1,2,2,2,2,2,2,2,3,
     $         1,1,1,1,1,1,1,1,2,3,
     $         1,1,1,1,1,1,1,1,2,3,
     $         1,1,1,1,1,1,1,2,2,3,
     $         1,1,1,1,1,1,2,2,3,3,
     $         2,2,2,2,2,2,2,3,3,0/),
     $         (/10,7/))

          bf_alignment(1,1) = nx-1
          bf_alignment(1,2) = bf_alignment(1,1) + size(grdpts_id,1) - (2*bc_size+1)
          bf_alignment(2,2) = ny-2
          bf_alignment(2,1) = bf_alignment(2,2) - size(grdpts_id,2) + (2*bc_size+1)

          if(detailled) then
             print '(''bf_alignment_E: '',4I4)', bf_alignment
          end if

          bf_layer_E_ptr => bf_interface_used%allocate_sublayer(
     $         E,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          call bf_layer_E_ptr%set_grdpts_id(grdpts_id)

          call bf_interface_used%resolve_bc_overlap_conflicts(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes)


          !test whether the N buffer layer has been correctly updated
          mainlayer_ptr  => bf_interface_used%get_mainlayer(N)
          bf_layer_N_ptr => mainlayer_ptr%get_head_sublayer()
          bf_alignment   = bf_layer_N_ptr%get_alignment_tab()

          if(detailled) then
             print '(''bf_alignment_N after: '',4I4)', bf_alignment
          end if

          test_validated = bf_alignment(1,1).eq.(nx-6)
          test_validated = test_validated.and.(bf_alignment(1,2).eq.(nx+3))
          test_validated = test_validated.and.(bf_alignment(2,1).eq.(ny-1))
          test_validated = test_validated.and.(bf_alignment(2,2).eq.(ny-1))


          !test whether the N buffer layer has been correctly updated
          mainlayer_ptr  => bf_interface_used%get_mainlayer(E)
          bf_layer_E_ptr => mainlayer_ptr%get_head_sublayer()
          bf_alignment   = bf_layer_E_ptr%get_alignment_tab()

          if(detailled) then
             print '(''bf_alignment_E after: '',4I4)', bf_alignment
          end if

          test_validated = test_validated.and.(bf_alignment(1,1).eq.(nx-1))
          test_validated = test_validated.and.(bf_alignment(1,2).eq.(nx+4))
          test_validated = test_validated.and.(bf_alignment(2,1).eq.(ny-4))
          test_validated = test_validated.and.(bf_alignment(2,2).eq.(ny-2))

        end function test_resolve_bc_overlap_conflicts1


        !test : resolve_bc_overlap_conflicts
        function test_resolve_bc_overlap_conflicts2(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_interface)                   :: bf_interface_used
          integer, dimension(:,:), allocatable :: grdpts_id
          integer(ikind), dimension(2,2)       :: bf_alignment
          type(bf_sublayer), pointer           :: bf_layer_S_ptr
          type(bf_sublayer), pointer           :: bf_layer_W_ptr
          real(rkind), dimension(nx)           :: interior_x_map
          real(rkind), dimension(ny)           :: interior_y_map
          real(rkind), dimension(nx,ny,ne)     :: interior_nodes
          type(bf_mainlayer), pointer          :: mainlayer_ptr


          !test case:
          !2 buffer layers: S and W at the interface SW
          !the South buffer layer should be increased to
          !be able to compute all the boundary grid points
          !                              
          !                    -------------------------------
          !                    |3 3|3 3 3 3|3 2|1 1| 
          !  West buffer layer |3 2|2 2 2 2|2 2|1 1| interior
          !                    |3 2|       | _ | _ | _ _ _ _ _
          !             4-     |3 2|       |   |   | 
          !             3-     |3 2|2_ _ _ | _ | _ | _ _ _ _ _
          !             2-     |3 3|2 2    |   |   |      2 2
          !             1-     |__3|3 2 2 2|2_2|2_2|2_2_2_2_3______
          !                        |3 3 3 3 3 3 3 3|3 3 3 3 3
          !                        --------------------------
          !                         | | | | | | | | 
          !                        -3 |-1 | 1 | 3 |
          !                          -2   0   2   4
          !
          !                             South buffer layer
          !------------------------------------------------
          call bf_interface_used%ini(interior_x_map,interior_y_map)


          !initialize the S buffer layer
          allocate(grdpts_id(13,5))

          grdpts_id = reshape((/
     $         3,3,3,3,3,3,3,3,3,3,3,3,3,
     $         3,2,2,2,2,2,2,2,2,2,2,2,3,
     $         2,2,1,1,1,1,1,1,1,1,1,2,2,
     $         2,1,1,1,1,1,1,1,1,1,1,1,1,
     $         1,1,1,1,1,1,1,1,1,1,1,1,1/),
     $         (/13,5/))

          bf_alignment(1,1) = -1
          bf_alignment(1,2) = bf_alignment(1,1) + size(grdpts_id,1) - (2*bc_size+1)
          bf_alignment(2,1) = bc_size
          bf_alignment(2,2) = bf_alignment(2,1) + size(grdpts_id,2) - (2*bc_size+1)
          
          if(detailled) then
             print '(''bf_alignment_S: '',4I4)', bf_alignment
          end if

          bf_layer_S_ptr => bf_interface_used%allocate_sublayer(
     $         S,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          call bf_layer_S_ptr%set_grdpts_id(grdpts_id)


          !initialize the W buffer layer
          allocate(grdpts_id(10,7))

          grdpts_id = reshape((/
     $         0,3,3,2,2,2,2,2,2,2,
     $         3,3,2,2,1,1,1,1,1,1,
     $         3,2,2,1,1,1,1,1,1,1,
     $         3,2,1,1,1,1,1,1,1,1,
     $         3,2,1,1,1,1,1,1,1,1,
     $         3,2,2,2,2,2,2,2,1,1,
     $         3,3,3,3,3,3,3,2,1,1/),
     $         (/10,7/))

          bf_alignment(1,2) = bc_size
          bf_alignment(1,1) = bf_alignment(1,2) - size(grdpts_id,1) + (2*bc_size+1)

          bf_alignment(2,1) = bc_size+1
          bf_alignment(2,2) = bf_alignment(2,1) + size(grdpts_id,2) - (2*bc_size+1)

          if(detailled) then
             print '(''bf_alignment_W: '',4I4)', bf_alignment
          end if

          bf_layer_W_ptr => bf_interface_used%allocate_sublayer(
     $         W,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          call bf_layer_W_ptr%set_grdpts_id(grdpts_id)

          call bf_interface_used%resolve_bc_overlap_conflicts(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes)


          !test whether the S buffer layer has been correctly updated
          mainlayer_ptr  => bf_interface_used%get_mainlayer(S)
          bf_layer_S_ptr => mainlayer_ptr%get_head_sublayer()
          bf_alignment   = bf_layer_S_ptr%get_alignment_tab()

          if(detailled) then
             print '(''bf_alignment_S after: '',4I4)', bf_alignment
          end if

          test_validated = bf_alignment(1,1).eq.(-2)
          test_validated = test_validated.and.(bf_alignment(1,2).eq.7)
          test_validated = test_validated.and.(bf_alignment(2,1).eq.2)
          test_validated = test_validated.and.(bf_alignment(2,2).eq.2)


          !test whether the W buffer layer has been correctly updated
          mainlayer_ptr  => bf_interface_used%get_mainlayer(W)
          bf_layer_W_ptr => mainlayer_ptr%get_head_sublayer()
          bf_alignment   = bf_layer_W_ptr%get_alignment_tab()

          if(detailled) then
             print '(''bf_alignment_W after: '',4I4)', bf_alignment
          end if

          test_validated = test_validated.and.(bf_alignment(1,1).eq.(bc_size+1))
          test_validated = test_validated.and.(bf_alignment(1,2).eq.(bc_size+3))
          test_validated = test_validated.and.(bf_alignment(2,1).eq.(-3))
          test_validated = test_validated.and.(bf_alignment(2,2).eq.2)

        end function test_resolve_bc_overlap_conflicts2

      end program test_bf_interface_prog
