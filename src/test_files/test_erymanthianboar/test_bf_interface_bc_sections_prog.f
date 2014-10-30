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


        implicit none

        type(bf_interface)                  :: interface_tested
        real(rkind)   , dimension(nx)       :: x_map
        real(rkind)   , dimension(ny)       :: y_map
        real(rkind)   , dimension(nx,ny,ne) :: nodes
        integer       , dimension(nx,ny)    :: grdpts_id
        integer       , dimension(2,2)      :: alignment
        integer                             :: mainlayer_id
        type(bf_sublayer), pointer          :: added_sublayer
        type(bf_sublayer), pointer          :: bf_sublayer_merged1
        type(bf_sublayer), pointer          :: bf_sublayer_merged2
        type(bf_sublayer), pointer          :: bf_sublayer_reallocated
        type(bf_sublayer), pointer          :: merged_sublayer

        integer(ikind) :: i,j
        integer        :: index

        index = 1

        !initialize the nodes and print them
        call ini_x_map(x_map)
        call ini_y_map(y_map)

        do j=1, ny
           do i=1, nx
              nodes(i,j,:)=1.0
           end do
        end do        

        call ini_grdpts_id(grdpts_id)
        call print_interior_data(
     $       x_map,
     $       y_map,
     $       nodes,
     $       grdpts_id, 
     $       'interior_x_map.dat',
     $       'interior_y_map.dat',
     $       'interior_nodes.dat',
     $       'interior_grdpts_id.dat',
     $       'interior_sizes.dat')
        

        !initialize the interface
        call interface_tested%ini()

        
        !add the N_E buffer layer
        call get_alignment(1, alignment, mainlayer_id)
        added_sublayer => interface_tested%allocate_sublayer(
     $       mainlayer_id, x_map, y_map, nodes, alignment)

        call test_bf_integration_borders(
     $       interface_tested,
     $       index)

        
        !add the W buffer layer under the N_W corner
        call get_alignment(4, alignment, mainlayer_id)
        added_sublayer => interface_tested%allocate_sublayer(
     $       mainlayer_id, x_map, y_map, nodes, alignment)

        call test_bf_integration_borders(
     $       interface_tested,
     $       index)
        

        !add the E buffer layer under the N_E corner
        call get_alignment(3, alignment, mainlayer_id)
        bf_sublayer_merged1 => interface_tested%allocate_sublayer(
     $       mainlayer_id, x_map, y_map, nodes, alignment)

        call test_bf_integration_borders(
     $       interface_tested,
     $       index)


        !add the N_W buffer layer
        call get_alignment(2, alignment, mainlayer_id)
        added_sublayer => interface_tested%allocate_sublayer(
     $       mainlayer_id, x_map, y_map, nodes, alignment)

        call test_bf_integration_borders(
     $       interface_tested,
     $       index)


        !add the E buffer layer
        call get_alignment(5, alignment, mainlayer_id)
        bf_sublayer_merged2 => interface_tested%allocate_sublayer(
     $       mainlayer_id, x_map, y_map, nodes, alignment)

        call test_bf_integration_borders(
     $       interface_tested,
     $       index)


        !add the W buffer layer
        call get_alignment(6, alignment, mainlayer_id)
        bf_sublayer_reallocated => interface_tested%allocate_sublayer(
     $       mainlayer_id, x_map, y_map, nodes, alignment)

        call test_bf_integration_borders(
     $       interface_tested,
     $       index)


        !add the S_E buffer layer
        call get_alignment(7, alignment, mainlayer_id)
        added_sublayer => interface_tested%allocate_sublayer(
     $       mainlayer_id, x_map, y_map, nodes, alignment)

        call test_bf_integration_borders(
     $       interface_tested,
     $       index)


        !add the S_W buffer layer
        call get_alignment(8, alignment, mainlayer_id)
        added_sublayer => interface_tested%allocate_sublayer(
     $       mainlayer_id, x_map, y_map, nodes, alignment)

        call test_bf_integration_borders(
     $       interface_tested,
     $       index)

        !merge the E buffer layers
        call get_alignment(9, alignment, mainlayer_id)
        merged_sublayer => interface_tested%merge_sublayers(
     $       bf_sublayer_merged1, bf_sublayer_merged2,
     $       x_map, y_map, nodes, alignment)

        call test_bf_integration_borders(
     $       interface_tested,
     $       index)
        

        !reallocate the W buffer layer
        call get_alignment(10, alignment, mainlayer_id)
        call interface_tested%reallocate_sublayer(
     $       bf_sublayer_reallocated,
     $       x_map, y_map, nodes, alignment)

        call test_bf_integration_borders(
     $       interface_tested,
     $       index)

        !add more sublayers on N,S
        do i=11,16
           call get_alignment(i, alignment, mainlayer_id)
           added_sublayer => interface_tested%allocate_sublayer(
     $          mainlayer_id, x_map, y_map, nodes, alignment)
        end do

        call test_bf_integration_borders(
     $       interface_tested,
     $       index)


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


        !test the determination of the integration borders
        subroutine test_bf_integration_borders(
     $     interface_used,
     $     index)
        
          implicit none
          
          class(bf_interface)             , intent(inout) :: interface_used
          integer                         , intent(inout) :: index


          !1) reinitialize the nodes of the boundary layers
          call reinitialize_nodes(interface_used)

          !2) colorize the integration borders in the buffer layers
          call colorize_integration_borders(interface_used)

          !3) print interface
          call print_output(interface_used, index)

          index = index+1

        end subroutine test_bf_integration_borders


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

      end program test_bf_interface_prog
