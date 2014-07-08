      !the idea with this test file is to recreate the
      !same results as presented in test_bf_layer_prog
      !but using the interface_abstract as the object
      !encapsulating the buffer layer objects      
      program test_bf_interface_prog

        use bf_mainlayer_class  , only : bf_mainlayer
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
        use sbf_list_class      , only : sbf_list


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
        type(bf_sublayer), pointer          :: merged_sublayer
        real(rkind)                         :: scale


        integer :: i, index


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
        merged_sublayer => interface_tested%merge_sublayers(
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

        
        !add more sublayers on N,S
        do i=11,16
           call get_alignment(i, alignment, mainlayer_id)
           added_sublayer => interface_tested%allocate_sublayer(
     $          mainlayer_id, nodes, alignment)
        end do

        call interface_tested%print_binary(
     $       'nodes6.dat',
     $       'grdpt_id6.dat',
     $       'sizes6.dat',
     $       '6.dat')

        index = 7


        !test bf_depends_on_neighbors
        !-------------------------------------------------------
        !for each sublayer in each mainlayer, test if the buffer layer
        !depends on neighbors
        call test_bf_layer_depends_on_neighbors(interface_tested)
        

        !test get_nbf_layers_sharing_grdpts_with
        !-------------------------------------------------------
        !for each sublayer in each mainlayer, test if the buffer layer
        !has neighbor dependencies and colorize the neighbors
        call test_get_nbf_layers_sharing_grdpts_with(interface_tested, index)

        
        !test if the neighboring buffer layers will be removed
        !-------------------------------------------------------
        call test_does_a_neighbor_remains(interface_tested)


        !test the removal of sublayer
        !----------------------------
        call test_remove_sublayer(interface_tested, merged_sublayer, index)


        !re-test if the neighboring buffer layers will be removed
        !--------------------------------------------------------
        call test_does_a_neighbor_remains(interface_tested)


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
                nodes(i,j,1)  = color
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
        subroutine test_remove_sublayer(interface_used, sublayer_ptr, index)
        
          implicit none

          class(bf_interface)       , intent(inout) :: interface_used
          type(bf_sublayer), pointer, intent(inout) :: sublayer_ptr
          integer                   , intent(inout) :: index

          print '(''test remove_sublayer'')'
          print '(''--------------------'')'

          call interface_used%remove_sublayer(sublayer_ptr)

          call print_output(interface_used, index)

          index = index+1

          print '()'

        end subroutine test_remove_sublayer

      
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
          write(filename_nodes, format_nodes) 'nodes', index, '.dat'
          write(filename_grdpt, format_grdpt) 'grdpt_id', index, '.dat'
          write(filename_sizes, format_nodes) 'sizes', index, '.dat'
          write(filename_nbsbl, format_nbsbl) index, '.dat'


          !write the outputs
          call interface_used%print_binary(
     $         filename_nodes,
     $         filename_grdpt,
     $         filename_sizes,
     $         filename_nbsbl)
          
        end subroutine print_output

      end program test_bf_interface_prog
