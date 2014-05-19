      !organization of the test cases for the bf_interface
      module test_cases_interface_module

        use bf_sublayer_class  , only : bf_sublayer
        use bf_interface_class , only : bf_interface
        use parameters_bf_layer, only : interior_pt, bc_interior_pt,
     $                                       bc_pt, no_pt
        use parameters_constant, only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input   , only : nx,ny,ne, bc_size
        use parameters_kind    , only : ikind, rkind


        implicit none

        private
        public :: ini_interface


        contains


        subroutine ini_interface(
     $       interface_used,
     $       test_case_nb,
     $       nodes,
     $       corner_id,
     $       corner_order,
     $       bf_corner_distance,
     $       over_allocated,
     $       mainlayer_id,
     $       nb_sublayers,
     $       relative_distances,
     $       relative_sizes,
     $       relative_distance,
     $       relative_size)

          implicit none 

          class(bf_interface)              , intent(inout) :: interface_used
          integer                          , intent(in)    :: test_case_nb
          real(rkind), dimension(nx,ny,ne) , intent(in)    :: nodes
          integer                , optional, intent(in)    :: corner_id
          integer                , optional, intent(in)    :: corner_order
          integer                , optional, intent(in)    :: bf_corner_distance
          integer                , optional, intent(in)    :: over_allocated
          integer                , optional, intent(in)    :: mainlayer_id
          integer                , optional, intent(in)    :: nb_sublayers
          integer, dimension(:)  , optional, intent(in)    :: relative_distances
          integer, dimension(:,:), optional, intent(in)    :: relative_sizes
          integer                , optional, intent(in)    :: relative_distance
          integer                , optional, intent(in)    :: relative_size

          
          integer                              :: corner_id_i
          integer                              :: corner_order_i
          integer                              :: bf_corner_distance_i
          integer                              :: over_allocated_i
          integer                              :: mainlayer_id_i
          integer                              :: nb_sublayers_i
          integer, dimension(:)  , allocatable :: relative_distances_i
          integer, dimension(:,:), allocatable :: relative_sizes_i
          integer                              :: relative_distance_i
          integer                              :: relative_size_i


          integer :: i

          !preprocess the optional arguments
          if(present(corner_id)) then
             corner_id_i = corner_id
          else
             corner_id_i = 1
          end if

          if(present(corner_order)) then
             corner_order_i = corner_order
          else
             corner_order_i = 1
          end if

          if(present(bf_corner_distance)) then
             bf_corner_distance_i = bf_corner_distance
          else
             bf_corner_distance_i = 0
          end if

          if(present(over_allocated)) then
             over_allocated_i = over_allocated
          else
             over_allocated_i = 0
          end if

          if(present(mainlayer_id)) then
             mainlayer_id_i = mainlayer_id
          else
             mainlayer_id_i = N
          end if          

          if(present(nb_sublayers)) then
             nb_sublayers_i = nb_sublayers
          else
             nb_sublayers_i = 2
          end if

          allocate(relative_distances_i(nb_sublayers_i))
          if(present(relative_distances)) then
             relative_distances_i = relative_distances
          else
             do i=1, nb_sublayers_i
                relative_distances_i(i) = 0
             end do
          end if

          allocate(relative_sizes_i(2,nb_sublayers_i))
          if(present(relative_sizes)) then
             relative_sizes_i = relative_sizes
          else
             do i=1, nb_sublayers_i
                relative_sizes_i(1,i) = 0
                relative_sizes_i(2,i) = 0
             end do
          end if

          if(present(relative_distance)) then
             relative_distance_i = relative_distance
          else
             relative_distance_i = 0
          end if

          if(present(relative_size)) then
             relative_size_i = relative_size
          else
             relative_size_i = 0
          end if


          !choose the test case
          select case(test_case_nb)
            case(0)
               call ini_interface_testcase0(interface_used)
            case(1)
               call ini_interface_testcase1(
     $              interface_used,
     $              nodes,
     $              corner_id_i,
     $              corner_order_i,
     $              bf_corner_distance_i,
     $              over_allocated_i)
            case(2)
               call ini_interface_testcase2(
     $              interface_used,
     $              nodes,
     $              corner_id_i,
     $              bf_corner_distance_i,
     $              over_allocated_i)

            case(21)
               call ini_interface_testcase21(
     $              interface_used,
     $              nodes,
     $              bf_corner_distance_i,
     $              over_allocated_i)

            case(3)
               call ini_interface_testcase3(
     $              interface_used,
     $              nodes,
     $              mainlayer_id_i,
     $              nb_sublayers_i,
     $              relative_distances_i,
     $              relative_sizes_i)

            case(4)
               call ini_interface_testcase4(
     $              interface_used,
     $              nodes,
     $              relative_distance_i,
     $              relative_size_i,
     $              over_allocated_i)

            case default
               print '(''test_cases_interface_module'')'
               print '(''ini_interface'')'
               print '(''test case not yet implemented'')'
               print '(''test case ID: '',I2)', test_case_nb
               stop 'the program will stop now'
          end select 

        end subroutine ini_interface


        subroutine ini_interface_testcase0(interface_used)

          implicit none

          class(bf_interface), intent(inout) :: interface_used

          print '(''********************************************'')'
          print '(''0: interface initialization with no sublayer'')'
          print '(''********************************************'')'

          call interface_used%ini()

        end subroutine ini_interface_testcase0


        subroutine ini_interface_testcase1(
     $     interface_used,
     $     nodes,
     $     corner_id,
     $     corner_order,
     $     bf_corner_distance,
     $     over_allocated)

          implicit none

          class(bf_interface)       , intent(inout) :: interface_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer                         , intent(in)    :: corner_id
          integer                         , intent(in)    :: corner_order
          integer                         , intent(in)    :: bf_corner_distance
          integer                         , intent(in)    :: over_allocated
          

          type(bf_sublayer), pointer     :: added_sublayer
          integer(ikind), dimension(2,2) :: alignment
          logical       , dimension(4)   :: neighbors

          if((corner_order.ne.1).and.(corner_order.ne.2)) then
             print '(''test_cases_interface_module'')'
             print '(''ini_interface_test_case1'')'
             stop 'chose corner_order correctly'
          end if

          print '(''************************************'')'
          print '(''1: interface with one sublayer      '')'
          print '('' : close to the corner chosen       '')'
          print '(''************************************'')'

          !initialization
          call interface_used%ini()

          neighbors = [.false.,.false.,.false.,.false.]

          select case(corner_id)

            case(N_E)

               !add north buffer layer close to the N_E corner
               if(corner_order.eq.1) then
                  alignment(1,1) = nx - bc_size - bf_corner_distance - 4
                  alignment(1,2) = nx - bc_size - bf_corner_distance
                  alignment(2,1) = ny-1
                  alignment(2,2) = ny-1

                  added_sublayer => interface_used%add_sublayer(
     $                 N, nodes, alignment)

                  call over_allocate_north(added_sublayer, over_allocated, nodes)
                  
               !add east buffer layer close to the N_E corner
               else
                  alignment(1,1) = nx - 1
                  alignment(1,2) = nx - 1
                  alignment(2,1) = ny - bc_size - bf_corner_distance - 4
                  alignment(2,2) = ny - bc_size - bf_corner_distance

                  added_sublayer => interface_used%add_sublayer(
     $                 E, nodes, alignment)
                  
                  call over_allocate_east(added_sublayer, over_allocated, nodes)

               end if

               
            case(N_W)

               !add north buffer layer close to the N_W corner
               if(corner_order.eq.2) then
                  alignment(1,1) = 1 + bc_size + bf_corner_distance
                  alignment(1,2) = 1 + bc_size + bf_corner_distance + 4
                  alignment(2,1) = ny-1
                  alignment(2,2) = ny-1

                  added_sublayer => interface_used%add_sublayer(
     $                 N, nodes, alignment)
                  
                  call over_allocate_north(added_sublayer, over_allocated, nodes)

               !add west buffer layer close to the N_W corner
               else
                  alignment(1,1) = nx - 1
                  alignment(1,2) = nx - 1
                  alignment(2,1) = ny - bc_size - bf_corner_distance - 4
                  alignment(2,2) = ny - bc_size - bf_corner_distance

                  added_sublayer => interface_used%add_sublayer(
     $                 W, nodes, alignment)
                  
                  call over_allocate_west(added_sublayer, over_allocated, nodes)

               end if


            case(S_E)

               !add south buffer layer close to the S_E corner
               if(corner_order.eq.2) then
                  alignment(1,1) = nx - bc_size - bf_corner_distance - 4
                  alignment(1,2) = nx - bc_size - bf_corner_distance
                  alignment(2,1) = 1
                  alignment(2,2) = 1

                  added_sublayer => interface_used%add_sublayer(
     $                 S, nodes, alignment)
                  
                  call over_allocate_south(added_sublayer, over_allocated, nodes)

               !add east buffer layer close to the S_E corner
               else
                  alignment(1,1) = nx - 1
                  alignment(1,2) = nx - 1
                  alignment(2,1) = 1 + bc_size + bf_corner_distance
                  alignment(2,2) = 1 + bc_size + bf_corner_distance + 4

                  added_sublayer => interface_used%add_sublayer(
     $                 E, nodes, alignment)
                  
                  call over_allocate_east(added_sublayer, over_allocated, nodes)

               end if


            case(S_W)

               !add south buffer layer close to the S_W corner
               if(corner_order.eq.1) then
                  alignment(1,1) = 1 + bc_size + bf_corner_distance
                  alignment(1,2) = 1 + bc_size + bf_corner_distance + 4
                  alignment(2,1) = 1
                  alignment(2,2) = 1

                  added_sublayer => interface_used%add_sublayer(
     $                 S, nodes, alignment)
                  
                  call over_allocate_south(added_sublayer, over_allocated, nodes)

               !add west buffer layer close to the S_W corner
               else
                  alignment(1,1) = nx - 1
                  alignment(1,2) = nx - 1
                  alignment(2,1) = 1 + bc_size + bf_corner_distance
                  alignment(2,2) = 1 + bc_size + bf_corner_distance + 4

                  added_sublayer => interface_used%add_sublayer(
     $                 W, nodes, alignment)
                  
                  call over_allocate_west(added_sublayer, over_allocated, nodes)

               end if

            end select

        end subroutine ini_interface_testcase1


        subroutine ini_interface_testcase2(
     $     interface_used,
     $     nodes,
     $     corner_id,
     $     bf_corner_distance,
     $     over_allocated)

          implicit none

          class(bf_interface)             , intent(inout) :: interface_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer                         , intent(in)    :: corner_id
          integer                         , intent(in)    :: bf_corner_distance
          integer                         , intent(in)    :: over_allocated

          print '(''************************************'')'
          print '(''2: interface with two sublayers     '')'
          print '('' : close to the corner chosen       '')'
          print '(''************************************'')'

          !initialization
          call interface_used%ini()

          !add two buffer layers around a corner
          call add_bf_layers_to_corner(
     $         interface_used,
     $         nodes,
     $         corner_id,
     $         bf_corner_distance,
     $         over_allocated)

        end subroutine ini_interface_testcase2


        subroutine ini_interface_testcase21(
     $     interface_used,
     $     nodes,
     $     bf_corner_distance,
     $     over_allocated)

          implicit none

          class(bf_interface)             , intent(inout) :: interface_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer                         , intent(in)    :: bf_corner_distance
          integer                         , intent(in)    :: over_allocated

          integer, dimension(4) :: corner_tab
          integer               :: i

          print '(''************************************'')'
          print '(''21: interface with two sublayers     '')'
          print '(''  : at each corner                   '')'
          print '(''************************************'')'

          corner_tab = [N_E,N_W,S_E,S_W]

          !initialization
          call interface_used%ini()

          do i=1, size(corner_tab,1)

             !add two buffer layers around each corner
             call add_bf_layers_to_corner(
     $            interface_used,
     $            nodes,
     $            corner_tab(i),
     $            bf_corner_distance,
     $            over_allocated)

          end do

        end subroutine ini_interface_testcase21


        subroutine add_bf_layers_to_corner(
     $     interface_used,
     $     nodes,
     $     corner_id,
     $     bf_corner_distance,
     $     over_allocated)

          implicit none

          class(bf_interface)             , intent(inout) :: interface_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer                         , intent(in)    :: corner_id
          integer                         , intent(in)    :: bf_corner_distance
          integer                         , intent(in)    :: over_allocated

          type(bf_sublayer), pointer     :: added_sublayer
          integer(ikind), dimension(2,2) :: alignment
          logical       , dimension(4)   :: neighbors

          neighbors = [.false.,.false.,.false.,.false.]

          select case(corner_id)

            case(N_E)               

               !add north buffer layer close to the N_E corner
               alignment(1,1) = nx - bc_size - bf_corner_distance - 4
               alignment(1,2) = nx - bc_size - bf_corner_distance
               alignment(2,1) = ny-1
               alignment(2,2) = ny-1
               
               added_sublayer => interface_used%add_sublayer(
     $              N, nodes, alignment)

               call over_allocate_north(added_sublayer, over_allocated, nodes)
               
               !add east buffer layer close to the N_E corner
               alignment(1,1) = nx - 1
               alignment(1,2) = nx - 1
               alignment(2,1) = ny - bc_size - bf_corner_distance - 4
               alignment(2,2) = ny - bc_size - bf_corner_distance
               
               added_sublayer => interface_used%add_sublayer(
     $              E, nodes, alignment)
               
               call over_allocate_east(added_sublayer, over_allocated, nodes)

            case(N_W)

               !add north buffer layer close to the N_W corner
               alignment(1,1) = 1 + bc_size + bf_corner_distance
               alignment(1,2) = 1 + bc_size + bf_corner_distance + 4
               alignment(2,1) = ny-1
               alignment(2,2) = ny-1
               
               added_sublayer => interface_used%add_sublayer(
     $              N, nodes, alignment)

               call over_allocate_north(added_sublayer, over_allocated, nodes)
               
               !add west buffer layer close to the N_W corner
               alignment(1,1) = nx - 1
               alignment(1,2) = nx - 1
               alignment(2,1) = ny - bc_size - bf_corner_distance - 4
               alignment(2,2) = ny - bc_size - bf_corner_distance
               
               added_sublayer => interface_used%add_sublayer(
     $              W, nodes, alignment)

               call over_allocate_west(added_sublayer, over_allocated, nodes)
               
            case(S_E)

               !add south buffer layer close to the S_E corner
               alignment(1,1) = nx - bc_size - bf_corner_distance - 4
               alignment(1,2) = nx - bc_size - bf_corner_distance
               alignment(2,1) = 1
               alignment(2,2) = 1
               
               added_sublayer => interface_used%add_sublayer(
     $              S, nodes, alignment)

               call over_allocate_south(added_sublayer, over_allocated, nodes)
               
               !add east buffer layer close to the S_E corner
               alignment(1,1) = nx - 1
               alignment(1,2) = nx - 1
               alignment(2,1) = 1 + bc_size + bf_corner_distance
               alignment(2,2) = 1 + bc_size + bf_corner_distance + 4
               
               added_sublayer => interface_used%add_sublayer(
     $              E, nodes, alignment)

               call over_allocate_east(added_sublayer, over_allocated, nodes)

            case(S_W)

               !add north buffer layer close to the S_W corner
               alignment(1,1) = 1 + bc_size + bf_corner_distance
               alignment(1,2) = 1 + bc_size + bf_corner_distance + 4
               alignment(2,1) = 1
               alignment(2,2) = 1
               
               added_sublayer => interface_used%add_sublayer(
     $              S, nodes, alignment)
               
               call over_allocate_south(added_sublayer, over_allocated, nodes)


               !add west buffer layer close to the N_W corner
               alignment(1,1) = nx - 1
               alignment(1,2) = nx - 1
               alignment(2,1) = 1 + bc_size + bf_corner_distance
               alignment(2,2) = 1 + bc_size + bf_corner_distance + 4
               
               added_sublayer => interface_used%add_sublayer(
     $              W, nodes, alignment)

               call over_allocate_west(added_sublayer, over_allocated, nodes)

            end select

        end subroutine add_bf_layers_to_corner


        subroutine ini_interface_testcase3(
     $     interface_used,
     $     nodes,
     $     mainlayer_id,
     $     nb_sublayers,
     $     relative_distance,
     $     relative_sizes)

          implicit none

          class(bf_interface)       , intent(inout) :: interface_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer                         , intent(in)    :: mainlayer_id
          integer                         , intent(in)    :: nb_sublayers
          integer, dimension(:)           , intent(in)    :: relative_distance
          integer, dimension(:,:)         , intent(in)    :: relative_sizes
          

          type(bf_sublayer), pointer     :: added_sublayer
          integer(ikind), dimension(2,2) :: alignment
          logical       , dimension(4)   :: neighbors

          integer :: i

          print '(''************************************'')'
          print '(''3: interface with two sublayers     '')'
          print '('' : on the same mainlayer            '')'
          print '(''************************************'')'

          !initialization
          call interface_used%ini()

          neighbors = [.false.,.false.,.false.,.false.]

          select case(mainlayer_id)

            case(N)

               !initialize the alignment
               alignment(1,2) = - bc_size
               alignment(2,1) = ny-1
               alignment(2,2) = ny

               !add the sublayers and print the content
               do i=1, nb_sublayers

                  alignment(1,1) = alignment(1,2)+2*bc_size+1
                  alignment(1,1) = alignment(1,1)+relative_distance(i)
                  alignment(1,2) = alignment(1,1)+relative_sizes(1,i)

                  added_sublayer => interface_used%add_sublayer(
     $                 mainlayer_id, nodes, alignment)

                  if(relative_sizes(2,i).gt.0) then
                     call over_allocate_north(
     $                    added_sublayer,
     $                    relative_sizes(2,i),
     $                    nodes)
                  end if

               end do               

            case(S)

               !initialize the alignment
               alignment(1,2) = - bc_size
               alignment(2,1) = 1
               alignment(2,2) = bc_size

               !add the sublayers and print the content
               do i=1, nb_sublayers

                  alignment(1,1) = alignment(1,2)+2*bc_size+1
                  alignment(1,1) = alignment(1,1)+relative_distance(i)
                  alignment(1,2) = alignment(1,1)+relative_sizes(1,i)

                  added_sublayer => interface_used%add_sublayer(
     $                 mainlayer_id, nodes, alignment)

                  if(relative_sizes(2,i).gt.0) then
                     call over_allocate_south(
     $                    added_sublayer,
     $                    relative_sizes(2,i),
     $                    nodes)
                  end if

               end do

            case(E)

               !initialize the alignment
               alignment(1,1) = nx-1
               alignment(1,2) = nx
               alignment(2,2) = - bc_size

               !add the sublayers and print the content
               do i=1, nb_sublayers

                  alignment(2,1) = alignment(2,2)+2*bc_size+1
                  alignment(2,1) = alignment(2,1)+relative_distance(i)
                  alignment(2,2) = alignment(2,1)+relative_sizes(1,i)

                  added_sublayer => interface_used%add_sublayer(
     $                 mainlayer_id, nodes, alignment)

                  if(relative_sizes(2,i).gt.0) then
                     call over_allocate_east(
     $                    added_sublayer,
     $                    relative_sizes(2,i),
     $                    nodes)
                  end if

               end do

            case(W)

               !initialize the alignment
               alignment(1,1) = 1
               alignment(1,2) = bc_size
               alignment(2,2) = - bc_size

               !add the sublayers and print the content
               do i=1, nb_sublayers

                  alignment(2,1) = alignment(2,2)+2*bc_size+1
                  alignment(2,1) = alignment(2,1)+relative_distance(i)
                  alignment(2,2) = alignment(2,1)+relative_sizes(1,i)

                  added_sublayer => interface_used%add_sublayer(
     $                 mainlayer_id, nodes, alignment)

                  if(relative_sizes(2,i).gt.0) then
                     call over_allocate_west(
     $                    added_sublayer,
     $                    relative_sizes(2,i), nodes)
                  end if

               end do

            end select

        end subroutine ini_interface_testcase3


        subroutine ini_interface_testcase4(
     $     interface_used,
     $     nodes,
     $     relative_distance,
     $     relative_size,
     $     over_allocated)

          implicit none

          class(bf_interface)             , intent(inout) :: interface_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes
          integer                         , intent(in)    :: relative_distance
          integer                         , intent(in)    :: relative_size
          integer                         , intent(in)    :: over_allocated

          type(bf_sublayer), pointer     :: added_sublayer
          integer(ikind), dimension(2,2) :: alignment
          logical       , dimension(4)   :: neighbors

          print '(''************************************'')'
          print '(''4: interface with one sublayers     '')'
          print '('' : for each mainlayer with an       '')'
          print '('' : over allcoation in each direction'')'
          print '(''************************************'')'

          !initialization
          call interface_used%ini()

          neighbors     = [.false.,.false.,.false.,.false.]


          !northern main layer
          !--------------------------------------------------
          !initialize the alignment
          alignment(1,1) = 1 + bc_size + relative_distance
          alignment(1,2) = alignment(1,1) + relative_size
          alignment(2,1) = ny+1
          alignment(2,2) = ny+1

          !add the sublayer
          added_sublayer => interface_used%add_sublayer(
     $         N, nodes, alignment)

          !over_allocated in all directions
          alignment(1,1) = alignment(1,1) - over_allocated
          alignment(1,2) = alignment(1,2) + over_allocated
          alignment(2,2) = alignment(2,2) + over_allocated
          call added_sublayer%reallocate_bf_layer(nodes, alignment)

          
          !southern main layer
          !--------------------------------------------------
          !initialize the alignment
          alignment(1,1) = 1 + bc_size + relative_distance
          alignment(1,2) = alignment(1,1) + relative_size
          alignment(2,1) = 0
          alignment(2,2) = 0

          !add the sublayer
          added_sublayer => interface_used%add_sublayer(
     $         S, nodes, alignment)

          !over_allocated in all directions
          alignment(1,1) = alignment(1,1) - over_allocated
          alignment(1,2) = alignment(1,2) + over_allocated
          alignment(2,1) = alignment(2,1) - over_allocated
          call added_sublayer%reallocate_bf_layer(nodes, alignment)


          !easthern main layer
          !--------------------------------------------------
          !initialize the alignment
          alignment(1,1) = nx+1
          alignment(1,2) = nx+1
          alignment(2,1) = 1 + bc_size + relative_distance
          alignment(2,2) = alignment(2,1) + relative_size

          !add the sublayer
          added_sublayer => interface_used%add_sublayer(
     $         E, nodes, alignment)

          !over_allocated in all directions
          alignment(1,2) = alignment(1,2) + over_allocated
          alignment(2,1) = alignment(2,1) - over_allocated
          alignment(2,2) = alignment(2,2) + over_allocated
          call added_sublayer%reallocate_bf_layer(nodes, alignment)
      

          !westhern main layer
          !--------------------------------------------------
          !initialize the alignment
          alignment(1,1) = 0
          alignment(1,2) = 0
          alignment(2,1) = 1 + bc_size + relative_distance
          alignment(2,2) = alignment(2,1) + relative_size

          !add the sublayer
          added_sublayer => interface_used%add_sublayer(
     $         W, nodes, alignment)

          !over_allocated in all directions
          alignment(1,1) = alignment(1,1) - over_allocated
          alignment(2,1) = alignment(2,1) - over_allocated
          alignment(2,2) = alignment(2,2) + over_allocated
          call added_sublayer%reallocate_bf_layer(nodes, alignment)

        end subroutine ini_interface_testcase4


        subroutine over_allocate_north(
     $     added_sublayer,
     $     over_allocated,
     $     nodes)

          implicit none

          type(bf_sublayer), pointer      , intent(inout) :: added_sublayer
          integer                         , intent(in)    :: over_allocated
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes

          integer(ikind), dimension(2)                :: new_sizes
          integer       , dimension(:,:), allocatable :: grdpts_id
          integer       , dimension(2,2)              :: alignment
          integer       , dimension(2,2)              :: new_alignment

          if(over_allocated.gt.0) then

             alignment          = added_sublayer%get_alignment_tab()
             new_alignment      = alignment
             new_alignment(2,2) = new_alignment(2,2) + over_allocated
             
             call added_sublayer%reallocate_bf_layer(
     $            nodes, new_alignment)

             new_sizes = added_sublayer%get_sizes()
             allocate(grdpts_id(new_sizes(1), new_sizes(2)))
             
             call ini_all_interior(grdpts_id)
             call make_east_layer(grdpts_id)
             call make_west_layer(grdpts_id)
             call make_north_layer(grdpts_id)
             call add_south_layer_exchange(grdpts_id)

             call added_sublayer%set_grdpts_id(grdpts_id)
             
          end if

        end subroutine over_allocate_north


        subroutine over_allocate_south(
     $     added_sublayer,
     $     over_allocated,
     $     nodes)

          implicit none

          type(bf_sublayer), pointer      , intent(in) :: added_sublayer
          integer                         , intent(in) :: over_allocated
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes

          integer(ikind), dimension(2)                :: new_sizes
          integer       , dimension(:,:), allocatable :: grdpts_id
          integer       , dimension(2,2)              :: alignment
          integer       , dimension(2,2)              :: new_alignment

          if(over_allocated.gt.0) then

             alignment          = added_sublayer%get_alignment_tab()
             new_alignment      = alignment
             new_alignment(2,1) = new_alignment(2,1) - over_allocated
             
             call added_sublayer%reallocate_bf_layer(
     $            nodes, new_alignment)

             new_sizes = added_sublayer%get_sizes()
             allocate(grdpts_id(new_sizes(1), new_sizes(2)))

             call ini_all_interior(grdpts_id)
             call make_east_layer(grdpts_id)
             call make_west_layer(grdpts_id)
             call make_south_layer(grdpts_id)
             call add_north_layer_exchange(grdpts_id)

             call added_sublayer%set_grdpts_id(grdpts_id)
             
          end if

        end subroutine over_allocate_south


        subroutine over_allocate_east(
     $     added_sublayer,
     $     over_allocated,
     $     nodes)

          implicit none

          type(bf_sublayer), pointer      , intent(in) :: added_sublayer
          integer                         , intent(in) :: over_allocated
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes

          integer(ikind), dimension(2)                :: new_sizes
          integer       , dimension(:,:), allocatable :: grdpts_id
          integer       , dimension(2,2)              :: alignment
          integer       , dimension(2,2)              :: new_alignment

          if(over_allocated.gt.0) then

             alignment          = added_sublayer%get_alignment_tab()
             new_alignment      = alignment
             new_alignment(1,2) = new_alignment(1,2) + over_allocated
             
             call added_sublayer%reallocate_bf_layer(
     $            nodes, new_alignment)
             
             new_sizes = added_sublayer%get_sizes()
             allocate(grdpts_id(new_sizes(1), new_sizes(2)))

             call ini_all_interior(grdpts_id)
             call make_east_layer(grdpts_id)
             call make_north_layer(grdpts_id)
             call make_south_layer(grdpts_id)
             call add_west_layer_exchange(grdpts_id)

             call added_sublayer%set_grdpts_id(grdpts_id)
             
          end if

        end subroutine over_allocate_east


        subroutine over_allocate_west(
     $     added_sublayer,
     $     over_allocated,
     $     nodes)

          implicit none

          type(bf_sublayer), pointer      , intent(in) :: added_sublayer
          integer                         , intent(in) :: over_allocated
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes

          integer(ikind), dimension(2)                :: new_sizes
          integer       , dimension(:,:), allocatable :: grdpts_id
          integer       , dimension(2,2)              :: alignment
          integer       , dimension(2,2)              :: new_alignment

          if(over_allocated.gt.0) then

             alignment          = added_sublayer%get_alignment_tab()
             new_alignment      = alignment
             new_alignment(1,1) = new_alignment(1,1) - over_allocated
             
             call added_sublayer%reallocate_bf_layer(
     $            nodes, new_alignment)

             new_sizes = added_sublayer%get_sizes()
             allocate(grdpts_id(new_sizes(1), new_sizes(2)))
             
             call ini_all_interior(grdpts_id)
             call make_west_layer(grdpts_id)
             call make_north_layer(grdpts_id)
             call make_south_layer(grdpts_id)
             call add_east_layer_exchange(grdpts_id)

             call added_sublayer%set_grdpts_id(grdpts_id)
             
          end if

        end subroutine over_allocate_west


        subroutine ini_all_interior(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i,j

          do j=1,size(grdpt_id,2)
             do i=1, size(grdpt_id,1)
                grdpt_id(i,j)=interior_pt
             end do
          end do

        end subroutine ini_all_interior


        subroutine make_west_layer(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: j

          do j=1,size(grdpt_id,2)
             grdpt_id(1,j)=bc_pt
             grdpt_id(2,j)=bc_interior_pt
          end do

        end subroutine make_west_layer


        subroutine make_east_layer(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: j

          do j=1,size(grdpt_id,2)
             grdpt_id(size(grdpt_id,1)-1,j)=bc_interior_pt
             grdpt_id(size(grdpt_id,1)  ,j)=bc_pt
          end do

        end subroutine make_east_layer


        subroutine make_north_layer(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i

          do i=bc_size,size(grdpt_id,1)-1
             grdpt_id(i,size(grdpt_id,2)-1) = bc_interior_pt
             grdpt_id(i,size(grdpt_id,2))   = bc_pt
          end do

        end subroutine make_north_layer


        subroutine make_south_layer(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i

          do i=bc_size,size(grdpt_id,1)-1
             grdpt_id(i,1)       = bc_pt
             grdpt_id(i,bc_size) = bc_interior_pt
          end do
        
        end subroutine make_south_layer

        
        subroutine add_NE_corner(grdpt_id)
        
          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          grdpt_id(size(grdpt_id,1),size(grdpt_id,2)-1) = bc_interior_pt
          grdpt_id(size(grdpt_id,1),size(grdpt_id,2))   = bc_pt

        end subroutine add_NE_corner


        subroutine add_NW_corner(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          grdpt_id(1,size(grdpt_id,2)-1) = bc_interior_pt
          grdpt_id(1,size(grdpt_id,2))   = bc_pt

        end subroutine add_NW_corner


        subroutine add_SW_corner(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          grdpt_id(1,1)   = bc_pt
          grdpt_id(1,bc_size) = bc_interior_pt

        end subroutine add_SW_corner

      
        subroutine add_SE_corner(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          grdpt_id(size(grdpt_id,1),1)       = bc_pt
          grdpt_id(size(grdpt_id,1),bc_size) = bc_interior_pt

        end subroutine add_SE_corner


        subroutine add_north_layer_exchange(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i,j
          integer :: exchange_pt
          
          exchange_pt = 4

          do j=size(grdpt_id,2)-bc_size+1, size(grdpt_id,2)
             do i=1, size(grdpt_id,1)
                grdpt_id(i,j)= exchange_pt
             end do
          end do

        end subroutine add_north_layer_exchange


        subroutine add_south_layer_exchange(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i,j
          integer :: exchange_pt
          
          exchange_pt = 4

          do j=1,bc_size
             do i=1, size(grdpt_id,1)
                grdpt_id(i,j)= exchange_pt
             end do
          end do

        end subroutine add_south_layer_exchange


        subroutine add_east_layer_exchange(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i,j
          integer :: exchange_pt
          
          exchange_pt = 4

          do j=1,size(grdpt_id,2)
             do i=size(grdpt_id,1)-bc_size+1, size(grdpt_id,1)
                grdpt_id(i,j)= exchange_pt
             end do
          end do

        end subroutine add_east_layer_exchange


        subroutine add_west_layer_exchange(grdpt_id)

          implicit none

          integer, dimension(:,:), intent(inout) :: grdpt_id

          integer :: i,j
          integer :: exchange_pt
          
          exchange_pt = 4

          do j=1,size(grdpt_id,2)
             do i=1,bc_size
                grdpt_id(i,j)=exchange_pt
             end do
          end do

        end subroutine add_west_layer_exchange

      end module test_cases_interface_module
