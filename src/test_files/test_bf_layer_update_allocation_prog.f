      program test_bf_layer_update_allocation_prog

        use bf_mainlayer_class               , only : bf_mainlayer
        use bf_layer_update_allocation_module, only : update_allocation_bf_layers
        use bf_layer_path_class              , only : bf_layer_path
        use bf_sublayer_class                , only : bf_sublayer
                                             
        use interface_abstract_class         , only : interface_abstract

        use parameters_constant              , only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input                 , only : nx,ny,ne,bc_size
        use parameters_kind                  , only : ikind, rkind
                                            
        use test_bf_layer_module             , only : print_nodes,
     $                                                print_grdpts_id,
     $                                                print_sizes,
     $                                                ini_grdpts_id,
     $                                                ini_nodes
        use test_cases_interface_module      , only : ini_interface
        use test_cases_path_module           , only : ini_path
        use interface_print_module           , only : print_interface
        use sublayers_ini_module             , only : ini_relative_sizes,
     $                                                ini_relative_distance,
     $                                                ini_alignment,
     $                                                ini_neighbors

        implicit none

        !parameters for the test
        integer, parameter :: interface_before = 0
        integer, parameter :: interface_after = 1
        integer, parameter :: test_case_id = 41
        integer, parameter :: mainlayer_id = 4
        integer, parameter :: random_seed = 86456


        !local variables
        real(rkind), dimension(nx,ny,ne) :: nodes
        integer    , dimension(nx,ny)    :: grdpts_id
        type(interface_abstract)         :: interface_used
        type(bf_layer_path)              :: current_path
        type(bf_sublayer), pointer       :: modified_sublayer


        !initialize the nodes for the test
        call ini_nodes(nodes)
        call ini_grdpts_id(grdpts_id)

        !initialize the test case: interface + path
        call ini_test_case(
     $       test_case_id,
     $       mainlayer_id,
     $       random_seed,
     $       interface_used,
     $       current_path)

        !modify the grdpts_id for the interior nodes to show the
        !alignment of the path
        call print_alignment_on_grdpts_id(
     $       mainlayer_id,
     $       current_path%alignment,
     $       grdpts_id)


        !print nodes before
        call print_sizes(nodes,'interior_sizes.dat')
        call print_grdpts_id(grdpts_id,'interior_grdpts_id.dat')
        call print_nodes(nodes,'interior_nodes.dat')

        !print interface before
        call print_interface(interface_used, interface_before)

        !update the buffer layer allocation
        modified_sublayer => update_allocation_bf_layers(
     $       interface_used,
     $       nodes,
     $       current_path)

        !print interface after
        call print_interface(interface_used, interface_after)


        contains

        subroutine ini_test_case(
     $       test_case_id,
     $       test_mainlayer_id,
     $       random_seed,
     $       interface_used,
     $       current_path)

          implicit none

          integer                  , intent(in)  :: test_case_id
          integer                  , intent(in)  :: test_mainlayer_id
          integer                  , intent(in)  :: random_seed
          class(interface_abstract), intent(out) :: interface_used
          type(bf_layer_path)      , intent(out) :: current_path

          integer                        :: corner_id, corner_order
          type(bf_mainlayer), pointer    :: mainlayer
          integer, dimension(2,2)        :: relative_sizes
          integer, dimension(3)          :: relative_distance
          integer(ikind), dimension(2,2) :: alignment

          select case(test_case_id)
            case(11,12,13,14)
               
               if(test_case_id.eq.11) then
                  print '(''************************************'')'
                  print '(''* test case 11                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with no sublayers      *'')'
                  print '(''* path with no neighbors           *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if
               if(test_case_id.eq.12) then
                  print '(''************************************'')'
                  print '(''* test case 12                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with no sublayers      *'')'
                  print '(''* path with one neighbor           *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if
               if(test_case_id.eq.13) then
                  print '(''************************************'')'
                  print '(''* test case 13                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with no sublayers      *'')'
                  print '(''* path with one neighbor           *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if
               if(test_case_id.eq.14) then
                  print '(''************************************'')'
                  print '(''* test case 14                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with no sublayers      *'')'
                  print '(''* path with two neighbors          *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if

               !no sublayers
               call ini_interface(
     $              interface_used,
     $              0,
     $              nodes)

               !path
               current_path%ends=.true.
               current_path%ends_with_corner=.false.
               current_path%corner_id=N_E
               current_path%mainlayer=test_mainlayer_id

               if(test_case_id.eq.11) then
                  current_path%neighbors=[.false.,.false.,.false.,.false.]
               end if
               if(test_case_id.eq.12) then
                  current_path%neighbors=[.false.,.true.,.false.,.true.]
               end if
               if(test_case_id.eq.13) then
                  current_path%neighbors=[.true.,.false.,.true.,.false.]
               end if
               if(test_case_id.eq.14) then
                  current_path%neighbors=[.true.,.true.,.true.,.true.]
               end if
               
               call ini_path_alignment(
     $              2, 2, test_mainlayer_id, 1, 2, random_seed,
     $              current_path%alignment)           

            case(21,22,23,24)
               if(test_case_id.eq.21) then
                  print '(''************************************'')'
                  print '(''* test case 21                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with one sublayer      *'')'
                  print '(''* that does not match the path     *'')'
                  print '(''* path with no neighbors           *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if
               if(test_case_id.eq.22) then
                  print '(''************************************'')'
                  print '(''* test case 22                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with one sublayer      *'')'
                  print '(''* that does not match the path     *'')'
                  print '(''* path with one neighbor           *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if
               if(test_case_id.eq.23) then
                  print '(''************************************'')'
                  print '(''* test case 23                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with one sublayer      *'')'
                  print '(''* that does not match the path     *'')'
                  print '(''* path with one neighbor           *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if
               if(test_case_id.eq.24) then
                  print '(''************************************'')'
                  print '(''* test case 24                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with one sublayer      *'')'
                  print '(''* that does not match the path     *'')'
                  print '(''* path with two neighbors          *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if

               !one sublayer that does not match the path
               select case(test_mainlayer_id)
                 case(N)
                    corner_id    = N_E
                    corner_order = 1
                 case(S)
                    corner_id    = S_E
                    corner_order = 2
                 case(E)
                    corner_id    = N_E
                    corner_order = 2
                 case(W)
                    corner_id    = N_W
                    corner_order = 1
                 case default
                    print '(''test_bf_layer_update_allocation_prog'')'
                    print '(''ini_test_case'')'
                    print '(''test mainlayer ID not recognized'')'
                    print '(''test_mainlayer_ID:'', I2)', test_mainlayer_id
                    stop 'choose another mainlayer'
               end select

               call ini_interface(
     $              interface_used,
     $              1,
     $              nodes,
     $              corner_id=corner_id,
     $              corner_order=corner_order,
     $              bf_corner_distance=2,
     $              over_allocated=1)

               !path
               current_path%ends=.true.
               current_path%ends_with_corner=.false.
               current_path%corner_id=N_E
               current_path%mainlayer=test_mainlayer_id

               if(test_case_id.eq.21) then
                  current_path%neighbors=[.false.,.false.,.false.,.false.]
               end if
               if(test_case_id.eq.22) then
                  current_path%neighbors=[.false.,.true.,.false.,.true.]
               end if
               if(test_case_id.eq.23) then
                  current_path%neighbors=[.true.,.false.,.true.,.false.]
               end if
               if(test_case_id.eq.24) then
                  current_path%neighbors=[.true.,.true.,.true.,.true.]
               end if
               
               call ini_path_alignment(
     $              2, 2, test_mainlayer_id, 1, 2, random_seed,
     $              current_path%alignment)

            case(31,32,33,34)
               if(test_case_id.eq.31) then
                  print '(''************************************'')'
                  print '(''* test case 31                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with one sublayer      *'')'
                  print '(''* that matches the path            *'')'
                  print '(''* path with no neighbors           *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if
               if(test_case_id.eq.32) then
                  print '(''************************************'')'
                  print '(''* test case 32                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with one sublayer      *'')'
                  print '(''* that matches the path            *'')'
                  print '(''* path with one neighbor           *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if
               if(test_case_id.eq.33) then
                  print '(''************************************'')'
                  print '(''* test case 33                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with one sublayer      *'')'
                  print '(''* that matches the path            *'')'
                  print '(''* path with one neighbor           *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if
               if(test_case_id.eq.34) then
                  print '(''************************************'')'
                  print '(''* test case 34                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with one sublayer      *'')'
                  print '(''* that matches the path            *'')'
                  print '(''* path with two neighbors          *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if

               !one sublayer that does not match the path
               select case(test_mainlayer_id)
                 case(N)
                    corner_id    = N_E
                    corner_order = 1
                 case(S)
                    corner_id    = S_W
                    corner_order = 1
                 case(E)
                    corner_id    = N_E
                    corner_order = 2
                 case(W)
                    corner_id    = N_W
                    corner_order = 1
                 case default
                    print '(''test_bf_layer_update_allocation_prog'')'
                    print '(''ini_test_case'')'
                    print '(''test mainlayer ID not recognized'')'
                    print '(''test_mainlayer_ID:'', I2)', test_mainlayer_id
                    stop 'choose another mainlayer'
               end select

               call ini_interface(
     $              interface_used,
     $              1,
     $              nodes,
     $              corner_id=corner_id,
     $              corner_order=corner_order,
     $              bf_corner_distance=6,
     $              over_allocated=1)

               !path
               current_path%ends=.true.
               current_path%ends_with_corner=.false.
               current_path%corner_id=N_E
               current_path%mainlayer=test_mainlayer_id

               mainlayer => interface_used%get_mainlayer(test_mainlayer_id)
               if(associated(mainlayer)) then
                  current_path%matching_sublayer => mainlayer%head_sublayer
               else
                  print '(''test_bf_layer_update_allocation_prog'')'
                  print '(''ini_test_case'')'
                  print '(''interface not correctly initialized'')'
                  print '(''no sublayer allocated for mainlayer_id'')'
                  print '(''mainlayer_id:'', I2)', test_mainlayer_id
                  stop 'verify the interface initialization'
               end if

               if(test_case_id.eq.31) then
                  current_path%neighbors=[.false.,.false.,.false.,.false.]
               end if
               if(test_case_id.eq.32) then
                  current_path%neighbors=[.false.,.true.,.false.,.true.]
               end if
               if(test_case_id.eq.33) then
                  current_path%neighbors=[.true.,.false.,.true.,.false.]
               end if
               if(test_case_id.eq.34) then
                  current_path%neighbors=[.true.,.true.,.true.,.true.]
               end if
               
               call ini_path_alignment(
     $              2, 2, test_mainlayer_id, 1, 2, random_seed,
     $              current_path%alignment)


            case(41,42,43,44)
               if(test_case_id.eq.41) then
                  print '(''************************************'')'
                  print '(''* test case 41                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with two sublayers     *'')'
                  print '(''* that matches the path            *'')'
                  print '(''* path with no neighbors           *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if
               if(test_case_id.eq.42) then
                  print '(''************************************'')'
                  print '(''* test case 42                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with two sublayers     *'')'
                  print '(''* that matches the path            *'')'
                  print '(''* path with one neighbor           *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if
               if(test_case_id.eq.43) then
                  print '(''************************************'')'
                  print '(''* test case 43                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with two sublayers     *'')'
                  print '(''* that matches the path            *'')'
                  print '(''* path with one neighbor           *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if
               if(test_case_id.eq.44) then
                  print '(''************************************'')'
                  print '(''* test case 44                     *'')'
                  print '(''************************************'')'
                  print '(''* interface with two sublayers     *'')'
                  print '(''* that matches the path            *'')'
                  print '(''* path with two neighbors          *'')'
                  print '(''************************************'')'
                  print '('''')'
               end if

               call ini_relative_sizes(
     $              2,
     $              random_seed,
     $              relative_sizes)
               
               call ini_relative_distance(
     $              2,
     $              random_seed,
     $              relative_distance)
               
               call ini_alignment(
     $              alignment,
     $              mainlayer_id,
     $              1,
     $              relative_sizes,
     $              relative_distance)

               !one sublayer that does not match the path
               call ini_interface(
     $              interface_used,
     $              3,
     $              nodes,
     $              mainlayer_id=test_mainlayer_id,
     $              nb_sublayers=2,
     $              relative_distance=relative_distance,
     $              relative_sizes=relative_sizes)

               !path
               current_path%ends=.true.
               current_path%ends_with_corner=.false.
               current_path%corner_id=N_E
               current_path%mainlayer=test_mainlayer_id

               mainlayer => interface_used%get_mainlayer(test_mainlayer_id)
               if(associated(mainlayer)) then
                  current_path%matching_sublayer => mainlayer%head_sublayer
               else
                  print '(''test_bf_layer_update_allocation_prog'')'
                  print '(''ini_test_case'')'
                  print '(''interface not correctly initialized'')'
                  print '(''no sublayer allocated for mainlayer_id'')'
                  print '(''mainlayer_id:'', I2)', test_mainlayer_id
                  stop 'verify the interface initialization'
               end if

               if(test_case_id.eq.41) then
                  current_path%neighbors=[.false.,.false.,.false.,.false.]
               end if
               if(test_case_id.eq.42) then
                  current_path%neighbors=[.false.,.true.,.false.,.true.]
               end if
               if(test_case_id.eq.43) then
                  current_path%neighbors=[.true.,.false.,.true.,.false.]
               end if
               if(test_case_id.eq.44) then
                  current_path%neighbors=[.true.,.true.,.true.,.true.]
               end if
               
               current_path%alignment=alignment

            case default
               print '(''test_bf_layer_update_allocation_prog'')'
               print '(''ini_test_case'')'
               print '(''test case ID not recognized'')'
               print '(''test_case_ID:'', I2)', test_case_ID
               stop 'choose another test case'
          end select
        end subroutine ini_test_case


        subroutine ini_path_alignment(
     $     size_case,
     $     distance_case,
     $     mainlayer_id,
     $     merge_id,
     $     nb_sublayers,
     $     seed,
     $     alignment)

          implicit none

          integer                , intent(in)  :: size_case
          integer                , intent(in)  :: distance_case
          integer                , intent(in)  :: mainlayer_id
          integer                , intent(in)  :: merge_id
          integer                , intent(in)  :: nb_sublayers
          integer                , intent(in)  :: seed
          integer, dimension(2,2), intent(out) :: alignment

          integer, dimension(:,:), allocatable :: relative_sizes
          integer, dimension(:)  , allocatable :: relative_distance


          allocate(relative_sizes(2,nb_sublayers))
          allocate(relative_distance(nb_sublayers+1))          

          call ini_relative_sizes(
     $         size_case,
     $         seed,
     $         relative_sizes)

          call ini_relative_distance(
     $         distance_case,
     $         seed,
     $         relative_distance)

          call ini_alignment(
     $         alignment,
     $         mainlayer_id,
     $         merge_id,
     $         relative_sizes,
     $         relative_distance)

          deallocate(relative_sizes)
          deallocate(relative_distance)

        end subroutine ini_path_alignment


        subroutine print_alignment_on_grdpts_id(
     $     mainlayer_id,
     $     alignment,
     $     grdpts_id)


          implicit none
          
          integer                         , intent(in)    :: mainlayer_id
          integer(ikind), dimension(2,2)  , intent(in)    :: alignment
          integer       , dimension(nx,ny), intent(inout) :: grdpts_id


          integer(ikind) :: i,j
          integer, parameter :: path_pt=5
          

          select case(mainlayer_id)
            case(N)
               j = ny-bc_size
               do i=alignment(1,1), alignment(1,2)
                  grdpts_id(i,j) = path_pt
               end do

            case(S)
               j = bc_size+1
               do i=alignment(1,1), alignment(1,2)
                  grdpts_id(i,j) = path_pt
               end do

            case(E)
               i = nx-bc_size
               do j=alignment(2,1), alignment(2,2)
                  grdpts_id(i,j) = path_pt
               end do

            case(W)
               i = bc_size+1
               do j=alignment(2,1), alignment(2,2)
                  grdpts_id(i,j) = path_pt
               end do

            case default
               print '(''test_bf_layer_update_allocation_prog'')'
               print '(''print_alignemnt_on_grdpts_id'')'
               print '(''mainlayer_id: '', I2)', mainlayer_id
               stop 'mainlayer_id not recognized'
          end select

        end subroutine print_alignment_on_grdpts_id

      end program test_bf_layer_update_allocation_prog
