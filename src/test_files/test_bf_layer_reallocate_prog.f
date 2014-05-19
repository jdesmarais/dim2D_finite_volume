      program test_bf_layer_reallocate_prog

        use bf_interface_class         , only : bf_interface
        use bf_mainlayer_class         , only : bf_mainlayer
        use bf_sublayer_class          , only : bf_sublayer

        use parameters_constant        , only : N,S,E,W
        use parameters_input           , only : nx,ny,ne,bc_size
        use parameters_kind            , only : ikind, rkind

        use test_bf_layer_module       , only : print_interior_data,
     $                                          ini_grdpts_id,
     $                                          ini_nodes
        use test_cases_interface_module, only : ini_interface


        implicit none


        integer, parameter :: test_case_id = 4
        integer, parameter :: test_relative_size = 1
        integer, parameter :: test_relative_distance = 1
        integer, parameter :: test_final_alignment = 1
        integer, parameter :: over_allocated = 0

        real(rkind), dimension(nx,ny,ne) :: nodes
        integer    , dimension(nx,ny)    :: grdpts_id
        type(bf_interface)               :: interface_used
        type(bf_mainlayer), pointer      :: mainlayer
        type(bf_sublayer) , pointer      :: sublayer

        integer :: i
        integer :: relative_size
        integer :: relative_distance
        integer(ikind), dimension(2,2) :: final_alignment


        !initialize the nodes for the test
        call ini_nodes(nodes)
        call ini_grdpts_id(grdpts_id)
        call print_interior_data(
     $       nodes,
     $       grdpts_id, 
     $       'interior_nodes.dat',
     $       'interior_grdpts_id.dat',
     $       'interior_sizes.dat')

        !initialize the relative_distance, relative_size
        !and over_allocated for the interface
        call ini_for_test(
     $       N,
     $       over_allocated,
     $       test_relative_size,
     $       relative_size,
     $       test_relative_distance,
     $       relative_distance,
     $       test_final_alignment,
     $       final_alignment)

        !initialize the interface
        call ini_interface(
     $       interface_used,
     $       test_case_id,
     $       nodes,
     $       relative_distance=relative_distance,
     $       relative_size=relative_size)

        !print interface before the test
        call interface_used%print_binary(
     $       'nodes0.dat',
     $       'grdpt_id0.dat',
     $       'sizes0.dat',
     $       '0.dat')


        !test the function
        do i=1,4

           !ini the data for the test
            call ini_for_test(
     $          i,
     $          over_allocated,
     $          test_relative_size,
     $          relative_size,
     $          test_relative_distance,
     $          relative_distance,
     $          test_final_alignment,
     $          final_alignment)

           !test
           mainlayer => interface_used%get_mainlayer(i)
           sublayer  => mainlayer%get_head_sublayer()

           call sublayer%reallocate_bf_layer(nodes, final_alignment)
           print *, 'final_alignment', final_alignment

        end do


        !print the interface after the test
        call interface_used%print_binary(
     $       'nodes1.dat',
     $       'grdpt_id1.dat',
     $       'sizes1.dat',
     $       '1.dat')


        contains


        subroutine ini_for_test(
     $       mainlayer_id,
     $       over_allocated,
     $       test_relative_size,
     $       relative_size,
     $       test_relative_distance,
     $       relative_distance,
     $       test_final_alignment,
     $       final_alignment)


          implicit none

          integer                , intent(in)  :: mainlayer_id
          integer                , intent(in)  :: over_allocated
          integer                , intent(in)  :: test_relative_size
          integer                , intent(out) :: relative_size
          integer                , intent(in)  :: test_relative_distance
          integer                , intent(out) :: relative_distance
          integer                , intent(in)  :: test_final_alignment
          integer, dimension(2,2), intent(out) :: final_alignment

          !test relative size
          select case(test_relative_size)
            case(1)
               relative_size = 0
            case(2)
               relative_size = 1
            case(3)
               relative_size = 5
            case default
               print '(''test_bf_layer_reallocate_prog'')'
               print '(''ini_for_test'')'
               print '(''test not recognized for relative_size'')'
               print '(''test_relative_size: '',I2)',
     $              test_relative_size
               stop 'change test_relative_size [1,3]'
          end select          

          !test relative distance
          select case(test_relative_distance)
            case(1)
               relative_distance = 0
            case(2)
               relative_distance = 1
            case(3)
               relative_distance = 5
            case(4)
               select case(mainlayer_id)
                 case(N,S)
                    relative_distance = nx - 2*bc_size - relative_size - 2
                 case(E,W)
                    relative_distance = ny - 2*bc_size - relative_size - 2
               end select
            case(5)
               select case(mainlayer_id)
                 case(N,S)
                    relative_distance = nx - 2*bc_size - relative_size - 1
                 case(E,W)
                    relative_distance = ny - 2*bc_size - relative_size - 1
               end select
            case default
               print '(''test_bf_layer_reallocate_prog'')'
               print '(''ini_for_test'')'
               print '(''test not recognized for relative_distance'')'
               print '(''test_relative_distance: '',I2)',
     $              test_relative_distance
               stop 'change test_relative_distance [1,3]'
          end select 


          !test final_alignment
          select case(mainlayer_id)
            case(N)
               final_alignment(2,1) = ny+1
               final_alignment(2,2) = ny+1+over_allocated

               select case(test_final_alignment)
                 case(1)
                    final_alignment(1,1) = bc_size+1
                    final_alignment(1,2) = nx-bc_size
                 case(2)
                    final_alignment(1,1) = bc_size+2
                    final_alignment(1,2) = nx-bc_size-1
                 case(3)
                    final_alignment(1,1) = bc_size+1+
     $                                     relative_distance-
     $                                     over_allocated
                    final_alignment(1,2) = final_alignment(1,1)+
     $                                     relative_size+
     $                                     2*over_allocated
                 case default
                    print '(''test_bf_layer_reallocate_prog'')'
                    print '(''ini_for_test'')'
                    print '(''test not recognized final_alignment'')'
                    print '(''test_final_alignment: '',I2)',
     $                   test_final_alignment
                    stop 'change test_final_alignment [1,3]'
               end select
 
            case(S)
               final_alignment(2,1) = 0-over_allocated
               final_alignment(2,2) = 0

               select case(test_final_alignment)
                 case(1)
                    final_alignment(1,1) = bc_size+1
                    final_alignment(1,2) = nx-bc_size
                 case(2)
                    final_alignment(1,1) = bc_size+2
                    final_alignment(1,2) = nx-bc_size-1
                 case(3)
                    final_alignment(1,1) = bc_size+1+
     $                                     relative_distance-
     $                                     over_allocated
                    final_alignment(1,2) = final_alignment(1,1)+
     $                                     relative_size+2*over_allocated
                 case default
                    print '(''test_bf_layer_reallocate_prog'')'
                    print '(''ini_for_test'')'
                    print '(''test not recognized final_alignment'')'
                    print '(''test_final_alignment: '',I2)',
     $                   test_final_alignment
                    stop 'change test_final_alignment [1,3]'
               end select
                 
            case(E)
               final_alignment(1,1) = nx+1
               final_alignment(1,2) = nx+1+over_allocated

               select case(test_final_alignment)
                 case(1)
                    final_alignment(2,1) = bc_size+1
                    final_alignment(2,2) = nx-bc_size
                 case(2)
                    final_alignment(2,1) = bc_size+2
                    final_alignment(2,2) = nx-bc_size-1
                 case(3)
                    final_alignment(2,1) = bc_size+1+
     $                                     relative_distance-
     $                                     over_allocated
                    final_alignment(2,2) = final_alignment(2,1)+
     $                                     relative_size+2*over_allocated
                 case default
                    print '(''test_bf_layer_reallocate_prog'')'
                    print '(''ini_for_test'')'
                    print '(''test not recognized final_alignment'')'
                    print '(''test_final_alignment: '',I2)',
     $                   test_final_alignment
                    stop 'change test_final_alignment [1,3]'
                 end select

              case(W)
                 final_alignment(1,1) = 0-over_allocated
                 final_alignment(1,2) = 0

                 select case(test_final_alignment)
                   case(1)
                      final_alignment(2,1) = bc_size+1
                      final_alignment(2,2) = nx-bc_size
                   case(2)
                      final_alignment(2,1) = bc_size+2
                      final_alignment(2,2) = nx-bc_size-1
                   case(3)
                      final_alignment(2,1) = bc_size+1+
     $                                       relative_distance-
     $                                       over_allocated
                      final_alignment(2,2) = final_alignment(2,1)+
     $                                       relative_size+2*over_allocated
                   case default
                      print '(''test_bf_layer_reallocate_prog'')'
                      print '(''ini_for_test'')'
                      print '(''test not recognized final_alignment'')'
                      print '(''test_final_alignment: '',I2)',
     $                     test_final_alignment
                      stop 'change test_final_alignment [1,3]'
                 end select

              case default
                 print '(''test_bf_layer_reallocate_prog'')'
                 print '(''ini_for_test'')'
                 print '(''mainlayer_id not recognized'')'
                 print '(''mainlayer_id: '',I2)',
     $                mainlayer_id
                 stop 'change mainlayer_id [1,4]'
            end select

        end subroutine ini_for_test

      end program test_bf_layer_reallocate_prog
