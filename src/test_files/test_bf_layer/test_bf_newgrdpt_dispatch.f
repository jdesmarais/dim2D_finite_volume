      program test_bf_newgrdpt_dispatch

        use bf_newgrdpt_dispatch_module, only :
     $     are_grdpts_available_to_get_newgrdpt_proc

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_kind, only :
     $       ikind

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled      = .true.
        test_validated = .true.


        test_loc = test_are_grdpts_available_to_get_newgrdpt_proc(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_are_grdpts_available_to_get_newgrdpt_proc: '',L1)', test_loc
        print '()'


        contains


        function test_are_grdpts_available_to_get_newgrdpt_proc(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer               :: k
          integer               :: bf_localization
          integer               :: size_x
          integer               :: size_y
          logical               :: bf_can_exchange_with_neighbor1
          logical               :: bf_can_exchange_with_neighbor2
          integer, dimension(2) :: bf_newgrdpt_coords
          logical               :: test_grdpts_available
          logical               :: grdpts_available
          logical               :: test_loc

          
          test_validated = .true.


          do k=1,210

             !input
             call get_param_test_are_grdpts_proc(
     $            k,
     $            bf_localization,
     $            size_x,size_y,
     $            bf_can_exchange_with_neighbor1,
     $            bf_can_exchange_with_neighbor2,
     $            bf_newgrdpt_coords,
     $            test_grdpts_available)

             !output
             grdpts_available = are_grdpts_available_to_get_newgrdpt_proc(
     $            bf_localization,
     $            size_x, size_y,
     $            bf_can_exchange_with_neighbor1,
     $            bf_can_exchange_with_neighbor2,
     $            bf_newgrdpt_coords)

             !validation
             test_loc = grdpts_available.eqv.test_grdpts_available
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''bf_localization      : '',I1)' , bf_localization
                print '(''size_x               : '',I1)' , size_x
                print '(''size_y               : '',I1)' , size_y
                print '(''neighbor1            : '',L1)' , bf_can_exchange_with_neighbor1
                print '(''neighbor2            : '',L1)' , bf_can_exchange_with_neighbor2
                print '(''bf_newgrdpt_coords   : '',2I2)', bf_newgrdpt_coords
                print '(''test_grdpts_available: '',L1)' , test_grdpts_available
                print '()'
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_are_grdpts_available_to_get_newgrdpt_proc


        subroutine get_param_test_are_grdpts_proc(
     $     test_id,
     $     bf_localization,
     $     size_x,size_y,
     $     bf_can_exchange_with_neighbor1,
     $     bf_can_exchange_with_neighbor2,
     $     bf_newgrdpt_coords,
     $     test_grdpts_available)

          implicit none

          integer               , intent(in)  :: test_id
          integer               , intent(out) :: bf_localization
          integer               , intent(out) :: size_x
          integer               , intent(out) :: size_y
          logical               , intent(out) :: bf_can_exchange_with_neighbor1
          logical               , intent(out) :: bf_can_exchange_with_neighbor2
          integer, dimension(2) , intent(out) :: bf_newgrdpt_coords
          logical               , intent(out) :: test_grdpts_available

          integer :: test_id_loc
          integer :: test_id_card

          !    _       _       _
          ! 8 |_|_____|_|_____|_|
          ! 7 |_|     |_|     |_|
          ! 6 |_|     |_|     |_|
          !    _|      _      |_
          ! 4 |_|     |_|     |_|   21 pts tested = 21 bf_newgrdpt_coords tested
          !    _|      _      |_        |
          ! 2 |_|     |_|     |_|       |___\  test_id_loc [0,20]
          ! 1 |_|_____|_|_____|_|           /
          ! 0 |_|     |_|     |_|
          !    0       3       6
          !
          !            X
          !
          ! 10 localizations tested:  _____\  test_id_card [1,10]
          !   - N                          /
          !   - S
          !   - E ( neighbor1=T, neighbor2=T)
          !   - E ( neighbor1=T, neighbor2=F)
          !   - E ( neighbor1=F, neighbor2=T)
          !   - E ( neighbor1=F, neighbor2=F)
          !   - W ( neighbor1=T, neighbor2=T)
          !   - W ( neighbor1=T, neighbor2=F)
          !   - W ( neighbor1=F, neighbor2=T)
          !   - W ( neighbor1=F, neighbor2=F)
          !------------------------------------------------------------

          size_x = 5
          size_y = 7

          test_id_loc  = mod(test_id-1,21)
          test_id_card = ((test_id - test_id_loc)/21)+1


          !bf_newgrdpt_coords(1)
          select case(mod(test_id_loc,3))
            case(0)
               bf_newgrdpt_coords(1) = 0
            case(1)
               bf_newgrdpt_coords(1) = 3
            case(2)
               bf_newgrdpt_coords(1) = 6
            case default
               print '(''get_param_test_are_grdpts_proc'')'
               print '(''wrong bf_newgrdpt_coords(1)'')'
               stop ''
          end select

          !bf_newgrdpt_coords(2)
          select case((test_id_loc-mod(test_id_loc,3))/3)
            case(0)
               bf_newgrdpt_coords(2) = 0
            case(1)
               bf_newgrdpt_coords(2) = 1
            case(2)
               bf_newgrdpt_coords(2) = 2
            case(3)
               bf_newgrdpt_coords(2) = 4
            case(4)
               bf_newgrdpt_coords(2) = 6
            case(5)
               bf_newgrdpt_coords(2) = 7
            case(6)
               bf_newgrdpt_coords(2) = 8
            case default
               print '(''get_param_test_are_grdpts_proc'')'
               print '(''wrong bf_newgrdpt_coords(2)'')'
               stop ''
          end select
            
          !bf_neighbors
          bf_can_exchange_with_neighbor1 = .true.
          bf_can_exchange_with_neighbor2 = .true.

          select case(test_id_card)
            case(1,2,3,7)
               bf_can_exchange_with_neighbor1 = .true.
               bf_can_exchange_with_neighbor2 = .true.
            case(4)
               bf_can_exchange_with_neighbor2 = .false.
            case(5)
               bf_can_exchange_with_neighbor1 = .false.
            case(6)
               bf_can_exchange_with_neighbor1 = .false.
               bf_can_exchange_with_neighbor2 = .false.
            case(8)
               bf_can_exchange_with_neighbor2 = .false.
            case(9)
               bf_can_exchange_with_neighbor1 = .false.
            case(10)
               bf_can_exchange_with_neighbor1 = .false.
               bf_can_exchange_with_neighbor2 = .false.
            case default
               print '(''get_param_test_are_grdpts_proc'')'
               print '(''wrong can_exchange_with_neighbor1'')'
               print '(''wrong can_exchange_with_neighbor2'')'
               stop ''
          end select

          !bf_localization
          select case(test_id_card)
            case(1)
               bf_localization = N
            case(2)
               bf_localization = S
            case(3,4,5,6)
               bf_localization = E
            case(7,8,9,10)
               bf_localization = W
            case default
               print '(''get_param_test_are_grdpts_proc'')'
               print '(''wrong bf_localization'')'
               stop
          end select

          select case(test_id_card)
            !N
            case(1)
               select case(test_id_loc+1)
                 case(1,2,3,4,5,6)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !S
            case(2)
               select case(test_id_loc+1)
                 case(16,17,18,19,20,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !E(T,T)
            case(3)
               select case(test_id_loc+1)
                 case(1,2,3,4,5,6,7,10,13,16,17,18,19,20,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !E(T,F)
            case(4)
               select case(test_id_loc+1)
                 case(1,2,3,4,5,6,7,10,13,16,19)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !E(F,T)
            case(5)
               select case(test_id_loc+1)
                 case(1,4,7,10,13,16,17,18,19,20,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !E(F,F)
            case(6)
               select case(test_id_loc+1)
                 case(1,4,7,10,13,16,19)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !W(T,T)
            case(7)
               select case(test_id_loc+1)
                 case(1,2,3,4,5,6,9,12,15,16,17,18,19,20,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !W(T,F)
            case(8)
               select case(test_id_loc+1)
                 case(1,2,3,4,5,6,9,12,15,18,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !W(F,T)
            case(9)
               select case(test_id_loc+1)
                 case(3,6,9,12,15,16,17,18,19,20,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            !W(F,F)
            case(10)
               select case(test_id_loc+1)
                 case(3,6,9,12,15,18,21)
                    test_grdpts_available = .false.
                 case default
                    test_grdpts_available = .true.
               end select

            case default
               print '(''get_param_test_are_grdpts_proc'')'
               print '(''wrong test_id_card'')'

          end select

        end subroutine get_param_test_are_grdpts_proc
        
      end program test_bf_newgrdpt_dispatch
