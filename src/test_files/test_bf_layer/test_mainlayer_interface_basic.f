      program test_mainlayer_interface_basic

        use mainlayer_interface_basic_class, only :
     $       mainlayer_interface_basic

        use bf_layer_errors_module, only :
     $       error_mainlayer_interface_type,
     $       error_mainlayer_interface_incompatible

        use bf_sublayer_class, only :
     $       bf_sublayer

        use parameters_bf_layer, only :
     $       NE_interface_type,
     $       NW_interface_type,
     $       SE_interface_type,
     $       SW_interface_type

        use parameters_constant, only :
     $       N,S,E,W


        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated
        

        detailled = .true.
        test_validated = .true.


        test_loc = test_ini(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_ini: '',L1)', test_loc
        print '()'


        test_loc = test_set_mainlayer_interface_bf_layer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_set_mainlayer_interface_bf_layer: '',L1)', test_loc
        print '()'


        test_loc = test_remove_mainlayer_interface_bf_layer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_remove_mainlayer_interface_bf_layer: '',L1)', test_loc
        print '()'


        contains


        function test_ini(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mainlayer_interface_basic) :: mainlayer_interface_used
          logical                         :: test_loc


          test_validated = .true.


          !output
          call mainlayer_interface_used%ini()


          !validation
          test_loc = .not.associated(mainlayer_interface_used%NE_interface_N_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NE_interface_N_ptr failed'')'
          end if

          test_loc = .not.associated(mainlayer_interface_used%NE_interface_E_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NE_interface_E_ptr failed'')'
          end if

          test_loc = .not.associated(mainlayer_interface_used%NW_interface_N_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NW_interface_N_ptr failed'')'
          end if

          test_loc = .not.associated(mainlayer_interface_used%NW_interface_W_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NW_interface_W_ptr failed'')'
          end if

          test_loc = .not.associated(mainlayer_interface_used%SE_interface_S_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SE_interface_S_ptr failed'')'
          end if

          test_loc = .not.associated(mainlayer_interface_used%SE_interface_E_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SE_interface_E_ptr failed'')'
          end if

          test_loc = .not.associated(mainlayer_interface_used%SW_interface_S_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SW_interface_S_ptr failed'')'
          end if

          test_loc = .not.associated(mainlayer_interface_used%SW_interface_W_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SW_interface_W_ptr failed'')'
          end if

        end function test_ini


        function test_set_mainlayer_interface_bf_layer(detailled)
     $     result(test_validated)

          implicit none
          
          logical, intent(in) :: detailled
          logical             :: test_validated

          
          integer :: k

          
          test_validated = .true.


          do k=1,12

             test_loc = perform_test_set_mainlayer_interface_bf_layer(
     $            k,detailled)
             test_validated = test_validated.and.test_loc
             
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_set_mainlayer_interface_bf_layer


        function perform_test_set_mainlayer_interface_bf_layer(
     $     test_id,
     $     detailled)
     $     result(test_validated)

          implicit none

          integer, intent(in) :: test_id
          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mainlayer_interface_basic) :: mainlayer_interface_used
          integer, dimension(5,2)         :: test_data
          integer                         :: k
          logical                         :: test_loc
          

          test_validated = .true.


          call mainlayer_interface_used%ini()


          !input+output
          call get_param_test_set_mainlayer_interface_bf_layer(
     $         test_id,
     $         mainlayer_interface_used,
     $         test_data)


          !validation
          do k=1,2

             select case(test_data(1,k))
               case(NE_interface_type)
                  select case(test_data(2,k))
                    case(N)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%NE_interface_N_ptr,
     $                      test_data(3:5,k),
     $                      detailled)
                       
                    case(E)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%NE_interface_E_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case default
                       call error_mainlayer_interface_incompatible(
     $                      'test_mainlayer_interface_basic',
     $                      'perform_test_set_mainlayer_interface_bf_layer',
     $                      test_data(1,k),
     $                      test_data(2,k))

                  end select

               case(NW_interface_type)
                  select case(test_data(2,k))
                    case(N)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%NW_interface_N_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case(W)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%NW_interface_W_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case default
                       call error_mainlayer_interface_incompatible(
     $                      'test_mainlayer_interface_basic',
     $                      'perform_test_set_mainlayer_interface_bf_layer',
     $                      test_data(1,k),
     $                      test_data(2,k))
                  end select

               case(SE_interface_type)
                  select case(test_data(2,k))
                    case(S)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%SE_interface_S_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case(E)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%SE_interface_E_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case default
                       call error_mainlayer_interface_incompatible(
     $                      'test_mainlayer_interface_basic',
     $                      'perform_test_set_mainlayer_interface_bf_layer',
     $                      test_data(1,k),
     $                      test_data(2,k))
                  end select

               case(SW_interface_type)
                  select case(test_data(2,k))
                    case(S)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%SW_interface_S_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case(W)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%SW_interface_W_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case default
                       call error_mainlayer_interface_incompatible(
     $                      'test_mainlayer_interface_basic',
     $                      'perform_test_set_mainlayer_interface_bf_layer',
     $                      test_data(1,k),
     $                      test_data(2,k))
                  end select

               case default
                  call error_mainlayer_interface_type(
     $                 'test_mainlayer_interface_basic',
     $                 'perform_test_set_mainlayer_interface_bf_layer',
     $                 test_data(1,k))
             end select

             test_validated = test_validated.and.test_loc

          end do


        end function perform_test_set_mainlayer_interface_bf_layer


        subroutine get_param_test_set_mainlayer_interface_bf_layer(
     $     test_id,
     $     mainlayer_interface_used,
     $     test_data)

          implicit none

          integer                        , intent(in)    :: test_id
          type(mainlayer_interface_basic), intent(inout) :: mainlayer_interface_used
          integer, dimension(5,2)        , intent(out)   :: test_data
          

          
          type(bf_sublayer), pointer :: bf_sublayer_ptr1
          type(bf_sublayer), pointer :: bf_sublayer_ptr2


          select case(test_id)
            case(1)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(N)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr1)

               test_data(:,1) = [NE_interface_type,N,1,2,0]
               test_data(:,2) = [NE_interface_type,E,0,2,0]

            case(2)
               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(E)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [NE_interface_type,N,0,2,0]
               test_data(:,2) = [NE_interface_type,E,1,2,0]

            case(3)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(N)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr1)

               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(E)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [NE_interface_type,N,1,2,1]
               test_data(:,2) = [NE_interface_type,E,1,2,1]

            case(4)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(N)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr1)

               test_data(:,1) = [NW_interface_type,N,1,1,0]
               test_data(:,2) = [NW_interface_type,W,0,2,0]

            case(5)
               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(W)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [NW_interface_type,N,0,1,0]
               test_data(:,2) = [NW_interface_type,W,1,2,0]

            case(6)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(N)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr1)

               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(W)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [NW_interface_type,N,1,1,1]
               test_data(:,2) = [NW_interface_type,W,1,2,1]

            case(7)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(S)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr1)

               test_data(:,1) = [SE_interface_type,S,1,2,0]
               test_data(:,2) = [SE_interface_type,E,0,1,0]

            case(8)
               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(E)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [SE_interface_type,S,0,2,0]
               test_data(:,2) = [SE_interface_type,E,1,1,0]

            case(9)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(S)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr1)

               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(E)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [SE_interface_type,S,1,2,1]
               test_data(:,2) = [SE_interface_type,E,1,1,1]

            case(10)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(S)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr1)

               test_data(:,1) = [SW_interface_type,S,1,1,0]
               test_data(:,2) = [SW_interface_type,W,0,1,0]

            case(11)
               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(W)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [SW_interface_type,S,0,1,0]
               test_data(:,2) = [SW_interface_type,W,1,1,0]

            case(12)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(S)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr1)

               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(W)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [SW_interface_type,S,1,1,1]
               test_data(:,2) = [SW_interface_type,W,1,1,1]

          end select

        end subroutine get_param_test_set_mainlayer_interface_bf_layer


        function verify_test_set_mainlayer_interface_bf_layer(
     $     bf_sublayer_ptr,
     $     test_data,
     $     detailled)
     $     result(test_validated)

          implicit none

          type(bf_sublayer), pointer, intent(in) :: bf_sublayer_ptr
          integer, dimension(3)     , intent(in) :: test_data
          logical                   , intent(in) :: detailled
          logical                                :: test_validated


          test_validated = .true.


          if(test_data(1).eq.1) then
             test_loc = associated(bf_sublayer_ptr)
             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''bf_sublayer_ptr failed'')'
             end if

             if(test_loc) then
                select case(test_data(2))
                  case(1)
                     test_loc = bf_sublayer_ptr%can_exchange_with_neighbor1().eqv.(test_data(3).eq.1)
                     test_validated = test_validated.and.test_loc

                     if(detailled.and.(.not.test_loc)) then
                        print '(''exchange_with_neighbor1 failed'')'
                     end if

                  case(2)
                     test_loc = bf_sublayer_ptr%can_exchange_with_neighbor2().eqv.(test_data(3).eq.1)
                     test_validated = test_validated.and.test_loc

                     if(detailled.and.(.not.test_loc)) then
                        print '(''exchange_with_neighbor2 failed'')'
                     end if

                end select
             end if

          else
             test_loc = .not.associated(bf_sublayer_ptr)
             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''bf_sublayer_ptr failed'')'
             end if
          end if

        end function verify_test_set_mainlayer_interface_bf_layer


        function test_remove_mainlayer_interface_bf_layer(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,16
            
             test_loc = perform_test_remove_mainlayer_interface_bf_layer(
     $            k,detailled)
             test_validated = test_validated.and.test_loc
             
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed '')', k
             end if
             
          end do

        end function test_remove_mainlayer_interface_bf_layer


        function perform_test_remove_mainlayer_interface_bf_layer(
     $     test_id, detailled)
     $     result(test_validated)

          implicit none

          integer, intent(in) :: test_id
          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mainlayer_interface_basic) :: mainlayer_interface_used
          integer, dimension(5,2)         :: test_data
          integer                         :: k
          logical                         :: test_loc

          test_validated = .true.

          
          call mainlayer_interface_used%ini()

          
          call get_param_test_remove_mainlayer_interface_bf_layer(
     $         test_id,
     $         mainlayer_interface_used,
     $         test_data)

          
          !validation
          do k=1,2

             select case(test_data(1,k))
               case(NE_interface_type)
                  select case(test_data(2,k))
                    case(N)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%NE_interface_N_ptr,
     $                      test_data(3:5,k),
     $                      detailled)
                       
                    case(E)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%NE_interface_E_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case default
                       call error_mainlayer_interface_incompatible(
     $                      'test_mainlayer_interface_basic',
     $                      'perform_test_set_mainlayer_interface_bf_layer',
     $                      test_data(1,k),
     $                      test_data(2,k))

                  end select

               case(NW_interface_type)
                  select case(test_data(2,k))
                    case(N)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%NW_interface_N_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case(W)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%NW_interface_W_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case default
                       call error_mainlayer_interface_incompatible(
     $                      'test_mainlayer_interface_basic',
     $                      'perform_test_set_mainlayer_interface_bf_layer',
     $                      test_data(1,k),
     $                      test_data(2,k))
                  end select

               case(SE_interface_type)
                  select case(test_data(2,k))
                    case(S)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%SE_interface_S_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case(E)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%SE_interface_E_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case default
                       call error_mainlayer_interface_incompatible(
     $                      'test_mainlayer_interface_basic',
     $                      'perform_test_set_mainlayer_interface_bf_layer',
     $                      test_data(1,k),
     $                      test_data(2,k))
                  end select

               case(SW_interface_type)
                  select case(test_data(2,k))
                    case(S)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%SW_interface_S_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case(W)
                       test_loc = verify_test_set_mainlayer_interface_bf_layer(
     $                      mainlayer_interface_used%SW_interface_W_ptr,
     $                      test_data(3:5,k),
     $                      detailled)

                    case default
                       call error_mainlayer_interface_incompatible(
     $                      'test_mainlayer_interface_basic',
     $                      'perform_test_set_mainlayer_interface_bf_layer',
     $                      test_data(1,k),
     $                      test_data(2,k))
                  end select

               case default
                  call error_mainlayer_interface_type(
     $                 'test_mainlayer_interface_basic',
     $                 'perform_test_set_mainlayer_interface_bf_layer',
     $                 test_data(1,k))
             end select

             test_validated = test_validated.and.test_loc

          end do

        end function perform_test_remove_mainlayer_interface_bf_layer


        subroutine get_param_test_remove_mainlayer_interface_bf_layer(
     $     test_id,
     $     mainlayer_interface_used,
     $     test_data)

          implicit none

          integer                        , intent(in)    :: test_id
          type(mainlayer_interface_basic), intent(inout) :: mainlayer_interface_used
          integer, dimension(5,2)        , intent(out)   :: test_data
          

          
          type(bf_sublayer), pointer :: bf_sublayer_ptr1
          type(bf_sublayer), pointer :: bf_sublayer_ptr2


          select case(test_id)

          !NE_interface
          !------------------------------------------------------------
            case(1)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(N)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr1)
               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr1)

               test_data(:,1) = [NE_interface_type,N,0,2,0]
               test_data(:,2) = [NE_interface_type,E,0,2,0]

            case(2)
               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(E)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr2)
               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [NE_interface_type,N,0,2,0]
               test_data(:,2) = [NE_interface_type,E,0,2,0]

            case(3)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(N)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr1)

               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(E)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr2)

               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr1)

               test_data(:,1) = [NE_interface_type,N,0,2,0]
               test_data(:,2) = [NE_interface_type,E,1,2,0]

            case(4)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(N)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr1)

               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(E)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr2)

               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [NE_interface_type,N,1,2,0]
               test_data(:,2) = [NE_interface_type,E,0,2,0]


          !NW_interface
          !------------------------------------------------------------
            case(5)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(N)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr1)
               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr1)

               test_data(:,1) = [NW_interface_type,N,0,1,0]
               test_data(:,2) = [NW_interface_type,W,0,2,0]

            case(6)
               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(W)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr2)
               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [NW_interface_type,N,0,1,0]
               test_data(:,2) = [NW_interface_type,W,0,2,0]

            case(7)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(N)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr1)

               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(W)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr2)

               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr1)

               test_data(:,1) = [NW_interface_type,N,0,1,0]
               test_data(:,2) = [NW_interface_type,W,1,2,0]

            case(8)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(N)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr1)

               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(W)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr2)

               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [NW_interface_type,N,1,1,0]
               test_data(:,2) = [NW_interface_type,W,0,2,0]


         !SE_interface
         !------------------------------------------------------------
            case(9)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(S)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr1)
               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr1)

               test_data(:,1) = [SE_interface_type,S,0,2,0]
               test_data(:,2) = [SE_interface_type,E,0,1,0]


            case(10)
               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(E)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr2)
               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [SE_interface_type,S,0,2,0]
               test_data(:,2) = [SE_interface_type,E,0,1,0]

            case(11)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(S)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr1)

               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(E)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr2)

               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr1)

               test_data(:,1) = [SE_interface_type,S,0,2,0]
               test_data(:,2) = [SE_interface_type,E,1,1,0]

            case(12)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(S)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr1)

               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(E)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr2)

               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [SE_interface_type,S,1,2,0]
               test_data(:,2) = [SE_interface_type,E,0,1,0]


         !SW_interface
         !------------------------------------------------------------
            case(13)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(S)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr1)
               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr1)

               test_data(:,1) = [SW_interface_type,S,0,1,0]
               test_data(:,2) = [SW_interface_type,W,0,1,0]

            case(14)
               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(W)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr2)
               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [SW_interface_type,S,0,1,0]
               test_data(:,2) = [SW_interface_type,W,0,1,0]

            case(15)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(S)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr1)

               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(W)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr2)

               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr1)

               test_data(:,1) = [SW_interface_type,S,0,1,0]
               test_data(:,2) = [SW_interface_type,W,1,1,0]

            case(16)
               allocate(bf_sublayer_ptr1)
               call bf_sublayer_ptr1%ini(S)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr1)

               allocate(bf_sublayer_ptr2)
               call bf_sublayer_ptr2%ini(W)
               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr2)

               call mainlayer_interface_used%remove_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr2)

               test_data(:,1) = [SW_interface_type,S,1,1,0]
               test_data(:,2) = [SW_interface_type,W,0,1,0]

          end select

        end subroutine get_param_test_remove_mainlayer_interface_bf_layer

      end program test_mainlayer_interface_basic
