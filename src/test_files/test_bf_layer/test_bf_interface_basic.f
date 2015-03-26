      program test_bf_interface_basic

        use bf_interface_basic_class, only :
     $     bf_interface_basic

        use bf_mainlayer_class, only :
     $       bf_mainlayer

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny

        use parameters_kind, only :
     $       rkind

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


        test_loc = test_get_mainlayer_ptr(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_mainlayer_ptr: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_ini(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_basic)   :: bf_interface_used
          real(rkind), dimension(nx) :: interior_x_map
          real(rkind), dimension(ny) :: interior_y_map
          integer                    :: i
          logical                    :: test_loc


          test_validated = .true.


          !output
          call bf_interface_used%ini(
     $         interior_x_map, interior_y_map)


          !validation of the initialization of the mainlayer pointers
          do i=1,4

             test_loc = .not.associated(
     $            bf_interface_used%mainlayer_pointers(i)%get_ptr())
             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''mainlayer_ptr('',I2,'') failed'')',i
             end if

          end do

          !validation of the initialization of the mainlayer interfaces
          test_loc = .not.associated(
     $         bf_interface_used%mainlayer_interfaces%NE_interface_N_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NE_interface_N failed'')'
          end if

          test_loc = .not.associated(
     $         bf_interface_used%mainlayer_interfaces%NE_interface_E_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NE_interface_E failed'')'
          end if

          test_loc = .not.associated(
     $         bf_interface_used%mainlayer_interfaces%NW_interface_N_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NW_interface_N failed'')'
          end if

          test_loc = .not.associated(
     $         bf_interface_used%mainlayer_interfaces%NW_interface_W_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NW_interface_W failed'')'
          end if

          test_loc = .not.associated(
     $         bf_interface_used%mainlayer_interfaces%SE_interface_S_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SE_interface_S failed'')'
          end if

          test_loc = .not.associated(
     $         bf_interface_used%mainlayer_interfaces%SE_interface_E_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SE_interface_E failed'')'
          end if

          test_loc = .not.associated(
     $         bf_interface_used%mainlayer_interfaces%SW_interface_S_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SW_interface_S failed'')'
          end if

          test_loc = .not.associated(
     $         bf_interface_used%mainlayer_interfaces%SW_interface_W_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SW_interface_W failed'')'
          end if          

        end function test_ini


        function test_get_mainlayer_ptr(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_basic)    :: bf_interface_used
          real(rkind), dimension(nx)  :: interior_x_map
          real(rkind), dimension(ny)  :: interior_y_map
          type(bf_mainlayer), pointer :: bf_mainlayer_ptr
          logical                     :: test_loc


          test_validated = .true.


          !input
          call bf_interface_used%ini(
     $         interior_x_map, interior_y_map)

          !N reference
          bf_mainlayer_ptr => bf_interface_used%get_mainlayer_ptr(N)

          test_loc = .not.associated(bf_mainlayer_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''mainlayer_pointer(N) failed'')'
          end if

          !S reference
          bf_mainlayer_ptr => bf_interface_used%get_mainlayer_ptr(S)

          test_loc = .not.associated(bf_mainlayer_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''mainlayer_pointer(S) failed'')'
          end if

          !E reference
          bf_mainlayer_ptr => bf_interface_used%get_mainlayer_ptr(E)

          test_loc = .not.associated(bf_mainlayer_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''mainlayer_pointer(E) failed'')'
          end if

          !W reference
          bf_mainlayer_ptr => bf_interface_used%get_mainlayer_ptr(W)

          test_loc = .not.associated(bf_mainlayer_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''mainlayer_pointer(W) failed'')'
          end if

        end function test_get_mainlayer_ptr

      end program test_bf_interface_basic
