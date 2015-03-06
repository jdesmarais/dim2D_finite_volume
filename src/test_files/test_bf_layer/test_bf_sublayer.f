      program test_bf_sublayer

        use bf_sublayer_class, only :
     $     bf_sublayer

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled = .true.
        test_validated  =.true.


        test_loc = test_ini(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_ini: '',L1)', test_loc
        print '()'


        contains

        function test_ini(detailled)
     $       result(test_validated)

          implicit none

          logical :: detailled
          logical :: test_validated


          type(bf_sublayer) :: bf_sublayer_used
          logical           :: test_loc

          type(bf_sublayer), pointer :: prev
          type(bf_sublayer), pointer :: next

          test_validated = .true.


          call bf_sublayer_used%ini(1)

          prev => bf_sublayer_used%get_prev()
          next => bf_sublayer_used%get_next()

          test_loc = .not.associated(prev)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test prev not allocated failed'')'
          end if

          test_loc = .not.associated(next)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test next not allocated failed'')'
          end if

        end function test_ini


      end program test_bf_sublayer
