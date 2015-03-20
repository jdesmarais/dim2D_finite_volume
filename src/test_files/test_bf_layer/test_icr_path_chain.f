      program test_icr_path_chain

        use icr_path_chain_class, only :
     $     icr_path_chain

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


          type(icr_path_chain) :: icr_path_chain_used
          logical              :: test_loc

          type(icr_path_chain), pointer :: prev
          type(icr_path_chain), pointer :: next

          test_validated = .true.


          call icr_path_chain_used%ini()

          prev => icr_path_chain_used%get_prev()
          next => icr_path_chain_used%get_next()

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


      end program test_icr_path_chain
