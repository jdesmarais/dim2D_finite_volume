      program test_dcr_chain

        use dcr_chain_class, only :
     $     dcr_chain

        use bf_sublayer_class, only :
     $       bf_sublayer


        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        test_loc = test_get_ptr(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_ptr: '',L1)', test_loc
        print '()'

        test_loc = test_set_ptr(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_set_ptr: '',L1)', test_loc
        print '()'

        test_loc = test_associated_ptr(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_associated_ptr: '',L1)', test_loc
        print '()'


        contains


        function test_get_ptr(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(dcr_chain)            :: dcr_chain_used
          type(bf_sublayer), pointer :: bf_sublayer1_ptr
          type(bf_sublayer), pointer :: bf_sublayer2_ptr
          logical                    :: test_loc


          test_validated = .true.


          bf_sublayer1_ptr => dcr_chain_used%get_ptr()
          test_loc = .not.associated(bf_sublayer1_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test get_ptr(1) failed'')'
          end if

          
          allocate(bf_sublayer2_ptr)
          dcr_chain_used%ptr => bf_sublayer2_ptr
          bf_sublayer1_ptr => dcr_chain_used%get_ptr()
          test_loc = associated(bf_sublayer1_ptr,bf_sublayer2_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test get_ptr(2) failed'')'
          end if

          deallocate(bf_sublayer2_ptr)

        end function test_get_ptr


        function test_set_ptr(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(dcr_chain)            :: dcr_chain_used
          type(bf_sublayer), pointer :: bf_sublayer1_ptr
          type(bf_sublayer), pointer :: bf_sublayer2_ptr
          logical                    :: test_loc


          test_validated = .true.


          call dcr_chain_used%set_ptr(bf_sublayer1_ptr)
          test_loc = .not.associated(dcr_chain_used%ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test set_ptr(1) failed'')'
          end if

          
          allocate(bf_sublayer2_ptr)
          call dcr_chain_used%set_ptr(bf_sublayer2_ptr)
          test_loc = associated(dcr_chain_used%ptr,bf_sublayer2_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test set_ptr(2) failed'')'
          end if

          deallocate(bf_sublayer2_ptr)

        end function test_set_ptr


        function test_associated_ptr(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(dcr_chain)            :: dcr_chain_used
          type(bf_sublayer), pointer :: bf_sublayer1_ptr
          type(bf_sublayer), pointer :: bf_sublayer2_ptr
          logical                    :: test_loc


          test_validated = .true.

          
          test_loc = .not.(dcr_chain_used%associated_ptr())
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test associated_ptr(1) failed'')'
          end if

          
          allocate(bf_sublayer1_ptr)
          dcr_chain_used%ptr => bf_sublayer1_ptr
          test_loc = dcr_chain_used%associated_ptr()
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test associated_ptr(2) failed'')'
          end if


          test_loc = dcr_chain_used%associated_ptr(bf_sublayer1_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test associated_ptr(3) failed'')'
          end if


          allocate(bf_sublayer2_ptr)
          test_loc = .not.(dcr_chain_used%associated_ptr(bf_sublayer2_ptr))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test associated_ptr(4) failed'')'
          end if

          deallocate(bf_sublayer1_ptr)
          deallocate(bf_sublayer2_ptr)

        end function test_associated_ptr

      end program test_dcr_chain
      
