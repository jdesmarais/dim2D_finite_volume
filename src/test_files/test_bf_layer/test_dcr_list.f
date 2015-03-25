      program test_dcr_list

        use dcr_list_class, only :
     $     dcr_list

        use bf_sublayer_class, only :
     $       bf_sublayer


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


        test_loc = test_add_ele(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_add_ele: '',L1)', test_loc
        print '()'

        
        test_loc = test_get_ele(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_ele: '',L1)', test_loc
        print '()'


        test_loc = test_does_not_contain(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_does_not_contain: '',L1)', test_loc
        print '()'


        test_loc = test_remove(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_remove: '',L1)', test_loc
        print '()'


        contains


        function test_ini(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(dcr_list) :: dcr_list_used
          logical        :: test_loc

          
          test_validated = .true.
          

          call dcr_list_used%ini(3)

          test_loc = dcr_list_used%nb_ele.eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_ele failed'')'
          end if

          
          test_loc = allocated(dcr_list_used%list)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''allocation failed'')'
          end if


          if(test_loc) then
             test_loc = size(dcr_list_used%list,1).eq.3
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''size failed'')'
             end if
          end if

        end function test_ini


        function test_add_ele(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(dcr_list)             :: dcr_list_used
          type(bf_sublayer), pointer :: bf_sublayer1_ptr
          type(bf_sublayer), pointer :: bf_sublayer2_ptr
          logical                    :: test_loc

          
          test_validated = .true.
          

          !add one element
          call dcr_list_used%ini(3)
          allocate(bf_sublayer1_ptr)
          call dcr_list_used%add_ele(bf_sublayer1_ptr)
          
          test_loc = dcr_list_used%nb_ele.eq.1
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_ele(1) failed'')'
          end if

          test_loc = dcr_list_used%list(1)%associated_ptr(bf_sublayer1_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''ele(1) failed'')'
          end if


          !add a second element
          allocate(bf_sublayer2_ptr)
          call dcr_list_used%add_ele(bf_sublayer2_ptr)
          
          test_loc = dcr_list_used%nb_ele.eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_ele(2) failed'')'
          end if

          test_loc = dcr_list_used%list(1)%associated_ptr(bf_sublayer1_ptr)
          test_validated = test_validated.and.test_loc
          test_loc = dcr_list_used%list(2)%associated_ptr(bf_sublayer2_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''ele(2) failed'')'
          end if
          

          !add a third element which is the same as the first one
          call dcr_list_used%add_ele(bf_sublayer1_ptr)
          
          test_loc = dcr_list_used%nb_ele.eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_ele(3) failed'')'
          end if

          test_loc = dcr_list_used%list(1)%associated_ptr(bf_sublayer1_ptr)
          test_validated = test_validated.and.test_loc
          test_loc = dcr_list_used%list(2)%associated_ptr(bf_sublayer2_ptr)
          test_validated = test_validated.and.test_loc
          test_loc = .not.(dcr_list_used%list(3)%associated_ptr())
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''ele(3) failed'')'
          end if

          call dcr_list_used%remove()

          deallocate(bf_sublayer1_ptr)
          deallocate(bf_sublayer2_ptr)

        end function test_add_ele


        function test_get_ele(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(dcr_list)             :: dcr_list_used
          type(bf_sublayer), pointer :: bf_sublayer1_ptr
          type(bf_sublayer), pointer :: bf_sublayer2_ptr
          type(bf_sublayer), pointer :: bf_sublayer_ptr
          logical                    :: test_loc

          
          test_validated = .true.
          

          !input
          call dcr_list_used%ini(3)
          allocate(bf_sublayer1_ptr)
          call dcr_list_used%add_ele(bf_sublayer1_ptr)
          allocate(bf_sublayer2_ptr)
          call dcr_list_used%add_ele(bf_sublayer2_ptr)
          

          !extract first element
          bf_sublayer_ptr => dcr_list_used%get_ele(1)

          test_loc = associated(bf_sublayer_ptr,bf_sublayer1_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''get_ele(1) failed'')'
          end if

          
          !extract second element
          bf_sublayer_ptr => dcr_list_used%get_ele(2)

          test_loc = associated(bf_sublayer_ptr,bf_sublayer2_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''get_ele(2) failed'')'
          end if

          call dcr_list_used%remove()

          deallocate(bf_sublayer1_ptr)
          deallocate(bf_sublayer2_ptr)

        end function test_get_ele


        function test_does_not_contain(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(dcr_list)             :: dcr_list_used
          type(bf_sublayer), pointer :: bf_sublayer1_ptr
          type(bf_sublayer), pointer :: bf_sublayer2_ptr
          type(bf_sublayer), pointer :: bf_sublayer3_ptr
          logical                    :: test_loc

          
          test_validated = .true.
          

          !input
          call dcr_list_used%ini(3)
          allocate(bf_sublayer1_ptr)
          call dcr_list_used%add_ele(bf_sublayer1_ptr)
          allocate(bf_sublayer2_ptr)
          call dcr_list_used%add_ele(bf_sublayer2_ptr)
          allocate(bf_sublayer3_ptr)
          

          !test look for first element
          test_loc = .not.(dcr_list_used%does_not_contain(bf_sublayer1_ptr))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test_does_not_contain(1) failed'')'
          end if


          !test look for second element
          test_loc = .not.(dcr_list_used%does_not_contain(bf_sublayer2_ptr))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test_does_not_contain(2) failed'')'
          end if


          !test look for an element which is not inside
          test_loc = dcr_list_used%does_not_contain(bf_sublayer3_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test_does_not_contain(3) failed'')'
          end if

          call dcr_list_used%remove()

          deallocate(bf_sublayer1_ptr)
          deallocate(bf_sublayer2_ptr)
          deallocate(bf_sublayer3_ptr)

        end function test_does_not_contain


        function test_remove(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(dcr_list)             :: dcr_list_used
          type(bf_sublayer), pointer :: bf_sublayer_ptr
          logical                    :: test_loc


          test_validated = .true.


          call dcr_list_used%remove()
          test_loc = .not.allocated(dcr_list_used%list)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test remove(1) failed'')'
          end if

          call dcr_list_used%ini(3)
          call dcr_list_used%remove()
          test_loc = .not.allocated(dcr_list_used%list)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test remove(2) failed'')'
          end if

          allocate(bf_sublayer_ptr)
          call dcr_list_used%ini(3)
          call dcr_list_used%add_ele(bf_sublayer_ptr)
          call dcr_list_used%remove()
          test_loc = .not.allocated(dcr_list_used%list)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test remove(3) failed'')'
          end if
          
        end function test_remove

      end program test_dcr_list
