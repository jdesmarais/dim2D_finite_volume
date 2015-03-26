      program test_icr_path_list

        use icr_path_list_class, only :
     $     icr_path_list

        use icr_path_chain_class, only :
     $     icr_path_chain

        use parameters_constant, only :
     $       N,S



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


        test_loc = test_add_path(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_add_path: '',L1)', test_loc
        print '()'


        test_loc = test_remove_path(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_remove_path: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_ini(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(icr_path_list) :: icr_path_list_used
          logical             :: test_loc

          test_validated = .true.


          !output
          call icr_path_list_used%ini()

          !validation
          test_loc = icr_path_list_used%nb_paths.eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_paths failed'')'
          end if

          test_loc = .not.associated(icr_path_list_used%head_path)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test head_path failed'')'
          end if
          
          test_loc = .not.associated(icr_path_list_used%tail_path)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test tail_path failed'')'
          end if

        end function test_ini


        function test_add_path(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(icr_path_list)           :: icr_path_list_used
          type(icr_path_chain), pointer :: added_path
          type(icr_path_chain), pointer :: current_path
          logical                       :: test_loc


          test_validated = .true.


          !input
          call icr_path_list_used%ini()

          
          !output: first path
          added_path => icr_path_list_used%add_path()
          added_path%mainlayer_id = N

          !validation: first path
          test_loc = icr_path_list_used%get_nb_paths().eq.1
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_paths(1) failed'')'
          end if

          current_path => icr_path_list_used%get_head_path()
          test_loc = current_path%mainlayer_id.eq.N
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test head_path(1) failed'')'
          end if

          current_path => icr_path_list_used%get_tail_path()
          test_loc = current_path%mainlayer_id.eq.N
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test tail_path(1) failed'')'
          end if


          !output: second path
          added_path => icr_path_list_used%add_path()
          added_path%mainlayer_id = S

          !validation: second path
          test_loc = icr_path_list_used%get_nb_paths().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_paths(2) failed'')'
          end if

          current_path => icr_path_list_used%get_head_path()
          test_loc = current_path%mainlayer_id.eq.N
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test head_path(2) failed'')'
          end if

          current_path => icr_path_list_used%get_tail_path()
          test_loc = current_path%mainlayer_id.eq.S
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test tail_path(2) failed'')'
          end if

        end function test_add_path


        function test_remove_path(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(icr_path_list)           :: icr_path_list_used
          type(icr_path_chain), pointer :: added_path1
          type(icr_path_chain), pointer :: added_path2
          type(icr_path_chain), pointer :: current_path
          logical                       :: test_loc


          test_validated = .true.


          !test 1: one path in the list
          !============================================================
          !we remove the only path in the list
          !------------------------------------------------------------

          !input
          call icr_path_list_used%ini()
          added_path1 => icr_path_list_used%add_path()
          added_path1%mainlayer_id = N

          !output
          call icr_path_list_used%remove_path(added_path1)

          !validation: first path
          test_loc = icr_path_list_used%get_nb_paths().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_paths(1) failed'')'
          end if

          current_path => icr_path_list_used%get_head_path()
          test_loc = .not.associated(current_path)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test head_path(1) failed'')'
          end if

          current_path => icr_path_list_used%get_tail_path()
          test_loc = .not.associated(current_path)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test tail_path(1) failed'')'
          end if


          !test 2: two paths in the list
          !============================================================
          !we remove the path2 and then path1
          !------------------------------------------------------------

          !input: second path
          call icr_path_list_used%ini()
          added_path1 => icr_path_list_used%add_path()
          added_path1%mainlayer_id = N
          added_path2 => icr_path_list_used%add_path()
          added_path2%mainlayer_id = S

          !output:1
          call icr_path_list_used%remove_path(added_path2)

          !validation:1
          test_loc = icr_path_list_used%get_nb_paths().eq.1
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_paths(2,1) failed'')'
          end if

          current_path => icr_path_list_used%get_head_path()
          test_loc = current_path%mainlayer_id.eq.N
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test head_path(2,1) failed'')'
          end if

          current_path => icr_path_list_used%get_tail_path()
          test_loc = current_path%mainlayer_id.eq.N
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test tail_path(2,1) failed'')'
          end if

          !output:2
          call icr_path_list_used%remove_path(added_path1)

          !validation:2
          test_loc = icr_path_list_used%get_nb_paths().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_paths(2,2) failed'')'
          end if

          current_path => icr_path_list_used%get_head_path()
          test_loc = .not.associated(current_path)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test head_path(2,2) failed'')'
          end if

          current_path => icr_path_list_used%get_tail_path()
          test_loc = .not.associated(current_path)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test tail_path(2,2) failed'')'
          end if


          !we remove the path1 and then path2
          !------------------------------------------------------------

          !input: second path
          call icr_path_list_used%ini()
          added_path1 => icr_path_list_used%add_path()
          added_path1%mainlayer_id = N
          added_path2 => icr_path_list_used%add_path()
          added_path2%mainlayer_id = S

          !output:1
          call icr_path_list_used%remove_path(added_path1)

          !validation:1
          test_loc = icr_path_list_used%get_nb_paths().eq.1
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_paths(3,1) failed'')'
          end if

          current_path => icr_path_list_used%get_head_path()
          test_loc = current_path%mainlayer_id.eq.S
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test head_path(3,1) failed'')'
          end if

          current_path => icr_path_list_used%get_tail_path()
          test_loc = current_path%mainlayer_id.eq.S
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test tail_path(3,1) failed'')'
          end if

          !output:2
          call icr_path_list_used%remove_path(added_path2)

          !validation:2
          test_loc = icr_path_list_used%get_nb_paths().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_paths(3,2) failed'')'
          end if

          current_path => icr_path_list_used%get_head_path()
          test_loc = .not.associated(current_path)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test head_path(3,2) failed'')'
          end if

          current_path => icr_path_list_used%get_tail_path()
          test_loc = .not.associated(current_path)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test tail_path(3,2) failed'')'
          end if

        end function test_remove_path

      end program test_icr_path_list
