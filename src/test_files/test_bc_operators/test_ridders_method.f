      program test_ridders_method

        use check_data_module, only :
     $       is_real_validated

        use ridders_method_fcts_module, only :
     $       root_fct1,
     $       root_fct2,
     $       root_fct3

        use ridders_method_module, only :
     $       root_fct_abstract,
     $       get_root_ridder_method

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        test_loc = test_get_root_ridder_method(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_root_ridder_method: '',L1)', test_loc
        print '()'
        
        print '(''test_validated: '',L1)', test_validated


        contains


        function test_get_root_ridder_method(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          logical         :: test_loc
          type(root_fct1) :: root_fct1_used
          type(root_fct2) :: root_fct2_used
          type(root_fct3) :: root_fct3_used

          test_validated = .true.


          !f(x) = 1+x [-2.0,0.0]
          test_loc = is_real_validated(
     $         get_root_ridder_method(root_fct1_used,-2.0d0,0.0d0,1.0e-12),
     $         -1.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''root_fct1, [-2,0]: failed'')'
          end if


          !f(x) = 1+x; [-2.0,1.0]
          test_loc = is_real_validated(
     $         get_root_ridder_method(root_fct1_used,-2.0d0,1.0d0,1.0e-12),
     $         -1.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''root_fct1, [-2,1]: failed'')'
          end if


          !f(x) = 1+x; [-2.0,1.0]
          test_loc = is_real_validated(
     $         get_root_ridder_method(root_fct1_used,-2.0d0,1.0d0,1.0e-12),
     $         -1.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''root_fct1, [-2,1]: failed'')'
          end if


          !f(x) = (x-1)(x-2); [-3.0,1.5]
          test_loc = is_real_validated(
     $         get_root_ridder_method(root_fct2_used,-3.0d0,1.5d0,1.0e-12),
     $         1.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''root_fct2, [-3,1.5]: failed'')'
          end if

          
          !f(x) = Tanh(x); [-3.0,2.0]
          test_loc = is_real_validated(
     $         get_root_ridder_method(root_fct3_used,-3.0d0,2.0d0,1.0e-12),
     $         0.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''root_fct3, [-3,2]: failed'')'
          end if

        end function test_get_root_ridder_method

      end program test_ridders_method
