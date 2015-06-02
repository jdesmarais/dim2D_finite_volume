      program test_bc_operators_wall_reflection

        use bc_operators_class, only :
     $     bc_operators

        use parameters_constant, only :
     $     

        use parameters_kind, only :
     $     ikind,
     $     rkind


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

        print '(''test_validated: '',L1)', test_validated

        contains

        function test_ini(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bc_operators) :: bc_used
          type(pmodel_eq)    :: pm_used
          logical :: test_loc

          call bc_used%ini(pm_used)

          test_loc = is_int_vector_validated(
     $         bc_used%bc_type,
     $         [bc_flux_and_node_choice, bc_flux_and_node_choice,
     $          bc_nodes_choice        , bc_flux_and_node_choice],
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_validated) then
             print '(test bc_type failed)'
          end if
          

          test_loc = is_int_vector_validated(
     $         bc_used%wall_prefactor,
     $         [],
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_validated) then
             print '(test bc_type failed)'
          end if


        end function test_ini


        subroutine test_inputs()

          implicit none

          logical :: test_loc
          logical :: test_inputs
          type(pmodel_eq) :: pm_used

          test_inputs = .true.

          ! verify that the physical model used is DIM
          test_loc = trim(pm_used%get_model_name()).eq."DIM2D"
          test_inputs = test_inputs.and.test_loc
          print '(''pm_model: '',L1)', test_loc


          if(.not.test_inputs) then
             print '(''inputs not correct for the test'')'
             stop ''
          end if

        end subroutine test_inputs


      end program test_bc_operators_wall_reflection
