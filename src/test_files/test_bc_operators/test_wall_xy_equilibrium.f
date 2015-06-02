      program test_wall_xy_equilibrium

        use check_data_module, only :
     $     is_real_validated,
     $     is_real_vector_validated

        use dim2d_parameters, only :
     $       cv_r,
     $       We,
     $       Pr

        use dim2d_state_eq_module, only :
     $       get_mass_density_vapor,
     $       get_mass_density_liquid

        use parameters_kind, only :
     $       rkind

        use ridders_method_fcts_module, only :
     $       root_fct2,
     $       root_fct4

        use wall_xy_equilibrium_module, only :
     $       dmddx,
     $       temperature,
     $       md_average,
     $       temperature_average,
     $       dwallInternalEnergy_dmd,
     $       wall_x_equilibrium_root_fct,
     $       wall_x_root_fct,
     $       get_wall_x_root_brackets

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled = .true.
        test_validated = .true.

        call test_inputs()

        test_loc = test_dmddx(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_dmddx: '',L1)', test_loc
        print '()'

        test_loc = test_temperature(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_temperature: '',L1)', test_loc
        print '()'

        test_loc = test_md_average(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_md_average: '',L1)', test_loc
        print '()'

        test_loc = test_temperature_average(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_temperature_average: '',L1)', test_loc
        print '()'

        test_loc = test_dwallInternalEnergy_dmd(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_dwallInternalEnergy_dmd: '',L1)', test_loc
        print '()'

        test_loc = test_wall_x_equilibrium_root_fct(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_wall_x_equilibrium_root_fct: '',L1)', test_loc
        print '()'

        test_loc = test_wall_x_root_fct(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_wall_x_root_fct: '',L1)', test_loc
        print '()'

        test_loc = test_get_wall_x_root_brackets(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_wall_x_root_brackets: '',L1)', test_loc
        print '()'

        print '(''test_validated: '',L1)', test_validated

        contains

        function test_get_wall_x_root_brackets(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          logical :: test_loc

          type(root_fct2) :: root_fct2_used
          type(root_fct4) :: root_fct4_used

          real(rkind) :: md_vap
          real(rkind) :: md_liq

          test_validated = .true.

          md_vap = get_mass_density_vapor(0.95d0)
          md_liq = get_mass_density_liquid(0.95d0)

          test_loc = is_real_vector_validated(
     $         get_wall_x_root_brackets(root_fct2_used,0.95,0.0),
     $         [md_vap,md_liq],
     $         detailled)
          test_validated = test_validated.and.test_loc


          test_loc = is_real_vector_validated(
     $         get_wall_x_root_brackets(root_fct4_used,0.95,1.5d0),
     $         [1.5d0,2.01d0],
     $         detailled)
          test_validated = test_validated.and.test_loc          

        end function test_get_wall_x_root_brackets
      

        function test_wall_x_root_fct(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(wall_x_root_fct) :: wall_x_root_fct_used

          call wall_x_root_fct_used%ini(
     $         [0.2d0,0.5d0],
     $         [1.2d0,0.9d0],
     $         0.2d0,
     $         0.3d0,
     $         1.6d0,
     $         0.5d0,
     $         0.8d0,
     $         45.0d0,
     $         0.2d0)
          
          test_validated = is_real_validated(
     $         wall_x_root_fct_used%f(1.5d0),
     $         wall_x_equilibrium_root_fct(
     $         1.5d0,
     $         [0.2d0,0.5d0],
     $         [1.2d0,0.9d0],
     $         0.2d0,
     $         0.3d0,
     $         1.6d0,
     $         0.5d0,
     $         0.8d0,
     $         45.0d0,
     $         0.2d0),
     $         detailled)

        end function test_wall_x_root_fct


        function test_wall_x_equilibrium_root_fct(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = is_real_validated(
     $         wall_x_equilibrium_root_fct(
     $         1.5d0,
     $         [0.2d0,0.5d0],
     $         [1.2d0,0.9d0],
     $         0.2d0,
     $         0.3d0,
     $         1.6d0,
     $         0.5d0,
     $         0.8d0,
     $         45.0d0,
     $         0.2d0),
     $         -0.375003199415867d0,
     $         detailled)

        end function test_wall_x_equilibrium_root_fct

        function test_dwallInternalEnergy_dmd(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = is_real_validated(
     $         dwallInternalEnergy_dmd(1.5d0,0.95d0,45.0d0),
     $         0.007151129d0,
     $         detailled)

        end function test_dwallInternalEnergy_dmd
        
        function test_temperature_average(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = is_real_validated(
     $         temperature_average(1.5d0,0.2d0,5.0d0),
     $         11.5d0,
     $         detailled)

        end function test_temperature_average

        function test_md_average(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = is_real_validated(
     $         md_average(1.5d0,0.5d0),
     $         1.0d0,
     $         detailled)

        end function test_md_average


        function test_temperature(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = is_real_validated(
     $         temperature(1.5d0,50.0d0,-0.1d0,0.2d0,6.3d0),
     $         -0.29875d0,
     $         detailled)

        end function test_temperature


        function test_dmddx(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          test_validated = is_real_validated(
     $         dmddx(0.5d0,-2.0d0,3.0d0),
     $         5.0d0,
     $         detailled)

        end function test_dmddx


        subroutine test_inputs()

          implicit none

          logical :: test_loc
          logical :: test_validated

          test_validated = .true.

          test_loc = is_real_validated(cv_r,2.5d0,detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''cv_r should equal 2.5d0'')'
          end if

          test_loc = is_real_validated(We,10.0d0,detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''we should equal 10.0d0'')'
          end if

          test_loc = is_real_validated(Pr,20.0d0,detailled)
          test_validated = test_validated.and.test_loc
          if(.not.test_loc) then
             print '(''Pr should equal 20.0d0'')'
          end if

          if(.not.test_validated) then
             stop ''
          end if

        end subroutine test_inputs

      end program test_wall_xy_equilibrium
