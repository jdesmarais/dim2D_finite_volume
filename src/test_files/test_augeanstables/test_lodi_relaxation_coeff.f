      program test_lodi_relaxation_coeff

        use lodi_relaxation_coeff_module, only:
     $       get_relaxation_normal_velocity,
     $       get_relaxation_trans_velocity,
     $       get_relaxation_temperature,
     $       get_relaxation_pressure

        use ns2d_parameters, only :
     $       gamma,
     $       mach_infty

        use parameters_constant, only :
     $       left

        use parameters_input, only :
     $       sigma_P

        use parameters_kind, only :
     $       rkind

        implicit none

        real(rkind), parameter :: l_domain_n = 2.5d0
        real(rkind), parameter :: M_un_infty = 0.1d0
        real(rkind), parameter :: M_local    = 0.5d0


        real(rkind), parameter :: r_normal_velocity = 9.9d0
        real(rkind), parameter :: r_trans_velocity  = 2.0d0
        real(rkind), parameter :: r_temperature     = 120.0d0
        real(rkind), parameter :: r_pressure        = 0.75d0
        
        
        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .false.
        if(
     $       (.not.is_test_validated(gamma,5.0d0/3.0d0,detailled)).or.
     $       (.not.is_test_validated(mach_infty,0.1d0,detailled)).or.
     $       (.not.is_test_validated(sigma_P,0.25d0,detailled))) then

           print '(''the test requires: '')'
           print '(''gamma=5/3'')'
           print '(''mach_infty=0.2'')'
           print '(''sigma_P=0.25'')'
           stop ''

        end if


        test_validated = .true.

        !test relaxation normal velocity
        test_loc = is_test_validated(
     $       get_relaxation_normal_velocity(l_domain_n,M_un_infty,left),
     $       r_normal_velocity,
     $       detailled)
        test_validated = test_validated.and.test_loc
        if(.not.detailled) then
           print '(''get_relaxation_normal_velocity: '',L2)', test_loc
        end if

        !test relaxation transverse velocity
        test_loc = is_test_validated(
     $       get_relaxation_trans_velocity(l_domain_n,M_local),
     $       r_trans_velocity,
     $       detailled)
        test_validated = test_validated.and.test_loc
        if(.not.detailled) then
           print '(''get_relaxation_trans_velocity: '',L2)', test_loc
        end if

        !test relaxation temperature
        test_loc = is_test_validated(
     $       get_relaxation_temperature(l_domain_n,M_local),
     $       r_temperature,
     $       detailled)
        test_validated = test_validated.and.test_loc
        if(.not.detailled) then
           print '(''get_relaxation_temperature: '',L2)', test_loc
        end if

        !test relaxation pressure
        test_loc = is_test_validated(
     $       get_relaxation_pressure(l_domain_n,M_local),
     $       r_pressure,
     $       detailled)
        test_validated = test_validated.and.test_loc
        if(.not.detailled) then
           print '(''get_relaxation_pressure: '',L2)', test_loc
        end if

        print '(''------------------'')'
        print '(''test_validated: '',L2)', test_validated

        contains


        function is_test_validated(var,cst,detailled) result(test_validated)

          implicit none
          
          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical    , intent(in) :: detailled
          logical                 :: test_validated
          
          if(detailled) then
             print *, int(var*1e5)
             print *, int(cst*1e5)
          end if
          
          test_validated=abs(
     $         nint(var*1e5)-
     $         nint(cst*1e5)).le.1
        
        end function is_test_validated

      end program test_lodi_relaxation_coeff
