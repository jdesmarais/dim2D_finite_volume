      program test_gaussian_perturbation

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated

        use field_class, only :
     $       field

        use gaussian_perturbation_module, only :
     $       add_gaussian_perturbation,
     $       compute_wave_numbers,
     $       compute_amplitudes,
     $       power_spectrum,
     $       compute_random_phases,
     $       compute_gaussian_intensity,
     $       smooth

        use parameters_constant, only :
     $       homogeneous_liquid

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       x_min,x_max,y_min,y_max,
     $       ntx,nty,npx,npy,
     $       T0,
     $       ic_choice,
     $       ic_perturbation_ac

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_smooth(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_smooth: '',L1)', test_loc
        print '()'


        test_loc = test_compute_gaussian_intensity(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_gaussian_intensity: '',L1)', test_loc
        print '()'


        test_loc = test_compute_random_phases(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_random_phases: '',L1)', test_loc
        print '()'


        test_loc = test_power_spectrum(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_power_spectrum: '',L1)', test_loc
        print '()'


        test_loc = test_compute_amplitudes(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_amplitudes: '',L1)', test_loc
        print '()'


        test_loc = test_compute_wave_numbers(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_wave_numbers: '',L1)', test_loc
        print '()'


        call test_add_gaussian_perturbation()


        print '(''test_validated: '',L1)', test_validated


        contains


        subroutine test_add_gaussian_perturbation()

          implicit none

          type(field)            :: field_used
          real(rkind), parameter :: amplitude=0.1d0

          ! input
          call field_used%ini()
          call field_used%write_data()

          ! output
          call add_gaussian_perturbation(
     $         field_used%x_map,
     $         field_used%y_map,
     $         field_used%nodes(:,:,1),
     $         amplitude)
          call field_used%write_data()

          print '(''check output files:'')'
          print '(''  - data0.nc: w/o perturbation'')'
          print '(''  - data1.nc: w/ perturbation'')'
          print '()'

        end subroutine test_add_gaussian_perturbation


        function test_compute_wave_numbers(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind)              , parameter :: x_min = 1.2d0
          real(rkind)              , parameter :: x_max = 5.1d0
          real(rkind), dimension(2), parameter :: test_kx = [1.61107315568707d0,3.2221463113742d0]

          real(rkind), dimension(2) :: kx


          !output
          call compute_wave_numbers(x_min,x_max,kx)

          !validation
          test_validated = is_real_vector_validated(
     $         kx,
     $         test_kx,
     $         detailled)

        end function test_compute_wave_numbers


        function test_compute_amplitudes(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(2), parameter :: kx = [1.3d0,2.7d0]
          real(rkind), dimension(2), parameter :: test_Ax = [0.88889160856286d0,0.5658279595543d0]

          real(rkind), dimension(2) :: Ax


          ! output
          call compute_amplitudes(kx,Ax)

          ! validation
          test_validated = is_real_vector_validated(
     $         Ax, test_Ax, detailled)
          

        end function test_compute_amplitudes


        function test_power_spectrum(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), parameter :: kx = 1.3d0
          real(rkind), parameter :: kx_max_spectrum = 5.0d0
          real(rkind), parameter :: test_power_spectrum_kx = 0.0027001414236d0

          real(rkind) :: power_spectrum_kx

          !output
          power_spectrum_kx = power_spectrum(kx,kx_max_spectrum)

          !validation
          test_validated = is_real_validated(
     $         power_spectrum_kx,
     $         test_power_spectrum_kx,
     $         detailled)

        end function test_power_spectrum


        function test_compute_random_phases(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), parameter    :: pi = ACOS(-1.0d0)
          integer                   :: seed_x
          real(rkind), dimension(2) :: Phix

          logical :: test_loc


          test_validated = .true.


          !input
          seed_x = 198738915

          !output
          call compute_random_phases(seed_x,Phix)

          !validation
          test_loc = (Phix(1).ge.0).and.(Phix(1).le.(2.0d0*pi))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''Phix(1) failed'')'
          end if

          test_loc = (Phix(2).ge.0).and.(Phix(2).le.(2.0d0*pi))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''Phix(2) failed'')'
          end if

          test_loc = .not.is_real_validated(Phix(1),Phix(2),.false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''Phix(1).ne.Phix(2) failed'')'
          end if

        end function test_compute_random_phases

        
        function test_compute_gaussian_intensity(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(2), parameter :: Ax = [1.0d0,2.0d0]
          real(rkind), dimension(2), parameter :: kx = [0.5d0,0.8d0]
          real(rkind)              , parameter :: x  = 1.2d0
          real(rkind), dimension(2), parameter :: Phix = [0.1d0, 3.5d0]
          real(rkind)              , parameter :: test_Px = 0.26540625963d0

          real(rkind) :: Px


          !output
          Px = compute_gaussian_intensity(Ax,kx,x,Phix)

          !validation
          test_validated = is_real_validated(
     $         Px,
     $         test_Px,
     $         detailled)

        end function test_compute_gaussian_intensity


        function test_smooth(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), parameter :: x_min = 1.0d0
          real(rkind), parameter :: x_max = 2.0d0
          real(rkind), parameter :: x = 1.2d0
          real(rkind) :: test_smooth_x = 0.64d0

          real(rkind) :: smooth_x


          !output
          smooth_x = smooth(x_min,x_max,x)

          !validation
          test_validated = is_real_validated(
     $         smooth_x,
     $         test_smooth_x,
     $         detailled)

        end function test_smooth


        subroutine check_inputs()

          implicit none

          if(.not.(
     $         (is_real_validated(x_min,-0.500000000d0,detailled)).and.
     $         (is_real_validated(x_max, 0.500000000d0,detailled)).and.
     $         (is_real_validated(y_min,-0.500000000d0,detailled)).and.
     $         (is_real_validated(y_max, 0.500000000d0,detailled)).and.
     $         (npx.eq.1).and.
     $         (npy.eq.1).and.
     $         (ntx.eq.105).and.
     $         (nty.eq.105).and.
     $         (is_real_validated(T0,0.99500000000d0,detailled)).and.
     $         (ic_choice.eq.homogeneous_liquid).and.
     $         (ic_perturbation_ac.eqv.(.false.)).and.
     $         (ne.eq.4))) then

             print '(''the test requires: '')'
             print '(''   - x_min = -0.5: '',L1)', is_real_validated(x_min,-0.500000000d0,.false.)
             print '(''   - x_max =  0.5: '',L1)', is_real_validated(x_max, 0.500000000d0,.false.)
             print '(''   - y_min = -0.5: '',L1)', is_real_validated(y_min,-0.500000000d0,.false.)
             print '(''   - y_max =  0.5: '',L1)', is_real_validated(y_max, 0.500000000d0,.false.)
             print '(''   - ntx   =  105: '',L1)', ntx.eq.105
             print '(''   - nty   =  105: '',L1)', nty.eq.105
             print '(''   - npx   =  1: '',L1)', npx.eq.1
             print '(''   - npy   =  1: '',L1)', npy.eq.1
             print '(''   - T0    =  0.995: '')', is_real_validated(T0,0.99500000000d0,.false.)
             print '(''   - ic_choice = homogeneous_liquid: '',L1)', ic_choice.eq.homogeneous_liquid
             print '(''   - ic_perturbation_ac=.false.: '', L1)', ic_perturbation_ac.eqv.(.false.)
             print '(''   - pm_model= dim2d: '',L1)', ne.eq.4
             print '()'
             stop ''

          end if

        end subroutine check_inputs

      end program test_gaussian_perturbation
