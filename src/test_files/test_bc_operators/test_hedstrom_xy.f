      program test_hedstrom_xy

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated

        use dim2d_parameters, only :
     $       cv_r

        use hedstrom_xy_module, only :
     $       compute_timedev_with_openbc

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right

        use parameters_constant, only :
     $       x_direction,
     $       obc_eigenqties_lin

        use parameters_input, only :
     $       ne,
     $       obc_eigenqties_strategy

        use parameters_kind, only :
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_R0

        implicit none

        logical :: test_validated
        logical :: test_loc
        logical :: detailled

        
        detailled = .true.
        test_validated = .true.


        test_loc       = test_compute_timedev_with_openbc(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_timedev_with_openbc: '',L1)', test_loc
        print '()'


        contains


        function test_compute_timedev_with_openbc(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind) :: t
          real(rkind) :: x
          real(rkind) :: y
          
          real(rkind), dimension(3,1,ne) :: nodes_L0
          real(rkind), dimension(3,1,ne) :: nodes_R0
          
          type(pmodel_eq) :: p_model
          
          real(rkind) :: dx
          
          real(rkind), dimension(ne) :: timedev_L0
          real(rkind), dimension(ne) :: timedev_R0
          
          call test_inputs()

        
          ! initialization of the nodes
          nodes_L0 = reshape((/
     $         1.46d0,  1.27d0,  1.47d0,
     $         0.146d0, 0.143d0, 0.145d0,
     $         0.01d0, 0.002d0,  0.05d0,
     $         4.89d0,  4.87d0,  4.86d0/),
     $         (/3,1,ne/))
          
          nodes_R0 = reshape((/
     $         1.47d0,  1.27d0,  1.46d0,
     $         -0.145d0,-0.143d0,-0.146d0,
     $         0.05d0, 0.002d0,  0.01d0,
     $         4.86d0,  4.87d0,  4.89d0/),
     $         (/3,1,ne/))
          
          t  = 0.0
          x  = 1.0
          y  = 1.5
          dx = 0.1
          
          timedev_L0 = compute_timedev_with_openbc(
     $         t,x,y,
     $         nodes_L0, 1,1,
     $         p_model,
     $         x_direction,
     $         gradient_x_x_oneside_L0, dx,
     $         incoming_left)
          
          timedev_R0 = compute_timedev_with_openbc(
     $         t,x,y,
     $         nodes_R0, 3,1,
     $         p_model,
     $         x_direction,
     $         gradient_x_x_oneside_R0, dx,
     $         incoming_right)

          test_validated = 
     $         is_real_vector_validated(
     $         timedev_L0,
     $         [-1.18764115727727d0,
     $           2.471895756d0,
     $           0.0d0,
     $          -2.70907495789801d0],
     $         detailled)

          test_validated = test_validated.and.
     $         is_real_vector_validated(
     $         timedev_R0,
     $         [-1.38610288419910d0,
     $          -2.97821061646847d0,
     $          -0.007727698696229735d0,
     $          -3.23478486279010d0],
     $         detailled)

       end function test_compute_timedev_with_openbc



       subroutine test_inputs()

         implicit none

         logical :: test_parameters

         test_parameters = ne.eq.4
         test_parameters = test_parameters.and.is_real_validated(cv_r,2.5d0,.false.)
         test_parameters = test_parameters.and.(obc_eigenqties_strategy.eq.obc_eigenqties_lin)

         if(.not.test_parameters) then
            print '(''ne=4              : '',L1)', ne.eq.4
            print '(''cv_r=2.5          : '',L1)', is_real_validated(cv_r,2.5d0,.false.)
            print '(''obc_eigenqties_lin: '',L1)', obc_eigenqties_strategy.eq.obc_eigenqties_lin
            print '()'
         end if

       end subroutine test_inputs

      end program test_hedstrom_xy
