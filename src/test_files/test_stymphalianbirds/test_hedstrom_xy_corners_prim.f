      program test_hedstrom_xy_corners_prim

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated

        use dim2d_parameters, only :
     $       cv_r

        use hedstrom_xy_corners_module, only :
     $       compute_timedev_with_openbc

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right

        use parameters_constant, only :
     $       n2_direction,
     $       obc_eigenqties_lin

        use parameters_input, only :
     $       ne,
     $       obc_eigenqties_strategy

        use parameters_kind, only :
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_fd_ncoords_module, only :
     $       gradient_n2_oneside_L0

        implicit none

        logical :: test_validated
        logical :: test_loc
        logical :: detailled

        
        detailled = .false.


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
          
          real(rkind), dimension(3,3,ne) :: nodes_L0
          
          type(pmodel_eq) :: p_model
          
          real(rkind) :: dn
          
          real(rkind), dimension(ne) :: timedev_L0
          
          call test_inputs()

        
          ! initialization of the nodes
          nodes_L0(1,1,:) = [ 1.46d0, 0.146d0, 0.010d0, 4.89d0]
          nodes_L0(2,2,:) = [ 1.27d0, 0.143d0, 0.002d0, 4.87d0]
          nodes_L0(3,3,:) = [ 1.47d0, 0.145d0, 0.050d0, 4.86d0]
          
          t  = 0.0
          x  = 1.0
          y  = 1.5
          dn = 0.1
          
          timedev_L0 = compute_timedev_with_openbc(
     $         t,x,y,
     $         nodes_L0, 1,1,
     $         p_model,
     $         n2_direction,
     $         gradient_n2_oneside_L0, dn,
     $         incoming_left)

          test_validated = 
     $         is_real_vector_validated(
     $         timedev_L0,
     $         [-1.151763244d0,
     $          -0.081441960d0,
     $           2.430955661d0,
     $          -2.700821882d0],
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

      end program test_hedstrom_xy_corners_prim
