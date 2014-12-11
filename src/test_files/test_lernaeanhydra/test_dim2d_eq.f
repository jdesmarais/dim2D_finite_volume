      !> @file
      !> test file for 'dim2d/pmodel_eq_class'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the procedures from pmodel_eq_class
      !
      !> @date
      ! 09_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_dim2d_eq

        use check_data_module, only :
     $       is_test_validated,
     $       is_vector_validated,
     $       is_matrix_validated

        use dim2d_parameters, only :
     $       viscous_r,
     $       re,pr,we,cv_r

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use sd_operators_fd_n_module, only :
     $       gradient_n1_xI_yI

        implicit none
        
        
        logical :: test_loc
        logical :: test_validated
        logical :: detailled        

        
        detailled = .true.


        !test the input parameters
        call test_inputs()

        !test_compute_flux
        test_loc = test_compute_flux(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_flux: '',L1)', test_loc
        print '()'

        !test_compute_x_eigenvalues
        test_loc = test_compute_x_eigenvalues(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_x_eigenvalues: '',L1)', test_loc
        print '()'

        !test_compute_y_eigenvalues
        test_loc = test_compute_y_eigenvalues(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_y_eigenvalues: '',L1)', test_loc
        print '()'

        !test_compute_x_lefteigenvector
        test_loc = test_compute_x_lefteigenvector(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_x_lefteigenvector: '',L1)', test_loc
        print '()'

        !test_compute_x_righteigenvector
        test_loc = test_compute_x_righteigenvector(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_x_righteigenvector: '',L1)', test_loc
        print '()'

        !test_compute_y_lefteigenvector
        test_loc = test_compute_y_lefteigenvector(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_y_lefteigenvector: '',L1)', test_loc
        print '()'

        !test_compute_y_righteigenvector
        test_loc = test_compute_y_righteigenvector(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_y_righteigenvector: '',L1)', test_loc
        print '()'

        !test_compute_x_transM
        test_loc = test_compute_x_transM(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_x_transM: '',L1)', test_loc
        print '()'

        !test_compute_y_transM
        test_loc = test_compute_y_transM(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_y_transM: '',L1)', test_loc
        print '()'

        !test_compute_x_leftConsLodiM
        test_loc = test_compute_x_leftConsLodiM(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_x_leftConsLodiM: '',L1)', test_loc
        print '()'

        !test_compute_y_leftConsLodiM
        test_loc = test_compute_y_leftConsLodiM(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_y_leftConsLodiM: '',L1)', test_loc
        print '()'
        
        !test_compute_n_gradient
        test_loc = test_compute_n_gradient(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_n_gradient: '',L1)', test_loc
        print '()'        


        contains

        !test the inputs
        subroutine test_inputs()
        
          implicit none

          logical :: test_parameters

          !<if nx<4, ny<4 then the test cannot be done
          if((nx.lt.5).or.(ny.lt.5).or.(ne.ne.4)) then
             stop 'nx and ny must be greater than 4 for the test'
          end if

          test_parameters=.true.
          test_parameters=test_parameters.and.(viscous_r.eq.-1.5d0)
          test_parameters=test_parameters.and.(re.eq.5d0)
          test_parameters=test_parameters.and.(pr.eq.20.0d0)
          test_parameters=test_parameters.and.(we.eq.10.0d0)
          test_parameters=test_parameters.and.(cv_r.eq.2.5d0)
          if(.not.test_parameters) then

             !< print the dim2d parameters used for the test
             print '(''WARNING: this test is designed for:'')'
             print '(''viscous_r: '', F16.6)', -1.5
             print '(''re:        '', F16.6)', 5.
             print '(''pr:        '', F16.6)', 20.
             print '(''we:        '', F16.6)', 10.
             print '(''cv_r:      '', F16.6)', 2.5
             print '(''it allows to see errors easily'')'
             print '('''')'

             stop 'dim2d_parameters not adapted for test'
          end if

        end subroutine test_inputs


        !test: compute_flux_x and compute_flux_y
        function test_compute_flux(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          logical :: test_loc

          real(rkind)                          :: dx
          real(rkind)                          :: dy
          real(rkind), dimension(nx,ny,ne)     :: nodes
          type(sd_operators)                   :: s
          type(pmodel_eq)                      :: pmodel_eq_tested
          integer(ikind)                       :: i,j
          real(rkind), dimension(nx+1,ny  ,ne) :: flux_x
          real(rkind), dimension(nx  ,ny+1,ne) :: flux_y
          real(rkind)                          :: prog_data
          real(rkind), dimension(8)            :: test_data
        

          test_validated = .true.

          !<initialize the tables for the field
          dx=0.5
          dy=0.6

          !<initialize the mass density
          nodes = reshape((/
     $          0.5d0,  0.2d0,  1.2d0,  5.0d0,  3.1d0,
     $          2.0d0,  4.2d0, 11.0d0, 10.6d0,  6.9d0,
     $        -14.2d0, 23.0d0,  9.8d0,  3.4d0,  2.3d0,
     $         2.45d0,  0.2d0,  9.0d0,  5.4d0, -1.0d0,
     $          6.3d0, -5.1d0,  4.2d0,  9.8d0, -0.3d0,
     $          9.5d0,  9.8d0,  8.8d0,  5.0d0,  6.9d0,
     $          8.0d0,  5.8d0, -1.0d0, -0.6d0,  3.1d0,
     $         24.2d0,-13.0d0,  0.2d0,  6.6d0,  7.7d0,
     $         7.55d0,  9.8d0,  1.0d0,  4.6d0, 11.0d0,
     $          3.7d0,  15.1d0, 5.8d0,  0.2d0, 10.3d0,
     $         -8.5d0,  9.4d0, -6.4d0,  5.0d0, -0.7d0,
     $         -4.0d0,  2.6d0, 23.0d0, 21.8d0, 10.7d0,
     $        -52.6d0, 59.0d0, 19.4d0, 0.20d0, -3.1d0,
     $        -2.65d0,-9.40d0, 17.0d0, 6.20d0,-13.0d0,
     $         8.9d0, -25.3d0,  2.6d0, 19.4d0,-10.9d0,
     $         -1.5d0, -1.8d0, -0.8d0,  3.0d0,  1.1d0,
     $          0.0d0,  2.2d0,  9.0d0,  8.6d0,  4.9d0,
     $        -16.2d0, 21.0d0,  7.8d0,  1.4d0,  0.3d0,
     $         0.45d0, -1.8d0,  7.0d0,  3.4d0, -3.0d0,
     $          4.3d0, -7.1d0,  2.2d0,  7.8d0, -2.3d0/),
     $         (/5,5,ne/))
        
          !<test the operators defined dim2d_fluxes
          i=3 !<index tested in the data along the x-axis
          j=3 !<index tested in the data along the y-axis
          
          !<test_data initialization        
          test_data(1) = -10.033333d0 !<flux_x_mass
          test_data(2) = -760.92652d0 !<flux_x_momentum_x
          test_data(3) =  -24.31047d0 !<flux_x_momentum_y
          test_data(4) = 1463.64782d0 !<flux_x_total_energy
          test_data(5) = 23.8500000d0 !<flux_y_mass
          test_data(6) =    6.16168d0 !<flux_y_momentum_x
          test_data(7) = -352.65586d0 !<flux_y_momentum_y
          test_data(8) = -840.36959d0!<flux_y_total_energy

          !< print the dim2d parameters used for the test
          if(detailled) then
             call print_variables_for_test(
     $            nodes,dx,dy)
          end if

        
          !<compute the flux_x and flux_y tables
          flux_x = pmodel_eq_tested%compute_flux_x(nodes,dx,dy,s)
          flux_y = pmodel_eq_tested%compute_flux_y(nodes,dx,dy,s)


          !<test of the operators
          if(detailled) then
             
             prog_data  = flux_x(i,j,1)
             test_loc = is_test_validated(prog_data, test_data(1),detailled)
             test_validated = test_validated.and.test_loc
             print '(''test flux_x_mass_density: '',1L1)', test_validated
             
             prog_data  = flux_x(i,j,2)
             test_loc = is_test_validated(prog_data, test_data(2),detailled)
             test_validated = test_validated.and.test_loc
             print '(''test flux_x_momentum_x: '',1L1)', test_validated
             
             prog_data  = flux_x(i,j,3)
             test_loc = is_test_validated(prog_data, test_data(3),detailled)
             test_validated = test_validated.and.test_loc
             print '(''test flux_x_momentum_y: '',1L1)', test_validated
             
             prog_data  = flux_x(i,j,4)
             test_loc = is_test_validated(prog_data, test_data(4),detailled)
             test_validated = test_validated.and.test_loc
             print '(''test flux_x_total_energy: '',1L1)', test_validated
             
             prog_data  = flux_y(i,j,1)
             test_loc = is_test_validated(prog_data, test_data(5),detailled)
             test_validated = test_validated.and.test_loc
             print '(''test flux_y_mass_density: '',1L1)', test_validated
             
             prog_data  = flux_y(i,j,2)
             test_loc = is_test_validated(prog_data, test_data(6),detailled)
             test_validated = test_validated.and.test_loc
             print '(''test flux_y_momentum_x: '',1L1)', test_validated
             
             prog_data  = flux_y(i,j,3)
             test_loc = is_test_validated(prog_data, test_data(7),detailled)
             test_validated = test_validated.and.test_loc
             print '(''test flux_y_momentum_y: '',1L1)', test_validated

             prog_data  = flux_y(i,j,4)
             test_loc = is_test_validated(prog_data, test_data(8),detailled)
             test_validated = test_validated.and.test_loc
             print '(''test flux_y_total_energy: '',1L1)', test_validated

          else

             test_loc = is_test_validated(flux_x(i,j,1), test_data(1),detailled)
             test_validated = test_validated.and.test_loc

             test_loc = is_test_validated(flux_x(i,j,2), test_data(2),detailled)
             test_validated = test_validated.and.test_loc

             test_loc = is_test_validated(flux_x(i,j,3), test_data(3),detailled)
             test_validated = test_validated.and.test_loc

             test_loc = is_test_validated(flux_x(i,j,4), test_data(4),detailled)
             test_validated = test_validated.and.test_loc

             test_loc = is_test_validated(flux_y(i,j,1), test_data(5),detailled)
             test_validated = test_validated.and.test_loc

             test_loc = is_test_validated(flux_y(i,j,2), test_data(6),detailled)
             test_validated = test_validated.and.test_loc

             test_loc = is_test_validated(flux_y(i,j,3), test_data(7),detailled)
             test_validated = test_validated.and.test_loc

             test_loc = is_test_validated(flux_y(i,j,4), test_data(8),detailled)
             test_validated = test_validated.and.test_loc

          end if

        end function test_compute_flux


        subroutine print_variables_for_test(nodes,dx,dy)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy


          print '(''mass:       '', 1X, 4F8.3)', nodes(1:5,5,1)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,4,1)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,3,1)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,2,1)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,1,1)
          print '('''')'
          
          print '(''momentum_x: '', 1X, 4F8.3)', nodes(1:5,5,2)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,4,2)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,3,2)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,2,2)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,1,2)
          print '('''')'
          
          print '(''momentum_y: '', 1X, 4F8.3)', nodes(1:5,5,3)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,4,3)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,3,3)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,2,3)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,1,3)
          print '('''')'

          print '(''energy:     '', 1X, 4F8.3)', nodes(1:5,5,4)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,4,4)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,3,4)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,2,4)
          print '(''            '', 1X, 4F8.3)', nodes(1:5,1,4)
          print '('''')'

          print '(''dx: '', 1X, F8.3)', dx
          print '(''dy: '', 1X, F8.3)', dy
          print '('''')'

          print '(''Re:        '', 1X, F12.3)', re
          print '(''We:        '', 1X, F12.3)', we
          print '(''Pr:        '', 1X, F12.3)', pr
          print '(''viscous_r: '', 1X, F12.3)', viscous_r
          print '(''cv_r:      '', 1X, F12.3)', cv_r
          print '('''')'
        
        end subroutine print_variables_for_test


        !test: compute_x_eigenvalues
        function test_compute_x_eigenvalues(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind), dimension(ne) :: nodes
          real(rkind), dimension(ne) :: eigenvalues
          real(rkind), dimension(ne) :: eigenvalues_test

          type(pmodel_eq) :: p_model

          nodes = [4.2d0,5.8d0,2.6d0,2.2d0]

          eigenvalues_test =
     $         [1.380952381d0,1.380952381d0,-2.708717144d0,5.470621906d0]

          eigenvalues = p_model%compute_x_eigenvalues(nodes)

          test_validated = is_vector_validated(
     $         eigenvalues,
     $         eigenvalues_test,
     $         detailled)

        end function test_compute_x_eigenvalues


        !test: compute_y_eigenvalues
        function test_compute_y_eigenvalues(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind), dimension(ne) :: nodes
          real(rkind), dimension(ne) :: eigenvalues
          real(rkind), dimension(ne) :: eigenvalues_test

          type(pmodel_eq) :: p_model

          nodes = [4.2d0,5.8d0,2.6d0,2.2d0]

          eigenvalues_test =
     $         [0.619047619d0,0.619047619d0,-3.470621906d0,4.708717144d0]

          eigenvalues = p_model%compute_y_eigenvalues(nodes)

          test_validated = is_vector_validated(
     $         eigenvalues,
     $         eigenvalues_test,
     $         detailled)

        end function test_compute_y_eigenvalues


        !test: compute_x_lefteigenvector
        function test_compute_x_lefteigenvector(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: eigenvect
          real(rkind), dimension(ne,ne) :: eigenvect_test

          type(pmodel_eq) :: p_model

          nodes = [4.2d0,5.8d0,2.6d0,2.2d0]

          eigenvect_test = reshape((/
     $         -0.14739229d0,	0.0d0	      ,  0.238095238d0,	 0.0d0        ,
     $          1.57515693d0,	-0.082566195d0,	-0.037012432d0,	 0.059789314d0,
     $        -1.986044512d0,	-1.354358572d0,	  0.30952381d0,	-0.5d0        ,
     $         -7.63368338d0,	 2.735310953d0,	  0.30952381d0,	-0.5d0        /),
     $         (/4,4/))

          eigenvect = p_model%compute_x_lefteigenvector(nodes)

          test_validated = is_matrix_validated(
     $         eigenvect,
     $         eigenvect_test,
     $         detailled)

        end function test_compute_x_lefteigenvector


        !test: compute_x_righteigenvector
        function test_compute_x_righteigenvector(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: eigenvect
          real(rkind), dimension(ne,ne) :: eigenvect_test

          type(pmodel_eq) :: p_model

          nodes = [4.2d0,5.8d0,2.6d0,2.2d0]

          eigenvect_test = reshape((/
     $         0.0d0,	1.0d0        ,   0.059789314d0,	 0.059789314d0,
     $         0.0d0,   1.380952381d0,	-0.161952339d0,	  0.32708473d0,
     $         4.2d0,   0.619047619d0,   0.037012432d0,	 0.037012432d0,
     $         2.6d0,  -7.329478458d0,	-1.775892941d0,	-1.100556035d0/),
     $         (/4,4/))

          eigenvect = p_model%compute_x_righteigenvector(nodes)

          test_validated = is_matrix_validated(
     $         eigenvect,
     $         eigenvect_test,
     $         detailled)

        end function test_compute_x_righteigenvector


        !test: compute_y_lefteigenvector
        function test_compute_y_lefteigenvector(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: eigenvect
          real(rkind), dimension(ne,ne) :: eigenvect_test

          type(pmodel_eq) :: p_model

          nodes = [4.2d0,5.8d0,2.6d0,2.2d0]

          eigenvect_test = reshape((/
     $         -0.328798186d0, 	0.238095238d0,  0.0d0        ,	0.0d0,
     $           1.57515693d0, -0.082566195d0, -0.037012432d0,	0.059789314d0,
     $         -3.544013854d0,   0.69047619d0, -1.735310953d0, -0.5d0,
     $         -6.075714037d0,   0.69047619d0, 	2.354358572d0, -0.5d0/),
     $         (/4,4/))

          eigenvect = p_model%compute_y_lefteigenvector(nodes)

          test_validated = is_matrix_validated(
     $         eigenvect,
     $         eigenvect_test,
     $         detailled)

        end function test_compute_y_lefteigenvector


        !test: compute_y_righteigenvector
        function test_compute_y_righteigenvector(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: eigenvect
          real(rkind), dimension(ne,ne) :: eigenvect_test

          type(pmodel_eq) :: p_model

          nodes = [4.2d0,5.8d0,2.6d0,2.2d0]

          eigenvect_test = reshape((/
     $         0.0d0,  1.0d0	    ,  0.059789314d0,  0.059789314d0,
     $         4.2d0,  1.380952381d0,  0.082566195d0,  0.082566195d0,
     $         0.0d0,  0.619047619d0, -0.207506102d0,  0.281530967d0,
     $         5.8d0, -7.329478458d0, -1.589593105d0, -1.286855871d0/),
     $         (/4,4/))

          eigenvect = p_model%compute_y_righteigenvector(nodes)

          test_validated = is_matrix_validated(
     $         eigenvect,
     $         eigenvect_test,
     $         detailled)

        end function test_compute_y_righteigenvector


        !test: compute_x_transM
        function test_compute_x_transM(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: transM
          real(rkind), dimension(ne,ne) :: transM_test

          type(pmodel_eq) :: p_model

          nodes = [4.2d0,5.8d0,2.6d0,2.2d0]

          transM_test = reshape((/
     $          0.0d0        ,	0.0d0        ,  1.0d0        , 0.0d0,
     $         -0.854875283d0,  0.619047619d0,  1.380952381d0, 0.0d0,
     $         -10.00294785d0,	1.380952381d0,  1.857142857d0,-1.0d0,
     $          8.936043624d0,	0.854875283d0, -23.67165533d0, 0.0d0/),
     $         (/4,4/))

          transM = p_model%compute_x_transM(nodes)

          test_validated = is_matrix_validated(
     $         transM,
     $         transM_test,
     $         detailled)

        end function test_compute_x_transM


        !test: compute_y_transM
        function test_compute_y_transM(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: transM
          real(rkind), dimension(ne,ne) :: transM_test

          type(pmodel_eq) :: p_model

          nodes = [4.2d0,5.8d0,2.6d0,2.2d0]

          transM_test = reshape((/
     $         0.0d0         ,	1.0d0        ,	0.0d0        , 0.0d0,
     $         -11.52675737d0,	4.142857143d0,	0.619047619d0,-1.0d0,
     $         -0.854875283d0,	0.619047619d0,	1.380952381d0, 0.0d0,
     $          19.93425116d0,	-22.1478458d0,	0.854875283d0, 0.0d0/),
     $         (/4,4/))

          transM = p_model%compute_y_transM(nodes)

          test_validated = is_matrix_validated(
     $         transM,
     $         transM_test,
     $         detailled)

        end function test_compute_y_transM   


        !test: compute_x_leftConsLodiM
        function test_compute_x_leftConsLodiM(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: leftConsLodiM
          real(rkind), dimension(ne,ne) :: leftConsLodiM_test

          type(pmodel_eq) :: p_model

          nodes = [4.2d0,5.8d0,2.6d0,2.2d0]

          leftConsLodiM_test = reshape((/
     $          -0.14739229d0, 0.0d0	     ,  0.238095238d0,  0.0d0,
     $          26.34512472d0, -1.380952381d0, -0.619047619d0,  1.0d0,
     $         -3.972089023d0, -2.708717144d0,	0.619047619d0, -1.0d0,
     $         -15.26736676d0, 	5.470621906d0,	0.619047619d0, -1.0d0/),
     $         (/4,4/))

          leftConsLodiM = p_model%compute_x_leftConslodiM(nodes)

          test_validated = is_matrix_validated(
     $         leftConsLodiM,
     $         leftConsLodiM_test,
     $         detailled)

        end function test_compute_x_leftConsLodiM


        !test: compute_x_leftConsLodiM
        function test_compute_y_leftConsLodiM(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: leftConsLodiM
          real(rkind), dimension(ne,ne) :: leftConsLodiM_test

          type(pmodel_eq) :: p_model

          nodes = [4.2d0,5.8d0,2.6d0,2.2d0]

          leftConsLodiM_test = reshape((/
     $         -0.328798186d0, 	0.238095238d0,  0.0d0	     ,  0.0d0,
     $          26.34512472d0, -1.380952381d0, -0.619047619d0,  1.0d0,
     $         -7.088027709d0, 	1.380952381d0, -3.470621906d0, -1.0d0,
     $         -12.15142807d0, 	1.380952381d0, 	4.708717144d0, -1.0d0/),
     $         (/4,4/))

          leftConsLodiM = p_model%compute_y_leftConslodiM(nodes)

          test_validated = is_matrix_validated(
     $         leftConsLodiM,
     $         leftConsLodiM_test,
     $         detailled)

        end function test_compute_y_leftConsLodiM


        !test: compute_n_gradient
        function test_compute_n_gradient(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind)                      :: dx
          real(rkind)                      :: dy
          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(ne)       :: gradient
          real(rkind), dimension(ne)       :: gradient_test

          type(pmodel_eq) :: p_model

          nodes = reshape((/
     $          0.5d0,  0.2d0,  1.2d0,  5.0d0,  3.1d0,
     $          2.0d0,  4.2d0, 11.0d0, 10.6d0,  6.9d0,
     $        -14.2d0, 23.0d0,  9.8d0,  3.4d0,  2.3d0,
     $         2.45d0,  0.2d0,  9.0d0,  5.4d0, -1.0d0,
     $          6.3d0, -5.1d0,  4.2d0,  9.8d0, -0.3d0,
     $          9.5d0,  9.8d0,  8.8d0,  5.0d0,  6.9d0,
     $          8.0d0,  5.8d0, -1.0d0, -0.6d0,  3.1d0,
     $         24.2d0,-13.0d0,  0.2d0,  6.6d0,  7.7d0,
     $         7.55d0,  9.8d0,  1.0d0,  4.6d0, 11.0d0,
     $          3.7d0,  15.1d0, 5.8d0,  0.2d0, 10.3d0,
     $         -8.5d0,  9.4d0, -6.4d0,  5.0d0, -0.7d0,
     $         -4.0d0,  2.6d0, 23.0d0, 21.8d0, 10.7d0,
     $        -52.6d0, 59.0d0, 19.4d0, 0.20d0, -3.1d0,
     $        -2.65d0,-9.40d0, 17.0d0, 6.20d0,-13.0d0,
     $         8.9d0, -25.3d0,  2.6d0, 19.4d0,-10.9d0,
     $         -1.5d0, -1.8d0, -0.8d0,  3.0d0,  1.1d0,
     $          0.0d0,  2.2d0,  9.0d0,  8.6d0,  4.9d0,
     $        -16.2d0, 21.0d0,  7.8d0,  1.4d0,  0.3d0,
     $         0.45d0, -1.8d0,  7.0d0,  3.4d0, -3.0d0,
     $          4.3d0, -7.1d0,  2.2d0,  7.8d0, -2.3d0/),
     $         (/5,5,ne/))

          dx = 0.1d0
          dy = 0.2d0

          gradient_test = [-65.76093065d0,65.76093065d0,-197.282792d0,-65.76093065d0]
          
          gradient = p_model%compute_n_gradient(nodes,3,3,gradient_n1_xI_yI,dx,dy)

          test_validated = is_vector_validated(
     $         gradient,
     $         gradient_test,
     $         detailled)

        end function test_compute_n_gradient

      end program test_dim2d_eq
