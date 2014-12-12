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
     $       re,pr,we,cv_r,
     $       gravity

        use parameters_bf_layer, only :
     $       interior_pt

        use parameters_constant, only :
     $       earth_gravity_choice

        use parameters_input, only :
     $       nx,ny,ne,
     $       gravity_choice

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

        !test_compute_flux_nopt
        test_loc = test_compute_flux_nopt(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_flux_nopt: '',L1)', test_loc
        print '()'

        !test_compute_flux_oneside
        test_loc = test_compute_flux_oneside(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_flux_oneside: '',L1)', test_loc
        print '()'

        !test_compute_flux_by_parts
        test_loc = test_compute_flux_by_parts(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_flux_by_parts: '',L1)', test_loc
        print '()'

        !test_compute_body_forces
        test_loc = test_compute_body_forces(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_body_forces: '',L1)', test_loc
        print '()'

        !test_get_viscous_coeff
        test_loc = test_get_viscous_coeff()
        test_validated = test_validated.and.test_loc
        print '(''test_get_viscous_coeff: '',L1)', test_loc
        print '()'

        !test_get_velocity
        test_loc = test_get_velocity(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_velocity: '',L1)', test_loc
        print '()'

        !test_are_open_bc_undermined
        test_loc = test_are_openbc_undermined()
        test_validated = test_validated.and.test_loc
        print '(''test_are_openbc_undermined: '',L1)', test_loc
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

          test_parameters=.true.

          test_parameters=test_parameters.and.(nx.eq.5)
          test_parameters=test_parameters.and.(ny.eq.5)
          test_parameters=test_parameters.and.(ne.eq.4)
          test_parameters=test_parameters.and.(gravity_choice.eq.earth_gravity_choice)

          test_parameters=test_parameters.and.is_test_validated(viscous_r,-1.5d0,.false.)
          test_parameters=test_parameters.and.is_test_validated(re,5.0d0,.false.)
          test_parameters=test_parameters.and.is_test_validated(pr,20.0d0,.false.)
          test_parameters=test_parameters.and.is_test_validated(we,10.0d0,.false.)
          test_parameters=test_parameters.and.is_test_validated(cv_r,2.5d0,.false.)
          test_parameters=test_parameters.and.is_test_validated(gravity,9.81d0,.false.)
          
          if(.not.test_parameters) then

             !< print the dim2d parameters used for the test
             print '(''WARNING: this test is designed for:'')'
             print '(''nx            : '', I2)', 5
             print '(''ny            : '', I2)', 5
             print '(''ne            : '', I2)', 4
             print '()'
             print '(''gravity_choice: '', I2)', earth_gravity_choice
             print '(''viscous_r     : '', F16.6)', -1.5
             print '(''re            : '', F16.6)', 5.
             print '(''pr            : '', F16.6)', 20.
             print '(''we            : '', F16.6)', 10.
             print '(''cv_r          : '', F16.6)', 2.5
             print '(''gravity       : '', F16.6)', 9.81
             print '()'
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

          real(rkind)                          :: dx
          real(rkind)                          :: dy
          real(rkind), dimension(nx,ny,ne)     :: nodes
          type(sd_operators)                   :: s
          type(pmodel_eq)                      :: pmodel_eq_tested
          integer(ikind)                       :: i,j
          real(rkind), dimension(nx+1,ny  ,ne) :: flux_x
          real(rkind), dimension(nx  ,ny+1,ne) :: flux_y
          real(rkind), dimension(8)            :: test_data
        

          !initialize the data for the test
          call initialize_data_for_test_flux(dx,dy,nodes,i,j,test_data)


          !print the dim2d parameters used for the test
          if(detailled) then
             call print_variables_for_test(
     $            nodes,dx,dy)
          end if

        
          !compute the flux_x and flux_y tables
          flux_x = pmodel_eq_tested%compute_flux_x(nodes,dx,dy,s)
          flux_y = pmodel_eq_tested%compute_flux_y(nodes,dx,dy,s)


          !test of the operators
          test_validated = test_flux_data(
     $         flux_x(i,j,:),
     $         flux_y(i,j,:),
     $         test_data,
     $         detailled)

        end function test_compute_flux


        !test: compute_flux_x_nopt and compute_flux_y_nopt
        function test_compute_flux_nopt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind)                          :: dx
          real(rkind)                          :: dy
          real(rkind), dimension(nx,ny,ne)     :: nodes
          type(sd_operators)                   :: s
          type(pmodel_eq)                      :: pmodel_eq_tested
          integer(ikind)                       :: i,j
          real(rkind), dimension(nx+1,ny  ,ne) :: flux_x
          real(rkind), dimension(nx  ,ny+1,ne) :: flux_y
          real(rkind), dimension(8)            :: test_data
        

          integer, dimension(2) :: x_borders
          integer, dimension(2) :: y_borders

          integer :: k,l
          integer, dimension(nx,ny) :: grdpts_id

          !initialize the data for the test
          call initialize_data_for_test_flux(dx,dy,nodes,i,j,test_data)
        
          !compute the flux_x and flux_y tables
          x_borders = [3,3]
          y_borders = [3,3]

          do l=1,5
             do k=1,5
                grdpts_id(k,l) = interior_pt
             end do
          end do

          call pmodel_eq_tested%compute_flux_x_nopt(
     $         nodes,dx,dy,s,
     $         grdpts_id,
     $         flux_x,
     $         x_borders,
     $         y_borders)

          call pmodel_eq_tested%compute_flux_y_nopt(
     $         nodes,dx,dy,s,
     $         grdpts_id,
     $         flux_y,
     $         x_borders,
     $         y_borders)

          !test of the operators
          test_validated = test_flux_data(
     $         flux_x(i,j,:),
     $         flux_y(i,j,:),
     $         test_data,
     $         detailled)

        end function test_compute_flux_nopt


        !test: compute_flux_x_oneside and compute_flux_y_oneside
        function test_compute_flux_oneside(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind)                          :: dx
          real(rkind)                          :: dy
          real(rkind), dimension(nx,ny,ne)     :: nodes
          type(sd_operators)                   :: s
          type(pmodel_eq)                      :: pmodel_eq_tested
          integer(ikind)                       :: i,j
          real(rkind), dimension(ne)           :: flux_x
          real(rkind), dimension(ne)           :: flux_y
          real(rkind), dimension(8)            :: test_data
        

          !initialize the data for the test
          call initialize_data_for_test_flux(dx,dy,nodes,i,j,test_data)

          !compute the flux_x and flux_y tables
          flux_x = pmodel_eq_tested%compute_flux_x_oneside(
     $         nodes,dx,dy,i,j,s)

          flux_y = pmodel_eq_tested%compute_flux_y_oneside(
     $         nodes,dx,dy,i,j,s)

          !test of the operators
          test_validated = test_flux_data(
     $         flux_x,
     $         flux_y,
     $         test_data,
     $         detailled)

        end function test_compute_flux_oneside


        !test: compute_flux_x_oneside and compute_flux_y_oneside
        function test_compute_flux_by_parts(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind)                      :: dx
          real(rkind)                      :: dy
          real(rkind), dimension(nx,ny,ne) :: nodes
          type(sd_operators)               :: s
          type(pmodel_eq)                  :: pmodel_eq_tested
          integer(ikind)                   :: i,j
          real(rkind), dimension(ne)       :: flux_x
          real(rkind), dimension(ne)       :: flux_y
          real(rkind), dimension(ne)       :: inviscid_flux_x
          real(rkind), dimension(ne)       :: inviscid_flux_y
          real(rkind), dimension(ne)       :: viscid_flux_x
          real(rkind), dimension(ne)       :: viscid_flux_y
          real(rkind), dimension(8)        :: test_data
          real(rkind), dimension(8)        :: inviscid_test_data
          real(rkind), dimension(8)        :: viscid_test_data
          logical                          :: test_loc

          !initialize the data for the test
          call initialize_data_for_test_flux(dx,dy,nodes,i,j,test_data)

          inviscid_test_data = [
     $         -6.4d0, -1004.159722d0, -16.47595386d0, 465.1271451d0,
     $         21.2d0, -0.847495362d0, -332.8217969d0,-750.4429049d0]

          viscid_test_data   = [
     $         0.0d0, 30.47379294d0, 18.75419382d0, 30.49323278d0,
     $         0.0d0,-0.348890097d0,-0.894392659d0,-2.090662951d0]


          !compute the flux_x and flux_y tables
          flux_x = pmodel_eq_tested%compute_flux_x_by_parts(
     $         nodes,dx,dy,i,j,s,
     $         inviscid_flux_x, viscid_flux_x)

          flux_y = pmodel_eq_tested%compute_flux_y_by_parts(
     $         nodes,dx,dy,i,j,s,
     $         inviscid_flux_y, viscid_flux_y)

          !test of the operators
          test_validated = .true.

          test_loc = test_flux_data(
     $         flux_x,
     $         flux_y,
     $         test_data,
     $         detailled)
          test_validated = test_loc

          test_loc = test_flux_data(
     $         inviscid_flux_x,
     $         inviscid_flux_y,
     $         inviscid_test_data,
     $         detailled)
          test_validated = test_loc

          test_loc = test_flux_data(
     $         viscid_flux_x,
     $         viscid_flux_y,
     $         viscid_test_data,
     $         detailled)
          test_validated = test_loc          

        end function test_compute_flux_by_parts


        subroutine initialize_data_for_test_flux(dx,dy,nodes,i,j,test_data)

          implicit none

          real(rkind)                     , intent(out) :: dx
          real(rkind)                     , intent(out) :: dy
          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
          integer(ikind)                  , intent(out) :: i
          integer(ikind)                  , intent(out) :: j
          real(rkind), dimension(8)       , intent(out) :: test_data

          
          !<initialize the tables for the field
          dx=0.5
          dy=0.6

          !<initialize the mass density
          nodes = reshape((/
     $         0.5d0,   0.2d0,  1.2d0,  5.0d0,  3.1d0,
     $         2.0d0,   4.2d0, 11.0d0, 10.6d0,  6.9d0,
     $        -14.2d0, 23.0d0,  9.8d0,  3.4d0,  2.3d0,
     $         2.45d0,  0.2d0,  9.0d0,  5.4d0, -1.0d0,
     $         6.3d0,  -5.1d0,  4.2d0,  9.8d0, -0.3d0,
     $         9.5d0,   9.8d0,  8.8d0,  5.0d0,  6.9d0,
     $         8.0d0,   5.8d0, -1.0d0, -0.6d0,  3.1d0,
     $         24.2d0,-13.0d0,  0.2d0,  6.6d0,  7.7d0,
     $         7.55d0,  9.8d0,  1.0d0,  4.6d0, 11.0d0,
     $         3.7d0,  15.1d0,  5.8d0,  0.2d0, 10.3d0,
     $        -8.5d0,   9.4d0, -6.4d0,  5.0d0, -0.7d0,
     $        -4.0d0,   2.6d0, 23.0d0, 21.8d0, 10.7d0,
     $        -52.6d0, 59.0d0, 19.4d0, 0.20d0, -3.1d0,
     $        -2.65d0,-9.40d0, 17.0d0, 6.20d0,-13.0d0,
     $         8.9d0, -25.3d0,  2.6d0, 19.4d0,-10.9d0,
     $        -1.5d0,  -1.8d0, -0.8d0,  3.0d0,  1.1d0,
     $         0.0d0,   2.2d0,  9.0d0,  8.6d0,  4.9d0,
     $        -16.2d0, 21.0d0,  7.8d0,  1.4d0,  0.3d0,
     $         0.45d0, -1.8d0,  7.0d0,  3.4d0, -3.0d0,
     $         4.3d0,  -7.1d0,  2.2d0,  7.8d0, -2.3d0/),
     $         (/5,5,ne/))
        
          !<test the operators defined dim2d_fluxes
          i=3 !<index tested in the data along the x-axis
          j=3 !<index tested in the data along the y-axis
          
          !<test_data initialization        
          test_data(1) =        -6.4d0 !<flux_x_mass
          test_data(2) =-734.7093239d0 !<flux_x_momentum_x
          test_data(3) =-13.62679262d0 !<flux_x_momentum_y
          test_data(4) = 1183.651164d0 !<flux_x_total_energy
          test_data(5) =        21.2d0 !<flux_y_mass
          test_data(6) = 0.542282658d0 !<flux_y_momentum_x
          test_data(7) =-318.0896337d0 !<flux_y_momentum_y
          test_data(8) =-721.1014212d0 !<flux_y_total_energy

        end subroutine initialize_data_for_test_flux

      
        function test_flux_data(flux_x, flux_y, test_data, detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(ne), intent(in) :: flux_x
          real(rkind), dimension(ne), intent(in) :: flux_y
          real(rkind), dimension(8) , intent(in) :: test_data
          logical                   , intent(in) :: detailled
          logical                                :: test_validated

          logical     :: test_loc
          real(rkind) :: prog_data

          test_validated = .true.

          prog_data  = flux_x(1)
          test_loc = is_test_validated(prog_data, test_data(1),.false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test flux_x_mass_density: '',F8.3,''->'',F8.3)',
     $            flux_x(1), test_data(1)
          end if
          
          prog_data  = flux_x(2)
          test_loc = is_test_validated(prog_data, test_data(2),.false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test flux_x_momentum_x: '',F8.3,''->'',F8.3)',
     $            flux_x(2), test_data(2)
          end if
          
          prog_data  = flux_x(3)
          test_loc = is_test_validated(prog_data, test_data(3),.false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test flux_x_momentum_y: '',F8.3,''->'',F8.3)',
     $            flux_x(3), test_data(3)
          end if
          
          prog_data  = flux_x(4)
          test_loc = is_test_validated(prog_data, test_data(4),.false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test flux_x_total_energy: '',F8.3,''->'',F8.3)',
     $            flux_x(4), test_data(4)
          end if
          
          prog_data  = flux_y(1)
          test_loc = is_test_validated(prog_data, test_data(5),.false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test flux_y_mass_density: '',F8.3,''->'',F8.3)',
     $            flux_y(1), test_data(5)
          end if
          
          prog_data  = flux_y(2)
          test_loc = is_test_validated(prog_data, test_data(6),.false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test flux_y_momentum_x: '',F8.3,''->'',F8.3)',
     $            flux_y(2), test_data(6)
          end if
          
          prog_data  = flux_y(3)
          test_loc = is_test_validated(prog_data, test_data(7),.false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test flux_y_momentum_y: '',F8.3,''->'',F8.3)',
     $            flux_y(3), test_data(7)
          end if

          prog_data  = flux_y(4)
          test_loc = is_test_validated(prog_data, test_data(8),.false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test flux_y_total_energy: '',F8.3,''->'',F8.3)',
     $            flux_y(4), test_data(8)
          end if

        end function test_flux_data


        subroutine print_variables_for_test(nodes,dx,dy)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy


          print '(''mass:       '', 1X, 5F8.3)', nodes(1:5,5,1)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,4,1)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,3,1)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,2,1)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,1,1)
          print '('''')'
          
          print '(''momentum_x: '', 1X, 5F8.3)', nodes(1:5,5,2)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,4,2)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,3,2)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,2,2)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,1,2)
          print '('''')'
          
          print '(''momentum_y: '', 1X, 5F8.3)', nodes(1:5,5,3)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,4,3)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,3,3)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,2,3)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,1,3)
          print '('''')'

          print '(''energy:     '', 1X, 5F8.3)', nodes(1:5,5,4)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,4,4)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,3,4)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,2,4)
          print '(''            '', 1X, 5F8.3)', nodes(1:5,1,4)
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


        !test: compute_body_forces
        function test_compute_body_forces(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(pmodel_eq)            :: p_model
          integer                    :: k
          real(rkind)                :: t,x,y
          real(rkind), dimension(ne) :: nodes
          real(rkind), dimension(ne) :: body_force
          real(rkind), dimension(ne) :: body_force_test
          
          nodes = [4.2d0,5.8d0,2.6d0,2.2d0]

          do k=1,ne
             body_force(k) =  p_model%compute_body_forces(t,x,y,nodes,k)
          end do

          body_force_test = [0.0d0, 0.0d0, -41.202d0, -25.506d0]

          test_validated = is_vector_validated(body_force, body_force_test,detailled)
          
        end function test_compute_body_forces


        !test: compute_body_forces
        function test_get_viscous_coeff()
     $       result(test_validated)

          implicit none

          logical             :: test_validated

          type(pmodel_eq) :: p_model
          real(rkind)     :: viscous_coeff
          real(rkind)     :: viscous_coeff_test
          
          viscous_coeff      = p_model%get_viscous_coeff()
          viscous_coeff_test = 0.2d0

          test_validated = is_test_validated(
     $         viscous_coeff,
     $         viscous_coeff_test,
     $         .false.)
          
        end function test_get_viscous_coeff


        !test: get_velocity
        function test_get_velocity(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(pmodel_eq)            :: p_model
          real(rkind), dimension(ne) :: nodes
          real(rkind), dimension(2)  :: velocity
          real(rkind), dimension(2)  :: velocity_test
          
          nodes = [4.2d0, 5.8d0, 2.6d0, 2.2d0]

          velocity      = p_model%get_velocity(nodes)
          velocity_test = [1.380952381d0, 0.619047619d0]

          test_validated = is_vector_validated(
     $         velocity,
     $         velocity_test,
     $         detailled)
          
        end function test_get_velocity


        !test: are_openbc_undermined
        function test_are_openbc_undermined()
     $       result(test_validated)

          implicit none

          logical :: test_validated

          type(pmodel_eq)            :: p_model
          real(rkind), dimension(ne) :: nodes
          logical                    :: undermined
          logical                    :: undermined_test
          
          nodes = [4.2d0, 5.8d0, 2.6d0, 2.2d0]

          undermined      = p_model%are_openbc_undermined(nodes)
          undermined_test = .true.

          test_validated = undermined.eqv.undermined_test
          
        end function test_are_openbc_undermined


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
