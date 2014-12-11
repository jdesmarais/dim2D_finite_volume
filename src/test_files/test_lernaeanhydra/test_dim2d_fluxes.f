      !test for 'dim2d_fluxes_module'
      program test_dim2d_fluxes

        use check_data_module, only :
     $       is_test_validated

        use dim2d_fluxes_module

        use dim2d_parameters, only :
     $       viscous_r,
     $       re,
     $       pr,
     $       we,
     $       cv_r

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use sd_operators_class, only :
     $       sd_operators

        implicit none

        
        real(rkind), dimension(nx,ny,ne) :: nodes
        real(rkind)                      :: dx
        real(rkind)                      :: dy
        type(sd_operators)               :: s

        logical                    :: detailled
        integer(ikind)             :: i,j
        real(rkind), dimension(26) :: test_data
        logical                    :: global, local
        
        detailled = .true.

        ! test the inputs
        call test_inputs()

        ! initialize data
        call initialize_data(
     $       nodes,dx,dy,i,j,
     $       test_data)

        ! dim2d_parameters
        if(detailled) then
           print '(''viscous_r: '', F16.6)', viscous_r
           print '(''re:        '', F16.6)', re
           print '(''pr:        '', F16.6)', pr
           print '(''we:        '', F16.6)', we
           print '(''cv_r:      '', F16.6)', cv_r
        end if

        detailled = .true.

        global = .true.

        local = is_test_validated(
     $          flux_x_mass_density(
     $          nodes,s,i,j),
     $          test_data(1),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_mass_density')

        local = is_test_validated(
     $          flux_y_mass_density(
     $          nodes,s,i,j),
     $          test_data(2),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_mass_density')

        local = is_test_validated(
     $          flux_x_inviscid_momentum_x(
     $          nodes,s,i,j),
     $          test_data(3),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_inv_momentum_x')

        local = is_test_validated(
     $          flux_x_viscid_momentum_x(
     $          nodes,s,i,j,dx,dy),
     $          test_data(4),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_vis_momentum_x')

        local = is_test_validated(
     $          flux_x_capillarity_momentum_x(
     $          nodes,s,i,j,dx,dy),
     $          test_data(5),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_cap_momentum_x')

        local = is_test_validated(
     $          flux_x_momentum_x(
     $          nodes,s,i,j,dx,dy),
     $          test_data(6),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_momentum_x')

        local = is_test_validated(
     $          flux_y_inviscid_momentum_x(
     $          nodes,s,i,j),
     $          test_data(7),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_inv_momentum_x')

        local = is_test_validated(
     $          flux_y_viscid_momentum_x(
     $          nodes,s,i,j,dx,dy),
     $          test_data(8),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_vis_momentum_x')

        local = is_test_validated(
     $          flux_y_capillarity_momentum_x(
     $          nodes,s,i,j,dx,dy),
     $          test_data(9),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_cap_momentum_x')

        local = is_test_validated(
     $          flux_y_momentum_x(
     $          nodes,s,i,j,dx,dy),
     $          test_data(10),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_momentum_x')

        local = is_test_validated(
     $          flux_x_inviscid_momentum_y(
     $          nodes,s,i,j),
     $          test_data(11),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_inv_momentum_y')

         local = is_test_validated(
     $          flux_x_viscid_momentum_y(
     $          nodes,s,i,j,dx,dy),
     $          test_data(12),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_vis_momentum_y')

         local = is_test_validated(
     $          flux_x_capillarity_momentum_y(
     $          nodes,s,i,j,dx,dy),
     $          test_data(13),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_cap_momentum_y')

        local = is_test_validated(
     $          flux_x_momentum_y(
     $          nodes,s,i,j,dx,dy),
     $          test_data(14),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_momentum_y')

        local = is_test_validated(
     $          flux_y_inviscid_momentum_y(
     $          nodes,s,i,j),
     $          test_data(15),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_inv_momentum_y')

        local = is_test_validated(
     $          flux_y_viscid_momentum_y(
     $          nodes,s,i,j,dx,dy),
     $          test_data(16),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_vis_momentum_y')

        local = is_test_validated(
     $          flux_y_capillarity_momentum_y(
     $          nodes,s,i,j,dx,dy),
     $          test_data(17),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_cap_momentum_y')

        local = is_test_validated(
     $          flux_y_momentum_y(
     $          nodes,s,i,j,dx,dy),
     $          test_data(18),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_momentum_y')

        local = is_test_validated(
     $          flux_x_inviscid_total_energy(
     $          nodes,s,i,j),
     $          test_data(19),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_inv_total_energy')

        local = is_test_validated(
     $          flux_x_viscid_total_energy(
     $          nodes,s,i,j,dx,dy),
     $          test_data(20),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_vis_total_energy')

        local = is_test_validated(
     $          flux_x_capillarity_total_energy(
     $          nodes,s,i,j,dx,dy),
     $          test_data(21),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_cap_total_energy')

        local = is_test_validated(
     $          flux_x_total_energy(
     $          nodes,s,i,j,dx,dy),
     $          test_data(22),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_x_total_energy')

        local = is_test_validated(
     $          flux_y_inviscid_total_energy(
     $          nodes,s,i,j),
     $          test_data(23),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_inv_total_energy')

        local = is_test_validated(
     $          flux_y_viscid_total_energy(
     $          nodes,s,i,j,dx,dy),
     $          test_data(24),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_vis_total_energy')

        local = is_test_validated(
     $          flux_y_capillarity_total_energy(
     $          nodes,s,i,j,dx,dy),
     $          test_data(25),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_cap_total_energy')

        local = is_test_validated(
     $          flux_y_total_energy(
     $          nodes,s,i,j,dx,dy),
     $          test_data(26),
     $          detailled)
        call print_screen(global,local,detailled,
     $       'flux_y_total_energy')

        contains

        subroutine print_screen(global,local,verbose,test_name)
          
          implicit none

          logical      , intent(inout) :: global
          logical      , intent(in)    :: local
          logical      , intent(in)    :: verbose
          character*(*), intent(in)    :: test_name

          if(verbose) then
             print '(A20,'': '', L1)', test_name, local
             print '()'
          end if

          global = global.and.local

        end subroutine print_screen


        subroutine initialize_data(
     $     nodes,dx,dy,i,j,
     $     test_data)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          real(rkind)                     , intent(out) :: dx
          real(rkind)                     , intent(out) :: dy
          integer                         , intent(out) :: i,j
          real(rkind), dimension(26)      , intent(out) :: test_data

          dx=0.5
          dy=0.6
          
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

          !<test the operators defined dim2d_prim
          i=3 !<index tested in the data along the x-axis
          j=3 !<index tested in the data along the y-axis
          
          !<test_data initialization
          test_data(1) = -6.4d0         !<mass flux_x
          test_data(2) = 21.2d0         !<mass flux_y

          test_data(3) = -1004.159722d0 !<momentum_x flux_x_inviscid   
          test_data(4) =  30.47379294d0 !<momentum_x flux_x_viscid
          test_data(5) = -2755.451565d0 !<momentum_x flux_x_capillarity
          test_data(6) = -734.7093239d0 !<momentum_x flux_x

          test_data(7) = -0.847495362d0 !<momentum_x flux_y_inviscid   
          test_data(8) = -0.348890097d0 !<momentum_x flux_y_viscid
          test_data(9) = -13.2d0        !<momentum_x flux_y_capillarity
          test_data(10)=  0.542282658d0 !<momentum_x flux_y

          test_data(11)= -16.47595386d0 !<momentum_y flux_x_inviscid   
          test_data(12)=  18.75419382d0 !<momentum_y flux_x_viscid
          test_data(13)= -66.0d0        !<momentum_y flux_x_capillarity
          test_data(14)= -13.62679262d0 !<momentum_y flux_x

          test_data(15)= -332.8217969d0 !<momentum_y flux_y_inviscid   
          test_data(16)= -0.894392659d0 !<momentum_y flux_y_viscid
          test_data(17)= -145.5328464d0 !<momentum_y flux_y_capillarity
          test_data(18)= -318.0896337d0 !<momentum_y flux_y

          test_data(19)=  465.1271451d0 !<total_energy flux_x_inviscid
          test_data(20)=  30.49323278d0 !<total_energy flux_x_viscid     
          test_data(21)= -7246.226655d0 !<total_energy flux_x_capillarity
          test_data(22)=  1183.651164d0 !<total_energy flux_x

          test_data(23)= -750.4429049d0 !<total_energy flux_y_inviscid   
          test_data(24)= -2.090662951d0 !<total_energy flux_y_viscid
          test_data(25)= -289.2335111d0 !<total_energy flux_y_capillarity
          test_data(26)= -721.1014212d0 !<total_energy flux_y

        end subroutine initialize_data


        subroutine test_inputs()
        
          implicit none

          logical :: test_parameters

          !<if nx<4, ny<4 then the test cannot be done
          if((nx.lt.4).or.(ny.lt.4).or.(ne.ne.4)) then
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

      end program test_dim2d_fluxes
