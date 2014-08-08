      program test_ns2d_eq_program

        use ns2d_parameters, only :
     $     mach_infty,
     $     gamma

        use ns2d_prim_module, only :
     $       mass_density,
     $       momentum_x,
     $       momentum_y,
     $       total_energy,
     $       velocity_x,
     $       velocity_y,
     $       pressure,
     $       temperature,
     $       qx_inviscid_x_flux,
     $       qy_inviscid_y_flux,
     $       qxy_transport,
     $       energy_inviscid_x_flux,
     $       energy_inviscid_y_flux

        use ns2d_fluxes_module, only :
     $       flux_x_mass_density, flux_y_mass_density,
     $       flux_x_momentum_x  , flux_y_momentum_x,
     $       flux_x_momentum_y  , flux_y_momentum_y,
     $       flux_x_total_energy, flux_y_total_energy

        use parameters_kind, only : ikind, rkind

        use sd_operators_class, only : sd_operators


        implicit none

        real(rkind), dimension(5,5,4) :: nodes
        real(rkind)                   :: dx
        real(rkind)                   :: dy
        type(sd_operators)            :: s
        logical                       :: detailled
        logical                       :: test_validated

        detailled = .false.
        if(
     $       (.not.is_test_validated(gamma,5.0d0/3.0d0,detailled)).or.
     $       (.not.is_test_validated(mach_infty,1.0d0,detailled))) then

           stop 'the test requires (gamma,mach)=(5/3,1)'

        end if


        !initialize the nodes
        call initialize_nodes(nodes,dx,dy)
        

        !compute the quantities for the NS equations
        !and compare with the test data
        detailled = .false.

        print '(''test ns2_prim'')'
        print '(''---------------------------'')'
        test_validated = test_ns2d_prim(nodes,detailled)
        if(detailled) print '(''---------------------------'')'
        print '(''test_validated: '',L3)', test_validated
        print '(''---------------------------'')'
        print '()'


        !compute the fluxes for the NS equations
        !and compare with the test data
        detailled = .false.

        print '(''test ns2_fluxes'')'
        print '(''---------------------------'')'
        test_validated = test_ns2d_fluxes(nodes,dx,dy,s,detailled)
        if(detailled) print '(''---------------------------'')'
        print '(''test_validated: '',L3)', test_validated
        print '(''---------------------------'')'
        print '()'


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
     $         int(var*10000.)-
     $         sign(int(abs(cst*10000.)),int(cst*10000.))).le.1
          
        end function is_test_validated


        subroutine initialize_nodes(nodes,dx,dy)

          implicit none

          real(rkind), dimension(:,:,:), intent(out) :: nodes
          real(rkind)                  , intent(out) :: dx
          real(rkind)                  , intent(out) :: dy

          integer(ikind) :: j


          !space steps
          dx = 0.5
          dy = 0.6


          !mass
          nodes(1,1,1) = 0.5
          nodes(2,1,1) = 0.2
          nodes(3,1,1) = 1.2
          nodes(4,1,1) = 5.0
          nodes(5,1,1) = 0.6

          nodes(1,2,1) = 3.0
          nodes(2,2,1) = 4.2
          nodes(3,2,1) = 11.0
          nodes(4,2,1) = 10.6
          nodes(5,2,1) = 5.2

          nodes(1,3,1) = -14.2
          nodes(2,3,1) = 23
          nodes(3,3,1) = 9.8
          nodes(4,3,1) = 3.4
          nodes(5,3,1) = 9.1

          nodes(1,4,1) = 2.45
          nodes(2,4,1) = 0.2
          nodes(3,4,1) = 9.0
          nodes(4,4,1) = 5.4
          nodes(5,4,1) =-2.3

          nodes(1,5,1) = 3.6
          nodes(2,5,1) = 0.1
          nodes(3,5,1) = 6.3
          nodes(4,5,1) = 8.9
          nodes(5,5,1) = -4.23


          !momentum-x
          nodes(1,1,2) = 7.012
          nodes(2,1,2) =-6.323
          nodes(3,1,2) = 3.012
          nodes(4,1,2) = 4.5
          nodes(5,1,2) = 9.6
                    
          nodes(1,2,2) = 4.26
          nodes(2,2,2) = 4.23
          nodes(3,2,2) = 4.5
          nodes(4,2,2) = 7.56
          nodes(5,2,2) = 7.21
                    
          nodes(1,3,2) = 0.23
          nodes(2,3,2) = 7.23
          nodes(3,3,2) = 3.1
          nodes(4,3,2) = 8.9
          nodes(5,3,2) = 9.3
                    
          nodes(1,4,2) = 8.23
          nodes(2,4,2) = -3.1
          nodes(3,4,2) = 6.03
          nodes(4,4,2) = 6.25
          nodes(5,4,2) = 5.12
                    
          nodes(1,5,2) = 3.2
          nodes(2,5,2) = 8.12
          nodes(3,5,2) = 8.9
          nodes(4,5,2) = 4.2
          nodes(5,5,2) = 7.8


          !momentum-y
          nodes(1,1,3) = 7.1
          nodes(2,1,3) = 1.052
          nodes(3,1,3) = 1.23
          nodes(4,1,3) = 7.89
          nodes(5,1,3) = 8.0
                    
          nodes(1,2,3) = 8.362
          nodes(2,2,3) = 4.56
          nodes(3,2,3) = 9.6
          nodes(4,2,3) = 8.96
          nodes(5,2,3) = -3.23
                    
          nodes(1,3,3) = 2.53
          nodes(2,3,3) = -3.23
          nodes(3,3,3) = 7.25
          nodes(4,3,3) = 1.02
          nodes(5,3,3) = 9.26
                    
          nodes(1,4,3) = 8.965
          nodes(2,4,3) = 4.789
          nodes(3,4,3) = 4.56
          nodes(4,4,3) = 3.012
          nodes(5,4,3) = -1.45
                    
          nodes(1,5,3) = 6.26
          nodes(2,5,3) = 5.201
          nodes(3,5,3) = 2.03
          nodes(4,5,3) = 7.89
          nodes(5,5,3) = 9.889


          !total energy
          nodes(1,1,4) = 6.23
          nodes(2,1,4) = 4.12
          nodes(3,1,4) = -3.6
          nodes(4,1,4) = -6.52
          nodes(5,1,4) = 9.57
                    
          nodes(1,2,4) = -0.12
          nodes(2,2,4) = 8.2
          nodes(3,2,4) = 1.2
          nodes(4,2,4) = 7.89
          nodes(5,2,4) = 5.62
                    
          nodes(1,3,4) = -6.23
          nodes(2,3,4) = 6.201
          nodes(3,3,4) = 6.7
          nodes(4,3,4) = 4.12
          nodes(5,3,4) = 1.29
                    
          nodes(1,4,4) = 1.2
          nodes(2,4,4) = 7.958
          nodes(3,4,4) = 1
          nodes(4,4,4) = -5.62
          nodes(5,4,4) = 0.36
                    
          nodes(1,5,4) = 9.6
          nodes(2,5,4) = 6.12
          nodes(3,5,4) = 8.9
          nodes(4,5,4) = 8.95
          nodes(5,5,4) = 6.3

          print '()'
          print '(''mass_density'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,1)
          end do
          print '()'

          print '()'
          print '(''momentum-x'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,2)
          end do
          print '()'

          print '()'
          print '(''momentum-y'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,3)
          end do
          print '()'

          print '()'
          print '(''total energy'')'
          do j=1,5
             print '(5F8.3)', nodes(1:5,6-j,4)
          end do
          print '()'

        end subroutine initialize_nodes


        function test_ns2d_prim(nodes,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated


          real(rkind), dimension(13) :: test_data
          logical                    :: loc
          integer(ikind)             :: i,j


          test_data(1)  =  9.8        !mass
          test_data(2)  =  3.1        !momentum-x
          test_data(3)  =  7.25       !momentum-y
          test_data(4)  =  6.7        !total energy
                        
          test_data(5)  =  0.31632653 !velocity-x
          test_data(6)  =  0.73979592 !velocity-y
          test_data(7)  =  2.35195578 !pressure
          test_data(8)  =  0.39999248 !temperature
          
          test_data(9)  =  3.33256802 !qx_inviscid_x_flux
          test_data(10) =  7.71547619 !qy_inviscid_y_flux
          test_data(11) =  2.29336734 !qxy_transport
          test_data(12) =  2.86337376 !energy_inviscid_x_flux
          test_data(13) =  6.69659994 !energy_inviscid_y_flux


          i = 3
          j = 3


          test_validated = .true.
          loc = is_test_validated(
     $         mass_density(nodes,i,j),
     $         test_data(1),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test mass: '', L3)', loc


          loc = is_test_validated(
     $         momentum_x(nodes,i,j),
     $         test_data(2),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test momentum-x: '', L3)', loc


          loc = is_test_validated(
     $         momentum_y(nodes,i,j),
     $         test_data(3),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test momentum-y: '', L3)', loc


          loc = is_test_validated(
     $         total_energy(nodes,i,j),
     $         test_data(4),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test total_energy: '', L3)', loc


          loc = is_test_validated(
     $         velocity_x(nodes,i,j),
     $         test_data(5),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test velocity_x: '', L3)', loc


          loc = is_test_validated(
     $         velocity_y(nodes,i,j),
     $         test_data(6),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test velocity_y: '', L3)', loc


          loc = is_test_validated(
     $         pressure(nodes,i,j),
     $         test_data(7),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test pressure: '', L3)', loc


          loc = is_test_validated(
     $         temperature(nodes,i,j),
     $         test_data(8),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test temperature: '', L3)', loc


          loc = is_test_validated(
     $         qx_inviscid_x_flux(nodes,i,j),
     $         test_data(9),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test qx_inviscid_x_flux: '', L3)', loc


          loc = is_test_validated(
     $         qy_inviscid_y_flux(nodes,i,j),
     $         test_data(10),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test qy_inviscid_y_flux: '', L3)', loc


          loc = is_test_validated(
     $         qxy_transport(nodes,i,j),
     $         test_data(11),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test qxy_transport: '', L3)', loc


          loc = is_test_validated(
     $         energy_inviscid_x_flux(nodes,i,j),
     $         test_data(12),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test energy_inviscid_x_flux: '', L3)', loc


          loc = is_test_validated(
     $         energy_inviscid_y_flux(nodes,i,j),
     $         test_data(13),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test energy_inviscid_y_flux: '', L3)', loc


        end function test_ns2d_prim


        function test_ns2d_fluxes(nodes,dx,dy,s,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          type(sd_operators)           , intent(in) :: s
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated


          real(rkind), dimension(8) :: test_data
          logical                   :: loc
          integer(ikind)            :: i,j


          test_data(1)  =  5.165      !flux_x_mass
          test_data(2)  =  5.039545704!flux_x_momentum_x
          test_data(3)  =  1.139892201!flux_x_momentum_y
          test_data(4)  =  3.210322639!flux_x_total_energy

          test_data(5)  =  8.425      !flux_y_mass
          test_data(6)  =  3.115780591!flux_y_momentum_x
          test_data(7)  =  6.840166285!flux_y_momentum_y
          test_data(8)  =  2.615690238!flux_y_total_energy

          i = 3
          j = 3


          !fluxes-x
          test_validated = .true.
          loc = is_test_validated(
     $         flux_x_mass_density(nodes,s,i,j),
     $         test_data(1),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_x_mass: '', L3)', loc


          test_validated = .true.
          loc = is_test_validated(
     $         flux_x_momentum_x(nodes,s,i,j,dx,dy),
     $         test_data(2),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_x_momentum_x: '', L3)', loc


          test_validated = .true.
          loc = is_test_validated(
     $         flux_x_momentum_y(nodes,s,i,j,dx,dy),
     $         test_data(3),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_x_momentum_y: '', L3)', loc


          test_validated = .true.
          loc = is_test_validated(
     $         flux_x_total_energy(nodes,s,i,j,dx,dy),
     $         test_data(4),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_x_total_energy: '', L3)', loc


          !fluxes-y
          test_validated = .true.
          loc = is_test_validated(
     $         flux_y_mass_density(nodes,s,i,j),
     $         test_data(5),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_y_mass: '', L3)', loc


          test_validated = .true.
          loc = is_test_validated(
     $         flux_y_momentum_x(nodes,s,i,j,dx,dy),
     $         test_data(6),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_y_momentum_x: '', L3)', loc


          test_validated = .true.
          loc = is_test_validated(
     $         flux_y_momentum_y(nodes,s,i,j,dx,dy),
     $         test_data(7),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_y_momentum_y: '', L3)', loc


          test_validated = .true.
          loc = is_test_validated(
     $         flux_y_total_energy(nodes,s,i,j,dx,dy),
     $         test_data(8),
     $         detailled)
          test_validated = test_validated.and.loc
          if(detailled) print '(''test flux_y_total_energy: '', L3)', loc

        end function test_ns2d_fluxes

      end program test_ns2d_eq_program
