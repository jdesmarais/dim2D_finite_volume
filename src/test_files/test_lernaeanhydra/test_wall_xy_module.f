      program test_wall_xy_module

        use sd_operators_class , only : sd_operators
        use dim2d_parameters   , only : re,we,pr,viscous_r,cv_r
        use parameters_constant, only : wall_xy_choice
        use parameters_input   , only : nx,ny,ne,bc_choice
        use parameters_kind    , only : ikind, rkind
        use wall_xy_module     , only : wall_fx_momentum_x,
     $                                  wall_fx_momentum_y,
     $                                  wall_fy_momentum_x,
     $                                  wall_fy_momentum_y


        implicit none

        
        !<operators tested
        real(rkind), dimension(nx,ny,ne) :: nodes
        real(rkind)                      :: dx
        real(rkind)                      :: dy
        type(sd_operators)               :: s

        !<test parameters
        logical, parameter :: detailled=.false.
        integer(ikind)     :: i,j
        logical            :: test_param
        real(rkind)        :: flux_x_momentum_x
        real(rkind)        :: flux_x_momentum_y
        real(rkind)        :: flux_y_momentum_x
        real(rkind)        :: flux_y_momentum_y
        real(rkind)        :: test_flux_x_momentum_x
        real(rkind)        :: test_flux_x_momentum_y
        real(rkind)        :: test_flux_y_momentum_x
        real(rkind)        :: test_flux_y_momentum_y
        logical            :: global,local


        !<test specifications
        if((nx.ne.10).or.(ny.ne.12).or.(ne.ne.4)) then
           stop 'the test requires (nx,ny,ne)=(10,12,4)'
        end if

        if(bc_choice.ne.wall_xy_choice) then
           stop 'the test is made for wall bc'
        end if


        !<initialize the data for the field
        !<mass data
        nodes(1,1,1)=0.5
        nodes(2,1,1)=0.2
        nodes(3,1,1)=1.2
        nodes(4,1,1)=5.0

        nodes(1,2,1)=2.0
        nodes(2,2,1)=4.2
        nodes(3,2,1)=11.0
        nodes(4,2,1)=10.6

        nodes(1,3,1)=-14.2
        nodes(2,3,1)=23.0
        nodes(3,3,1)=9.8
        nodes(4,3,1)=3.4
      
        nodes(1,4,1)=2.45
        nodes(2,4,1)=0.2
        nodes(3,4,1)=9.0
        nodes(4,4,1)=5.4


        !<momentum_x data
        nodes(1,1,2)=9.5
        nodes(2,1,2)=9.8
        nodes(3,1,2)=8.8
        nodes(4,1,2)=5.0

        nodes(1,2,2)=8.0
        nodes(2,2,2)=5.8
        nodes(3,2,2)=-1.0
        nodes(4,2,2)=-0.6

        nodes(1,3,2)=24.2
        nodes(2,3,2)=-13.0
        nodes(3,3,2)=0.2
        nodes(4,3,2)=6.6
      
        nodes(1,4,2)=7.55
        nodes(2,4,2)=9.8
        nodes(3,4,2)=1.0
        nodes(4,4,2)=4.6


        !<momentum_y data
        nodes(1,1,3)=-8.5
        nodes(2,1,3)=-9.4
        nodes(3,1,3)=-6.4
        nodes(4,1,3)=5.0
                             
        nodes(1,2,3)=-4.0
        nodes(2,2,3)=2.6
        nodes(3,2,3)=23.0
        nodes(4,2,3)=21.8
                             
        nodes(1,3,3)=-52.6
        nodes(2,3,3)=59.0
        nodes(3,3,3)=19.4
        nodes(4,3,3)=0.2
                             
        nodes(1,4,3)=-2.65
        nodes(2,4,3)=-9.4
        nodes(3,4,3)=17.0
        nodes(4,4,3)=6.2

        
        !<total energy data
        nodes(1,1,4)=-1.5
        nodes(2,1,4)=-1.8
        nodes(3,1,4)=-0.8
        nodes(4,1,4)=3.0
                             
        nodes(1,2,4)=0.0
        nodes(2,2,4)=2.2
        nodes(3,2,4)=9.0
        nodes(4,2,4)=8.6
                             
        nodes(1,3,4)=-16.2
        nodes(2,3,4)=21.0
        nodes(3,3,4)=7.8
        nodes(4,3,4)=1.4
                             
        nodes(1,4,4)=0.45
        nodes(2,4,4)=-1.8
        nodes(3,4,4)=7.0
        nodes(4,4,4)=3.4


        !< initialize the dx and dy data
        dx=0.5
        dy=0.6


        !<check if the DIM2d parameters
        !>are the same
        test_param=(re.eq.5)
        test_param=test_param.and.(we.eq.10.0)
        test_param=test_param.and.(pr.eq.20.0)
        test_param=test_param.and.(viscous_r.eq.-1.5)
        test_param=test_param.and.(cv_r.eq.2.5)
        if(.not.test_param) then
           stop 'the dim2d parameters are not correct for the test'
        end if

        if(detailled) then
           !< print the data
           print *, 'mass_density'
           do j=1,4
              print '(4F7.2)', nodes(1:4,5-j,1)
           end do
           print *, ''
              
           print *, 'momentum_x'
           do j=1,4
              print '(4F7.2)', nodes(1:4,5-j,2)
           end do
           print *, ''
           
           print *, 'momentum_y'
           do j=1,4
              print '(4F7.2)', nodes(1:4,5-j,3)
           end do
           print *, ''
           
           print *, 'energy'
           do j=1,4
              print '(4F7.2)', nodes(1:4,5-j,4)
           end do
           print *, ''
           
           print *, 'dim2d_parameters'
           print '(''dx        '',F7.2)', dx
           print '(''dy        '',F7.2)', dy
           print '(''Re        '',F7.2)', re
           print '(''We        '',F7.2)', we
           print '(''Pr        '',F7.2)', pr
           print '(''viscous_r '',F7.2)', viscous_r
           print '(''cv_r      '',F7.2)', cv_r
           print *,''
        end if

        !< initialized the test data for the fluxes
        test_flux_x_momentum_x=-284.1147599
        test_flux_x_momentum_y=-0.588744589
        test_flux_y_momentum_x=0.648723257
        test_flux_y_momentum_y=-815.9718776



        !<compute the wall fluxes at (2,2)
        i=2
        j=2

        flux_x_momentum_x=wall_fx_momentum_x(nodes,dx,dy,s,i,j)
        flux_x_momentum_y=wall_fx_momentum_y(nodes,dx,s,i,j)

        flux_y_momentum_x=wall_fy_momentum_x(nodes,dy,s,i,j)
        flux_y_momentum_y=wall_fy_momentum_y(nodes,dx,dy,s,i,j)


        !<compare the data
        if(detailled) then
           print '(''program        | excel         '')'
           print '(''-------------------------------'')'
           print '(F14.7,'' | '', F14.7)',
     $          flux_x_momentum_x, test_flux_x_momentum_x
           print '(F14.7,'' | '', F14.7)',
     $          flux_x_momentum_y, test_flux_x_momentum_y
           print '(F14.7,'' | '', F14.7)',
     $          flux_y_momentum_x, test_flux_y_momentum_x
           print '(F14.7,'' | '', F14.7)',
     $          flux_y_momentum_y, test_flux_y_momentum_y
        end if

        global = .true.
        local = is_test_validated(
     $       flux_x_momentum_x,
     $       test_flux_x_momentum_x)
        global = global.and.local

        local = is_test_validated(
     $       flux_x_momentum_y,
     $       test_flux_x_momentum_y)
        global = global.and.local

        local = is_test_validated(
     $       flux_y_momentum_x,
     $       test_flux_y_momentum_x)
        global = global.and.local

        local = is_test_validated(
     $       flux_y_momentum_y,
     $       test_flux_y_momentum_y)
        global = global.and.local

        print '(''test validated :'', 1L)', global


        contains

        function is_test_validated(var,cst) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical                 :: test_validated

          test_validated=(
     $         int(var*10000.)-
     $         sign(int(abs(cst*10000.)),int(cst*10000.))).eq.0
          
        end function is_test_validated

      end program test_wall_xy_module
