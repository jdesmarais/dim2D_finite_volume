      !> @file
      !> test file for the object 'rk3tvd'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the integration of simple test equations
      !> using Runge-Kutta 3rd order TVD time integration
      !> scheme with periodic boundary conditions
      !
      !> @date
      ! 14_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_rk3tvd

        use field_abstract_class, only : field_abstract
        use parameters_constant , only : periodic_xy_choice
        use parameters_input    , only : nx,ny,ne,bc_choice
        use parameters_kind     , only : ikind, rkind
        use td_integrator_class , only : td_integrator

        implicit none

        
        !<operators tested
        type(field_abstract)             :: field_tested
        type(td_integrator)              :: ti
        real(rkind), dimension(nx,ny,ne) :: nodes
        real(rkind), parameter           :: dt=1.0
        real(rkind)                      :: dx
        real(rkind)                      :: dy


        !<test parameters
        logical, parameter         :: detailled=.false.
        integer(ikind)             :: i,j
        real(rkind), dimension(6,2):: test_data
        logical                    :: test_validated


        !<if nx<4, ny<4 then the test cannot be done
        if((nx.ne.10).or.(ny.ne.6).or.(ne.ne.1)) then
           stop 'the test needs: (nx,ny,ne)=(10,6,1)'
        end if

        if(bc_choice.ne.periodic_xy_choice) then
           stop 'the test needs periodic bc'
        end if
        

        !<initialize the tables for the field
        dx=1.0
        dy=1.0

        do j=1, ny
           do i=1, nx
              nodes(i,j,1) = i + (j-1)*nx
           end do
        end do

        call field_tested%ini()
        call field_tested%set_nodes(nodes)
        call field_tested%set_dx(dx)
        call field_tested%set_dy(dy)


        !<integrate the field for dt
        call ti%integrate(field_tested,dt)


        !<check the field after integration
        test_data(1,1) =-1409.028d0
        test_data(2,1) =-1471.806d0
        test_data(3,1) =-1496.111d0
        test_data(4,1) =-1472.889d0
        test_data(5,1) =-1538.861d0
        test_data(6,1) =-1478.306d0
        
        test_data(1,2) = 1600.972d0
        test_data(2,2) = 1538.194d0
        test_data(3,2) = 1513.889d0
        test_data(4,2) = 1537.111d0
        test_data(5,2) = 1471.139d0
        test_data(6,2) = 1531.694d0


        call field_tested%get_nodes(nodes)

        j=1
        i=1
        test_validated=.true.
        do while (test_validated.and.(j.le.(size(test_data,2))))
           do while(test_validated.and.(i.le.(size(test_data,1))))

              test_validated=is_test_validated(nodes(i+2,j+2,1),
     $                                         test_data(i,j))
              i=i+1

           end do
           j=j+1
        end do

        print *, nodes
        

        !<display the result of the test
        if(test_validated) then
           print '(''test_validated: '', 1L)', test_validated
        else
           print '(''WARNING: periodic b.c. necessary'')'
           print '(''test_failed at ('',I2,'','',I2,'')'')', i,j
        end if
        
        contains

        function is_test_validated(var,cst) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical                 :: test_validated

          if(detailled) then
             print *, int(var*1e2)
             print *, int(cst*1e2)
          end if
          
          test_validated=(
     $         int(var*100.)-
     $         sign(int(abs(cst*100.)),int(cst*100.))).eq.0
          
        end function is_test_validated

      end program test_rk3tvd
