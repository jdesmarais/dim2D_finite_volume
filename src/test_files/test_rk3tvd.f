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

        use cg_operators_class , only : cg_operators
        use field_class        , only : field
        use fv_operators_class , only : fv_operators
        use parameters_kind    , only : ikind, rkind
        use simpletest_eq_class, only : simpletest_eq
        use rk3tvd_class       , only : rk3tvd

        implicit none

        
        !<operators tested
        type(field) :: field_tested

        integer(ikind), parameter :: nx=10
        integer(ikind), parameter :: ny=6
        integer       , parameter :: ne=1

        type(cg_operators)     :: sd
        type(simpletest_eq)    :: p_model
        type(fv_operators)     :: td
        type(rk3tvd)           :: ti
        real(rkind), parameter :: dt=1.0


        !<CPU recorded times
        real    :: time1, time2

        !<test parameters
        logical, parameter         :: detailled=.false.
        integer(ikind)             :: i,j
        real(rkind), dimension(6,2):: test_data
        logical                    :: test_validated


        !<get the initial CPU time
        call CPU_TIME(time1)


        !<allocate the tables for the field
        call field_tested%allocate_tables(nx,ny,ne)


        !<initialize the tables for the field
        field_tested%dx=1.0
        field_tested%dy=1.0

        do j=1, ny
           do i=1, nx
              field_tested%nodes(i,j,1) = i + (j-1)*nx
           end do
        end do


        !<integrate the field for dt
        call ti%integrate(field_tested,sd,p_model,td,dt)


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

        j=1
        i=1
        test_validated=.true.
        do while (test_validated.and.(j.le.(size(test_data,2))))
           do while(test_validated.and.(i.le.(size(test_data,1))))

              test_validated=is_test_validated(
     $             field_tested%nodes(i+2,j+2,1),
     $             test_data(i,j))
              i=i+1

           end do
           j=j+1
        end do
        

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
