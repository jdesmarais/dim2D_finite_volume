      !> @file
      !> test file for the module 'dim2d_prim'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the primary variable operators by comparing the
      !> results with previous computations
      !
      !> @date
      ! 09_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_dim2d_prim

        use dim2d_prim_module
        use field_class      , only : field
        use dim2d_parameters , only : viscous_r,re,pr,we,cv_r
        use parameters_kind  , only : ikind, rkind

        implicit none
        
        
        !<operators tested
        type(field)               :: field_tested
        integer(ikind), parameter :: nx=4
        integer(ikind), parameter :: ny=4
        integer       , parameter :: ne=4
        
        !<CPU recorded times
        real    :: time1, time2

        !<test parameters
        logical, parameter         :: detailled=.true.
        integer(ikind)             :: i,j
        real(rkind)                :: computed_data
        real(rkind), dimension(19) :: test_data
        logical                    :: test_validated


        !<get the initial CPU time
        call CPU_TIME(time1)

        !<allocate the tables for the field
        call field_tested%allocate_tables(nx,ny,ne)


        !<initialize the tables for the field
        field_tested%dx=0.5
        field_tested%dy=0.6

        !<initialize the mass density
        field_tested%nodes(1,1,1)=0.5d0
        field_tested%nodes(2,1,1)=0.2d0
        field_tested%nodes(3,1,1)=1.2d0
        field_tested%nodes(4,1,1)=5.0d0

        field_tested%nodes(1,2,1)=2.0d0
        field_tested%nodes(2,2,1)=4.2d0
        field_tested%nodes(3,2,1)=11.0d0
        field_tested%nodes(4,2,1)=10.6d0

        field_tested%nodes(1,3,1)=-14.2d0
        field_tested%nodes(2,3,1)= 23.0d0
        field_tested%nodes(3,3,1)=  9.8d0
        field_tested%nodes(4,3,1)=  3.4d0

        field_tested%nodes(1,4,1)= 2.45d0
        field_tested%nodes(2,4,1)= 0.2d0
        field_tested%nodes(3,4,1)= 9.0d0
        field_tested%nodes(4,4,1)= 5.4d0

        !<initialize the momentum_x
        field_tested%nodes(1,1,2)= 9.5d0
        field_tested%nodes(2,1,2)= 9.8d0
        field_tested%nodes(3,1,2)= 8.8d0
        field_tested%nodes(4,1,2)= 5.0d0

        field_tested%nodes(1,2,2)= 8.0d0
        field_tested%nodes(2,2,2)= 5.8d0
        field_tested%nodes(3,2,2)=-1.0d0
        field_tested%nodes(4,2,2)=-0.6d0

        field_tested%nodes(1,3,2)=-24.2d0
        field_tested%nodes(2,3,2)=-13.0d0
        field_tested%nodes(3,3,2)= 0.2d0
        field_tested%nodes(4,3,2)= 6.6d0

        field_tested%nodes(1,4,2)= 7.55d0
        field_tested%nodes(2,4,2)= 9.8d0
        field_tested%nodes(3,4,2)= 1.0d0
        field_tested%nodes(4,4,2)= 4.6d0

        !<initialize the momentum_y
        field_tested%nodes(1,1,3)=-8.5d0
        field_tested%nodes(2,1,3)=-9.4d0
        field_tested%nodes(3,1,3)=-6.4d0
        field_tested%nodes(4,1,3)= 5.0d0

        field_tested%nodes(1,2,3)=-4.0d0
        field_tested%nodes(2,2,3)= 2.6d0
        field_tested%nodes(3,2,3)= 23.0d0
        field_tested%nodes(4,2,3)= 21.8d0

        field_tested%nodes(1,3,3)=-52.6d0
        field_tested%nodes(2,3,3)= 59.0d0
        field_tested%nodes(3,3,3)= 19.4d0
        field_tested%nodes(4,3,3)= 0.20d0

        field_tested%nodes(1,4,3)=-2.65d0
        field_tested%nodes(2,4,3)=-9.40d0
        field_tested%nodes(3,4,3)= 17.0d0
        field_tested%nodes(4,4,3)= 6.20d0

        !<initialize the total energy
        field_tested%nodes(1,1,4)=-1.5d0
        field_tested%nodes(2,1,4)=-1.8d0
        field_tested%nodes(3,1,4)=-0.8d0
        field_tested%nodes(4,1,4)= 3.0d0

        field_tested%nodes(1,2,4)= 0.0d0
        field_tested%nodes(2,2,4)= 2.2d0
        field_tested%nodes(3,2,4)= 9.0d0
        field_tested%nodes(4,2,4)= 8.6d0

        field_tested%nodes(1,3,4)=-16.2d0
        field_tested%nodes(2,3,4)= 21.0d0
        field_tested%nodes(3,3,4)= 7.8d0
        field_tested%nodes(4,3,4)= 1.4d0

        field_tested%nodes(1,4,4)= 0.45d0
        field_tested%nodes(2,4,4)=-1.8d0
        field_tested%nodes(3,4,4)= 7.0d0
        field_tested%nodes(4,4,4)= 3.4d0
        
        !<test the operators defined dim2d_prim
        i=2 !<index tested in the data along the x-axis
        j=2 !<index tested in the data along the y-axis

        !<test_data initialization
        test_data(1) =  4.2d0      !<mass
        test_data(2) =  5.8d0      !<momentum_x
        test_data(3) =  2.6d0      !<momentum_y
        test_data(4) =  2.2d0      !<total_energy
        test_data(5) =  1.380952d0 !<velocity_x
        test_data(6) =  0.619048d0 !<velocity_y
        test_data(7) = -103.23047d0!<classical pressure
        test_data(8) =  6.71678d0  !<temperature_eff
        test_data(9) =  8.009524d0 !<qx_transport_x
        test_data(10)=  3.590476d0 !<qy_transport_y
        test_data(11)=  3.590476d0 !<qx_transport_y
        test_data(12)=  1.609524d0 !<qy_transport_y
        test_data(13)=  3.038095d0 !<energy transport_x
        test_data(14)=  1.361905d0 !<energy transport_y
        test_data(15)= -0.833333d0 !<capillarity_pressure
        test_data(16)= -1.15079d0  !<capillarity_pressure_xwork
        test_data(17)= -0.51587d0  !<capillarity_pressure_ywork
        test_data(18)= -142.55637d0!<classic pressure work along x-axis> 
        test_data(19)= -63.90458d0 !<classic pressure work along y-axis>


        !< print the dim2d parameters used for the test
        if(detailled) then
           print '(''viscous_r: '', F16.6)', viscous_r
           print '(''re:        '', F16.6)', re
           print '(''pr:        '', F16.6)', pr
           print '(''we:        '', F16.6)', we
           print '(''cv_r:      '', F16.6)', cv_r
        else
           print '(''WARNING: this test is designed for:'')'
           print '(''viscous_r: '', F16.6)', -1.5
           print '(''re:        '', F16.6)', 5
           print '(''pr:        '', F16.6)', 20
           print '(''we:        '', F16.6)', 10
           print '(''cv_r:      '', F16.6)', 2.5
           print '(''it allows to see errors easily'')'
        end if

        !<test of the operators
        if(detailled) then
           !DEC$ FORCEINLINE RECURSIVE
           print '(''test mass_density: '',1L)',
     $          is_test_validated(
     $          mass_density(
     $          field_tested,
     $          i,j),
     $          test_data(1))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test momentum_x: '',1L)',
     $          is_test_validated(
     $          momentum_x(
     $          field_tested,
     $          i,j),
     $          test_data(2))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test momentum_y: '',1L)',
     $          is_test_validated(
     $          momentum_y(
     $          field_tested,
     $          i,j),
     $          test_data(3))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test energy: '',1L)',
     $          is_test_validated(
     $          total_energy(
     $          field_tested,
     $          i,j),
     $          test_data(4))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test velocity_x: '',1L)',
     $          is_test_validated(
     $          velocity_x(
     $          field_tested,
     $          i,j),
     $          test_data(5))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test velocity_y: '',1L)',
     $          is_test_validated(
     $          velocity_y(
     $          field_tested,
     $          i,j),
     $          test_data(6))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test classical_pressure: '',1L)',
     $          is_test_validated(
     $          classical_pressure(
     $          field_tested,
     $          i,j),
     $          test_data(7))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test temperature_eff: '',1L)',
     $          is_test_validated(
     $          temperature_eff(
     $          field_tested,
     $          i,j),
     $          test_data(8))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test qx_transport_x: '',1L)',
     $          is_test_validated(
     $          qx_transport_x(
     $          field_tested,
     $          i,j),
     $          test_data(9))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test qy_transport_x: '',1L)',
     $          is_test_validated(
     $          qy_transport_x(
     $          field_tested,
     $          i,j),
     $          test_data(10))

           !DEC$ FORCEINLINE RECURSIVE
           computed_data = qx_transport_y(field_tested,i,j)
           print '(''test qx_transport_y: '',1L)',
     $          is_test_validated(
     $          computed_data,
     $          test_data(11))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test qy_transport_y: '',1L)',
     $          is_test_validated(
     $          qy_transport_y(
     $          field_tested,
     $          i,j),
     $          test_data(12))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test energy_transport_x: '',1L)',
     $          is_test_validated(
     $          energy_transport_x(
     $          field_tested,
     $          i,j),
     $          test_data(13))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test energy_transport_y: '',1L)',
     $          is_test_validated(
     $          energy_transport_y(
     $          field_tested,
     $          i,j),
     $          test_data(14))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test capillarity_pressure: '',1L)',
     $          is_test_validated(
     $          capillarity_pressure(
     $          field_tested,
     $          i,j),
     $          test_data(15))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test capillarity_pressure_xwork: '',1L)',
     $          is_test_validated(
     $          capillarity_pressure_xwork(
     $          field_tested,
     $          i,j),
     $          test_data(16))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test capillarity_pressure_ywork: '',1L)',
     $          is_test_validated(
     $          capillarity_pressure_ywork(
     $          field_tested,
     $          i,j),
     $          test_data(17))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test classical_pressure_xwork: '',1L)',
     $          is_test_validated(
     $          classical_pressure_xwork(
     $          field_tested,
     $          i,j),
     $          test_data(18))

           !DEC$ FORCEINLINE RECURSIVE
           print '(''test classical_pressure_ywork: '',1L)',
     $          is_test_validated(
     $          classical_pressure_ywork(
     $          field_tested,
     $          i,j),
     $          test_data(19))

        else

           test_validated=
     $          is_test_validated(
     $          mass_density(
     $          field_tested,
     $          i,j),
     $          test_data(1))
           
           test_validated=test_validated.and.
     $          is_test_validated(
     $          momentum_x(
     $          field_tested,
     $          i,j),
     $          test_data(2))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          momentum_y(
     $          field_tested,
     $          i,j),
     $          test_data(3))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          total_energy(
     $          field_tested,
     $          i,j),
     $          test_data(4))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          velocity_x(
     $          field_tested,
     $          i,j),
     $          test_data(5))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          velocity_y(
     $          field_tested,
     $          i,j),
     $          test_data(6))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          classical_pressure(
     $          field_tested,
     $          i,j),
     $          test_data(7))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          temperature_eff(
     $          field_tested,
     $          i,j),
     $          test_data(8))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          qx_transport_x(
     $          field_tested,
     $          i,j),
     $          test_data(9))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          qy_transport_x(
     $          field_tested,
     $          i,j),
     $          test_data(10))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          qx_transport_y(
     $          field_tested,
     $          i,j),
     $          test_data(11))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          qy_transport_y(
     $          field_tested,
     $          i,j),
     $          test_data(12))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          energy_transport_x(
     $          field_tested,
     $          i,j),
     $          test_data(13))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          energy_transport_y(
     $          field_tested,
     $          i,j),
     $          test_data(14))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          capillarity_pressure(
     $          field_tested,
     $          i,j),
     $          test_data(15))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          capillarity_pressure_xwork(
     $          field_tested,
     $          i,j),
     $          test_data(16))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          capillarity_pressure_ywork(
     $          field_tested,
     $          i,j),
     $          test_data(17))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          classical_pressure_xwork(
     $          field_tested,
     $          i,j),
     $          test_data(18))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          capillarity_pressure_ywork(
     $          field_tested,
     $          i,j),
     $          test_data(19))           

           print '(''test_validated: '', 1L)', test_validated

        end if


        !<get the last CPU time
        call CPU_TIME(time2)
        print *, 'time elapsed: ', time2-time1


        contains

        function is_test_validated(var,cst) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical                 :: test_validated

          if(detailled) then
             print *, int(var*1e5)
             print *, int(cst*1e5)
          end if
          
          test_validated=(
     $         int(var*10000.)-
     $         sign(int(abs(cst*10000.)),int(cst*10000.))).eq.0
          
        end function is_test_validated


      end program test_dim2d_prim
