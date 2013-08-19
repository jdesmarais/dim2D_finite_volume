      !> @file
      !> test file for the object 'cg_operators'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the space discretization operators by comparing the
      !> results with previous computations
      !
      !> @date
      ! 07_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_cg_operators

        use cg_operators_class, only : cg_operators
        use dim2d_prim_module , only : mass_density
        use field_class       , only : field
        use parameters_input  , only : nx,ny
        use parameters_kind   , only : ikind, rkind

        
        implicit none


        !<operators tested
        type(field)        :: field_tested
        type(cg_operators) :: sd_operators_tested

        !<CPU recorded times
        real :: time1, time2

        !<test parameters
        logical, parameter         :: detailled=.true.
        integer(ikind)             :: i,j
        real(rkind), dimension(12) :: test_data
        logical                    :: test_validated


        !<if nx<4, ny<4 then the test cannot be done
        if((nx.lt.4).or.(ny.lt.4)) then
           stop 'nx and ny must be greater than 4 for the test'
        end if


        !<get the initial CPU time
        call CPU_TIME(time1)


        !<initialize the tables for the field
        field_tested%dx=0.5
        field_tested%dy=0.6

        field_tested%nodes(1,1,1)=0.5
        field_tested%nodes(2,1,1)=0.2
        field_tested%nodes(3,1,1)=1.2
        field_tested%nodes(4,1,1)=5.0

        field_tested%nodes(1,2,1)=3.0
        field_tested%nodes(2,2,1)=4.2
        field_tested%nodes(3,2,1)=11.0
        field_tested%nodes(4,2,1)=10.6

        field_tested%nodes(1,3,1)=-14.2
        field_tested%nodes(2,3,1)=23
        field_tested%nodes(3,3,1)=9.8
        field_tested%nodes(4,3,1)=3.4

        field_tested%nodes(1,4,1)=2.45
        field_tested%nodes(2,4,1)=0.2
        field_tested%nodes(3,4,1)=9.0
        field_tested%nodes(4,4,1)=5.4

        
        !<test the operators defined in cg_operators
        i=2 !<index tested in the data along the x-axis
        j=2 !<index tested in the data along the y-axis

        test_data(1) =  7.7333333d0 !<test f
        test_data(2) =  13.6d0      !<test dfdx
        test_data(3) =  16.395833d0 !<test dfdy
        test_data(4) = -3.199999d0  !<test d2fdx2
        test_data(5) =  13.680555d0 !<test d2fdy2
        test_data(6) = -23.666667d0 !<test d2fdxdy

        test_data(7) =  15.833333d0 !<test g
        test_data(8) =  18.0625d0   !<test dgdx
        test_data(9) =  31.333333d0 !<test dgdy
        test_data(10)= -108.64999d0 !<test d2gdx2
        test_data(11)= -37.222210d0 !<test d2gdy2
        test_data(12)=  26.666666d0 !<test d2gdxdy


        if(detailled) then
           
           !TAG INLINE
           print '(''test %f: '',1L)',
     $          is_test_validated(
     $          sd_operators_tested%f(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(1))

           !TAG INLINE
           print '(''test %dfdx: '',1L)',
     $          is_test_validated(
     $          sd_operators_tested%dfdx(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(2))
           !TAG INLINE
           print '(''test %dfdy: '',1L)',
     $          is_test_validated(
     $          sd_operators_tested%dfdy(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(3))
           !TAG INLINE
           print '(''test %d2fdx2: '',1L)',
     $          is_test_validated(
     $          sd_operators_tested%d2fdx2(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(4))

           !TAG INLINE
           print '(''test %d2fdy2: '',1L)',
     $          is_test_validated(
     $          sd_operators_tested%d2fdy2(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(5))

           !TAG INLINE
           print '(''test %d2fdxdy: '',1L)',
     $          is_test_validated(
     $          sd_operators_tested%d2fdxdy(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(6))

           !TAG INLINE
           print '(''test %g: '',1L)',
     $          is_test_validated(
     $          sd_operators_tested%g(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(7))

           !TAG INLINE
           print '(''test %dgdx: '',1L)',
     $          is_test_validated(
     $          sd_operators_tested%dgdx(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(8))

           !TAG INLINE
           print '(''test %dgdy: '',1L)',
     $          is_test_validated(
     $          sd_operators_tested%dgdy(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(9))

           !TAG INLINE
           print '(''test %d2gdx2: '',1L)',
     $          is_test_validated(
     $          sd_operators_tested%d2gdx2(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(10))

           !TAG INLINE
           print '(''test %d2gdy2: '',1L)',
     $          is_test_validated(
     $          sd_operators_tested%d2gdy2(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(11))

           !TAG INLINE
           print '(''test %d2gdxdy: '',1L)',
     $          is_test_validated(
     $          sd_operators_tested%d2gdxdy(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(12))
        else
           test_validated=is_test_validated(
     $          sd_operators_tested%f(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(1))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          sd_operators_tested%dfdx(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(2))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          sd_operators_tested%dfdy(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(3))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          sd_operators_tested%d2fdx2(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(3))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          sd_operators_tested%d2fdx2(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(4))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          sd_operators_tested%d2fdy2(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(5))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          sd_operators_tested%d2fdxdy(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(6))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          sd_operators_tested%g(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(7))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          sd_operators_tested%dgdx(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(8))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          sd_operators_tested%dgdy(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(9))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          sd_operators_tested%d2gdx2(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(10))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          sd_operators_tested%d2gdy2(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(11))

           test_validated=test_validated.and.
     $          is_test_validated(
     $          sd_operators_tested%d2gdxdy(
     $          field_tested,
     $          i,j,
     $          mass_density),
     $          test_data(12))

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


      end program test_cg_operators
