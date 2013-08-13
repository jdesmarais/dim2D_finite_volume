      program test_bc_operators

        use bc_operators_class , only : bc_operators
        use cg_operators_class , only : cg_operators
        use field_class        , only : field
        use parameters_constant, only : periodic_xy_choice
        use parameters_input   , only : bc_choice
        use parameters_kind    , only : ikind, rkind

        implicit none

        
        !<operators tested
        type(field) :: field_tested

        integer(ikind), parameter :: nx=10
        integer(ikind), parameter :: ny=12
        integer       , parameter :: ne=1

        type(cg_operators) :: s
        type(bc_operators) :: bc_used 


        !<CPU recorded times
        real    :: time1, time2

        !<test parameters
        logical, parameter          :: detailled=.true.
        integer(ikind)              :: i,j
        real(rkind), dimension(nx,2):: test_north
        real(rkind), dimension(nx,2):: test_south
        real(rkind), dimension(2,8) :: test_east
        real(rkind), dimension(2,8) :: test_west
        logical                     :: test_validated


        !<get the initial CPU time
        call CPU_TIME(time1)


        !<allocate the tables for the field
        call field_tested%allocate_tables(nx,ny,ne)


        !<initialize the tables for the field
        do j=1, ny
           do i=1, nx
              field_tested%nodes(i,j,1) = i + (j-1)*nx
           end do
        end do


        !<initialize the test data
        select case(bc_choice)
          case(periodic_xy_choice)

             do j=1,2
                do i=3,nx-2
                   test_north(i,j) = i + (j-1+2)*nx
                   test_south(i,j) = i + (j-1+8)*nx
                end do
             end do

             test_north(1,1)=27
             test_north(2,1)=28
             test_north(1,2)=37
             test_north(2,2)=38

             test_north(9,1)=23
             test_north(10,1)=24
             test_north(9,2)=33
             test_north(10,2)=34

             test_south(1,1)=87
             test_south(2,1)=88
             test_south(1,2)=97
             test_south(2,2)=98

             test_south(9,1)=83
             test_south(10,1)=84
             test_south(9,2)=93
             test_south(10,2)=94

             do j=1,8
                do i=1,2
                   test_east(i,j) = i+2 + (j-1+2)*nx
                   test_west(i,j)=  i+6 + (j-1+2)*nx
                end do
             end do             

          case default
             print '(''test_bc'')'
             stop 'test for this boundary condition not implemented'
        end select 


        !<perform the test
        select case(bc_choice)
          case(periodic_xy_choice)
             call bc_used%apply_bc_on_nodes(field_tested,s)
             
             !<check if the boundary conditions are applied correctly
             !north and south tests
             test_validated=.true.
             i=1
             j=1
             do while(test_validated.and.(j.le.2))
                do while(test_validated.and.(i.le.nx))
                   test_validated=
     $                  field_tested%nodes(i,j+10,1).eq.test_north(i,j)
                   test_validated=test_validated.and.
     $                  field_tested%nodes(i,j,1).eq.test_south(i,j)
                   i=i+1           
                end do
                j=j+1
             end do
             
             
             !east and west tests
             if(test_validated) then
                i=1
                j=1
             end if
             do while(test_validated.and.(j.le.8))
                do while(test_validated.and.(i.le.2))
                   test_validated=
     $                  field_tested%nodes(i,j+2,1).eq.test_west(i,j)
                   test_validated=test_validated.and.
     $                  field_tested%nodes(i+8,j+2,1).eq.test_east(i,j)
                   i=i+1           
                end do
                j=j+1
             end do

             if(.not.test_validated) then
                print *, 'test_failed at: ', i,j
             end if
 
          case default
             print '(''test_bc'')'
             stop 'test for this boundary condition not implemented'
        end select


        !<print if the test is validated
        print *, 'test_validated: ', test_validated        

      end program test_bc_operators
