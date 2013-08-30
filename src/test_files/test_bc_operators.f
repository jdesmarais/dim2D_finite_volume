      !> @file
      !> test file for the object 'bc_operators'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the application of the boundary conditions
      !> on the gridpoints and compare the results with the
      !> expected data
      !
      !> @date
      ! 13_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_bc_operators

        use bc_operators_class , only : bc_operators
        use cg_operators_class , only : cg_operators
        use dim2d_eq_class     , only : dim2d_eq
        use field_class        , only : field
        use parameters_constant, only : periodic_xy_choice
        use parameters_input   , only : nx,ny,ne,bc_choice
        use parameters_kind    , only : ikind, rkind

        implicit none

        
        !<operators tested
        type(field)        :: field_tested
        type(dim2d_eq)     :: p_model
        type(cg_operators) :: s
        type(bc_operators) :: bc_used


        !<CPU recorded times
        real    :: time1, time2

        !<test parameters
        logical, parameter             :: detailled=.true.
        integer(ikind)                 :: i,j
        integer                        :: k
        real(rkind), dimension(10,2,ne):: test_north
        real(rkind), dimension(10,2,ne):: test_south
        real(rkind), dimension(2,8,ne) :: test_east
        real(rkind), dimension(2,8,ne) :: test_west
        logical                        :: test_validated


        !<test specifications
        if((nx.ne.10).or.(ny.ne.12).or.(ne.ne.4)) then
           stop 'the test requires (nx,ny,ne)=(10,12,4)'
        end if


        !<get the initial CPU time
        call CPU_TIME(time1)


        !<initialize the tables for the field
        do k=1, ne
           do j=1, ny
              do i=1, nx
                 field_tested%nodes(i,j,k) = 100*(k-1) + i + (j-1)*nx
              end do
           end do
        end do


        !<initialize the test data
        select case(bc_choice)
          case(periodic_xy_choice)

             do k=1, ne
                do j=1,2
                   do i=3,nx-2
                      test_north(i,j,k) = 100*(k-1) + i + (j-1+2)*nx
                      test_south(i,j,k) = 100*(k-1) + i + (j-1+8)*nx
                   end do
                end do
             end do

             test_north(1,1,1)=27
             test_north(2,1,1)=28
             test_north(1,2,1)=37
             test_north(2,2,1)=38

             test_north(9 ,1,1)=23
             test_north(10,1,1)=24
             test_north(9 ,2,1)=33
             test_north(10,2,1)=34

             test_south(1,1,1)=87
             test_south(2,1,1)=88
             test_south(1,2,1)=97
             test_south(2,2,1)=98

             test_south(9 ,1,1)=83
             test_south(10,1,1)=84
             test_south(9 ,2,1)=93
             test_south(10,2,1)=94

             do k=2,ne

                test_north(1 ,1,k) = test_north(1 ,1,k-1) + 100
                test_north(2 ,1,k) = test_north(2 ,1,k-1) + 100
                test_north(1 ,2,k) = test_north(1 ,2,k-1) + 100
                test_north(2 ,2,k) = test_north(2 ,2,k-1) + 100

                test_north(9 ,1,k) = test_north(9 ,1,k-1) + 100
                test_north(10,1,k) = test_north(10,1,k-1) + 100
                test_north(9 ,2,k) = test_north(9 ,2,k-1) + 100
                test_north(10,2,k) = test_north(10,2,k-1) + 100

                test_south(1 ,1,k) = test_south(1 ,1,k-1) + 100
                test_south(2 ,1,k) = test_south(2 ,1,k-1) + 100
                test_south(1 ,2,k) = test_south(1 ,2,k-1) + 100
                test_south(2 ,2,k) = test_south(2 ,2,k-1) + 100

                test_south(9 ,1,k) = test_south(9 ,1,k-1) + 100
                test_south(10,1,k) = test_south(10,1,k-1) + 100
                test_south(9 ,2,k) = test_south(9 ,2,k-1) + 100
                test_south(10,2,k) = test_south(10,2,k-1) + 100

             end do


             do k=1,ne
                do j=1,8
                   do i=1,2
                      test_east(i,j,k) = 100*(k-1) + i+2 + (j-1+2)*nx
                      test_west(i,j,k)=  100*(k-1) + i+6 + (j-1+2)*nx
                   end do
                end do
             end do    

c$$$          case default
c$$$             print '(''test_bc'')'
c$$$             stop 'test for this boundary condition not implemented'
        end select 


        !< apply the boundary conditions
        call bc_used%apply_bc_on_nodes(field_tested,p_model,s)


        !< perform the test
        select case(bc_choice)
          case(periodic_xy_choice)

             
             !<check if the boundary conditions are applied correctly
             !north and south tests
             test_validated=.true.
             k=1
             do while(test_validated.and.(k.le.ne))
                j=1
                do while(test_validated.and.(j.le.2))
                   i=1
                   do while(test_validated.and.(i.le.nx))

                      test_validated=
     $                  field_tested%nodes(i,j+10,k).eq.test_north(i,j,k)

                      test_validated=test_validated.and.
     $                  field_tested%nodes(i,j,k).eq.test_south(i,j,k)

                      i=i+1           
                   end do
                   j=j+1
                end do
                k=k+1
             end do
             
             
             !east and west tests
             if(test_validated) then
                i=1
                j=1
                k=1
             end if

             k=1
             do while(test_validated.and.(k.le.ne))
                j=1
                do while(test_validated.and.(j.le.8))
                   i=1
                   do while(test_validated.and.(i.le.2))

                      test_validated=
     $                  field_tested%nodes(i,j+2,k).eq.test_west(i,j,k)

                      test_validated=test_validated.and.
     $                  field_tested%nodes(i+8,j+2,k).eq.test_east(i,j,k)

                      i=i+1           
                   end do
                   j=j+1
                end do
                k=k+1
             end do

             if(.not.test_validated) then
                print *, 'test_failed at: ', i,j,k
             end if
 
c$$$          case default
c$$$             print '(''test_bc'')'
c$$$             stop 'test for this boundary condition not implemented'
        end select


        !<print if the test is validated
        print *, 'test_validated: ', test_validated        

      end program test_bc_operators
