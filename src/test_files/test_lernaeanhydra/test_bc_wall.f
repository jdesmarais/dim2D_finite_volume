      !> @file
      !> test file for the object 'bc_operators'
      !> for wall boundary conditions
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the application of the wall boundary
      !> conditions on the gridpoints and compare the
      !> results with the expected data
      !
      !> @date
      ! 25_09_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      !in this test, the application of the boundary conditions on the 
      !nodes is checked by comparing the saved data with the data
      !computed by the wall b.c.
      !the application of the boundary conditions on the fluxes is only
      !checked by looking if the fluxes modified are the correct ones,
      !if the fluxes for the mass density are equal to zero, and if the
      !energy fluxes are equal to zero unless the wall is heated
      !-----------------------------------------------------------------
      program test_bc_wall

        use bc_operators_class , only : bc_operators
        use sd_operators_class , only : sd_operators
        use pmodel_eq_class    , only : pmodel_eq
        use parameters_constant, only : wall_xy_choice
        use parameters_input   , only : nx,ny,ne,bc_choice,bc_size
        use parameters_kind    , only : ikind, rkind

        implicit none

        
        !<operators tested
        real(rkind), dimension(nx,ny,ne) :: nodes
        real(rkind)                      :: dx
        real(rkind)                      :: dy
        type(pmodel_eq)                  :: p_model
        type(sd_operators)               :: s
        type(bc_operators)               :: bc_used


        !<CPU recorded times
        real    :: time1, time2

        !<test parameters
        logical, parameter                 :: detailled=.true.
        integer(ikind)                     :: i,j
        integer                            :: k
        integer    , dimension(4)          :: prefactor_x
        integer    , dimension(4)          :: prefactor_y
        real(rkind), dimension(nx+1,ny,ne) :: flux_x
        real(rkind), dimension(nx,ny+1,ne) :: flux_y
        real(rkind), dimension(10,2,ne)    :: test_north
        real(rkind), dimension(10,2,ne)    :: test_south
        real(rkind), dimension(2,8,ne)     :: test_east
        real(rkind), dimension(2,8,ne)     :: test_west
        logical                            :: test_validated, wall_test


        !<test specifications
        if((nx.ne.10).or.(ny.ne.12).or.(ne.ne.4)) then
           stop 'the test requires (nx,ny,ne)=(10,12,4)'
        end if

        if(bc_choice.ne.wall_xy_choice) then
           stop 'the test is made for wall bc'
        end if

        !<get the initial CPU time
        call CPU_TIME(time1)


        !<initialize the tables for the field
        do k=1, ne
           do j=1, ny
              do i=1, nx
                 nodes(i,j,k) = 100*(k-1) + i + (j-1)*nx
              end do
           end do
        end do

        !dx=0.5
        !dy=0.6

        !<initialize the flux tables
        do k=1, ne
           do j=1, ny
              do i=1,nx+1
                 flux_x(i,j,k) = k
              end do
           end do
        end do

        do k=1, ne
           do j=1, ny+1
              do i=1,nx
                 flux_y(i,j,k) = k
              end do
           end do
        end do


        !<initialize the test data
        prefactor_x = [1,-1,-1,1]
        prefactor_y = [1,-1,-1,1]

        do k=1, ne
           do j=1,2
              do i=3,nx-2
                 test_north(i,j,k) = prefactor_y(k)*(100*(k-1) + i + (11-j-1)*nx)
                 test_south(i,j,k) = prefactor_y(k)*(100*(k-1) + i + (5-j-1)*nx)
              end do
           end do
        end do

        test_north(1,1,1)=94
        test_north(2,1,1)=93
        test_north(1,2,1)=84
        test_north(2,2,1)=83

        test_north(9 ,1,1)=98
        test_north(10,1,1)=97
        test_north(9 ,2,1)=88
        test_north(10,2,1)=87

        test_south(1,1,1)=34
        test_south(2,1,1)=33
        test_south(1,2,1)=24
        test_south(2,2,1)=23

        test_south(9 ,1,1)=38
        test_south(10,1,1)=37
        test_south(9 ,2,1)=28
        test_south(10,2,1)=27

        do k=2,ne

           test_north(1 ,1,k) = prefactor_x(k)*prefactor_y(k)*(test_north(1 ,1,1) + (k-1)*100)
           test_north(2 ,1,k) = prefactor_x(k)*prefactor_y(k)*(test_north(2 ,1,1) + (k-1)*100)
           test_north(1 ,2,k) = prefactor_x(k)*prefactor_y(k)*(test_north(1 ,2,1) + (k-1)*100)
           test_north(2 ,2,k) = prefactor_x(k)*prefactor_y(k)*(test_north(2 ,2,1) + (k-1)*100)

           test_north(9 ,1,k) = prefactor_x(k)*prefactor_y(k)*(test_north(9 ,1,1) + (k-1)*100)
           test_north(10,1,k) = prefactor_x(k)*prefactor_y(k)*(test_north(10,1,1) + (k-1)*100)
           test_north(9 ,2,k) = prefactor_x(k)*prefactor_y(k)*(test_north(9 ,2,1) + (k-1)*100)
           test_north(10,2,k) = prefactor_x(k)*prefactor_y(k)*(test_north(10,2,1) + (k-1)*100)

           test_south(1 ,1,k) = prefactor_x(k)*prefactor_y(k)*(test_south(1 ,1,1) + (k-1)*100)
           test_south(2 ,1,k) = prefactor_x(k)*prefactor_y(k)*(test_south(2 ,1,1) + (k-1)*100)
           test_south(1 ,2,k) = prefactor_x(k)*prefactor_y(k)*(test_south(1 ,2,1) + (k-1)*100)
           test_south(2 ,2,k) = prefactor_x(k)*prefactor_y(k)*(test_south(2 ,2,1) + (k-1)*100)

           test_south(9 ,1,k) = prefactor_x(k)*prefactor_y(k)*(test_south(9 ,1,1) + (k-1)*100)
           test_south(10,1,k) = prefactor_x(k)*prefactor_y(k)*(test_south(10,1,1) + (k-1)*100)
           test_south(9 ,2,k) = prefactor_x(k)*prefactor_y(k)*(test_south(9 ,2,1) + (k-1)*100)
           test_south(10,2,k) = prefactor_x(k)*prefactor_y(k)*(test_south(10,2,1) + (k-1)*100)

        end do


        do k=1,ne
           do j=1,8
              do i=1,2
                 test_east(i,j,k) = prefactor_x(k)*(100*(k-1) + 9-i + (j-1+2)*nx)
                 test_west(i,j,k)=  prefactor_x(k)*(100*(k-1) + 5-i + (j-1+2)*nx)
              end do
           end do
        end do

        !< initialize the boundary conditions
        call bc_used%ini(p_model)

        !< apply the boundary conditions
        call bc_used%apply_bc_on_nodes(nodes)


        !< perform the test
             
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
     $             nodes(i,j+10,k).eq.test_north(i,j,k)

                 test_validated=test_validated.and.
     $             nodes(i,j,k).eq.test_south(i,j,k)

                 i=i+1           
              end do
              j=j+1
           end do
           k=k+1
        end do
        
        if(.not.test_validated) then
           print *, 'test_failed at: ', i,j,k, nodes(i,j,k)
           stop 'nodes test failed for north and south'
        end if
        

        !east and west tests
        test_validated=.true.
        k=1
        do while(test_validated.and.(k.le.ne))
           j=1
           do while(test_validated.and.(j.le.8))
              i=1
              do while(test_validated.and.(i.le.2))

                 test_validated=
     $             nodes(i,j+2,k).eq.test_west(i,j,k)

                 test_validated=test_validated.and.
     $             nodes(i+8,j+2,k).eq.test_east(i,j,k)

                 i=i+1           
              end do
              j=j+1
           end do
           k=k+1
        end do

        if(.not.test_validated) then
           print *, 'test_failed at: ', i,j,k, nodes(i,j,k)
           stop 'nodes test failed for east and west'
        end if

        
        !<test the application of the wall b.c. on the fluxes
        call bc_used%apply_bc_on_fluxes(nodes,dx,dy,s,flux_x,flux_y)

        !< check if the points modified by the b.c. are only the 
        !> points that should be modified
        i=1
        j=1
        k=1
        test_validated=.true.
        do while(test_validated.and.(k.le.ne))
           j=1
           do while(test_validated.and.(j.le.ny))
              i=1
              do while(test_validated.and.(i.le.nx+1))
                 wall_test = (i.eq.(bc_size+1)).or.(i.eq.(nx-bc_size+1))
                 wall_test = wall_test.and.(j.ge.(bc_size+1))
                 wall_test = wall_test.and.(j.le.(ny-bc_size))
                 if(wall_test) then
                    test_validated=int(flux_x(i,j,k)).ne.k
                    test_validated=test_validated.and.(flux_x(i,j,1).eq.0)
                    test_validated=test_validated.and.(flux_x(i,j,4).eq.0)
                 else
                    test_validated=int(flux_x(i,j,k)).eq.k
                 end if
                 i=i+1
              end do
              j=j+1
           end do
           k=k+1
        end do
              
        if(.not.test_validated) then
           print *, 'test_failed at: ', i,j,k, flux_x(i,j,k), wall_test
           stop 'fluxes test failed for modified gridpts flux_x'
        end if   

 
        i=1
        j=1
        k=1
        test_validated=.true.
        do while(test_validated.and.(k.le.ne))
           j=1
           do while(test_validated.and.(j.le.ny+1))
              i=1
              do while(test_validated.and.(i.le.nx))
                 wall_test = (j.eq.(bc_size+1)).or.(j.eq.(ny-bc_size+1))
                 wall_test = wall_test.and.(i.ge.(bc_size+1))
                 wall_test = wall_test.and.(i.le.(nx-bc_size))
                 if(wall_test) then
                    test_validated=int(flux_y(i,j,k)).ne.k
                    test_validated=test_validated.and.(flux_y(i,j,1).eq.0)
                    test_validated=test_validated.and.(flux_y(i,j,4).eq.0)
                 else
                    test_validated=int(flux_y(i,j,k)).eq.k
                 end if
                 i=i+1
              end do
              j=j+1
           end do
           k=k+1
        end do

        if(.not.test_validated) then
           print *, 'test_failed at: ', i,j,k, flux_y(i,j,k), wall_test
           stop 'fluxes test failed for modified gridpts flux_y'
        end if


        !<print if the test is validated
        print *, 'test_validated: ', test_validated  

        call CPU_TIME(time2)
        print '(''time elapsed:'', F6.2)', time2-time1

        !print *, flux_x
        !print *, flux_y

      end program test_bc_wall
