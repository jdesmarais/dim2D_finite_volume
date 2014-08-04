      program test_hedstrom_xy

        use hedstrom_xy_module, only :
     $     compute_x_timedev_with_openbc

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right

        use parameters_kind , only : ikind, rkind
        use parameters_input, only : nx,ny,ne,bc_size
        use pmodel_eq_class , only : pmodel_eq

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_L1,
     $       gradient_x_x_oneside_R1,
     $       gradient_x_x_oneside_R0

        implicit none

        real(rkind), dimension(nx,ny,ne) :: nodes
        real(rkind)                      :: dx, dy
        type(pmodel_eq)                  :: p_model
        real(rkind), dimension(nx,ny,ne) :: test_data
        real(rkind), dimension(nx,ny,ne) :: time_dev

        integer(ikind) :: i,j,k
        logical        :: detailled
        logical        :: loc
        logical        :: test_validated


        if((nx.ne.7).or.(ny.ne.5)) then
           stop 'the test requires (nx,ny)=(7,5)'
        end if


        !initialization of the nodes
        call initialize_nodes(nodes,dx,dy)

        !initialization of the test data
        call initialize_test_data(test_data)


        !computation of the time derivatives using open b.c.
        do j=1, ny

           i=1
           time_dev(i,j,:) = compute_x_timedev_with_openbc(
     $          nodes, i, j, p_model, dx,
     $          gradient_x_x_oneside_L0, incoming_left)

           i=bc_size
           time_dev(i,j,:) = compute_x_timedev_with_openbc(
     $          nodes, i, j, p_model, dx,
     $          gradient_x_x_oneside_L1, incoming_left)

           i=nx-1
           time_dev(i,j,:) = compute_x_timedev_with_openbc(
     $          nodes, i, j, p_model, dx,
     $          gradient_x_x_oneside_R1, incoming_right)

           i=nx
           time_dev(i,j,:) = compute_x_timedev_with_openbc(
     $          nodes, i, j, p_model, dx,
     $          gradient_x_x_oneside_R0, incoming_right)

        end do


        !checking the data
        detailled = .false.
        if(detailled) then
           do k=1,2
              do j=1,5
                 do i=1,2
                    loc = is_test_validated(time_dev(i,j,1), test_data(i,j,1), detailled)
                    print '(''test('',I2,'','',I2,'','',I2,''): '', L3)', i,j,k,loc
                 end do
                 
                 do i=5,6
                    loc = is_test_validated(time_dev(i,j,1), test_data(i,j,1), detailled)
                    print '(''test('',I2,'','',I2,'','',I2,''): '', L3)', i,j,k,loc
                 end do
              end do
           end do

        else
           test_validated=.true.
           do k=1,2
              do j=1,5
                 do i=1,2
                    loc = is_test_validated(time_dev(i,j,1), test_data(i,j,1), detailled)
                    test_validated=test_validated.and.loc 
                 end do
                 
                 do i=5,6
                    loc = is_test_validated(time_dev(i,j,1), test_data(i,j,1), detailled)
                    test_validated=test_validated.and.loc 
                 end do
              end do
           end do
           print '(''test_validated: '',L1)', test_validated
        end if

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

          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
          real(rkind)                     , intent(out) :: dx
          real(rkind)                     , intent(out) :: dy
          
          dx=0.5
          dy=0.6

          !position
          nodes(1,1,1)=0.5
          nodes(2,1,1)=0.2
          nodes(3,1,1)=1.2
          nodes(4,1,1)=5.0
          nodes(5,1,1)=0.6
          nodes(6,1,1)=-3.6
          nodes(7,1,1)=-6.52

          nodes(1,2,1)=3.0
          nodes(2,2,1)=4.2
          nodes(3,2,1)=11.0
          nodes(4,2,1)=10.6
          nodes(5,2,1)=5.2
          nodes(6,2,1)=1.2
          nodes(7,2,1)=7.89

          nodes(1,3,1)=-14.2
          nodes(2,3,1)=23
          nodes(3,3,1)=9.8
          nodes(4,3,1)=3.4
          nodes(5,3,1)=9.1
          nodes(6,3,1)=6.7
          nodes(7,3,1)=4.12

          nodes(1,4,1)=2.45
          nodes(2,4,1)=0.2
          nodes(3,4,1)=9.0
          nodes(4,4,1)=5.4
          nodes(5,4,1)=-2.3
          nodes(6,4,1)=1.0
          nodes(7,4,1)=-5.62

          nodes(1,5,1)=3.6
          nodes(2,5,1)=0.1
          nodes(3,5,1)=6.3
          nodes(4,5,1)=8.9
          nodes(5,5,1)=-4.23
          nodes(6,5,1)=8.9
          nodes(7,5,1)=8.95


          !velocity_x
          nodes(1,1,2)=7.012
          nodes(2,1,2)=-6.323
          nodes(3,1,2)=3.012
          nodes(4,1,2)=4.5
          nodes(5,1,2)=9.6
          nodes(6,1,2)=9.57
          nodes(7,1,2)=6.012
                    
          nodes(1,2,2)=4.26
          nodes(2,2,2)=4.23
          nodes(3,2,2)=4.5
          nodes(4,2,2)=7.56
          nodes(5,2,2)=7.21
          nodes(6,2,2)=5.62
          nodes(7,2,2)=3.62
                    
          nodes(1,3,2)=0.23
          nodes(2,3,2)=7.23
          nodes(3,3,2)=3.1
          nodes(4,3,2)=8.9
          nodes(5,3,2)=9.3
          nodes(6,3,2)=1.29
          nodes(7,3,2)=7.26
                    
          nodes(1,4,2)=8.23
          nodes(2,4,2)=-3.1
          nodes(3,4,2)=6.03
          nodes(4,4,2)=6.25
          nodes(5,4,2)=5.12
          nodes(6,4,2)=0.36
          nodes(7,4,2)=6.89
                    
          nodes(1,5,2)=3.2
          nodes(2,5,2)=8.12
          nodes(3,5,2)=8.9
          nodes(4,5,2)=4.2
          nodes(5,5,2)=7.8
          nodes(6,5,2)=6.3
          nodes(7,5,2)=1.2          

        end subroutine initialize_nodes


        subroutine initialize_test_data(test_data)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          !position
          test_data(1,1,1)=-13.635
          test_data(2,1,1)=-1.65
          test_data(6,1,1)= 1.766
          test_data(7,1,1)=-0.638

          test_data(1,2,1)=1.17
          test_data(2,2,1)=4.12
          test_data(6,2,1)=-3.14
          test_data(7,2,1)=-8.69

          test_data(1,3,1)=44.2
          test_data(2,3,1)=13.435
          test_data(6,3,1)= 1.47
          test_data(7,3,1)= 8.55

          test_data(1,4,1)=-13.58
          test_data(2,4,1)=2.175
          test_data(6,4,1)= 2.545
          test_data(7,4,1)=13.15

          test_data(1,5,1)=1.42
          test_data(2,5,1)=4.2
          test_data(6,5,1)=-9.89
          test_data(7,5,1)=-5.15

          !velocity_x
          test_data(1,1,2)=-13.635
          test_data(2,1,2)=-1.65
          test_data(6,1,2)=-1.766
          test_data(7,1,2)=0.638
                    
          test_data(1,2,2)=1.17
          test_data(2,2,2)=4.12
          test_data(6,2,2)=3.14
          test_data(7,2,2)=8.69
                    
          test_data(1,3,2)=44.2
          test_data(2,3,2)=13.435
          test_data(6,3,2)=-1.47
          test_data(7,3,2)=-8.55
                    
          test_data(1,4,2)=-13.58
          test_data(2,4,2)=2.175
          test_data(6,4,2)=-2.545
          test_data(7,4,2)=-13.15
                    
          test_data(1,5,2)=1.42
          test_data(2,5,2)=4.2
          test_data(6,5,2)=9.89
          test_data(7,5,2)=5.15         

        end subroutine initialize_test_data

      end program test_hedstrom_xy
