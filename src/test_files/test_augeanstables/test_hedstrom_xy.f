      program test_hedstrom_xy

        use bc_operators_class, only :
     $     bc_operators

        use hedstrom_xy_module, only :
     $     compute_x_timedev_with_openbc,
     $     compute_timedev_xlayer

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
        type(bc_operators)               :: bc_used

        logical :: detailled


        if((nx.ne.7).or.(ny.ne.5)) then
           stop 'the test requires (nx,ny)=(7,5)'
        end if


        print '()'
        print '(''*************************'')'
        print '(''WARNING: use wave1d model'')'
        print '(''*************************'')'
        print '()'


        !initialization of the nodes
        call initialize_nodes(nodes,dx,dy)


        !test compute_x_timedev_with_openbc
        detailled = .false.
        print '(''----------------------------------'')'
        print '(''test_compute_x_timedev_with_openbc'')'
        print '(''----------------------------------'')'
        call test_compute_x_timedev_with_openbc(
     $       nodes,dx,dy,
     $       p_model,
     $       detailled)
        print '()'


        !test compute_timedev_x_layer
        detailled = .false.
        print '(''----------------------------------'')'
        print '(''test_compute_timedev_xlayer'')'
        print '(''----------------------------------'')'
        call test_compute_timedev_xlayer(
     $       nodes,dx,dy,
     $       p_model,
     $       detailled)
        print '()'


        !test apply_bc_on_timedev_2ndorder
        detailled=.false.
        print '(''----------------------------------'')'
        print '(''test_apply_bc_on_timedev'')'
        print '(''----------------------------------'')'
        call test_apply_bc_on_timedev_2ndorder(
     $       nodes,dx,dy,
     $       p_model,
     $       bc_used,
     $       detailled)
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


        subroutine test_compute_x_timedev_with_openbc(
     $     nodes,dx,dy,
     $     p_model,
     $     detailled)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind)                     , intent(in) :: dx
          real(rkind)                     , intent(in) :: dy
          type(pmodel_eq)                 , intent(in) :: p_model
          logical                         , intent(in) :: detailled

          
          real(rkind), dimension(nx,ny,ne) :: test_data
          real(rkind), dimension(nx,ny,ne) :: time_dev
          integer(ikind)    :: i,j
          integer           :: k
          logical           :: test_validated
          logical           :: loc
          character(len=45) :: fmt
          real(rkind)       :: dy_s

          dy_s = dy

          !initialize the test data
          call initialize_test_data(test_data)


          !computation of the time derivatives using open b.c.
          do j=1, ny
          
             i=1
             time_dev(i,j,:) = compute_x_timedev_with_openbc(
     $            nodes, i, j, p_model, dx,
     $            gradient_x_x_oneside_L0, incoming_left)
          
             i=bc_size
             time_dev(i,j,:) = compute_x_timedev_with_openbc(
     $            nodes, i, j, p_model, dx,
     $            gradient_x_x_oneside_L1, incoming_left)
          
             i=nx-1
             time_dev(i,j,:) = compute_x_timedev_with_openbc(
     $            nodes, i, j, p_model, dx,
     $            gradient_x_x_oneside_R1, incoming_right)
          
             i=nx
             time_dev(i,j,:) = compute_x_timedev_with_openbc(
     $            nodes, i, j, p_model, dx,
     $            gradient_x_x_oneside_R0, incoming_right)
          
          end do


          !checking the data
          if(detailled) then
             do k=1,ne
                do j=1,ny

                   do i=1,bc_size
                      loc = is_test_validated(
     $                     time_dev(i,j,1), test_data(i,j,1), detailled)
                      fmt ='(''test('',I2,'','',I2,'','',I2,''): '',L3)'
                      print fmt, i,j,k,loc
                   end do
                   
                   do i=nx-bc_size+1,nx
                      loc = is_test_validated(
     $                     time_dev(i,j,1), test_data(i,j,1), detailled)
                      fmt ='(''test('',I2,'','',I2,'','',I2,''): '',L3)'
                      print fmt, i,j,k,loc

                   end do
                end do
             end do
          
          else
             test_validated=.true.
             do k=1,ne
                do j=1,ny
                   do i=1,bc_size
                      loc = is_test_validated(time_dev(i,j,1), test_data(i,j,1), detailled)
                      test_validated=test_validated.and.loc 
                   end do
                   
                   do i=nx-bc_size+1,nx
                      loc = is_test_validated(time_dev(i,j,1), test_data(i,j,1), detailled)
                      test_validated=test_validated.and.loc 
                   end do
                end do
             end do
             print '(''test_validated: '',L1)', test_validated
          end if

        end subroutine test_compute_x_timedev_with_openbc


        subroutine test_compute_timedev_xlayer(
     $     nodes,dx,dy,
     $     p_model,
     $     detailled)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind)                     , intent(in) :: dx
          real(rkind)                     , intent(in) :: dy
          type(pmodel_eq)                 , intent(in) :: p_model
          logical                         , intent(in) :: detailled

          
          real(rkind), dimension(nx,ny,ne)   :: test_data
          real(rkind), dimension(nx,ny,ne)   :: time_dev
          real(rkind), dimension(nx,ny+1,ne) :: flux_y
          integer(ikind)    :: i,j
          integer           :: k
          logical           :: test_validated
          logical           :: loc
          character(len=45) :: fmt


          !initialize the test data
          call initialize_test_data(test_data)


          !initialize the y flux at zero
          do k=1,ne
             do j=1,ny+1
                do i=1,nx
                   flux_y(i,j,k)=0.0d0
                end do
             end do
          end do
          

          !computation of the time derivatives using open b.c.
          do j=1, ny
          
             i=1
             call compute_timedev_xlayer(
     $            nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_L0, incoming_left,
     $            time_dev)
          
             i=bc_size
             call compute_timedev_xlayer(
     $            nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_L1, incoming_left,
     $            time_dev)
          
             i=nx-1
             call compute_timedev_xlayer(
     $            nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_R1, incoming_right,
     $            time_dev)
          
             i=nx
             call compute_timedev_xlayer(
     $            nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_R0, incoming_right,
     $            time_dev)
          
          end do


          !checking the data
          if(detailled) then
             do k=1,ne
                do j=1,ny

                   do i=1,bc_size
                      loc = is_test_validated(
     $                     time_dev(i,j,1), test_data(i,j,1), detailled)
                      fmt ='(''test('',I2,'','',I2,'','',I2,''): '',L3)'
                      print fmt, i,j,k,loc
                   end do
                   
                   do i=nx-bc_size+1,nx
                      loc = is_test_validated(
     $                     time_dev(i,j,1), test_data(i,j,1), detailled)
                      fmt ='(''test('',I2,'','',I2,'','',I2,''): '',L3)'
                      print fmt, i,j,k,loc

                   end do
                end do
             end do
          
          else
             test_validated=.true.
             do k=1,ne
                do j=1,ny
                   do i=1,bc_size
                      loc = is_test_validated(time_dev(i,j,1), test_data(i,j,1), detailled)
                      test_validated=test_validated.and.loc 
                   end do
                   
                   do i=nx-bc_size+1,nx
                      loc = is_test_validated(time_dev(i,j,1), test_data(i,j,1), detailled)
                      test_validated=test_validated.and.loc 
                   end do
                end do
             end do
             print '(''test_validated: '',L1)', test_validated
          end if

        end subroutine test_compute_timedev_xlayer


        subroutine test_apply_bc_on_timedev_2ndorder(
     $     nodes,dx,dy,
     $     p_model,
     $     bc_used,
     $     detailled)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind)                     , intent(in) :: dx
          real(rkind)                     , intent(in) :: dy
          type(pmodel_eq)                 , intent(in) :: p_model
          type(bc_operators)              , intent(in) :: bc_used
          logical                         , intent(in) :: detailled

          
          real(rkind), dimension(nx,ny,ne)   :: test_data
          real(rkind), dimension(nx,ny,ne)   :: time_dev
          real(rkind), dimension(nx+1,ny,ne) :: flux_x
          real(rkind), dimension(nx,ny+1,ne) :: flux_y
          integer(ikind)    :: i,j
          integer           :: k
          logical           :: test_validated
          logical           :: loc
          character(len=45) :: fmt


          !initialize the test data
          call initialize_test_data(test_data)


          !initialize the y flux at zero
          do k=1,ne
             do j=1,ny+1
                do i=1,nx
                   flux_y(i,j,k)=0.0d0
                end do
             end do
          end do


          !initialize the flux_x randomly
          call srand(10)

          do k=1,ne
             do j=1,ny
                do i=1,nx+1
                   flux_x(i,j,k) = RAND()
                end do
             end do
          end do
          

          !computation of the time derivatives using open b.c.
          call bc_used%apply_bc_on_timedev(
     $         nodes,dx,dy,
     $         p_model,
     $         flux_x,flux_y,
     $         time_dev)


          !checking the data
          if(detailled) then
             do k=1,ne
                do j=bc_size+1,ny-bc_size

                   do i=1,bc_size
                      loc = is_test_validated(
     $                     time_dev(i,j,1), test_data(i,j,1), detailled)
                      fmt ='(''test('',I2,'','',I2,'','',I2,''): '',L3)'
                      print fmt, i,j,k,loc
                   end do
                   
                   do i=nx-bc_size+1,nx
                      loc = is_test_validated(
     $                     time_dev(i,j,1), test_data(i,j,1), detailled)
                      fmt ='(''test('',I2,'','',I2,'','',I2,''): '',L3)'
                      print fmt, i,j,k,loc

                   end do
                end do
             end do
          
          else
             test_validated=.true.
             do k=1,ne
                do j=bc_size+1,ny-bc_size
                   do i=1,bc_size
                      loc = is_test_validated(time_dev(i,j,1), test_data(i,j,1), detailled)
                      test_validated=test_validated.and.loc 
                   end do
                   
                   do i=nx-bc_size+1,nx
                      loc = is_test_validated(time_dev(i,j,1), test_data(i,j,1), detailled)
                      test_validated=test_validated.and.loc 
                   end do
                end do
             end do
             print '(''test_validated: '',L1)', test_validated
          end if

        end subroutine test_apply_bc_on_timedev_2ndorder
          

        subroutine initialize_test_data(test_data)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          !position
          test_data(1,1,1)=-6.8175
          test_data(2,1,1)=-0.825
          test_data(6,1,1)= 0.883
          test_data(7,1,1)=-0.319

          test_data(1,2,1)=0.585
          test_data(2,2,1)=2.06
          test_data(6,2,1)=-1.57
          test_data(7,2,1)=-4.345

          test_data(1,3,1)=22.1
          test_data(2,3,1)=6.7175
          test_data(6,3,1)= 0.735
          test_data(7,3,1)= 4.275

          test_data(1,4,1)=-6.79
          test_data(2,4,1)=1.0875
          test_data(6,4,1)= 1.2725
          test_data(7,4,1)=6.575

          test_data(1,5,1)=0.71
          test_data(2,5,1)=2.1
          test_data(6,5,1)=-4.945
          test_data(7,5,1)=-2.575

          !velocity_x
          test_data(1,1,2)=-6.8175
          test_data(2,1,2)=-0.825
          test_data(6,1,2)=-0.883
          test_data(7,1,2)= 0.319
                        
          test_data(1,2,2)=0.585
          test_data(2,2,2)=2.06
          test_data(6,2,2)= 1.57
          test_data(7,2,2)= 4.345
                        
          test_data(1,3,2)=22.1
          test_data(2,3,2)=6.7175
          test_data(6,3,2)=-0.735
          test_data(7,3,2)=-4.275
                        
          test_data(1,4,2)=-6.79
          test_data(2,4,2)=1.0875
          test_data(6,4,2)=-1.2725
          test_data(7,4,2)=-6.575
                        
          test_data(1,5,2)=0.71
          test_data(2,5,2)=2.1
          test_data(6,5,2)=4.945
          test_data(7,5,2)=2.575       

        end subroutine initialize_test_data

      end program test_hedstrom_xy
