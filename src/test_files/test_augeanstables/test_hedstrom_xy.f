      program test_hedstrom_xy

        use ifport

        use bc_operators_class, only :
     $       bc_operators

        use hedstrom_xy_module, only :
     $       compute_x_timedev_with_openbc,
     $       compute_timedev_xlayer

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       bc_size

        use pmodel_eq_class , only :
     $       pmodel_eq

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_L1,
     $       gradient_x_x_oneside_R1,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_L1,
     $       gradient_y_y_oneside_R1,
     $       gradient_y_y_oneside_R0


        implicit none

        real(rkind)                      :: t
        real(rkind), dimension(nx)       :: x_map
        real(rkind), dimension(ny)       :: y_map
        real(rkind), dimension(nx,ny,ne) :: nodes
        real(rkind)                      :: dx, dy
        type(pmodel_eq)                  :: p_model
        type(bc_operators)               :: bc_used

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        if((nx.ne.7).or.(ny.ne.5)) then
           stop 'the test requires (nx,ny)=(7,5)'
        end if


        print '()'
        print '(''*************************'')'
        print '(''WARNING: use wave1d model'')'
        print '(''*************************'')'
        print '()'

        test_validated = .true.


        !initialization of the nodes
        call initialize_nodes(x_map,y_map,nodes,dx,dy)


        !test compute_x_timedev_with_openbc
        detailled = .false.
        test_loc = test_compute_x_timedev_with_openbc(
     $       nodes,dx,dy,
     $       p_model,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_x_timedev_with_openbc: '', L1)', test_loc
        print '()'


        !test compute_timedev_x_layer
        detailled = .false.
        test_loc = test_compute_timedev_xlayer(
     $       t, x_map, y_map,
     $       nodes,dx,dy,
     $       p_model,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_timedev_xlayer: '', L1)', test_loc
        print '()'


        !test apply_bc_on_timedev_2ndorder
        detailled=.false.
        test_loc = test_apply_bc_on_timedev_2ndorder(
     $       nodes,dx,dy,
     $       p_model,
     $       bc_used,
     $       detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_apply_bc_on_timedev: '',L1)', test_loc
        print '()'

        print '(''test_hedstrom_xy: '',L1)', test_validated
        print '(''------------------------'')'
        print '()'

        contains


        function is_test_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical    , intent(in) :: detailled
          logical                 :: test_validated

          integer(ikind) :: diff

          if(detailled) then
             print '(I12)', int(var*1e5)
             print '(I12)', int(cst*1e5)
             print '(I12)', abs(nint(var*1e5)-nint(cst*1e5))
          end if
          
          diff = abs(nint(var*1e5)-nint(cst*1e5))

          test_validated=diff.le.1
          
        end function is_test_validated        


        subroutine initialize_nodes(x_map,y_map,nodes,dx,dy)

          implicit none

          real(rkind), dimension(nx)      , intent(out) :: x_map
          real(rkind), dimension(ny)      , intent(out) :: y_map
          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
          real(rkind)                     , intent(out) :: dx
          real(rkind)                     , intent(out) :: dy
          
          integer(ikind) :: i,j,k

          dx=0.5
          dy=0.6

          do i=1, size(x_map,1)
             x_map(i) = (i-1)*dx
          end do

          do i=1, size(y_map,1)
             y_map(i) = (i-1)*dy
          end do

          !position
          nodes(1,1,1)=  0.5d0
          nodes(2,1,1)=  0.2d0
          nodes(3,1,1)=  1.2d0
          nodes(4,1,1)=  5.0d0
          nodes(5,1,1)=  0.6d0
          nodes(6,1,1)= -3.6d0
          nodes(7,1,1)=-6.52d0

          nodes(1,2,1)=  3.0d0
          nodes(2,2,1)=  4.2d0
          nodes(3,2,1)= 11.0d0
          nodes(4,2,1)= 10.6d0
          nodes(5,2,1)=  5.2d0
          nodes(6,2,1)=  1.2d0
          nodes(7,2,1)= 7.89d0

          nodes(1,3,1)=-14.2d0
          nodes(2,3,1)=   23d0
          nodes(3,3,1)=  9.8d0
          nodes(4,3,1)=  3.4d0
          nodes(5,3,1)=  9.1d0
          nodes(6,3,1)=  6.7d0
          nodes(7,3,1)= 4.12d0

          nodes(1,4,1)= 2.45d0
          nodes(2,4,1)=  0.2d0
          nodes(3,4,1)=  9.0d0
          nodes(4,4,1)=  5.4d0
          nodes(5,4,1)= -2.3d0
          nodes(6,4,1)=  1.0d0
          nodes(7,4,1)=-5.62d0

          nodes(1,5,1)=  3.6d0
          nodes(2,5,1)=  0.1d0
          nodes(3,5,1)=  6.3d0
          nodes(4,5,1)=  8.9d0
          nodes(5,5,1)=-4.23d0
          nodes(6,5,1)=  8.9d0
          nodes(7,5,1)= 8.95d0


          !velocity_x
          nodes(1,1,2)= 7.012d0
          nodes(2,1,2)=-6.323d0
          nodes(3,1,2)= 3.012d0
          nodes(4,1,2)=   4.5d0
          nodes(5,1,2)=   9.6d0
          nodes(6,1,2)=  9.57d0
          nodes(7,1,2)= 6.012d0
                    
          nodes(1,2,2)=4.26d0
          nodes(2,2,2)=4.23d0
          nodes(3,2,2)= 4.5d0
          nodes(4,2,2)=7.56d0
          nodes(5,2,2)=7.21d0
          nodes(6,2,2)=5.62d0
          nodes(7,2,2)=3.62d0
                    
          nodes(1,3,2)=0.23d0
          nodes(2,3,2)=7.23d0
          nodes(3,3,2)= 3.1d0
          nodes(4,3,2)= 8.9d0
          nodes(5,3,2)= 9.3d0
          nodes(6,3,2)=1.29d0
          nodes(7,3,2)=7.26d0
                    
          nodes(1,4,2)=8.23d0
          nodes(2,4,2)=-3.1d0
          nodes(3,4,2)=6.03d0
          nodes(4,4,2)=6.25d0
          nodes(5,4,2)=5.12d0
          nodes(6,4,2)=0.36d0
          nodes(7,4,2)=6.89d0
                    
          nodes(1,5,2)= 3.2d0
          nodes(2,5,2)=8.12d0
          nodes(3,5,2)= 8.9d0
          nodes(4,5,2)= 4.2d0
          nodes(5,5,2)= 7.8d0
          nodes(6,5,2)= 6.3d0
          nodes(7,5,2)= 1.2d0


          do k=3,3
             do j=1,ny
                do i=1,nx
                   nodes(i,j,k) = 0.0d0
                end do
             end do
          end do

        end subroutine initialize_nodes


        function test_compute_x_timedev_with_openbc(
     $     nodes,dx,dy,
     $     p_model,
     $     detailled)
     $     result(test_validated)

          implicit none
          
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind)                     , intent(in) :: dx
          real(rkind)                     , intent(in) :: dy
          type(pmodel_eq)                 , intent(in) :: p_model
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated

          
          real(rkind), dimension(nx,ny,ne) :: test_data
          real(rkind), dimension(nx,ny,ne) :: time_dev
          integer(ikind)    :: i,j
          integer           :: k
          logical           :: test_loc
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

          test_validated = .true.


          !checking the data
          do k=1,ne
             do j=1,ny
                do i=1,bc_size

                   test_loc = is_test_validated(
     $                  time_dev(i,j,k), test_data(i,j,k), detailled)
                   if((.not.test_loc).and.detailled) then
                      print '(''['',3I2,'']: '',F9.4,'' -> '',F9.4)',
     $                     i,j,k,time_dev(i,j,k), test_data(i,j,k)
                   end if
                   test_validated = test_validated.and.test_loc
                end do
                
                do i=nx-bc_size+1,nx
                   test_loc = is_test_validated(
     $                  time_dev(i,j,k), test_data(i,j,k), detailled)
                   if((.not.test_loc).and.detailled) then
                      print '(''['',3I2,'']: '',F9.4,'' -> '',F9.4)',
     $                     i,j,k,time_dev(i,j,k), test_data(i,j,k)
                   end if
                   test_validated = test_validated.and.test_loc

                end do
             end do
          end do

        end function test_compute_x_timedev_with_openbc


        function test_compute_timedev_xlayer(
     $     t,x_map,y_map,
     $     nodes,dx,dy,
     $     p_model,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind)                     , intent(in) :: t
          real(rkind), dimension(nx)      , intent(in) :: x_map
          real(rkind), dimension(ny)      , intent(in) :: y_map
          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind)                     , intent(in) :: dx
          real(rkind)                     , intent(in) :: dy
          type(pmodel_eq)                 , intent(in) :: p_model
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated

          
          real(rkind), dimension(nx,ny,ne)   :: test_data
          real(rkind), dimension(nx,ny,ne)   :: time_dev
          real(rkind), dimension(nx,ny+1,ne) :: flux_y
          integer(ikind)    :: i,j
          integer           :: k
          logical           :: test_loc


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
     $            t,x_map,y_map,nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_L0, incoming_left,
     $            time_dev)
          
             i=bc_size
             call compute_timedev_xlayer(
     $            t,x_map,y_map,nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_L1, incoming_left,
     $            time_dev)
          
             i=nx-1
             call compute_timedev_xlayer(
     $            t,x_map,y_map,nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_R1, incoming_right,
     $            time_dev)
          
             i=nx
             call compute_timedev_xlayer(
     $            t,x_map,y_map,nodes, i,j, dx,dy, p_model, flux_y,
     $            gradient_x_x_oneside_R0, incoming_right,
     $            time_dev)
          
          end do


          !checking the data
          test_validated=.true.

          do k=1,ne
             do j=1,ny
                do i=1,bc_size

                   test_loc = is_test_validated(
     $                  time_dev(i,j,k), test_data(i,j,k), detailled)
                   if((.not.test_loc).and.detailled) then
                      print '(''['',3I2,'']: '',F9.4,'' -> '',F9.4)',
     $                     i,j,k,time_dev(i,j,k), test_data(i,j,k)
                   end if
                   test_validated = test_validated.and.test_loc

                end do
                
                do i=nx-bc_size+1,nx

                   test_loc = is_test_validated(
     $                  time_dev(i,j,k), test_data(i,j,k), detailled)
                   if((.not.test_loc).and.detailled) then
                      print '(''['',3I2,'']: '',F9.4,'' -> '',F9.4)',
     $                     i,j,k,time_dev(i,j,k), test_data(i,j,k)
                   end if
                   test_validated = test_validated.and.test_loc
                   
                end do
             end do
          end do

        end function test_compute_timedev_xlayer


        function test_apply_bc_on_timedev_2ndorder(
     $     nodes,dx,dy,
     $     p_model,
     $     bc_used,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: nodes
          real(rkind)                     , intent(in) :: dx
          real(rkind)                     , intent(in) :: dy
          type(pmodel_eq)                 , intent(in) :: p_model
          type(bc_operators)              , intent(in) :: bc_used
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated
          
          real(rkind), dimension(nx,ny,ne)   :: test_data
          real(rkind), dimension(nx,ny,ne)   :: time_dev
          real(rkind), dimension(nx+1,ny,ne) :: flux_x
          real(rkind), dimension(nx,ny+1,ne) :: flux_y
          integer(ikind)    :: i,j
          integer           :: k
          logical           :: test_loc

          real(rkind), dimension(nx) :: x_map
          real(rkind), dimension(ny) :: y_map
          real(rkind)                :: t


          do i=1,nx
             x_map(i) = (i-1)*dx
          end do

          do j=1,ny
             y_map(j) = (j-1)*dy
          end do


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
     $         p_model,
     $         t,nodes,x_map,y_map,
     $         flux_x,flux_y,
     $         time_dev)


          !checking the data
          test_validated = .true.

          do k=1,ne
             do j=bc_size+1,ny-bc_size
                
                do i=1,bc_size
                   
                   test_loc = is_test_validated(
     $                  time_dev(i,j,k), test_data(i,j,k), detailled)
                   if((.not.test_loc).and.detailled) then
                      print '(''['',3I2,'']: '',F9.4,'' -> '',F9.4)',
     $                     i,j,k,time_dev(i,j,k), test_data(i,j,k)
                   end if
                   test_validated = test_validated.and.test_loc
                   
                end do
                
                do i=nx-bc_size+1,nx
                   
                   test_loc = is_test_validated(
     $                  time_dev(i,j,k), test_data(i,j,k), detailled)
                   if((.not.test_loc).and.detailled) then
                      print '(''['',3I2,'']: '',F9.4,'' -> '',F9.4)',
     $                     i,j,k,time_dev(i,j,k), test_data(i,j,k)
                   end if
                   test_validated = test_validated.and.test_loc
                   
                end do
             end do
          end do

        end function test_apply_bc_on_timedev_2ndorder
          

        subroutine initialize_test_data(test_data)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: test_data

          !position
          test_data(1,1,1)=-6.8175d0
          test_data(2,1,1)= -0.825d0
          test_data(6,1,1)=  0.883d0
          test_data(7,1,1)= -0.319d0

          test_data(1,2,1)= 0.585d0
          test_data(2,2,1)=  2.06d0
          test_data(6,2,1)= -1.57d0
          test_data(7,2,1)=-4.345d0

          test_data(1,3,1)=  22.1d0
          test_data(2,3,1)=6.7175d0
          test_data(6,3,1)= 0.735d0
          test_data(7,3,1)= 4.275d0

          test_data(1,4,1)=  -6.79d0
          test_data(2,4,1)= 1.0875d0
          test_data(6,4,1)= 1.2725d0
          test_data(7,4,1)=  6.575d0

          test_data(1,5,1)=  0.71d0
          test_data(2,5,1)=   2.1d0
          test_data(6,5,1)=-4.945d0
          test_data(7,5,1)=-2.575d0

          !velocity_x
          test_data(1,1,2)=-6.8175d0
          test_data(2,1,2)= -0.825d0
          test_data(6,1,2)= -0.883d0
          test_data(7,1,2)=  0.319d0
                        
          test_data(1,2,2)= 0.585d0
          test_data(2,2,2)=  2.06d0
          test_data(6,2,2)=  1.57d0
          test_data(7,2,2)= 4.345d0
                        
          test_data(1,3,2)=  22.1d0
          test_data(2,3,2)=6.7175d0
          test_data(6,3,2)=-0.735d0
          test_data(7,3,2)=-4.275d0
                        
          test_data(1,4,2)=  -6.79d0
          test_data(2,4,2)= 1.0875d0
          test_data(6,4,2)=-1.2725d0
          test_data(7,4,2)= -6.575d0
                        
          test_data(1,5,2)= 0.71d0
          test_data(2,5,2)=  2.1d0
          test_data(6,5,2)=4.945d0
          test_data(7,5,2)=2.575d0

        end subroutine initialize_test_data

      end program test_hedstrom_xy
