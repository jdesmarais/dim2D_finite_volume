      program test_bc_operators_hedstrom_xy_corners

        use bc_operators_class, only :
     $     bc_operators

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right

        use parameters_constant, only :
     $       left,right

        use parameters_input, only :
     $       nx,ny,ne,bc_size
        
        use parameters_kind, only :
     $       ikind, rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_fd_module, only :
     $       gradient_x_interior,
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_interior,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_R0

        use wave2d_parameters, only :
     $       c

        implicit none


        logical :: test_validated
        logical :: test_loc
        logical :: detailled


        test_validated = .true.
        
        if(
     $       (nx.ne.6).or.
     $       (ny.ne.6).or.
     $       (ne.ne.3).or.
     $       (.not.is_test_validated(c,0.5d0,.false.))
     $       ) then
           print '(''nx.eq.6: '',L1)', nx.eq.6
           print '(''ny.eq.6: '',L1)', ny.eq.6
           print '(''ne.eq.3: '',L1)', ne.eq.3
           print '(''c=0.5  : '',L1)', is_test_validated(c,0.5d0,.false.)
           stop ''
        end if


        print '()'
        print '(''*************************'')'
        print '(''WARNING: use wave2d model'')'
        print '(''*************************'')'
        print '()'


        detailled = .true.
        test_loc = test_apply_bc_on_timedev(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_apply_bc_on_timedev: '',L1)', test_loc
        print '()'


        detailled = .true.
        test_loc = test_apply_bc_on_timedev_xy_corner(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_apply_bc_on_timedev_xy_corner: '',L1)', test_loc
        print '()'


        contains


        function test_apply_bc_on_timedev(detailled)
     $       result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind), dimension(nx,ny,ne)   :: nodes
          real(rkind), dimension(nx)         :: x_map
          real(rkind), dimension(ny)         :: y_map
          real(rkind)                        :: t
          real(rkind), dimension(nx+1,ny,ne) :: flux_x
          real(rkind), dimension(nx,ny+1,ne) :: flux_y
          real(rkind), dimension(nx,ny,ne)   :: timedev_test
          type(bc_operators)                 :: bc_operators_used
          type(pmodel_eq)                    :: p_model
          real(rkind), dimension(nx,ny,ne)   :: timedev

          
          logical :: test_loc
          integer :: i,j,k


          test_validated = .true.


          !initialize the nodes and data for test
          call initialize_nodes(x_map,y_map,nodes,timedev_test)


          !compute the time derivatives at the edge of
          !the computational domain resulting from the
          !application of the hedstrom_xy_corner b.c.
          call bc_operators_used%apply_bc_on_timedev(
     $         p_model,
     $         t,nodes,x_map,y_map,
     $         flux_x,flux_y,
     $         timedev)


          !test the corner timedev
          do k=1,ne
             do j=1, bc_size
                do i=1,bc_size
                   test_loc = is_test_validated(timedev(i,j,k),timedev_test(i,j,k),.false.)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print '(''['',3I2,'']: '',F10.4,'' -> '',F10.4)',
     $                     i,j,k,
     $                     timedev(i,j,k), timedev_test(i,j,k)
                   end if
                end do

                do i=nx-1,nx
                   test_loc = is_test_validated(timedev(i,j,k),timedev_test(i,j,k),.false.)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print '(''['',3I2,'']: '',F10.4,'' -> '',F10.4)',
     $                     i,j,k,
     $                     timedev(i,j,k), timedev_test(i,j,k)
                   end if
                end do
             end do

             do j=ny-1,ny
                do i=1,bc_size
                   test_loc = is_test_validated(timedev(i,j,k),timedev_test(i,j,k),.false.)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                     print '(''['',3I2,'']: '',F10.4,'' -> '',F10.4)',
     $                     i,j,k,
     $                     timedev(i,j,k), timedev_test(i,j,k)
                  end if
               end do
               
               do i=nx-1,nx
                  test_loc = is_test_validated(timedev(i,j,k),timedev_test(i,j,k),.false.)
                  test_validated = test_validated.and.test_loc
                  if(detailled.and.(.not.test_loc)) then
                     print '(''['',3I2,'']: '',F10.4,'' -> '',F10.4)',
     $                    i,j,k,
     $                    timedev(i,j,k), timedev_test(i,j,k)
                  end if
               end do
            end do
         end do                   
         
        end function test_apply_bc_on_timedev


        function test_apply_bc_on_timedev_xy_corner(detailled)
     $       result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          real(rkind), dimension(nx,ny,ne)   :: nodes
          real(rkind), dimension(nx)         :: x_map
          real(rkind), dimension(ny)         :: y_map
          real(rkind)                        :: t
          real(rkind), dimension(nx,ny,ne)   :: timedev_test
          type(bc_operators)                 :: bc_operators_used
          type(pmodel_eq)                    :: p_model
          real(rkind), dimension(nx,ny,ne)   :: timedev

          
          logical :: test_loc
          integer :: i,j,k
          logical :: side_x
          logical :: side_y


          test_validated = .true.


          !initialize the nodes and data for test
          call initialize_nodes(x_map,y_map,nodes,timedev_test)


          !compute the time derivatives at the edge of
          !the computational domain resulting from the
          !application of the hedstrom_xy_corner b.c.

          !SW corner
          side_x = left
          side_y = left

          j=1
          i=1
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_L0)

          i=2
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_L0)

          j=2
          i=1
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_interior)

          i=2
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_interior,
     $         gradient_y_interior)


          !SE corner
          side_x = right
          side_y = left

          j=1
          i=5
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_L0)

          i=6
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_L0)

          j=2
          i=5
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_interior,
     $         gradient_y_interior)

          i=6
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_interior)

          !NW corner
          side_x = left
          side_y = right

          j=ny-1
          i=1
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_interior)

          i=2
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_interior,
     $         gradient_y_interior)

          j=ny
          i=1
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_L0,
     $         gradient_y_y_oneside_R0)

          i=2
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_R0)


          !NE corner
          side_x = right
          side_y = right

          j=ny-1
          i=nx-1
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_interior,
     $         gradient_y_interior)

          i=nx
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_interior)

          j=ny
          i=nx-1
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_R0)

          i=nx
          timedev(i,j,:) = bc_operators_used%apply_bc_on_timedev_xy_corner(
     $         p_model,
     $         t,nodes,x_map,y_map,i,j,
     $         side_x,side_y,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_R0)


          !test the corner timedev
          do k=1,ne
             do j=1, bc_size
                do i=1,bc_size
                   test_loc = is_test_validated(timedev(i,j,k),timedev_test(i,j,k),.false.)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print '(''['',3I2,'']: '',F10.4,'' -> '',F10.4)',
     $                     i,j,k,
     $                     timedev(i,j,k), timedev_test(i,j,k)
                   end if
                end do

                do i=nx-1,nx
                   test_loc = is_test_validated(timedev(i,j,k),timedev_test(i,j,k),.false.)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print '(''['',3I2,'']: '',F10.4,'' -> '',F10.4)',
     $                     i,j,k,
     $                     timedev(i,j,k), timedev_test(i,j,k)
                   end if
                end do
             end do

             do j=ny-1,ny
                do i=1,bc_size
                   test_loc = is_test_validated(timedev(i,j,k),timedev_test(i,j,k),.false.)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                     print '(''['',3I2,'']: '',F10.4,'' -> '',F10.4)',
     $                     i,j,k,
     $                     timedev(i,j,k), timedev_test(i,j,k)
                  end if
               end do
               
               do i=nx-1,nx
                  test_loc = is_test_validated(timedev(i,j,k),timedev_test(i,j,k),.false.)
                  test_validated = test_validated.and.test_loc
                  if(detailled.and.(.not.test_loc)) then
                     print '(''['',3I2,'']: '',F10.4,'' -> '',F10.4)',
     $                    i,j,k,
     $                    timedev(i,j,k), timedev_test(i,j,k)
                  end if
               end do
            end do
         end do                   
         
        end function test_apply_bc_on_timedev_xy_corner


        subroutine initialize_nodes(x_map,y_map,nodes,timedev_test)

          implicit none

          real(rkind), dimension(nx)      , intent(out) :: x_map
          real(rkind), dimension(ny)      , intent(out) :: y_map
          real(rkind), dimension(nx,ny,ne), intent(out) :: nodes
          real(rkind), dimension(nx,ny,ne), intent(out) :: timedev_test

          real(rkind)    :: dx,dy
          integer(ikind) :: i


          !space steps
          dx = 0.1d0
          dy = 0.2d0

          do i=1, nx
             x_map(i) = (i-1)*dx
          end do

          do i=1, ny
             y_map(i) = (i-1)*dy
          end do


          !SW corner
          nodes(1:3,1,1) = [ 1.2d0, 2.3d0, 8.9d0]
          nodes(1:3,2,1) = [-2.6d0,7.45d0,2.36d0]
          nodes(1:3,3,1) = [1.23d0,8.52d0,1.02d0]

          nodes(1:3,1,2) = [  0.5d0, -9.8d0, 5.63d0]
          nodes(1:3,2,2) = [ 1.23d0, 4.52d0, -2.5d0]
          nodes(1:3,3,2) = [ -8.9d0, 4.53d0, -7.1d0]

          nodes(1:3,1,3) = [  3.6d0, -1.2d0,-6.25d0]
          nodes(1:3,2,3) = [-7.56d0,-3.25d0,  8.9d0]
          nodes(1:3,3,3) = [-3.16d0,-4.15d0,-7.45d0]

          timedev_test(1,1,:) = [-25.76023178d0,-5.728046412d0,-5.728046412d0]
          timedev_test(2,1,:) = [-2.580423664d0, 7.277705062d0, 7.277705062d0]
          timedev_test(1,2,:) = [25.66465806d0 , 20.41237509d0, 20.41237509d0]
          timedev_test(2,2,:) = [13.72271983d0 , 7.748408409d0, 7.748408409d0]


          !SE corner
          nodes(4:6,1,1) = [ 8.9d0, 2.3d0, 1.2d0]
          nodes(4:6,2,1) = [2.36d0,7.45d0,-2.6d0]
          nodes(4:6,3,1) = [1.02d0,8.52d0,1.23d0]

          nodes(4:6,1,2) = [-5.63d0,  9.8d0, -0.5d0]
          nodes(4:6,2,2) = [  2.5d0,-4.52d0,-1.23d0]
          nodes(4:6,3,2) = [  7.1d0,-4.53d0,  8.9d0]

          nodes(4:6,1,3) = [-6.25d0, -1.2d0,  3.6d0]
          nodes(4:6,2,3) = [  8.9d0,-3.25d0,-7.56d0]
          nodes(4:6,3,3) = [-7.45d0,-4.15d0,-3.16d0]

          timedev_test(5,1,:) = [-2.580423664d0,-7.277705062d0, 7.277705062d0]
          timedev_test(6,1,:) = [-25.76023178d0, 5.728046412d0,-5.728046412d0]
          timedev_test(5,2,:) = [13.72271983d0 ,-7.748408409d0, 7.748408409d0]
          timedev_test(6,2,:) = [25.66465806d0 ,-20.41237509d0, 20.41237509d0]


          !NW corner
          nodes(1:3,4,1) = [1.23d0,8.52d0,1.02d0]
          nodes(1:3,5,1) = [-2.6d0,7.45d0,2.36d0]
          nodes(1:3,6,1) = [1.2d0,2.3d0,8.9d0]

          nodes(1:3,4,2) = [-8.9d0,4.53d0,-7.1d0]
          nodes(1:3,5,2) = [1.23d0,4.52d0,-2.5d0]
          nodes(1:3,6,2) = [0.5d0,-9.8d0,5.63d0]

          nodes(1:3,4,3) = [3.16d0,4.15d0,7.45d0]
          nodes(1:3,5,3) = [7.56d0,3.25d0,-8.9d0]
          nodes(1:3,6,3) = [-3.6d0,1.2d0,6.25d0]

          timedev_test(1,5,:) = [25.66465806d0 , 20.41237509d0,-20.41237509d0]
          timedev_test(2,5,:) = [13.72271983d0 , 7.748408409d0,-7.748408409d0]
          timedev_test(1,6,:) = [-25.76023178d0,-5.728046412d0, 5.728046412d0]
          timedev_test(2,6,:) = [-2.580423664d0, 7.277705062d0,-7.277705062d0]
 

          !NE corner
          nodes(4:6,4,1) = [1.02d0,8.52d0,1.23d0]
          nodes(4:6,5,1) = [2.36d0,7.45d0,-2.6d0]
          nodes(4:6,6,1) = [8.9d0 ,2.3d0 ,1.2d0 ]

          nodes(4:6,4,2) = [  7.1d0,-4.53d0, 8.9d0 ]
          nodes(4:6,5,2) = [  2.5d0,-4.52d0,-1.23d0]
          nodes(4:6,6,2) = [-5.63d0,  9.8d0,-0.5d0 ]

          nodes(4:6,4,3) = [7.45d0,4.15d0,3.16d0]
          nodes(4:6,5,3) = [-8.9d0,3.25d0,7.56d0]
          nodes(4:6,6,3) = [6.25d0,1.2d0 ,-3.6d0]

          timedev_test(5,5,:) = [13.72271983d0 ,-7.748408409d0,-7.748408409d0]
          timedev_test(6,5,:) = [25.66465806d0 ,-20.41237509d0,-20.41237509d0]
          timedev_test(5,6,:) = [-2.580423664d0,-7.277705062d0,-7.277705062d0]
          timedev_test(6,6,:) = [-25.76023178d0, 5.728046412d0, 5.728046412d0]

        end subroutine initialize_nodes


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
     $         nint(var*1e5)-
     $         nint(cst*1e5)).le.1
          
        end function is_test_validated


      end program test_bc_operators_hedstrom_xy_corners
