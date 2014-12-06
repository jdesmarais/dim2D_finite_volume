      !test file for module hedstrom_xy_corners_module.f
      program test_hedstrom_xy_corners

        use hedstrom_xy_corners_module, only :
     $       compute_n_timedev_with_openbc

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right

        use parameters_constant, only :
     $       n1_direction,
     $       n2_direction

        use parameters_input, only :
     $       nx,ny,ne
        
        use parameters_kind, only :
     $       ikind, rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_fd_n_module, only :
     $       gradient_n1_xL0_yI,
     $       gradient_n1_xI_yI,
     $       gradient_n1_xL0_yR0,
     $       gradient_n1_xI_yR0,
     $       gradient_n2_xL0_yI,
     $       gradient_n2_xI_yI,
     $       gradient_n2_xL0_yR0,
     $       gradient_n2_xI_yR0,
     $       
     $       gradient_n1_xR0_yI,
     $       gradient_n1_xI_yR0,
     $       gradient_n1_xR0_yR0,
     $       gradient_n2_xR0_yI,
     $       gradient_n2_xI_yR0,
     $       gradient_n2_xR0_yR0,
     $       
     $       gradient_n1_xI_yL0,
     $       gradient_n1_xR0_yL0,
     $       gradient_n1_xR0_yI,
     $       gradient_n2_xI_yL0,
     $       gradient_n2_xR0_yL0,
     $       gradient_n2_xR0_yI,
     $       
     $       gradient_n1_xL0_yL0,
     $       gradient_n1_xI_yL0,
     $       gradient_n1_xL0_yI,
     $       gradient_n2_xL0_yL0,
     $       gradient_n2_xI_yL0,
     $       gradient_n2_xL0_yI



        use wave2d_parameters, only :
     $       c


        implicit none

        logical :: test_loc
        logical :: test_validated
        logical :: detailled

        test_validated = .true.
        
        if(
     $       (nx.ne.5).or.
     $       (ny.ne.5).or.
     $       (ne.ne.3).or.
     $       (.not.is_test_validated(c,0.5d0,.false.))
     $       ) then
           print '(''nx.eq.5: '',L1)', nx.eq.5
           print '(''ny.eq.5: '',L1)', ny.eq.5
           print '(''ne.eq.3: '',L1)', ne.eq.3
           print '(''c=0.5  : '',L1)', is_test_validated(c,0.5d0,.false.)
           stop ''
        end if


        print '()'
        print '(''*************************'')'
        print '(''WARNING: use wave2d model'')'
        print '(''*************************'')'
        print '()'


        !test the computation of the NW corner
        detailled = .true.
        test_loc = test_compute_n_timedev_with_openbc_NW(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test compute_n_timedev_with_openbc_NW: '',L1)', test_validated
        print '()'


        !test the computation of the NE corner
        detailled = .true.
        test_loc = test_compute_n_timedev_with_openbc_NE(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test compute_n_timedev_with_openbc_NE: '',L1)', test_validated
        print '()'


        !test the computation of the SE corner
        detailled = .true.
        test_loc = test_compute_n_timedev_with_openbc_SE(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test compute_n_timedev_with_openbc_SE: '',L1)', test_validated
        print '()'


        !test the computation of the SW corner
        detailled = .true.
        test_loc = test_compute_n_timedev_with_openbc_SW(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test compute_n_timedev_with_openbc_SW: '',L1)', test_validated
        print '()'


        print '(''test_hedstrom_xy_corners: '',L1)', test_validated
        print '(''--------------------------------'')'
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
     $         nint(var*1e5)-
     $         nint(cst*1e5)).le.1
          
        end function is_test_validated


        function test_compute_n_timedev_with_openbc_NW(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind)   , dimension(nx,ny,ne) :: nodes
          real(rkind)                         :: dx
          real(rkind)                         :: dy
          real(rkind)   , dimension(ne,4)     :: timedev_test
          real(rkind)   , dimension(ne,4)     :: timedev_comp
          logical                             :: test_loc
          integer(ikind), dimension(2,4)      :: indices_tested
          integer                             :: k

          type(pmodel_eq) :: p_model


          test_validated = .true.


          !initialization of the nodes
          nodes(1:3,3,1) = [1.23d0,8.52d0,1.02d0]
          nodes(1:3,4,1) = [-2.6d0,7.45d0,2.36d0]
          nodes(1:3,5,1) = [1.2d0,2.3d0,8.9d0]

          nodes(1:3,3,2) = [-8.9d0,4.53d0,-7.1d0]
          nodes(1:3,4,2) = [1.23d0,4.52d0,-2.5d0]
          nodes(1:3,5,2) = [0.5d0,-9.8d0,5.63d0]

          nodes(1:3,3,3) = [3.16d0,4.15d0,7.45d0]
          nodes(1:3,4,3) = [7.56d0,3.25d0,-8.9d0]
          nodes(1:3,5,3) = [-3.6d0,1.2d0,6.25d0]

          dx = 0.1d0
          dy = 0.2d0


          !initialization of the test data
          timedev_test(:,1) = [25.66465806d0 , 20.41237509d0,-20.41237509d0]
          timedev_test(:,2) = [13.72271983d0 , 7.748408409d0,-7.748408409d0]
          timedev_test(:,3) = [-25.76023178d0,-5.728046412d0, 5.728046412d0]
          timedev_test(:,4) = [-2.580423664d0, 7.277705062d0,-7.277705062d0]
          
          indices_tested(:,1) = [1,4]
          indices_tested(:,2) = [2,4]
          indices_tested(:,3) = [1,5]
          indices_tested(:,4) = [2,5]
          

          !computation of the NW corner
          timedev_comp(:,1) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,1),
     $         indices_tested(2,1),
     $         p_model, dx, dy,
     $         gradient_n1_xL0_yI, gradient_n2_xL0_yI,
     $         incoming_left,
     $         n1_direction)

          timedev_comp(:,2) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,2),
     $         indices_tested(2,2),
     $         p_model, dx, dy,
     $         gradient_n1_xI_yI, gradient_n2_xI_yI,
     $         incoming_left,
     $         n1_direction)

          timedev_comp(:,3) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,3),
     $         indices_tested(2,3),
     $         p_model, dx, dy,
     $         gradient_n1_xL0_yR0, gradient_n2_xL0_yR0,
     $         incoming_left,
     $         n1_direction)

          timedev_comp(:,4) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,4),
     $         indices_tested(2,4),
     $         p_model, dx, dy,
     $         gradient_n1_xI_yR0, gradient_n2_xI_yR0,
     $         incoming_left,
     $         n1_direction)

          do k=1,4

             test_loc = compare_timedev(
     $            timedev_comp(:,k),
     $            timedev_test(:,k),
     $            detailled)

             test_validated =
     $            test_validated.and.test_loc

          end do

        end function test_compute_n_timedev_with_openbc_NW


        function test_compute_n_timedev_with_openbc_NE(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind)   , dimension(nx,ny,ne) :: nodes
          real(rkind)                         :: dx
          real(rkind)                         :: dy
          real(rkind)   , dimension(ne,4)     :: timedev_test
          real(rkind)   , dimension(ne,4)     :: timedev_comp
          logical                             :: test_loc
          integer(ikind), dimension(2,4)      :: indices_tested
          integer                             :: k

          type(pmodel_eq) :: p_model


          test_validated = .true.


          !initialization of the nodes
          nodes(3:5,3,1) = [1.02d0,8.52d0,1.23d0]
          nodes(3:5,4,1) = [2.36d0,7.45d0,-2.6d0]
          nodes(3:5,5,1) = [8.9d0 ,2.3d0 ,1.2d0 ]

          nodes(3:5,3,2) = [  7.1d0,-4.53d0, 8.9d0 ]
          nodes(3:5,4,2) = [  2.5d0,-4.52d0,-1.23d0]
          nodes(3:5,5,2) = [-5.63d0,  9.8d0,-0.5d0 ]

          nodes(3:5,3,3) = [7.45d0,4.15d0,3.16d0]
          nodes(3:5,4,3) = [-8.9d0,3.25d0,7.56d0]
          nodes(3:5,5,3) = [6.25d0,1.2d0 ,-3.6d0]

          dx = 0.1d0
          dy = 0.2d0


          !initialization of the test data
          timedev_test(:,1) = [13.72271983d0 ,-7.748408409d0,-7.748408409d0]
          timedev_test(:,2) = [25.66465806d0 ,-20.41237509d0,-20.41237509d0]
          timedev_test(:,3) = [-2.580423664d0,-7.277705062d0,-7.277705062d0]
          timedev_test(:,4) = [-25.76023178d0, 5.728046412d0, 5.728046412d0]

          
          indices_tested(:,1) = [4,4]
          indices_tested(:,2) = [5,4]
          indices_tested(:,3) = [4,5]
          indices_tested(:,4) = [5,5]
          

          !computation of the NE corner
          timedev_comp(:,1) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,1),
     $         indices_tested(2,1),
     $         p_model, dx, dy,
     $         gradient_n2_xI_yI, gradient_n1_xI_yI,
     $         incoming_right,
     $         n2_direction)

          timedev_comp(:,2) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,2),
     $         indices_tested(2,2),
     $         p_model, dx, dy,
     $         gradient_n2_xR0_yI, gradient_n1_xR0_yI,
     $         incoming_right,
     $         n2_direction)

          timedev_comp(:,3) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,3),
     $         indices_tested(2,3),
     $         p_model, dx, dy,
     $         gradient_n2_xI_yR0, gradient_n1_xI_yR0,
     $         incoming_right,
     $         n2_direction)

          timedev_comp(:,4) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,4),
     $         indices_tested(2,4),
     $         p_model, dx, dy,
     $         gradient_n2_xR0_yR0, gradient_n1_xR0_yR0,
     $         incoming_right,
     $         n2_direction)

          do k=1,4

             test_loc = compare_timedev(
     $            timedev_comp(:,k),
     $            timedev_test(:,k),
     $            detailled)

             test_validated =
     $            test_validated.and.test_loc

          end do

        end function test_compute_n_timedev_with_openbc_NE      


        function test_compute_n_timedev_with_openbc_SE(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind)   , dimension(nx,ny,ne) :: nodes
          real(rkind)                         :: dx
          real(rkind)                         :: dy
          real(rkind)   , dimension(ne,4)     :: timedev_test
          real(rkind)   , dimension(ne,4)     :: timedev_comp
          logical                             :: test_loc
          integer(ikind), dimension(2,4)      :: indices_tested
          integer                             :: k

          type(pmodel_eq) :: p_model


          test_validated = .true.


          !initialization of the nodes
          nodes(3:5,1,1) = [ 8.9d0, 2.3d0, 1.2d0]
          nodes(3:5,2,1) = [2.36d0,7.45d0,-2.6d0]
          nodes(3:5,3,1) = [1.02d0,8.52d0,1.23d0]

          nodes(3:5,1,2) = [-5.63d0,  9.8d0, -0.5d0]
          nodes(3:5,2,2) = [  2.5d0,-4.52d0,-1.23d0]
          nodes(3:5,3,2) = [  7.1d0,-4.53d0,  8.9d0]

          nodes(3:5,1,3) = [-6.25d0, -1.2d0,  3.6d0]
          nodes(3:5,2,3) = [  8.9d0,-3.25d0,-7.56d0]
          nodes(3:5,3,3) = [-7.45d0,-4.15d0,-3.16d0]

          dx = 0.1d0
          dy = 0.2d0


          !initialization of the test data
          timedev_test(:,1) = [-2.580423664d0,-7.277705062d0, 7.277705062d0]
          timedev_test(:,2) = [-25.76023178d0, 5.728046412d0,-5.728046412d0]
          timedev_test(:,3) = [13.72271983d0 ,-7.748408409d0, 7.748408409d0]
          timedev_test(:,4) = [25.66465806d0 ,-20.41237509d0, 20.41237509d0]

          
          indices_tested(:,1) = [4,1]
          indices_tested(:,2) = [5,1]
          indices_tested(:,3) = [4,2]
          indices_tested(:,4) = [5,2]
          

          !computation of the SE corner
          timedev_comp(:,1) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,1),
     $         indices_tested(2,1),
     $         p_model, dx, dy,
     $         gradient_n1_xI_yL0, gradient_n2_xI_yL0,
     $         incoming_right,
     $         n1_direction)

          timedev_comp(:,2) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,2),
     $         indices_tested(2,2),
     $         p_model, dx, dy,
     $         gradient_n1_xR0_yL0, gradient_n2_xR0_yL0,
     $         incoming_right,
     $         n1_direction)

          timedev_comp(:,3) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,3),
     $         indices_tested(2,3),
     $         p_model, dx, dy,
     $         gradient_n1_xI_yI, gradient_n2_xI_yI,
     $         incoming_right,
     $         n1_direction)

          timedev_comp(:,4) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,4),
     $         indices_tested(2,4),
     $         p_model, dx, dy,
     $         gradient_n1_xR0_yI, gradient_n2_xR0_yI,
     $         incoming_right,
     $         n1_direction)

          do k=1,4

             test_loc = compare_timedev(
     $            timedev_comp(:,k),
     $            timedev_test(:,k),
     $            detailled)

             test_validated =
     $            test_validated.and.test_loc

          end do

        end function test_compute_n_timedev_with_openbc_SE


        function test_compute_n_timedev_with_openbc_SW(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind)   , dimension(nx,ny,ne) :: nodes
          real(rkind)                         :: dx
          real(rkind)                         :: dy
          real(rkind)   , dimension(ne,4)     :: timedev_test
          real(rkind)   , dimension(ne,4)     :: timedev_comp
          logical                             :: test_loc
          integer(ikind), dimension(2,4)      :: indices_tested
          integer                             :: k

          type(pmodel_eq) :: p_model


          test_validated = .true.


          !initialization of the nodes
          nodes(1:3,1,1) = [ 1.2d0, 2.3d0, 8.9d0]
          nodes(1:3,2,1) = [-2.6d0,7.45d0,2.36d0]
          nodes(1:3,3,1) = [1.23d0,8.52d0,1.02d0]

          nodes(1:3,1,2) = [  0.5d0, -9.8d0, 5.63d0]
          nodes(1:3,2,2) = [ 1.23d0, 4.52d0, -2.5d0]
          nodes(1:3,3,2) = [ -8.9d0, 4.53d0, -7.1d0]

          nodes(1:3,1,3) = [  3.6d0, -1.2d0,-6.25d0]
          nodes(1:3,2,3) = [-7.56d0,-3.25d0,  8.9d0]
          nodes(1:3,3,3) = [-3.16d0,-4.15d0,-7.45d0]

          dx = 0.1d0
          dy = 0.2d0


          !initialization of the test data
          timedev_test(:,1) = [-25.76023178d0,-5.728046412d0,-5.728046412d0]
          timedev_test(:,2) = [-2.580423664d0, 7.277705062d0, 7.277705062d0]
          timedev_test(:,3) = [25.66465806d0 , 20.41237509d0, 20.41237509d0]
          timedev_test(:,4) = [13.72271983d0 , 7.748408409d0, 7.748408409d0]

          
          indices_tested(:,1) = [1,1]
          indices_tested(:,2) = [2,1]
          indices_tested(:,3) = [1,2]
          indices_tested(:,4) = [2,2]
          

          !computation of the SW corner
          timedev_comp(:,1) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,1),
     $         indices_tested(2,1),
     $         p_model, dx, dy,
     $         gradient_n2_xL0_yL0, gradient_n1_xL0_yL0,
     $         incoming_left,
     $         n2_direction)

          timedev_comp(:,2) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,2),
     $         indices_tested(2,2),
     $         p_model, dx, dy,
     $         gradient_n2_xI_yL0, gradient_n1_xI_yL0,
     $         incoming_left,
     $         n2_direction)

          timedev_comp(:,3) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,3),
     $         indices_tested(2,3),
     $         p_model, dx, dy,
     $         gradient_n2_xL0_yI, gradient_n1_xL0_yI,
     $         incoming_left,
     $         n2_direction)

          timedev_comp(:,4) = compute_n_timedev_with_openbc(
     $         nodes,
     $         indices_tested(1,4),
     $         indices_tested(2,4),
     $         p_model, dx, dy,
     $         gradient_n2_xI_yI, gradient_n1_xI_yI,
     $         incoming_left,
     $         n2_direction)

          do k=1,4

             test_loc = compare_timedev(
     $            timedev_comp(:,k),
     $            timedev_test(:,k),
     $            detailled)

             test_validated =
     $            test_validated.and.test_loc

          end do

        end function test_compute_n_timedev_with_openbc_SW


        function compare_timedev(timedev_comp,timedev_test,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(ne), intent(in) :: timedev_comp
          real(rkind), dimension(ne), intent(in) :: timedev_test
          logical                   , intent(in) :: detailled
          logical                                :: test_validated

          integer :: k
          logical :: test_loc


          test_validated = .true.

          do k=1, ne

             test_loc = is_test_validated(
     $            timedev_comp(k),
     $            timedev_test(k),
     $            .false.)

             if(detailled.and.(.not.test_loc)) then
                print '(''['',I2,'']: '',F8.3,''->'',F8.3)',
     $               k,
     $               timedev_comp(k),
     $               timedev_test(k)
             end if

             test_validated = test_validated.and.test_loc

          end do          

        end function compare_timedev
       
      end program test_hedstrom_xy_corners
