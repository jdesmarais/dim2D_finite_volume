      program test_bf_newgrdpt

        use bf_newgrdpt_class, only :
     $       bf_newgrdpt

        use parameters_constant, only :
     $       right,
     $       n2_direction

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_R0

        use sd_operators_fd_n_module, only :
     $       gradient_n1_xR0_yR1,
     $       gradient_n1_xR1_yR0,
     $       gradient_n1_xR0_yR0

        use wave2d_parameters, only :
     $       c

        implicit none


        type(bf_newgrdpt) :: bf_newgrdpt_used
        logical           :: test_validated
        logical           :: detailled


        !test requirements
        if(ne.ne.3) then
           stop 'the test requires ne=3'
        end if

        if(.not.is_test_validated(c,0.5d0,.false.)) then
           stop 'the test requires c=0.5'
        end if

        !test config
        detailled = .true.


        !test of get_interpolation_coeff_1D
        test_validated = test_get_interpolation_coeff_1D(bf_newgrdpt_used,detailled)
        print '(''test_get_interpolation_coeff_1D: '', L1)', test_validated


        !test of interpolate_1D
        test_validated = test_interpolate_1D(bf_newgrdpt_used,detailled)
        print '(''test_interpolate_1D: '', L1)', test_validated


        !test of compute_NewtonCotes_integration
        test_validated = test_compute_NewtonCotes_integration(bf_newgrdpt_used,detailled)
        print '(''test_compute_NewtonCotes_integration: '', L1)', test_validated


        !test of compute_newgrdpt_x
        test_validated = test_compute_newgrdpt_x(bf_newgrdpt_used,detailled)
        print '(''test_compute_newgrdpt_x: '', L1)', test_validated


        !test of compute_newgrdpt_y
        test_validated = test_compute_newgrdpt_y(bf_newgrdpt_used,detailled)
        print '(''test_compute_newgrdpt_y: '', L1)', test_validated


        !test of get_interpolation_coeff_2D
        test_validated = test_get_interpolation_coeff_2D(bf_newgrdpt_used,detailled)
        print '(''test_get_interpolation_coeff_2D: '', L1)', test_validated


        !test of interpolate_2D
        test_validated = test_interpolate_2D(bf_newgrdpt_used,detailled)
        print '(''test_interpolate_2D: '', L1)', test_validated


        !test of compute_newgrdpt_xy
        test_validated = test_compute_newgrdpt_xy(bf_newgrdpt_used,detailled)
        print '(''test_compute_newgrdpt_xy: '', L1)', test_validated


        contains

        function test_get_interpolation_coeff_1D(bf_newgrdpt_used,detailled)
     $       result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
          logical           , intent(in) :: detailled
          logical                        :: test_validated

          real(rkind), dimension(2)    :: x_map
          real(rkind), dimension(2,ne) :: nodes
          real(rkind), dimension(2,ne) :: nodes_inter
          real(rkind), dimension(2,ne) :: nodes_inter_data
          integer                      :: i,k
          logical                      :: test_loc

          test_validated = .true.

          x_map = [1.0d0,2.0d0]

          nodes = reshape( (/
     $         1.0d0,2.0d0,
     $         0.5d0,-0.5d0,
     $         0.25d0,-0.75d0
     $         /),
     $         (/2,ne/))

          nodes_inter_data = reshape( (/
     $          1.0d0, 0.0d0,
     $         -1.0d0, 1.5d0,
     $         -1.0d0, 1.25d0
     $         /),
     $         (/2,ne/))

          nodes_inter = bf_newgrdpt_used%get_interpolation_coeff_1D(
     $         x_map, nodes)

          do k=1, ne
             do i=1,2
                test_loc = is_test_validated(
     $               nodes_inter(i,k),
     $               nodes_inter_data(i,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_get_interpolation_coeff_1D


        function test_interpolate_1D(bf_newgrdpt_used,detailled)
     $       result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
          logical           , intent(in) :: detailled
          logical                        :: test_validated

          real(rkind)                  :: x
          real(rkind), dimension(2,ne) :: inter_coeff
          real(rkind), dimension(ne)   :: nodes_inter_data
          real(rkind), dimension(ne)   :: nodes_inter
          logical                      :: test_loc

          integer :: k
          
          test_validated = .true.

          x=1.0
          inter_coeff = reshape( (/
     $          1.0d0, 0.0d0,
     $         -1.0d0, 1.5d0,
     $         -1.0d0, 1.25d0
     $         /),
     $         (/2,ne/))
          
          nodes_inter_data = [1.0d0,0.5d0,0.25d0]

          nodes_inter = bf_newgrdpt_used%interpolate_1D(
     $         x,
     $         inter_coeff)

          do k=1, ne
             test_loc = is_test_validated(
     $            nodes_inter(k),
     $            nodes_inter_data(k),
     $            detailled)
             test_validated = test_validated.and.test_loc
          end do

        end function test_interpolate_1D


        function test_compute_NewtonCotes_integration(bf_newgrdpt_used,detailled)
     $       result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
          logical           , intent(in) :: detailled
          logical                        :: test_validated

          real(rkind), dimension(ne) :: data0
          real(rkind), dimension(ne) :: data1
          real(rkind)                :: dt
          real(rkind), dimension(ne) :: data_test
          real(rkind), dimension(ne) :: data_integrated
          integer                    :: k
          logical                    :: test_loc

          
          test_validated = .true.
                    
          data0 = [1.0d0 , 1.25d0, 0.5d0]
          data1 = [2.0d0, -0.25d0, 1.5d0]

          dt = 0.5d0

          data_test = [0.75d0, 0.25d0,0.50d0]

          data_integrated = bf_newgrdpt_used%compute_NewtonCotes_integration(
     $         data0,data1,dt)

          do k=1,ne

             test_loc = is_test_validated(
     $            data_test(k),
     $            data_integrated(k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

        end function test_compute_NewtonCotes_integration


        function test_compute_newgrdpt_x(bf_newgrdpt_used, detailled)
     $     result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
          logical           , intent(in) :: detailled
          logical                        :: test_validated

          type(pmodel_eq)                :: p_model
          real(rkind)                    :: t
          real(rkind)                    :: dt
          
          integer(ikind), dimension(2,2) :: bf_align0
          real(rkind), dimension(3)      :: bf_x_map0
          real(rkind), dimension(2)      :: bf_y_map0
          real(rkind), dimension(3,2,ne) :: bf_nodes0

          integer(ikind), dimension(2,2) :: bf_align1
          real(rkind), dimension(3)      :: bf_x_map1
          real(rkind), dimension(2)      :: bf_y_map1
          real(rkind), dimension(3,2,ne) :: bf_nodes1

          real(rkind), dimension(ne)     :: newgrdpt_data
          integer(ikind)                 :: i1
          integer(ikind)                 :: j1
          logical                        :: side_x

          integer                        :: k
          logical                        :: test_loc
          

          test_validated = .true.

          !initialization of the inputs
          t=0.0d0
          dt=0.25d0
          
          bf_align0(1,1) = 0
          bf_x_map0 = [0.5d0,1.5d0,2.5d0]
          bf_y_map0 = [0.25d0,0.5d0]
          bf_nodes0 = reshape((/
     $         1.0d0,2.0d0, 3.0d0,
     $         0.5d0,-0.5d0,1.25d0,
     $         0.25d0,-0.75d0,3.26d0,
     $         0.1d0 ,-0.45d0,6.15d0,
     $         2.05d0,-8.25d0,3.26d0,
     $         9.26d0, 7.85d0,9.23d0
     $         /),
     $         (/3,2,ne/))

          bf_align1(1,1) = 1
          bf_x_map1 = [1.5d0,2.5d0,3.5d0]
          bf_y_map1 = [0.25d0,0.5d0]
          bf_nodes1 = reshape((/
     $         1.0d0,2.45d0, 3.0d0,
     $         0.5d0,-0.26d0,2.25d0,
     $         0.25d0,-0.75d0,3.26d0,
     $         0.1d0 ,-8.52d0,7.15d0,
     $         2.05d0,-2.15d0,3.26d0,
     $         9.26d0, 7.85d0,6.23d0
     $         /),
     $         (/3,2,ne/))

          i1 = 3
          j1 = 2

          side_x = right

          !tested data
          newgrdpt_data = [-3.345117188d0,4.943867188d0,9.87d0]

          !test
          call bf_newgrdpt_used%compute_newgrdpt_x(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1, side_x, gradient_y_y_oneside_R0)


          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            bf_nodes1(i1,j1,k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

        end function test_compute_newgrdpt_x


        function test_compute_newgrdpt_y(bf_newgrdpt_used, detailled)
     $     result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
          logical           , intent(in) :: detailled
          logical                        :: test_validated

          type(pmodel_eq)                :: p_model
          real(rkind)                    :: t
          real(rkind)                    :: dt
          
          integer(ikind), dimension(2,2) :: bf_align0
          real(rkind), dimension(2)      :: bf_x_map0
          real(rkind), dimension(3)      :: bf_y_map0
          real(rkind), dimension(2,3,ne) :: bf_nodes0

          integer(ikind), dimension(2,2) :: bf_align1
          real(rkind), dimension(2)      :: bf_x_map1
          real(rkind), dimension(3)      :: bf_y_map1
          real(rkind), dimension(2,3,ne) :: bf_nodes1

          real(rkind), dimension(ne)     :: newgrdpt_data
          integer(ikind)                 :: i1
          integer(ikind)                 :: j1
          logical                        :: side_y

          integer                        :: k
          logical                        :: test_loc
          

          test_validated = .true.

          !initialization of the inputs
          t=0.0d0
          dt=0.25d0
          
          bf_align0(2,1) = 0
          bf_x_map0 = [0.25d0,0.5d0]
          bf_y_map0 = [0.5d0,1.5d0,2.5d0]
          bf_nodes0 = reshape((/
     $         1.0d0, 0.5d0,
     $         2.0d0,-0.5d0,
     $         3.0d0, 1.25d0,
     $         2.05d0,9.26d0,
     $         -8.25d0,7.85d0,
     $         3.26d0,9.23d0,
     $         0.25d0,0.1d0,
     $         -0.75d0,-0.45d0,
     $         3.26d0,6.15d0
     $         /),
     $         (/2,3,ne/))

          bf_align1(2,1) = 1
          bf_x_map1 = [0.25d0,0.5d0]
          bf_y_map1 = [1.5d0,2.5d0,3.5d0]
          bf_nodes1 = reshape((/
     $         1.0d0,0.5d0,
     $         2.45d0,-0.26d0,
     $         3.0d0,2.25d0,
     $         2.05d0,9.26d0,
     $         -2.15d0,7.85d0,
     $         3.26d0,6.23d0,
     $         0.25d0,0.1d0,
     $         -0.75d0,-8.52d0,
     $         3.26d0,7.15d0
     $         /),
     $         (/2,3,ne/))

          i1 = 2
          j1 = 3

          side_y = right

          !tested data
          newgrdpt_data = [-3.345117188d0,9.87d0,4.943867188d0]

          !test
          call bf_newgrdpt_used%compute_newgrdpt_y(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1, side_y, gradient_x_x_oneside_R0)


          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            bf_nodes1(i1,j1,k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

        end function test_compute_newgrdpt_y


        function test_compute_newgrdpt_xy(bf_newgrdpt_used, detailled)
     $     result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
          logical           , intent(in) :: detailled
          logical                        :: test_validated

          type(pmodel_eq)                :: p_model
          real(rkind)                    :: t
          real(rkind)                    :: dt
          
          integer(ikind), dimension(2,2) :: bf_align0
          real(rkind), dimension(3)      :: bf_x_map0
          real(rkind), dimension(3)      :: bf_y_map0
          real(rkind), dimension(3,3,ne) :: bf_nodes0

          integer(ikind), dimension(2,2) :: bf_align1
          real(rkind), dimension(4)      :: bf_x_map1
          real(rkind), dimension(4)      :: bf_y_map1
          real(rkind), dimension(4,4,ne) :: bf_nodes1

          real(rkind), dimension(ne)     :: newgrdpt_data
          integer(ikind)                 :: i1
          integer(ikind)                 :: j1
          integer                        :: n_direction
          logical                        :: side_n
          integer, dimension(2)          :: eigen_indices
          integer, dimension(2,3)        :: inter_indices1

          integer                        :: k
          logical                        :: test_loc
          

          test_validated = .true.

          !initialization of the inputs
          t=0.0d0
          dt=0.25d0
          
          bf_align0(1,1) = 0
          bf_align0(2,1) = 0
          bf_x_map0 = [0.5d0, 1.5d0 , 2.5d0]
          bf_y_map0 = [0.0d0, 0.25d0, 0.5d0]
          bf_nodes0 = reshape((/
     $          0.6d0,  0.2d0, 0.8d0,
     $          1.0d0,  2.0d0, 3.0d0,
     $          0.5d0, -0.5d0, 1.25d0,
     $        -3.25d0, 6.12d0, 7.15d0,
     $         0.25d0,-0.75d0, 3.26d0,
     $          0.1d0,-0.45d0, 6.15d0,
     $         9.26d0, 3.25d0, 4.15d0,
     $         2.05d0,-8.25d0, 3.26d0,
     $         9.26d0, 7.85d0, 9.23d0/),
     $         (/3,3,ne/))

          bf_align1(1,1) = 0
          bf_align0(2,1) = 0
          bf_x_map1 = [0.5d0, 1.5d0,2.5d0, 3.5d0]
          bf_y_map1 = [0.0d0,0.25d0,0.5d0,0.75d0]
          bf_nodes1 = reshape((/
     $         2.3d0,  7.8d0,   1.2d0, 1.5d0,
     $         8.9d0,  1.0d0,  2.45d0, 3.0d0,
     $         0.2d0,  0.5d0, -0.26d0,2.25d0,
     $        6.23d0,-5.15d0,  2.36d0, 0.0d0,
     $        -5.2d0, 1.23d0,  7.15d0, 6.2d0,
     $        9.26d0, 0.25d0, -0.75d0,3.26d0,
     $         2.3d0, 0.1d0, -8.52d0, 7.15d0,
     $         0.0d0, 0.0d0,   0.0d0,  0.0d0,
     $        9.63d0, 1.2d0,  7.32d0, 1.52d0,
     $        1.25d0, 2.05d0,-2.15d0, 3.26d0,
     $        7.26d0, 9.26d0, 7.85d0, 6.23d0,
     $         0.0d0,  0.0d0,  0.0d0,  0.0d0/),
     $         (/4,4,ne/))

          i1 = 4
          j1 = 4

          n_direction         = n2_direction
          side_n              = right
          eigen_indices       = [3,3]
          inter_indices1(:,1) = [3,2]
          inter_indices1(:,2) = [2,3]
          inter_indices1(:,3) = [3,3]

          !tested data
          newgrdpt_data = [-11.06536693d0,7.194275495d0,8.134275495d0]

          !test
          call bf_newgrdpt_used%compute_newgrdpt_xy(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         n_direction,
     $         side_n,
     $         gradient_n1_xR0_yR1,
     $         gradient_n1_xR1_yR0,
     $         gradient_n1_xR0_yR0,
     $         eigen_indices,
     $         inter_indices1)

          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            bf_nodes1(i1,j1,k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

        end function test_compute_newgrdpt_xy


        function test_get_interpolation_coeff_2D(bf_newgrdpt_used,detailled)
     $       result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
          logical           , intent(in) :: detailled
          logical                        :: test_validated

          real(rkind), dimension(3)    :: x_map
          real(rkind), dimension(3)    :: y_map
          real(rkind), dimension(3,ne) :: nodes
          real(rkind), dimension(3,ne) :: nodes_inter
          real(rkind), dimension(3,ne) :: nodes_inter_data
          integer                      :: i,k
          logical                      :: test_loc

          test_validated = .true.

          x_map = [1.0d0,1.0d0,2.0d0]
          y_map = [1.0d0,2.0d0,2.0d0]

          nodes = reshape( (/
     $         1.0d0,1.0d0,1.0d0,
     $         1.0d0,1.0d0,0.0d0,
     $         0.0d0,1.0d0,1.0d0
     $         /),
     $         (/3,ne/))

          nodes_inter_data = reshape( (/
     $           0.0d0, 0.0d0, 1.0d0,
     $          -1.0d0, 0.0d0, 2.0d0,
     $           0.0d0, 1.0d0,-1.0d0
     $         /),
     $         (/3,ne/))

          nodes_inter = bf_newgrdpt_used%get_interpolation_coeff_2D(
     $         x_map, y_map,nodes)

          do k=1, ne
             do i=1,3
                test_loc = is_test_validated(
     $               nodes_inter(i,k),
     $               nodes_inter_data(i,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_get_interpolation_coeff_2D


        function test_interpolate_2D(bf_newgrdpt_used,detailled)
     $       result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
          logical           , intent(in) :: detailled
          logical                        :: test_validated

          real(rkind)                  :: x
          real(rkind)                  :: y
          real(rkind), dimension(3,ne) :: inter_coeff
          real(rkind), dimension(ne)   :: nodes_inter_data
          real(rkind), dimension(ne)   :: nodes_inter
          logical                      :: test_loc

          integer :: k
          
          test_validated = .true.

          x=0.25d0
          y=0.5d0
          inter_coeff = reshape( (/
     $          0.0d0, 0.0d0, 1.0d0,
     $         -1.0d0, 0.0d0, 1.0d0,
     $          0.0d0,-1.0d0, 1.0d0
     $         /),
     $         (/3,ne/))
          
          nodes_inter_data = [1.0d0,0.75d0,0.5d0]

          nodes_inter = bf_newgrdpt_used%interpolate_2D(
     $         x,y,
     $         inter_coeff)

          do k=1, ne
             test_loc = is_test_validated(
     $            nodes_inter(k),
     $            nodes_inter_data(k),
     $            detailled)
             test_validated = test_validated.and.test_loc
          end do

        end function test_interpolate_2D

        
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
     $         int(var*1e5)-
     $         int(cst*1e5)).le.1
          
        end function is_test_validated

      end program test_bf_newgrdpt
