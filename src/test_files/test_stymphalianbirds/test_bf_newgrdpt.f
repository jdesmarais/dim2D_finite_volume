      program test_bf_newgrdpt

        use bf_compute_class, only :
     $       bf_compute

        use bf_layer_bc_procedure_module, only : 
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type

        use bf_layer_class, only :
     $       bf_layer

        use bf_layer_newgrdpt_procedure_module, only :
     $       no_gradient_type,
     $       gradient_R0_type

        use bf_newgrdpt_class, only :
     $       bf_newgrdpt

        use parameters_constant, only :
     $       left,right,
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


        !test of compute_newgrdpt
        test_validated = test_compute_newgrdpt(bf_newgrdpt_used,detailled)
        print '(''test_compute_newgrdpt: '', L1)', test_validated


        !test of bf_compute/compute_newgrdpt
        test_validated = test_bf_compute_compute_newgrdpt(detailled)
        print '(''test_bf_compute_compute_newgrdpt: '',L1)', test_validated


        !test of bf_layer/compute_newgrdpt
        test_validated = test_bf_layer_compute_newgrdpt(detailled)
        print '(''test_bf_layer_compute_newgrdpt: '',L1)', test_validated

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
          real(rkind), dimension(ne)     :: newgrdpt
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
          newgrdpt = bf_newgrdpt_used%compute_newgrdpt_x(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1, side_x, gradient_y_y_oneside_R0)


          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            newgrdpt(k),
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
          real(rkind), dimension(ne)     :: newgrdpt
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
          newgrdpt = bf_newgrdpt_used%compute_newgrdpt_y(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1, side_y, gradient_x_x_oneside_R0)


          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            newgrdpt(k),
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
          real(rkind), dimension(ne)     :: newgrdpt
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
          newgrdpt = bf_newgrdpt_used%compute_newgrdpt_xy(
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
     $            newgrdpt(k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

        end function test_compute_newgrdpt_xy


        function test_compute_newgrdpt(bf_newgrdpt_used, detailled)
     $     result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
          logical           , intent(in) :: detailled
          logical                        :: test_validated

          type(pmodel_eq)                :: p_model
          real(rkind)                    :: t
          real(rkind)                    :: dt
          
          integer(ikind), dimension(2,2)             :: bf_align0
          real(rkind), dimension(:)    , allocatable :: bf_x_map0
          real(rkind), dimension(:)    , allocatable :: bf_y_map0
          real(rkind), dimension(:,:,:), allocatable :: bf_nodes0

          integer(ikind), dimension(2,2)             :: bf_align1
          real(rkind), dimension(:)    , allocatable :: bf_x_map1
          real(rkind), dimension(:)    , allocatable :: bf_y_map1
          real(rkind), dimension(:,:,:), allocatable :: bf_nodes1

          real(rkind), dimension(ne)     :: newgrdpt_data
          real(rkind), dimension(ne)     :: newgrdpt
          integer(ikind)                 :: i1
          integer(ikind)                 :: j1

          integer                        :: k
          logical                        :: test_loc
          logical                        :: test_NE_corner
          logical                        :: test_NW_corner
          logical                        :: test_SE_corner
          logical                        :: test_SW_corner
          logical                        :: test_E_edge
          logical                        :: test_W_edge
          logical                        :: test_N_edge
          logical                        :: test_S_edge
          

          allocate(bf_x_map0(3))
          allocate(bf_y_map0(3))
          allocate(bf_nodes0(3,3,ne))
          
          allocate(bf_x_map1(4))
          allocate(bf_y_map1(4))
          allocate(bf_nodes1(4,4,ne))

          !------------------------------------------------
          !test of the corners
          !------------------------------------------------

          test_validated = .true.

          !test of the NE corner
          !------------------------------------------------
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

          !tested data
          newgrdpt_data = [-11.06536693d0,7.194275495d0,8.134275495d0]

          !test
          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
     $         p_model, t, dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         NE_corner_type,no_gradient_type)

          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            newgrdpt(k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

          test_NE_corner = test_validated
          if(detailled) then
             print '(''test_NE_corner: '',L1)', test_NE_corner
          end if

          
          test_validated = .true.

          !test of the NW corner
          !------------------------------------------------
          !initialization of the inputs
          t=0.0d0
          dt=0.25d0
          
          bf_align0(1,1) = 1
          bf_align0(2,1) = 0
          bf_x_map0 = [ 0.5d0, 1.5d0 , 2.5d0]
          bf_y_map0 = [ 0.0d0, 0.25d0, 0.5d0]
          call reflect_x(bf_nodes0)

          bf_align1(1,1) = 0
          bf_align0(2,1) = 0
          bf_x_map1 = [-0.5d0, 0.5d0, 1.5d0, 2.5d0]
          bf_y_map1 = [ 0.0d0,0.25d0, 0.5d0,0.75d0]
          call reflect_x(bf_nodes1)

          i1 = 1
          j1 = 4

          !tested data
          newgrdpt_data = [-11.06536693d0,-7.194275495d0,8.134275495d0]

          !test
          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         NW_corner_type,no_gradient_type)

          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            newgrdpt(k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

          test_NW_corner = test_validated
          if(detailled) then
             print '(''test_NW_corner: '',L1)', test_NW_corner
          end if


          test_validated = .true.

          !test of the SW corner
          !------------------------------------------------
          !initialization of the inputs
          t=0.0d0
          dt=0.25d0
          
          bf_align0(1,1) = 1
          bf_align0(2,1) = 1
          bf_x_map0 = [ 0.5d0, 1.5d0 , 2.5d0]
          bf_y_map0 = [ 0.0d0, 0.25d0, 0.5d0]
          call reflect_y(bf_nodes0)

          bf_align1(1,1) = 0
          bf_align1(2,1) = 0
          bf_x_map1 = [ -0.5d0, 0.5d0, 1.5d0, 2.5d0]
          bf_y_map1 = [-0.25d0, 0.0d0,0.25d0, 0.5d0]
          call reflect_y(bf_nodes1)

          i1 = 1
          j1 = 1

          !tested data
          newgrdpt_data = [-11.06536693d0,-7.194275495d0,-8.134275495d0]

          !test
          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         SW_corner_type,no_gradient_type)

          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            newgrdpt(k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

          test_SW_corner = test_validated
          if(detailled) then
             print '(''test_SW_corner: '',L1)', test_SW_corner
          end if


           test_validated = .true.

          !test of the SE corner
          !------------------------------------------------
          !initialization of the inputs
          t=0.0d0
          dt=0.25d0
          
          bf_align0(1,1) = 0
          bf_align0(2,1) = 1
          bf_x_map0 = [ 0.5d0, 1.5d0 , 2.5d0]
          bf_y_map0 = [ 0.0d0, 0.25d0, 0.5d0]
          call reflect_x(bf_nodes0)

          bf_align1(1,1) = 0
          bf_align1(2,1) = 0
          bf_x_map1 = [  0.5d0, 1.5d0,  2.5d0, 3.5d0]
          bf_y_map1 = [-0.25d0, 0.0d0, 0.25d0, 0.5d0]
          call reflect_x(bf_nodes1)

          i1 = 4
          j1 = 1

          !tested data
          newgrdpt_data = [-11.06536693d0, 7.194275495d0,-8.134275495d0]

          !test
          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         SE_corner_type,no_gradient_type)

          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            newgrdpt(k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

          test_SE_corner = test_validated
          if(detailled) then
             print '(''test_SE_corner: '',L1)', test_SE_corner
          end if

          deallocate(bf_x_map0)
          deallocate(bf_y_map0)
          deallocate(bf_nodes0)
          
          deallocate(bf_x_map1)
          deallocate(bf_y_map1)
          deallocate(bf_nodes1)


          !------------------------------------------------
          !test of the x_edge
          !------------------------------------------------

          allocate(bf_x_map0(3))
          allocate(bf_y_map0(2))
          allocate(bf_nodes0(3,2,ne))

          allocate(bf_x_map1(3))
          allocate(bf_y_map1(2))
          allocate(bf_nodes1(3,2,ne))

          test_validated = .true.

          !test E_edge
          t=0.0d0
          dt=0.25d0
          
          bf_align0(1,1) = 0
          bf_align0(2,1) = 0
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
          bf_align1(2,1) = 0
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

          !tested data
          newgrdpt_data = [-3.345117188d0,4.943867188d0,9.87d0]

          !test
          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         E_edge_type,gradient_R0_type)


          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            newgrdpt(k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

          test_E_edge = test_validated

          if(detailled) then
             print '(''test_E_edge: '',L1)', test_E_edge
          end if

           test_validated = .true.

          !test W_edge
          t=0.0d0
          dt=0.25d0
          
          bf_align0(1,1) = 1
          bf_align0(2,1) = 0
          bf_x_map0 = [0.5d0,1.5d0,2.5d0]
          bf_y_map0 = [0.25d0,0.5d0]
          call reflect_x(bf_nodes0)

          bf_align1(1,1) = 0
          bf_align1(2,1) = 0
          bf_x_map1 = [-0.5d0,0.5d0,1.5d0]
          bf_y_map1 = [0.25d0,0.5d0]
          call reflect_x(bf_nodes1)

          i1 = 1
          j1 = 2

          !tested data
          newgrdpt_data = [-3.345117188d0,-4.943867188d0,9.87d0]

          !test
          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         W_edge_type,gradient_R0_type)


          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            newgrdpt(k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

          test_W_edge = test_validated

          if(detailled) then
             print '(''test_W_edge: '',L1)', test_W_edge
          end if

          deallocate(bf_x_map0)
          deallocate(bf_y_map0)
          deallocate(bf_nodes0)
          
          deallocate(bf_x_map1)
          deallocate(bf_y_map1)
          deallocate(bf_nodes1)


          !------------------------------------------------
          !test of the y_edge
          !------------------------------------------------

          allocate(bf_x_map0(2))
          allocate(bf_y_map0(3))
          allocate(bf_nodes0(2,3,ne))

          allocate(bf_x_map1(2))
          allocate(bf_y_map1(3))
          allocate(bf_nodes1(2,3,ne))

          test_validated = .true.

          !test N_edge
          t=0.0d0
          dt=0.25d0
          
          bf_align0(1,1) = 0
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

          bf_align1(1,1) = 0
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

          !tested data
          newgrdpt_data = [-3.345117188d0,9.87d0,4.943867188d0]

          !test
          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         N_edge_type,
     $         gradient_R0_type)


          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            newgrdpt(k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

          test_N_edge = test_validated

          if(detailled) then
             print '(''test_N_edge: '',L1)', test_N_edge
          end if

          
          !test S_edge
          t=0.0d0
          dt=0.25d0
          
          bf_align0(1,1) = 0
          bf_align0(2,1) = 1
          bf_x_map0 = [0.25d0,0.5d0]
          bf_y_map0 = [0.5d0,1.5d0,2.5d0]
          call reflect_y(bf_nodes0)

          bf_align1(1,1) = 0
          bf_align1(2,1) = 0
          bf_x_map1 = [0.25d0,0.5d0]
          bf_y_map1 = [-0.5d0,0.5d0,1.5d0]
          call reflect_y(bf_nodes1)

          i1 = 2
          j1 = 1

          !tested data
          newgrdpt_data = [-3.345117188d0,9.87d0,-4.943867188d0]

          !test
          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         S_edge_type,
     $         gradient_R0_type)


          !comparison
          do k=1,ne

             test_loc = is_test_validated(
     $            newgrdpt_data(k),
     $            newgrdpt(k),
     $            detailled)
             test_validated = test_validated.and.test_loc

          end do

          test_S_edge = test_validated

          if(detailled) then
             print '(''test_S_edge: '',L1)', test_S_edge
          end if

          deallocate(bf_x_map0)
          deallocate(bf_y_map0)
          deallocate(bf_nodes0)
          
          deallocate(bf_x_map1)
          deallocate(bf_y_map1)
          deallocate(bf_nodes1)


          test_validated =
     $         test_NE_corner.and.
     $         test_NW_corner.and.
     $         test_SW_corner.and.
     $         test_SE_corner.and.
     $         test_E_edge.and.
     $         test_W_edge.and.
     $         test_N_edge.and.
     $         test_S_edge

        end function test_compute_newgrdpt


        function test_bf_compute_compute_newgrdpt(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer    , dimension(:,:)  , allocatable :: align0
          integer    , dimension(:,:)  , allocatable :: grdpts_id0
          real(rkind), dimension(:)    , allocatable :: x_map0
          real(rkind), dimension(:)    , allocatable :: y_map0
          real(rkind), dimension(:,:,:), allocatable :: nodes0

          integer    , dimension(:,:)  , allocatable :: align1
          integer    , dimension(:,:)  , allocatable :: grdpts_id1
          real(rkind), dimension(:)    , allocatable :: x_map1
          real(rkind), dimension(:)    , allocatable :: y_map1
          real(rkind), dimension(:,:,:), allocatable :: nodes1

          integer(ikind)                             :: i1,j1
          integer                                    :: k,l
          logical                                    :: test_loc

          integer    , dimension(:,:)  , allocatable :: compute_index
          real(rkind), dimension(:,:)  , allocatable :: newgrdpt_data
          real(rkind), dimension(ne)                 :: newgrdpt

          type(bf_compute) :: bf_compute_used
          type(pmodel_eq)  :: p_model

          real(rkind) :: t
          real(rkind) :: dt

          !t and dt initialization
          t  = 0.0d0
          dt = 0.25d0


          !allocations
          allocate(align0(2,2))
          allocate(grdpts_id0(14,12))
          allocate(x_map0(14))
          allocate(y_map0(12))
          allocate(nodes0(14,12,ne))

          allocate(align1(2,2))
          allocate(grdpts_id1(16,14))
          allocate(x_map1(16))
          allocate(y_map1(14))
          allocate(nodes1(16,14,ne))


          !initializations

          !initialization of the alignments
          !--------------------------------------
          align0(1,1) = 4
          align0(1,2) = 13
          align0(2,1) = 3
          align0(2,2) = 10

          align1(1,1) = 3
          align1(1,2) = 14
          align1(2,1) = 2
          align1(2,2) = 11
          

          !initialization of the grdpts
          !--------------------------------------
          grdpts_id0 = reshape((/
     $         0,0,0,0,3,3,3,3,3,0,0,0,0,0,
     $         3,3,3,3,3,2,2,2,3,3,3,3,3,3,
     $         3,2,2,2,2,2,1,2,2,2,2,2,2,3,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,2,2,1,1,1,1,1,1,2,2,2,3,
     $         3,3,3,2,1,1,1,1,1,1,2,3,3,3,
     $         0,0,3,2,1,1,1,1,1,1,2,3,0,0,
     $         0,0,3,2,2,2,1,2,2,2,2,3,0,0,
     $         0,0,3,3,3,2,1,2,3,3,3,3,0,0,
     $         0,0,0,0,3,2,2,2,3,0,0,0,0,0,
     $         0,0,0,0,3,3,3,3,3,0,0,0,0,0/),
     $         (/14,12/))

          grdpts_id1 = reshape((/
     $         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     $         0,0,0,0,0,3,3,3,3,3,0,0,0,0,0,0,
     $         0,3,3,3,3,3,2,2,2,3,3,3,3,3,3,0,
     $         0,3,2,2,2,2,2,1,2,2,2,2,2,2,3,0,
     $         0,3,2,1,1,1,1,1,1,1,1,1,1,2,3,0,
     $         0,3,2,1,1,1,1,1,1,1,1,1,1,2,3,0,
     $         0,3,2,2,2,1,1,1,1,1,1,2,2,2,3,0,
     $         0,3,3,3,2,1,1,1,1,1,1,2,3,3,3,0,
     $         0,0,0,3,2,1,1,1,1,1,1,2,3,0,0,0,
     $         0,0,0,3,2,2,2,1,2,2,2,2,3,0,0,0,
     $         0,0,0,3,3,3,2,1,2,3,3,3,3,0,0,0,
     $         0,0,0,0,0,3,2,2,2,3,0,0,0,0,0,0,
     $         0,0,0,0,0,3,3,3,3,3,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/),
     $         (/16,14/))


          !initialization of the maps
          !--------------------------------------
          do k=1, size(x_map0,1)
             x_map0(k) = 1.0d0 + 0.25d0*(k-12)
          end do

          do k=1, size(y_map0,1)
             y_map0(k) = 0.25d0 + 0.25d0*(k-9) !1.0
          end do

          do k=1, size(x_map1,1)
             x_map1(k) = 0.75d0 + 0.25d0*(k-12)
          end do

          do k=1, size(y_map1,1)
             y_map1(k) = 0.0d0 + 0.25d0*(k-9) !1.0
          end do


          !initialization of the nodes
          !--------------------------------------
          !S
          nodes0(8:9,1:2,:) = reshape((/
     $         3.0d0,  1.25d0,
     $         2.0d0,  -0.5d0,
     $         3.26d0, 9.23d0,
     $        -8.25d0, 7.85d0,
     $        -3.26d0,-6.15d0,
     $         0.75d0, 0.45d0/),
     $         (/2,2,3/))

          nodes1(9:10,2:3,:) = reshape((/
     $         2.45d0, -0.26d0,
     $          1.0d0,   0.5d0,
     $        -2.15d0,  7.85d0,
     $         2.05d0,  9.26d0,
     $         0.75d0,  8.52d0,
     $        -0.25d0,  -0.1d0/),
     $         (/2,2,3/))

          !SW_1
          nodes0(1:3,2:4,:) = reshape((/
     $          1.25d0,  -0.5d0,   0.5d0,
     $           3.0d0,   2.0d0,   1.0d0,
     $           0.8d0,   0.2d0,   0.6d0,
     $         -6.15d0,  0.45d0,  -0.1d0,
     $         -3.26d0,  0.75d0, -0.25d0,
     $         -7.15d0, -6.12d0,  3.25d0,
     $         -9.23d0, -7.85d0, -9.26d0,
     $         -3.26d0,  8.25d0, -2.05d0,
     $         -4.15d0, -3.25d0, -9.26d0/),
     $         (/3,3,3/))

          nodes1(2:4,3:5,:) = reshape((/
     $        -0.26d0,   0.5d0,   0.2d0,
     $         2.45d0,   1.0d0,   8.9d0,
     $          1.2d0,   7.8d0,   2.3d0,
     $         8.52d0,  -0.1d0,  -2.3d0,
     $         0.75d0, -0.25d0, -9.26d0,
     $        -7.15d0, -1.23d0,   5.2d0,
     $        -7.85d0, -9.26d0, -7.26d0,
     $         2.15d0, -2.05d0, -1.25d0,
     $        -7.32d0,  -1.2d0, -9.63d0/),
     $         (/3,3,3/))

          !NE_1
          nodes0(10:12,8:10,:) = reshape((/
     $          0.6d0,  0.2d0, 0.8d0,
     $          1.0d0,  2.0d0, 3.0d0,
     $          0.5d0, -0.5d0,1.25d0,
     $        -3.25d0, 6.12d0,7.15d0,
     $         0.25d0,-0.75d0,3.26d0,
     $          0.1d0,-0.45d0,6.15d0,
     $         9.26d0, 3.25d0,4.15d0,
     $         2.05d0,-8.25d0,3.26d0,
     $         9.26d0, 7.85d0,9.23d0/),
     $         (/3,3,3/))

          nodes1(11:13,9:11,:) = reshape((/
     $         2.3d0,  7.8d0,  1.2d0,
     $         8.9d0,  1.0d0, 2.45d0,
     $         0.2d0,  0.5d0,-0.26d0,
     $        -5.2d0, 1.23d0, 7.15d0,
     $        9.26d0, 0.25d0,-0.75d0,
     $         2.3d0,  0.1d0,-8.52d0,
     $        9.63d0,  1.2d0, 7.32d0,
     $        1.25d0, 2.05d0,-2.15d0,
     $        7.26d0, 9.26d0, 7.85d0/),
     $         (/3,3,3/))


          !set the nodes in bf_compute object
          call bf_compute_used%set_alignment(align0)
          call bf_compute_used%set_grdpts_id(grdpts_id0)
          call bf_compute_used%set_x_map(x_map0)
          call bf_compute_used%set_y_map(y_map0)
          call bf_compute_used%set_nodes(nodes0)

          
          !test the computation of the new grid point
          allocate(compute_index(2,3))
          allocate(newgrdpt_data(3,3))

          compute_index(:,1) = [10,1]
          newgrdpt_data(:,1) = [-2.77171875d0,9.87d0,-4.37046875d0]

          compute_index(:,2) = [1,2]
          newgrdpt_data(:,2) = [-11.75157856d0,-6.859983049d0,-7.799983049d0]
          

          compute_index(:,3) = [14,12]
          newgrdpt_data(:,3) = [-11.75157856d0, 6.859983049d0, 7.799983049d0]


          !test the computation of the new grid points
          do k=1, size(compute_index,2)
             
             !indices of the new grid point
             i1 = compute_index(1,k)
             j1 = compute_index(2,k)

             !get the new grid point
             newgrdpt = bf_compute_used%compute_newgrdpt(
     $            p_model,t,dt,
     $            align1, x_map1, y_map1, nodes1,
     $            i1, j1)

             !compare with the data
             test_validated = .true.
             do l=1,ne
                test_loc = is_test_validated(
     $               newgrdpt(l),
     $               newgrdpt_data(l,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
             end do

             !print the comparison results
             if(detailled) then
                if(test_loc) then
                   print '(''test '',I2, '' validated'')', k
                else
                   print '(''**test '',I2, '' failed**'')', k
                end if
             end if

          end do

        end function test_bf_compute_compute_newgrdpt


        function test_bf_layer_compute_newgrdpt(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer    , dimension(:,:)  , allocatable :: align0
          integer    , dimension(:,:)  , allocatable :: grdpts_id0
          real(rkind), dimension(:)    , allocatable :: x_map0
          real(rkind), dimension(:)    , allocatable :: y_map0
          real(rkind), dimension(:,:,:), allocatable :: nodes0

          integer    , dimension(:,:)  , allocatable :: align1
          integer    , dimension(:,:)  , allocatable :: grdpts_id1
          real(rkind), dimension(:)    , allocatable :: x_map1
          real(rkind), dimension(:)    , allocatable :: y_map1
          real(rkind), dimension(:,:,:), allocatable :: nodes1

          integer(ikind)                             :: i1,j1
          integer                                    :: k,l
          logical                                    :: test_loc

          integer    , dimension(:,:)  , allocatable :: compute_index
          real(rkind), dimension(:,:)  , allocatable :: newgrdpt_data
          real(rkind), dimension(ne)                 :: newgrdpt

          type(bf_layer)  :: bf_layer_used
          type(pmodel_eq) :: p_model

          real(rkind) :: t
          real(rkind) :: dt

          !t and dt initialization
          t  = 0.0d0
          dt = 0.25d0


          !allocations
          allocate(align0(2,2))
          allocate(grdpts_id0(14,12))
          allocate(x_map0(14))
          allocate(y_map0(12))
          allocate(nodes0(14,12,ne))

          allocate(align1(2,2))
          allocate(grdpts_id1(16,14))
          allocate(x_map1(16))
          allocate(y_map1(14))
          allocate(nodes1(16,14,ne))


          !initializations

          !initialization of the alignments
          !--------------------------------------
          align0(1,1) = 4
          align0(1,2) = 13
          align0(2,1) = 3
          align0(2,2) = 10

          align1(1,1) = 3
          align1(1,2) = 14
          align1(2,1) = 2
          align1(2,2) = 11
          

          !initialization of the grdpts
          !--------------------------------------
          grdpts_id0 = reshape((/
     $         0,0,0,0,3,3,3,3,3,0,0,0,0,0,
     $         3,3,3,3,3,2,2,2,3,3,3,3,3,3,
     $         3,2,2,2,2,2,1,2,2,2,2,2,2,3,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
     $         3,2,2,2,1,1,1,1,1,1,2,2,2,3,
     $         3,3,3,2,1,1,1,1,1,1,2,3,3,3,
     $         0,0,3,2,1,1,1,1,1,1,2,3,0,0,
     $         0,0,3,2,2,2,1,2,2,2,2,3,0,0,
     $         0,0,3,3,3,2,1,2,3,3,3,3,0,0,
     $         0,0,0,0,3,2,2,2,3,0,0,0,0,0,
     $         0,0,0,0,3,3,3,3,3,0,0,0,0,0/),
     $         (/14,12/))

          grdpts_id1 = reshape((/
     $         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
     $         0,0,0,0,0,3,3,3,3,3,0,0,0,0,0,0,
     $         0,3,3,3,3,3,2,2,2,3,3,3,3,3,3,0,
     $         0,3,2,2,2,2,2,1,2,2,2,2,2,2,3,0,
     $         0,3,2,1,1,1,1,1,1,1,1,1,1,2,3,0,
     $         0,3,2,1,1,1,1,1,1,1,1,1,1,2,3,0,
     $         0,3,2,2,2,1,1,1,1,1,1,2,2,2,3,0,
     $         0,3,3,3,2,1,1,1,1,1,1,2,3,3,3,0,
     $         0,0,0,3,2,1,1,1,1,1,1,2,3,0,0,0,
     $         0,0,0,3,2,2,2,1,2,2,2,2,3,0,0,0,
     $         0,0,0,3,3,3,2,1,2,3,3,3,3,0,0,0,
     $         0,0,0,0,0,3,2,2,2,3,0,0,0,0,0,0,
     $         0,0,0,0,0,3,3,3,3,3,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/),
     $         (/16,14/))


          !initialization of the maps
          !--------------------------------------
          do k=1, size(x_map0,1)
             x_map0(k) = 1.0d0 + 0.25d0*(k-12)
          end do

          do k=1, size(y_map0,1)
             y_map0(k) = 0.25d0 + 0.25d0*(k-9) !1.0
          end do

          do k=1, size(x_map1,1)
             x_map1(k) = 0.75d0 + 0.25d0*(k-12)
          end do

          do k=1, size(y_map1,1)
             y_map1(k) = 0.0d0 + 0.25d0*(k-9) !1.0
          end do


          !initialization of the nodes
          !--------------------------------------
          !S
          nodes0(8:9,1:2,:) = reshape((/
     $         3.0d0,  1.25d0,
     $         2.0d0,  -0.5d0,
     $         3.26d0, 9.23d0,
     $        -8.25d0, 7.85d0,
     $        -3.26d0,-6.15d0,
     $         0.75d0, 0.45d0/),
     $         (/2,2,3/))

          nodes1(9:10,2:3,:) = reshape((/
     $         2.45d0, -0.26d0,
     $          1.0d0,   0.5d0,
     $        -2.15d0,  7.85d0,
     $         2.05d0,  9.26d0,
     $         0.75d0,  8.52d0,
     $        -0.25d0,  -0.1d0/),
     $         (/2,2,3/))

          !SW_1
          nodes0(1:3,2:4,:) = reshape((/
     $          1.25d0,  -0.5d0,   0.5d0,
     $           3.0d0,   2.0d0,   1.0d0,
     $           0.8d0,   0.2d0,   0.6d0,
     $         -6.15d0,  0.45d0,  -0.1d0,
     $         -3.26d0,  0.75d0, -0.25d0,
     $         -7.15d0, -6.12d0,  3.25d0,
     $         -9.23d0, -7.85d0, -9.26d0,
     $         -3.26d0,  8.25d0, -2.05d0,
     $         -4.15d0, -3.25d0, -9.26d0/),
     $         (/3,3,3/))

          nodes1(2:4,3:5,:) = reshape((/
     $        -0.26d0,   0.5d0,   0.2d0,
     $         2.45d0,   1.0d0,   8.9d0,
     $          1.2d0,   7.8d0,   2.3d0,
     $         8.52d0,  -0.1d0,  -2.3d0,
     $         0.75d0, -0.25d0, -9.26d0,
     $        -7.15d0, -1.23d0,   5.2d0,
     $        -7.85d0, -9.26d0, -7.26d0,
     $         2.15d0, -2.05d0, -1.25d0,
     $        -7.32d0,  -1.2d0, -9.63d0/),
     $         (/3,3,3/))

          !NE_1
          nodes0(10:12,8:10,:) = reshape((/
     $          0.6d0,  0.2d0, 0.8d0,
     $          1.0d0,  2.0d0, 3.0d0,
     $          0.5d0, -0.5d0,1.25d0,
     $        -3.25d0, 6.12d0,7.15d0,
     $         0.25d0,-0.75d0,3.26d0,
     $          0.1d0,-0.45d0,6.15d0,
     $         9.26d0, 3.25d0,4.15d0,
     $         2.05d0,-8.25d0,3.26d0,
     $         9.26d0, 7.85d0,9.23d0/),
     $         (/3,3,3/))

          nodes1(11:13,9:11,:) = reshape((/
     $         2.3d0,  7.8d0,  1.2d0,
     $         8.9d0,  1.0d0, 2.45d0,
     $         0.2d0,  0.5d0,-0.26d0,
     $        -5.2d0, 1.23d0, 7.15d0,
     $        9.26d0, 0.25d0,-0.75d0,
     $         2.3d0,  0.1d0,-8.52d0,
     $        9.63d0,  1.2d0, 7.32d0,
     $        1.25d0, 2.05d0,-2.15d0,
     $        7.26d0, 9.26d0, 7.85d0/),
     $         (/3,3,3/))


          !set the nodes in bf_compute object
          call bf_layer_used%bf_compute_used%set_alignment(align0)
          call bf_layer_used%bf_compute_used%set_grdpts_id(grdpts_id0)
          call bf_layer_used%bf_compute_used%set_x_map(x_map0)
          call bf_layer_used%bf_compute_used%set_y_map(y_map0)
          call bf_layer_used%bf_compute_used%set_nodes(nodes0)

          !set the nodes in the bf_layer_object
          call bf_layer_used%set_alignment_tab(align1)
          call bf_layer_used%set_x_map(x_map1)
          call bf_layer_used%set_y_map(y_map1)
          call bf_layer_used%set_nodes(nodes1)

          
          !test the computation of the new grid point
          allocate(compute_index(2,3))
          allocate(newgrdpt_data(3,3))

          compute_index(:,1) = [10,1]
          newgrdpt_data(:,1) = [-2.77171875d0,9.87d0,-4.37046875d0]

          compute_index(:,2) = [1,2]
          newgrdpt_data(:,2) = [-11.75157856d0,-6.859983049d0,-7.799983049d0]
          
          compute_index(:,3) = [14,12]
          newgrdpt_data(:,3) = [-11.75157856d0, 6.859983049d0, 7.799983049d0]


          !test the computation of the new grid points
          do k=1, size(compute_index,2)
             
             !indices of the new grid point
             i1 = compute_index(1,k)
             j1 = compute_index(2,k)

             !get the new grid point
             call bf_layer_used%compute_newgrdpt(
     $            p_model,t,dt,
     $            i1, j1)

             !compare with the data
             test_validated = .true.
             newgrdpt = bf_layer_used%get_nodes([i1,j1])
             do l=1,ne
                test_loc = is_test_validated(
     $               newgrdpt(l),
     $               newgrdpt_data(l,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
             end do

             !print the comparison results
             if(detailled) then
                if(test_loc) then
                   print '(''test '',I2, '' validated'')', k
                else
                   print '(''**test '',I2, '' failed**'')', k
                end if
             end if

          end do

        end function test_bf_layer_compute_newgrdpt


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


        subroutine reflect_x(nodes)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes

          real(rkind), dimension(:,:,:), allocatable :: nodes_copy
          integer :: nx_nodes, ny_nodes
          integer :: i,j,k

          nx_nodes = size(nodes,1)
          ny_nodes = size(nodes,2)
          allocate(nodes_copy(nx_nodes,ny_nodes,ne))
          nodes_copy(:,:,:) = nodes(:,:,:)

          do k=1,ne
             do j=1,ny_nodes
                do i=1,nx_nodes

                   if(k.ne.2) then
                      nodes(i,j,k) = nodes_copy(nx_nodes-i+1,j,k)
                   else
                      nodes(i,j,k) =-nodes_copy(nx_nodes-i+1,j,k)
                   end if

                end do
             end do
          end do

          deallocate(nodes_copy)

        end subroutine reflect_x


        subroutine reflect_y(nodes)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes

          real(rkind), dimension(:,:,:), allocatable :: nodes_copy
          integer :: nx_nodes, ny_nodes
          integer :: i,j,k

          nx_nodes = size(nodes,1)
          ny_nodes = size(nodes,2)
          allocate(nodes_copy(nx_nodes,ny_nodes,ne))
          nodes_copy(:,:,:) = nodes(:,:,:)

          do k=1,ne
             do j=1,ny_nodes
                do i=1,nx_nodes

                   if(k.ne.3) then
                      nodes(i,j,k) = nodes_copy(i,ny_nodes-j+1,k)
                   else
                      nodes(i,j,k) =-nodes_copy(i,ny_nodes-j+1,k)
                   end if

                end do
             end do
          end do

          deallocate(nodes_copy)

        end subroutine reflect_y

      end program test_bf_newgrdpt
