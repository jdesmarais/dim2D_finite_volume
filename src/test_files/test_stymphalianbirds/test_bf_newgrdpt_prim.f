      program test_bf_newgrdpt_prim

        use bf_compute_class, only :
     $       bf_compute

        use bf_layer_class, only :
     $       bf_layer

        use bf_newgrdpt_class, only :
     $       bf_newgrdpt

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated,
     $       is_real_matrix_validated

        use parameters_bf_layer, only :
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
     $       SW_corner_type,
     $       no_gradient_type,
     $       gradient_L0_type,
     $       gradient_R0_type,
     $       gradient_I_type,
     $       gradient_xLR0_yI_type,
     $       gradient_xI_yLR0_type,
     $       gradient_xLR0_yLR0_type

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

        use dim2d_parameters, only :
     $       cv_r

        use sd_operators_fd_module, only :
     $       gradient_x_interior,
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_interior,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_R0

        implicit none


        type(bf_newgrdpt) :: bf_newgrdpt_used
        logical           :: test_validated
        logical           :: detailled

        type(pmodel_eq)            :: p_model
        real(rkind), dimension(ne) :: far_field

        !test requirements
        if(ne.ne.4) then
           stop 'the test requires ne=4: DIM2D model'
        end if

        if(.not.is_real_validated(cv_r,2.5d0,.false.)) then
           stop 'the test requires c_v_r=2.5'
        end if

        far_field = p_model%get_far_field(0.0d0,1.0d0,1.0d0)
        if(.not.is_real_vector_validated(
     $       far_field,
     $       [1.465102139d0,0.146510214d0,0.0d0,2.84673289d0],
     $       .true.)) then
           print '(''the test requires p_model%get_far_field(t,x,y)='')'
           print '(''[1.465102139d0,0.146510214d0,0.0d0,2.84673289d0]'')'
           stop ''
           
        end if

        !tes config
        detailled = .false.


        !test of get_interpolation_coeff_1D
        test_validated = test_get_interpolation_coeff_1D(bf_newgrdpt_used,detailled)
        print '(''test_get_interpolation_coeff_1D: '', L1)', test_validated


        !test of interpolate_1D
        test_validated = test_interpolate_1D(bf_newgrdpt_used,detailled)
        print '(''test_interpolate_1D: '', L1)', test_validated


        !test of compute_NewtonCotes_integration
        test_validated = test_compute_NewtonCotes_integration(bf_newgrdpt_used,detailled)
        print '(''test_compute_NewtonCotes_integration: '', L1)', test_validated
c$$$
c$$$
c$$$        !test of compute_newgrdpt_x
c$$$        test_validated = test_compute_newgrdpt_x(bf_newgrdpt_used,detailled)
c$$$        print '(''test_compute_newgrdpt_x: '', L1)', test_validated
c$$$
c$$$
c$$$        !test of compute_newgrdpt_y
c$$$        test_validated = test_compute_newgrdpt_y(bf_newgrdpt_used,detailled)
c$$$        print '(''test_compute_newgrdpt_y: '', L1)', test_validated


        !test of get_interpolation_coeff_2D
        test_validated = test_get_interpolation_coeff_2D(bf_newgrdpt_used,detailled)
        print '(''test_get_interpolation_coeff_2D: '', L1)', test_validated


        !test of interpolate_2D
        test_validated = test_interpolate_2D(bf_newgrdpt_used,detailled)
        print '(''test_interpolate_2D: '', L1)', test_validated


        detailled = .true.


        !test of compute_newgrdpt_xy
        test_validated = test_compute_newgrdpt_xy(bf_newgrdpt_used,detailled)
        print '(''test_compute_newgrdpt_xy: '', L1)', test_validated
c$$$
c$$$
c$$$        !test of compute_newgrdpt
c$$$        test_validated = test_compute_newgrdpt(bf_newgrdpt_used,detailled)
c$$$        print '(''test_compute_newgrdpt: '', L1)', test_validated
c$$$
c$$$
c$$$        !test of bf_compute/compute_newgrdpt
c$$$        test_validated = test_bf_compute_compute_newgrdpt(detailled)
c$$$        print '(''test_bf_compute_compute_newgrdpt: '',L1)', test_validated
c$$$
c$$$
c$$$        !test of bf_layer/compute_newgrdpt
c$$$        test_validated = test_bf_layer_compute_newgrdpt(detailled)
c$$$        print '(''test_bf_layer_compute_newgrdpt: '',L1)', test_validated

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


          x_map = [1.0d0,2.0d0]

          nodes = reshape( (/
     $         1.0d0,2.0d0,
     $         0.5d0,-0.5d0,
     $         0.25d0,-0.75d0,
     $         0.25d0,-0.75d0
     $         /),
     $         (/2,ne/))

          nodes_inter_data = reshape( (/
     $          1.0d0, 0.0d0,
     $         -1.0d0, 1.5d0,
     $         -1.0d0, 1.25d0,
     $         -1.0d0, 1.25d0
     $         /),
     $         (/2,ne/))

          nodes_inter = bf_newgrdpt_used%get_interpolation_coeff_1D(
     $         x_map, nodes)

          test_validated = is_real_matrix_validated(
     $         nodes_inter,
     $         nodes_inter_data,
     $         detailled)

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


          x=1.0
          inter_coeff = reshape( (/
     $          1.0d0, 0.0d0,
     $         -1.0d0, 1.5d0,
     $         -1.0d0, 1.25d0,
     $         -1.0d0, 1.25d0
     $         /),
     $         (/2,ne/))
          
          nodes_inter_data = [1.0d0,0.5d0,0.25d0,0.25d0]

          nodes_inter = bf_newgrdpt_used%interpolate_1D(
     $         x,
     $         inter_coeff)

          test_validated = is_real_vector_validated(
     $            nodes_inter,
     $            nodes_inter_data,
     $            detailled)

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

          
          data0 = [1.0d0 , 1.25d0, 0.5d0, 0.5d0]
          data1 = [2.0d0, -0.25d0, 1.5d0, 1.5d0]

          dt = 0.5d0

          data_test = [0.75d0, 0.25d0,0.50d0,0.5d0]

          data_integrated = bf_newgrdpt_used%compute_NewtonCotes_integration(
     $         data0,data1,dt)

          test_validated = is_real_vector_validated(
     $         data_test,
     $         data_integrated,
     $         detailled)

        end function test_compute_NewtonCotes_integration


c$$$        function test_compute_newgrdpt_x(bf_newgrdpt_used, detailled)
c$$$     $     result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
c$$$          logical           , intent(in) :: detailled
c$$$          logical                        :: test_validated
c$$$
c$$$          type(pmodel_eq)                :: p_model
c$$$          real(rkind)                    :: t
c$$$          real(rkind)                    :: dt
c$$$          
c$$$          integer(ikind), dimension(2,2) :: bf_align0
c$$$          real(rkind), dimension(3)      :: bf_x_map0
c$$$          real(rkind), dimension(2)      :: bf_y_map0
c$$$          real(rkind), dimension(3,2,ne) :: bf_nodes0
c$$$
c$$$          integer(ikind), dimension(2,2) :: bf_align1
c$$$          real(rkind), dimension(3)      :: bf_x_map1
c$$$          real(rkind), dimension(2)      :: bf_y_map1
c$$$          real(rkind), dimension(3,2,ne) :: bf_nodes1
c$$$
c$$$          real(rkind), dimension(ne)     :: newgrdpt_data
c$$$          real(rkind), dimension(ne)     :: newgrdpt
c$$$          integer(ikind)                 :: i1
c$$$          integer(ikind)                 :: j1
c$$$          logical                        :: side_x
c$$$
c$$$          integer                        :: k
c$$$          logical                        :: test_loc
c$$$          
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          !initialization of the inputs
c$$$          t=0.0d0
c$$$          dt=0.25d0
c$$$          
c$$$          bf_align0(1,1) = 0
c$$$          bf_x_map0 = [0.5d0,1.5d0,2.5d0]
c$$$          bf_y_map0 = [0.25d0,0.5d0]
c$$$          bf_nodes0 = reshape((/
c$$$     $         1.0d0,2.0d0, 3.0d0,
c$$$     $         0.5d0,-0.5d0,1.25d0,
c$$$     $         0.25d0,-0.75d0,3.26d0,
c$$$     $         0.1d0 ,-0.45d0,6.15d0,
c$$$     $         2.05d0,-8.25d0,3.26d0,
c$$$     $         9.26d0, 7.85d0,9.23d0
c$$$     $         /),
c$$$     $         (/3,2,ne/))
c$$$
c$$$          bf_align1(1,1) = 1
c$$$          bf_x_map1 = [1.5d0,2.5d0,3.5d0]
c$$$          bf_y_map1 = [0.25d0,0.5d0]
c$$$          bf_nodes1 = reshape((/
c$$$     $         1.0d0,2.45d0, 3.0d0,
c$$$     $         0.5d0,-0.26d0,2.25d0,
c$$$     $         0.25d0,-0.75d0,3.26d0,
c$$$     $         0.1d0 ,-8.52d0,7.15d0,
c$$$     $         2.05d0,-2.15d0,3.26d0,
c$$$     $         9.26d0, 7.85d0,6.23d0
c$$$     $         /),
c$$$     $         (/3,2,ne/))
c$$$
c$$$          i1 = 3
c$$$          j1 = 2
c$$$
c$$$          side_x = right
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-3.345117188d0,4.943867188d0,9.87d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt_x(
c$$$     $         p_model, t,dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1, side_x, gradient_y_y_oneside_R0)
c$$$
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$        end function test_compute_newgrdpt_x
c$$$
c$$$
c$$$        function test_compute_newgrdpt_y(bf_newgrdpt_used, detailled)
c$$$     $     result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
c$$$          logical           , intent(in) :: detailled
c$$$          logical                        :: test_validated
c$$$
c$$$          type(pmodel_eq)                :: p_model
c$$$          real(rkind)                    :: t
c$$$          real(rkind)                    :: dt
c$$$          
c$$$          integer(ikind), dimension(2,2) :: bf_align0
c$$$          real(rkind), dimension(2)      :: bf_x_map0
c$$$          real(rkind), dimension(3)      :: bf_y_map0
c$$$          real(rkind), dimension(2,3,ne) :: bf_nodes0
c$$$
c$$$          integer(ikind), dimension(2,2) :: bf_align1
c$$$          real(rkind), dimension(2)      :: bf_x_map1
c$$$          real(rkind), dimension(3)      :: bf_y_map1
c$$$          real(rkind), dimension(2,3,ne) :: bf_nodes1
c$$$
c$$$          real(rkind), dimension(ne)     :: newgrdpt_data
c$$$          real(rkind), dimension(ne)     :: newgrdpt
c$$$          integer(ikind)                 :: i1
c$$$          integer(ikind)                 :: j1
c$$$          logical                        :: side_y
c$$$
c$$$          integer                        :: k
c$$$          logical                        :: test_loc
c$$$          
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          !initialization of the inputs
c$$$          t=0.0d0
c$$$          dt=0.25d0
c$$$          
c$$$          bf_align0(2,1) = 0
c$$$          bf_x_map0 = [0.25d0,0.5d0]
c$$$          bf_y_map0 = [0.5d0,1.5d0,2.5d0]
c$$$          bf_nodes0 = reshape((/
c$$$     $         1.0d0, 0.5d0,
c$$$     $         2.0d0,-0.5d0,
c$$$     $         3.0d0, 1.25d0,
c$$$     $         2.05d0,9.26d0,
c$$$     $         -8.25d0,7.85d0,
c$$$     $         3.26d0,9.23d0,
c$$$     $         0.25d0,0.1d0,
c$$$     $         -0.75d0,-0.45d0,
c$$$     $         3.26d0,6.15d0
c$$$     $         /),
c$$$     $         (/2,3,ne/))
c$$$
c$$$          bf_align1(2,1) = 1
c$$$          bf_x_map1 = [0.25d0,0.5d0]
c$$$          bf_y_map1 = [1.5d0,2.5d0,3.5d0]
c$$$          bf_nodes1 = reshape((/
c$$$     $         1.0d0,0.5d0,
c$$$     $         2.45d0,-0.26d0,
c$$$     $         3.0d0,2.25d0,
c$$$     $         2.05d0,9.26d0,
c$$$     $         -2.15d0,7.85d0,
c$$$     $         3.26d0,6.23d0,
c$$$     $         0.25d0,0.1d0,
c$$$     $         -0.75d0,-8.52d0,
c$$$     $         3.26d0,7.15d0
c$$$     $         /),
c$$$     $         (/2,3,ne/))
c$$$
c$$$          i1 = 2
c$$$          j1 = 3
c$$$
c$$$          side_y = right
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-3.345117188d0,9.87d0,4.943867188d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt_y(
c$$$     $         p_model, t,dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1, side_y, gradient_x_x_oneside_R0)
c$$$
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$        end function test_compute_newgrdpt_y
c$$$
c$$$

        function test_sym_compute_newgrdpt_xy(bf_newgrdpt_used, detailled)
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
          

          !initialization of the inputs
          t=0.0d0
          dt=0.25d0


          !computation of the new grdpt          
          bf_align0(1,1) = 0
          bf_align0(2,1) = 0
          bf_x_map0 = [0.5d0, 1.5d0 , 2.5d0]
          bf_y_map0 = [0.0d0, 0.25d0, 0.5d0]
          bf_nodes0 = reshape((/
     $         1.48d0, 1.30d0, 1.35d0,
     $         1.26d0, 1.45d0, 1.4d0,
     $         1.46d0, 1.27d0, 1.47d0,
     $         
     $         0.128d0, 0.127d0, 0.142d0,
     $         1.138d0, 0.148d0, 0.132d0,
     $         0.146d0, 0.143d0, 0.145d0,
     $         
     $         0.0050d0, 0.020d0, 0.060d0,
     $         0.0025d0, 0.001d0, 0.015d0,
     $         0.0100d0, 0.002d0, 0.050d0,
     $         
     $         4.88d0, 4.870d0,	4.855d0,
     $         4.85d0, 4.865d0, 4.845d0,
     $         4.89d0, 4.870d0, 4.860d0/),
     $         (/3,3,ne/))

          bf_align1(1,1) = 0
          bf_align0(2,1) = 0
          bf_x_map1 = [0.5d0, 1.5d0,2.5d0, 3.5d0]
          bf_y_map1 = [0.0d0,0.25d0,0.5d0,0.75d0]
          bf_nodes1 = reshape((/
     $         1.50d0, 1.455d0, 1.48d0, 0.0d0,
     $         1.20d0, 1.350d0, 1.25d0, 0.0d0,
     $         1.49d0, 1.250d0, 1.40d0, 0.0d0,
     $         0.00d0, 0.000d0, 0.00d0, 0.0d0,
     $         
     $         0.128d0, 0.450d0, 0.135d0, 0.0d0,
     $         0.148d0, 0.150d0, 0.122d0, 0.0d0,
     $         0.142d0, 1.152d0, 0.236d0, 0.0d0,
     $         0.000d0, 0.000d0, 0.000d0, 0.0d0,
     $         
     $         0.006d0,	0.0600d0, 0.020d0, 0.0d0,
     $         0.000d0,	0.0028d0, 0.035d0, 0.0d0,
     $         0.020d0,	0.0030d0, 0.040d0, 0.0d0,
     $         0.000d0, 0.0000d0, 0.000d0, 0.0d0,
     $         
     $         4.876d0, 4.825d0, 4.862d0, 0.0d0,
     $         4.890d0, 4.871d0, 4.892d0, 0.0d0,
     $         4.865d0, 4.757d0, 4.895d0, 0.0d0,
     $         0.000d0, 0.0000d0, 0.000d0, 0.0d0
     $         /),
     $         (/4,4,ne/))

          i1 = 4
          j1 = 4

          n_direction         = n2_direction
          side_n              = right
          eigen_indices       = [3,3]
          inter_indices1(:,1) = [3,2]
          inter_indices1(:,2) = [2,3]
          inter_indices1(:,3) = [3,3]

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


          !compute the symmetrized newgrdpt
          bf_align0(1,1) = 0
          bf_align0(2,1) = 1
          bf_x_map0 = [0.5d0, 1.5d0 , 2.5d0]
          bf_y_map0 = [0.0d0, 0.25d0, 0.5d0]
          bf_nodes0 = reshape((/
     $         1.48d0, 1.30d0, 1.35d0,
     $         1.26d0, 1.45d0, 1.4d0,
     $         1.46d0, 1.27d0, 1.47d0,
     $         
     $         0.128d0, 0.127d0, 0.142d0,
     $         1.138d0, 0.148d0, 0.132d0,
     $         0.146d0, 0.143d0, 0.145d0,
     $         
     $         0.0050d0, 0.020d0, 0.060d0,
     $         0.0025d0, 0.001d0, 0.015d0,
     $         0.0100d0, 0.002d0, 0.050d0,
     $         
     $         4.88d0, 4.870d0,	4.855d0,
     $         4.85d0, 4.865d0, 4.845d0,
     $         4.89d0, 4.870d0, 4.860d0/),
     $         (/3,3,ne/))

          bf_align1(1,1) = 0
          bf_align0(2,1) = 0
          bf_x_map1 = [  0.5d0, 1.5d0, 2.5d0, 3.5d0]
          bf_y_map1 = [-0.25d0, 0.0d0, 0.25d0,0.5d0]
          bf_nodes1 = reshape((/
     $         0.00d0, 0.000d0, 0.00d0, 0.0d0,
     $         1.49d0, 1.250d0, 1.40d0, 0.0d0,
     $         1.20d0, 1.350d0, 1.25d0, 0.0d0,
     $         1.50d0, 1.455d0, 1.48d0, 0.0d0,
     $         
     $         0.000d0, 0.000d0, 0.000d0, 0.0d0,
     $         0.142d0, 1.152d0, 0.236d0, 0.0d0,
     $         0.148d0, 0.150d0, 0.122d0, 0.0d0,
     $         0.128d0, 0.450d0, 0.135d0, 0.0d0,
     $         
     $        -0.000d0, -0.0000d0, -0.000d0, -0.0d0,
     $        -0.020d0,	-0.0030d0, -0.040d0, -0.0d0,
     $        -0.000d0,	-0.0028d0, -0.035d0, -0.0d0,
     $        -0.006d0,	-0.0600d0, -0.020d0, -0.0d0,
     $         
     $         0.000d0, 0.0000d0, 0.000d0, 0.0d0
     $         4.865d0, 4.757d0, 4.895d0, 0.0d0,
     $         4.890d0, 4.871d0, 4.892d0, 0.0d0,
     $         4.876d0, 4.825d0, 4.862d0, 0.0d0
     $         /),
     $         (/4,4,ne/))

          i1 = 1
          j1 = 1

          n_direction         = n1_direction
          side_n              = right
          eigen_indices       = [3,2]
          inter_indices1(:,1) = [2,2]
          inter_indices1(:,2) = [3,2]
          inter_indices1(:,3) = [3,3]

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
          


        end function test_sym_compute_newgrdpt_xy


c$$$        function test_compute_newgrdpt_xy(bf_newgrdpt_used, detailled)
c$$$     $     result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
c$$$          logical           , intent(in) :: detailled
c$$$          logical                        :: test_validated
c$$$
c$$$          type(pmodel_eq)                :: p_model
c$$$          real(rkind)                    :: t
c$$$          real(rkind)                    :: dt
c$$$          
c$$$          integer(ikind), dimension(2,2) :: bf_align0
c$$$          real(rkind), dimension(3)      :: bf_x_map0
c$$$          real(rkind), dimension(3)      :: bf_y_map0
c$$$          real(rkind), dimension(3,3,ne) :: bf_nodes0
c$$$
c$$$          integer(ikind), dimension(2,2) :: bf_align1
c$$$          real(rkind), dimension(4)      :: bf_x_map1
c$$$          real(rkind), dimension(4)      :: bf_y_map1
c$$$          real(rkind), dimension(4,4,ne) :: bf_nodes1
c$$$
c$$$          real(rkind), dimension(ne)     :: newgrdpt_data
c$$$          real(rkind), dimension(ne)     :: newgrdpt
c$$$          integer(ikind)                 :: i1
c$$$          integer(ikind)                 :: j1
c$$$          integer                        :: n_direction
c$$$          logical                        :: side_n
c$$$          integer, dimension(2)          :: eigen_indices
c$$$          integer, dimension(2,3)        :: inter_indices1
c$$$
c$$$          integer                        :: k
c$$$          logical                        :: test_loc
c$$$          
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          !initialization of the inputs
c$$$          t=0.0d0
c$$$          dt=0.25d0
c$$$          
c$$$          bf_align0(1,1) = 0
c$$$          bf_align0(2,1) = 0
c$$$          bf_x_map0 = [0.5d0, 1.5d0 , 2.5d0]
c$$$          bf_y_map0 = [0.0d0, 0.25d0, 0.5d0]
c$$$          bf_nodes0 = reshape((/
c$$$     $          0.6d0,  0.2d0, 0.8d0,
c$$$     $          1.0d0,  2.0d0, 3.0d0,
c$$$     $          0.5d0, -0.5d0, 1.25d0,
c$$$     $        -3.25d0, 6.12d0, 7.15d0,
c$$$     $         0.25d0,-0.75d0, 3.26d0,
c$$$     $          0.1d0,-0.45d0, 6.15d0,
c$$$     $         9.26d0, 3.25d0, 4.15d0,
c$$$     $         2.05d0,-8.25d0, 3.26d0,
c$$$     $         9.26d0, 7.85d0, 9.23d0/),
c$$$     $         (/3,3,ne/))
c$$$
c$$$          bf_align1(1,1) = 0
c$$$          bf_align0(2,1) = 0
c$$$          bf_x_map1 = [0.5d0, 1.5d0,2.5d0, 3.5d0]
c$$$          bf_y_map1 = [0.0d0,0.25d0,0.5d0,0.75d0]
c$$$          bf_nodes1 = reshape((/
c$$$     $         2.3d0,  7.8d0,   1.2d0, 1.5d0,
c$$$     $         8.9d0,  1.0d0,  2.45d0, 3.0d0,
c$$$     $         0.2d0,  0.5d0, -0.26d0,2.25d0,
c$$$     $        6.23d0,-5.15d0,  2.36d0, 0.0d0,
c$$$     $        -5.2d0, 1.23d0,  7.15d0, 6.2d0,
c$$$     $        9.26d0, 0.25d0, -0.75d0,3.26d0,
c$$$     $         2.3d0, 0.1d0, -8.52d0, 7.15d0,
c$$$     $         0.0d0, 0.0d0,   0.0d0,  0.0d0,
c$$$     $        9.63d0, 1.2d0,  7.32d0, 1.52d0,
c$$$     $        1.25d0, 2.05d0,-2.15d0, 3.26d0,
c$$$     $        7.26d0, 9.26d0, 7.85d0, 6.23d0,
c$$$     $         0.0d0,  0.0d0,  0.0d0,  0.0d0/),
c$$$     $         (/4,4,ne/))
c$$$
c$$$          i1 = 4
c$$$          j1 = 4
c$$$
c$$$          n_direction         = n2_direction
c$$$          side_n              = right
c$$$          eigen_indices       = [3,3]
c$$$          inter_indices1(:,1) = [3,2]
c$$$          inter_indices1(:,2) = [2,3]
c$$$          inter_indices1(:,3) = [3,3]
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-11.06536693d0,7.194275495d0,8.134275495d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt_xy(
c$$$     $         p_model, t,dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         n_direction,
c$$$     $         side_n,
c$$$     $         gradient_n1_xR0_yR1,
c$$$     $         gradient_n1_xR1_yR0,
c$$$     $         gradient_n1_xR0_yR0,
c$$$     $         eigen_indices,
c$$$     $         inter_indices1)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$        end function test_compute_newgrdpt_xy


c$$$        function test_compute_newgrdpt(bf_newgrdpt_used, detailled)
c$$$     $     result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
c$$$          logical           , intent(in) :: detailled
c$$$          logical                        :: test_validated
c$$$
c$$$          type(pmodel_eq)                :: p_model
c$$$          real(rkind)                    :: t
c$$$          real(rkind)                    :: dt
c$$$          
c$$$          integer(ikind), dimension(2,2)             :: bf_align0
c$$$          real(rkind), dimension(:)    , allocatable :: bf_x_map0
c$$$          real(rkind), dimension(:)    , allocatable :: bf_y_map0
c$$$          real(rkind), dimension(:,:,:), allocatable :: bf_nodes0
c$$$
c$$$          integer(ikind), dimension(2,2)             :: bf_align1
c$$$          real(rkind), dimension(:)    , allocatable :: bf_x_map1
c$$$          real(rkind), dimension(:)    , allocatable :: bf_y_map1
c$$$          real(rkind), dimension(:,:,:), allocatable :: bf_nodes1
c$$$
c$$$          real(rkind), dimension(ne)     :: newgrdpt_data
c$$$          real(rkind), dimension(ne)     :: newgrdpt
c$$$          integer(ikind)                 :: i1
c$$$          integer(ikind)                 :: j1
c$$$
c$$$          integer                        :: k
c$$$          logical                        :: test_loc
c$$$          logical                        :: test_NE_corner
c$$$          logical                        :: test_NW_corner
c$$$          logical                        :: test_SE_corner
c$$$          logical                        :: test_SW_corner
c$$$          logical                        :: test_E_edge
c$$$          logical                        :: test_W_edge
c$$$          logical                        :: test_N_edge
c$$$          logical                        :: test_S_edge
c$$$          logical                        :: test_NE_edge_I
c$$$          logical                        :: test_NE_edge_xLR0          
c$$$          logical                        :: test_NE_edge_yLR0
c$$$          logical                        :: test_NE_edge_xyLR0
c$$$          logical                        :: test_NW_edge_I
c$$$          logical                        :: test_NW_edge_xLR0
c$$$          logical                        :: test_NW_edge_yLR0
c$$$          logical                        :: test_NW_edge_xyLR0
c$$$          logical                        :: test_SW_edge_I
c$$$          logical                        :: test_SW_edge_xLR0
c$$$          logical                        :: test_SW_edge_yLR0
c$$$          logical                        :: test_SW_edge_xyLR0
c$$$          logical                        :: test_SE_edge_I
c$$$          logical                        :: test_SE_edge_xLR0
c$$$          logical                        :: test_SE_edge_yLR0
c$$$          logical                        :: test_SE_edge_xyLR0
c$$$
c$$$          allocate(bf_x_map0(3))
c$$$          allocate(bf_y_map0(3))
c$$$          allocate(bf_nodes0(3,3,ne))
c$$$          
c$$$          allocate(bf_x_map1(4))
c$$$          allocate(bf_y_map1(4))
c$$$          allocate(bf_nodes1(4,4,ne))
c$$$
c$$$          !------------------------------------------------
c$$$          !test of the corners
c$$$          !------------------------------------------------
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          !test of the NE corner
c$$$          !------------------------------------------------
c$$$          !initialization of the inputs
c$$$          t=0.0d0
c$$$          dt=0.25d0
c$$$          
c$$$          bf_align0(1,1) = 0
c$$$          bf_align0(2,1) = 0
c$$$          bf_x_map0 = [0.5d0, 1.5d0 , 2.5d0]
c$$$          bf_y_map0 = [0.0d0, 0.25d0, 0.5d0]
c$$$          bf_nodes0 = reshape((/
c$$$     $          0.6d0,  0.2d0, 0.8d0,
c$$$     $          1.0d0,  2.0d0, 3.0d0,
c$$$     $          0.5d0, -0.5d0, 1.25d0,
c$$$     $        -3.25d0, 6.12d0, 7.15d0,
c$$$     $         0.25d0,-0.75d0, 3.26d0,
c$$$     $          0.1d0,-0.45d0, 6.15d0,
c$$$     $         9.26d0, 3.25d0, 4.15d0,
c$$$     $         2.05d0,-8.25d0, 3.26d0,
c$$$     $         9.26d0, 7.85d0, 9.23d0/),
c$$$     $         (/3,3,ne/))
c$$$
c$$$          bf_align1(1,1) = 0
c$$$          bf_align0(2,1) = 0
c$$$          bf_x_map1 = [0.5d0, 1.5d0,2.5d0, 3.5d0]
c$$$          bf_y_map1 = [0.0d0,0.25d0,0.5d0,0.75d0]
c$$$          bf_nodes1 = reshape((/
c$$$     $         2.3d0,  7.8d0,   1.2d0, 1.5d0,
c$$$     $         8.9d0,  1.0d0,  2.45d0, 3.0d0,
c$$$     $         0.2d0,  0.5d0, -0.26d0,2.25d0,
c$$$     $        6.23d0,-5.15d0,  2.36d0, 0.0d0,
c$$$     $        -5.2d0, 1.23d0,  7.15d0, 6.2d0,
c$$$     $        9.26d0, 0.25d0, -0.75d0,3.26d0,
c$$$     $         2.3d0, 0.1d0, -8.52d0, 7.15d0,
c$$$     $         0.0d0, 0.0d0,   0.0d0,  0.0d0,
c$$$     $        9.63d0, 1.2d0,  7.32d0, 1.52d0,
c$$$     $        1.25d0, 2.05d0,-2.15d0, 3.26d0,
c$$$     $        7.26d0, 9.26d0, 7.85d0, 6.23d0,
c$$$     $         0.0d0,  0.0d0,  0.0d0,  0.0d0/),
c$$$     $         (/4,4,ne/))
c$$$
c$$$          i1 = 4
c$$$          j1 = 4
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-11.06536693d0,7.194275495d0,8.134275495d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         NE_corner_type,no_gradient_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_NE_corner = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_NE_corner: '',L1)', test_NE_corner
c$$$          end if
c$$$
c$$$          
c$$$          test_validated = .true.
c$$$
c$$$          !test of the NW corner
c$$$          !------------------------------------------------
c$$$          !initialization of the inputs
c$$$          t=0.0d0
c$$$          dt=0.25d0
c$$$          
c$$$          bf_align0(1,1) = 1
c$$$          bf_align0(2,1) = 0
c$$$          bf_x_map0 = [ 0.5d0, 1.5d0 , 2.5d0]
c$$$          bf_y_map0 = [ 0.0d0, 0.25d0, 0.5d0]
c$$$          call reflect_x(bf_nodes0)
c$$$
c$$$          bf_align1(1,1) = 0
c$$$          bf_align0(2,1) = 0
c$$$          bf_x_map1 = [-0.5d0, 0.5d0, 1.5d0, 2.5d0]
c$$$          bf_y_map1 = [ 0.0d0,0.25d0, 0.5d0,0.75d0]
c$$$          call reflect_x(bf_nodes1)
c$$$
c$$$          i1 = 1
c$$$          j1 = 4
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-11.06536693d0,-7.194275495d0,8.134275495d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t,dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         NW_corner_type,no_gradient_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_NW_corner = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_NW_corner: '',L1)', test_NW_corner
c$$$          end if
c$$$
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          !test of the SW corner
c$$$          !------------------------------------------------
c$$$          !initialization of the inputs
c$$$          t=0.0d0
c$$$          dt=0.25d0
c$$$          
c$$$          bf_align0(1,1) = 1
c$$$          bf_align0(2,1) = 1
c$$$          bf_x_map0 = [ 0.5d0, 1.5d0 , 2.5d0]
c$$$          bf_y_map0 = [ 0.0d0, 0.25d0, 0.5d0]
c$$$          call reflect_y(bf_nodes0)
c$$$
c$$$          bf_align1(1,1) = 0
c$$$          bf_align1(2,1) = 0
c$$$          bf_x_map1 = [ -0.5d0, 0.5d0, 1.5d0, 2.5d0]
c$$$          bf_y_map1 = [-0.25d0, 0.0d0,0.25d0, 0.5d0]
c$$$          call reflect_y(bf_nodes1)
c$$$
c$$$          i1 = 1
c$$$          j1 = 1
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-11.06536693d0,-7.194275495d0,-8.134275495d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t,dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         SW_corner_type,no_gradient_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_SW_corner = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_SW_corner: '',L1)', test_SW_corner
c$$$          end if
c$$$
c$$$
c$$$           test_validated = .true.
c$$$
c$$$          !test of the SE corner
c$$$          !------------------------------------------------
c$$$          !initialization of the inputs
c$$$          t=0.0d0
c$$$          dt=0.25d0
c$$$          
c$$$          bf_align0(1,1) = 0
c$$$          bf_align0(2,1) = 1
c$$$          bf_x_map0 = [ 0.5d0, 1.5d0 , 2.5d0]
c$$$          bf_y_map0 = [ 0.0d0, 0.25d0, 0.5d0]
c$$$          call reflect_x(bf_nodes0)
c$$$
c$$$          bf_align1(1,1) = 0
c$$$          bf_align1(2,1) = 0
c$$$          bf_x_map1 = [  0.5d0, 1.5d0,  2.5d0, 3.5d0]
c$$$          bf_y_map1 = [-0.25d0, 0.0d0, 0.25d0, 0.5d0]
c$$$          call reflect_x(bf_nodes1)
c$$$
c$$$          i1 = 4
c$$$          j1 = 1
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-11.06536693d0, 7.194275495d0,-8.134275495d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t,dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         SE_corner_type,no_gradient_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_SE_corner = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_SE_corner: '',L1)', test_SE_corner
c$$$          end if
c$$$
c$$$          deallocate(bf_x_map0)
c$$$          deallocate(bf_y_map0)
c$$$          deallocate(bf_nodes0)
c$$$          
c$$$          deallocate(bf_x_map1)
c$$$          deallocate(bf_y_map1)
c$$$          deallocate(bf_nodes1)
c$$$
c$$$
c$$$          !------------------------------------------------
c$$$          !test of the x_edge
c$$$          !------------------------------------------------
c$$$
c$$$          allocate(bf_x_map0(3))
c$$$          allocate(bf_y_map0(2))
c$$$          allocate(bf_nodes0(3,2,ne))
c$$$
c$$$          allocate(bf_x_map1(3))
c$$$          allocate(bf_y_map1(2))
c$$$          allocate(bf_nodes1(3,2,ne))
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          !test E_edge
c$$$          t=0.0d0
c$$$          dt=0.25d0
c$$$          
c$$$          bf_align0(1,1) = 0
c$$$          bf_align0(2,1) = 0
c$$$          bf_x_map0 = [0.5d0,1.5d0,2.5d0]
c$$$          bf_y_map0 = [0.25d0,0.5d0]
c$$$          bf_nodes0 = reshape((/
c$$$     $         1.0d0,2.0d0, 3.0d0,
c$$$     $         0.5d0,-0.5d0,1.25d0,
c$$$     $         0.25d0,-0.75d0,3.26d0,
c$$$     $         0.1d0 ,-0.45d0,6.15d0,
c$$$     $         2.05d0,-8.25d0,3.26d0,
c$$$     $         9.26d0, 7.85d0,9.23d0
c$$$     $         /),
c$$$     $         (/3,2,ne/))
c$$$
c$$$          bf_align1(1,1) = 1
c$$$          bf_align1(2,1) = 0
c$$$          bf_x_map1 = [1.5d0,2.5d0,3.5d0]
c$$$          bf_y_map1 = [0.25d0,0.5d0]
c$$$          bf_nodes1 = reshape((/
c$$$     $         1.0d0,2.45d0, 3.0d0,
c$$$     $         0.5d0,-0.26d0,2.25d0,
c$$$     $         0.25d0,-0.75d0,3.26d0,
c$$$     $         0.1d0 ,-8.52d0,7.15d0,
c$$$     $         2.05d0,-2.15d0,3.26d0,
c$$$     $         9.26d0, 7.85d0,6.23d0
c$$$     $         /),
c$$$     $         (/3,2,ne/))
c$$$
c$$$          i1 = 3
c$$$          j1 = 2
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-3.345117188d0,4.943867188d0,9.87d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t,dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         E_edge_type,gradient_R0_type)
c$$$
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_E_edge = test_validated
c$$$
c$$$          if(detailled) then
c$$$             print '(''test_E_edge: '',L1)', test_E_edge
c$$$          end if
c$$$
c$$$           test_validated = .true.
c$$$
c$$$          !test W_edge
c$$$          t=0.0d0
c$$$          dt=0.25d0
c$$$          
c$$$          bf_align0(1,1) = 1
c$$$          bf_align0(2,1) = 0
c$$$          bf_x_map0 = [0.5d0,1.5d0,2.5d0]
c$$$          bf_y_map0 = [0.25d0,0.5d0]
c$$$          call reflect_x(bf_nodes0)
c$$$
c$$$          bf_align1(1,1) = 0
c$$$          bf_align1(2,1) = 0
c$$$          bf_x_map1 = [-0.5d0,0.5d0,1.5d0]
c$$$          bf_y_map1 = [0.25d0,0.5d0]
c$$$          call reflect_x(bf_nodes1)
c$$$
c$$$          i1 = 1
c$$$          j1 = 2
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-3.345117188d0,-4.943867188d0,9.87d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t,dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         W_edge_type,gradient_R0_type)
c$$$
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_W_edge = test_validated
c$$$
c$$$          if(detailled) then
c$$$             print '(''test_W_edge: '',L1)', test_W_edge
c$$$          end if
c$$$
c$$$          deallocate(bf_x_map0)
c$$$          deallocate(bf_y_map0)
c$$$          deallocate(bf_nodes0)
c$$$          
c$$$          deallocate(bf_x_map1)
c$$$          deallocate(bf_y_map1)
c$$$          deallocate(bf_nodes1)
c$$$
c$$$
c$$$          !------------------------------------------------
c$$$          !test of the y_edge
c$$$          !------------------------------------------------
c$$$
c$$$          allocate(bf_x_map0(2))
c$$$          allocate(bf_y_map0(3))
c$$$          allocate(bf_nodes0(2,3,ne))
c$$$
c$$$          allocate(bf_x_map1(2))
c$$$          allocate(bf_y_map1(3))
c$$$          allocate(bf_nodes1(2,3,ne))
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          !test N_edge
c$$$          t=0.0d0
c$$$          dt=0.25d0
c$$$          
c$$$          bf_align0(1,1) = 0
c$$$          bf_align0(2,1) = 0
c$$$          bf_x_map0 = [0.25d0,0.5d0]
c$$$          bf_y_map0 = [0.5d0,1.5d0,2.5d0]
c$$$          bf_nodes0 = reshape((/
c$$$     $         1.0d0, 0.5d0,
c$$$     $         2.0d0,-0.5d0,
c$$$     $         3.0d0, 1.25d0,
c$$$     $         2.05d0,9.26d0,
c$$$     $         -8.25d0,7.85d0,
c$$$     $         3.26d0,9.23d0,
c$$$     $         0.25d0,0.1d0,
c$$$     $         -0.75d0,-0.45d0,
c$$$     $         3.26d0,6.15d0
c$$$     $         /),
c$$$     $         (/2,3,ne/))
c$$$
c$$$          bf_align1(1,1) = 0
c$$$          bf_align1(2,1) = 1
c$$$          bf_x_map1 = [0.25d0,0.5d0]
c$$$          bf_y_map1 = [1.5d0,2.5d0,3.5d0]
c$$$          bf_nodes1 = reshape((/
c$$$     $         1.0d0,0.5d0,
c$$$     $         2.45d0,-0.26d0,
c$$$     $         3.0d0,2.25d0,
c$$$     $         2.05d0,9.26d0,
c$$$     $         -2.15d0,7.85d0,
c$$$     $         3.26d0,6.23d0,
c$$$     $         0.25d0,0.1d0,
c$$$     $         -0.75d0,-8.52d0,
c$$$     $         3.26d0,7.15d0
c$$$     $         /),
c$$$     $         (/2,3,ne/))
c$$$
c$$$          i1 = 2
c$$$          j1 = 3
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-3.345117188d0,9.87d0,4.943867188d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t,dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         N_edge_type,
c$$$     $         gradient_R0_type)
c$$$
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_N_edge = test_validated
c$$$
c$$$          if(detailled) then
c$$$             print '(''test_N_edge: '',L1)', test_N_edge
c$$$          end if
c$$$
c$$$          
c$$$          !test S_edge
c$$$          t=0.0d0
c$$$          dt=0.25d0
c$$$          
c$$$          bf_align0(1,1) = 0
c$$$          bf_align0(2,1) = 1
c$$$          bf_x_map0 = [0.25d0,0.5d0]
c$$$          bf_y_map0 = [0.5d0,1.5d0,2.5d0]
c$$$          call reflect_y(bf_nodes0)
c$$$
c$$$          bf_align1(1,1) = 0
c$$$          bf_align1(2,1) = 0
c$$$          bf_x_map1 = [0.25d0,0.5d0]
c$$$          bf_y_map1 = [-0.5d0,0.5d0,1.5d0]
c$$$          call reflect_y(bf_nodes1)
c$$$
c$$$          i1 = 2
c$$$          j1 = 1
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-3.345117188d0,9.87d0,-4.943867188d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t,dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         S_edge_type,
c$$$     $         gradient_R0_type)
c$$$
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_S_edge = test_validated
c$$$
c$$$          if(detailled) then
c$$$             print '(''test_S_edge: '',L1)', test_S_edge
c$$$          end if
c$$$
c$$$          deallocate(bf_x_map0)
c$$$          deallocate(bf_y_map0)
c$$$          deallocate(bf_nodes0)
c$$$          
c$$$          deallocate(bf_x_map1)
c$$$          deallocate(bf_y_map1)
c$$$          deallocate(bf_nodes1)
c$$$
c$$$
c$$$          !------------------------------------------------
c$$$          !test of the edge
c$$$          !------------------------------------------------
c$$$
c$$$          test_validated = .true.
c$$$
c$$$          allocate(bf_x_map0(4))
c$$$          allocate(bf_y_map0(4))
c$$$          allocate(bf_nodes0(4,4,ne))
c$$$          
c$$$          allocate(bf_x_map1(4))
c$$$          allocate(bf_y_map1(4))
c$$$          allocate(bf_nodes1(4,4,ne))
c$$$
c$$$
c$$$          !test of the NE edge: I gradient
c$$$          !------------------------------------------------
c$$$          !initialization of the inputs
c$$$          t=0.0d0
c$$$          dt=0.25d0
c$$$          
c$$$          bf_align0(1,1) = 0
c$$$          bf_align0(2,1) = 0
c$$$          bf_x_map0 = [0.5d0, 1.5d0 , 2.5d0,  3.5d0]
c$$$          bf_y_map0 = [0.0d0, 0.25d0, 0.5d0, 0.75d0]
c$$$          bf_nodes0 = reshape((/
c$$$     $          0.6d0,  0.2d0, 0.8d0, -2.3d0,
c$$$     $          1.0d0,  2.0d0, 3.0d0, 1.23d0,
c$$$     $          0.5d0, -0.5d0, 0.0d0,  0.0d0,
c$$$     $          2.3d0,-6.23d0, 0.0d0,  0.0d0,
c$$$     $         
c$$$     $        -3.25d0, 6.12d0, 7.15d0,1.23d0,
c$$$     $         0.25d0,-0.75d0, 3.26d0,5.86d0,
c$$$     $          0.1d0,-0.45d0,  0.0d0, 0.0d0,
c$$$     $        -2.56d0, 1.23d0,  0.0d0, 0.0d0,
c$$$     $         
c$$$     $         9.26d0, 3.25d0, 4.15d0,9.56d0,
c$$$     $         2.05d0,-8.25d0, 3.26d0,1.23d0,
c$$$     $         9.26d0, 7.85d0,  0.0d0, 0.0d0,
c$$$     $        -6.23d0, 2.36d0,  0.0d0, 0.0d0/),
c$$$     $         (/4,4,ne/))
c$$$
c$$$          bf_align1(1,1) = 0
c$$$          bf_align1(2,1) = 0
c$$$          bf_x_map1 = [0.5d0, 1.5d0,2.5d0, 3.5d0]
c$$$          bf_y_map1 = [0.0d0,0.25d0,0.5d0,0.75d0]
c$$$          bf_nodes1 = reshape((/
c$$$     $         2.3d0,  7.8d0,   1.2d0, 1.5d0,
c$$$     $         8.9d0,  1.0d0,  2.45d0, 3.0d0,
c$$$     $         0.2d0,  0.5d0,   0.0d0, 0.0d0,
c$$$     $        6.23d0,-5.15d0,   0.0d0, 0.0d0,
c$$$     $         
c$$$     $        -5.2d0, 1.23d0,  7.15d0, 6.2d0,
c$$$     $        9.26d0, 0.25d0, -0.75d0,3.26d0,
c$$$     $         2.3d0,  0.1d0,   0.0d0, 0.0d0,
c$$$     $       7.154d0,  1.2d0,   0.0d0, 0.0d0,
c$$$     $         
c$$$     $        9.63d0, 1.2d0,  7.32d0, 1.52d0,
c$$$     $        1.25d0, 2.05d0,-2.15d0, 3.26d0,
c$$$     $        7.26d0, 9.26d0,  0.0d0,  0.0d0,
c$$$     $        5.26d0, 1.45d0,  0.0d0,  0.0d0/),
c$$$     $         (/4,4,ne/))
c$$$
c$$$          i1 = 3
c$$$          j1 = 3
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-4.992695567d0,-3.40475565d0,12.39524435d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         NE_edge_type,gradient_I_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_NE_edge_I = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_NE_edge_I: '',L1)', test_NE_edge_I
c$$$          end if
c$$$
c$$$
c$$$          !test of the NE edge: xLR0_yI gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$          !tested data
c$$$          newgrdpt_data = [-5.065722713d0,-3.426851262d0,12.37314874d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         NE_edge_type,gradient_xLR0_yI_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_NE_edge_xLR0 = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_NE_edge_xLR0_yI: '',L1)', test_NE_edge_xLR0
c$$$          end if
c$$$
c$$$
c$$$          !test of the NE edge: xI_yLR0 gradient
c$$$          !------------------------------------------------
c$$$          !tested data
c$$$          newgrdpt_data = [-5.624210148d0,-3.423421794d0,12.37657821d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         NE_edge_type,gradient_xI_yLR0_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_NE_edge_yLR0 = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_NE_edge_xI_yLR0: '',L1)', test_NE_edge_yLR0
c$$$          end if
c$$$
c$$$
c$$$          !test of the NE edge: xLR0_yLR0 gradient
c$$$          !------------------------------------------------
c$$$          !tested data
c$$$          newgrdpt_data = [-5.697237294d0,-3.445517406d0,12.35448259d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         NE_edge_type,gradient_xLR0_yLR0_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_NE_edge_xyLR0 = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_NE_edge_xLR0_yLR0: '',L1)', test_NE_edge_xyLR0
c$$$          end if
c$$$
c$$$
c$$$          !test of the NW edge: I gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$          call reflect_x(bf_nodes0)
c$$$          call reflect_x(bf_nodes1)
c$$$
c$$$          i1 = 2
c$$$          j1 = 3
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-4.992695567d0,3.40475565d0,12.39524435d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         NW_edge_type,gradient_I_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_NW_edge_I = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_NW_edge_I: '',L1)', test_NW_edge_I
c$$$          end if
c$$$
c$$$
c$$$          !test of the NW edge: xLR0 gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$          !tested data
c$$$          newgrdpt_data = [-5.065722713d0,3.426851262d0,12.37314874d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         NW_edge_type,gradient_xLR0_yI_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_NW_edge_xLR0 = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_NW_edge_xLR0_yI: '',L1)', test_NW_edge_xLR0
c$$$          end if
c$$$
c$$$          !test of the NW edge: yLR0 gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$          !tested data
c$$$          newgrdpt_data = [-5.624210148d0, 3.423421794d0,12.37657821d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         NW_edge_type,gradient_xI_yLR0_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_NW_edge_yLR0 = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_NW_edge_xI_yLR0: '',L1)', test_NW_edge_yLR0
c$$$          end if
c$$$
c$$$
c$$$          !test of the NW edge: xLR0_yLR0 gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$          !tested data
c$$$          newgrdpt_data = [-5.697237294d0, 3.445517406d0,12.35448259d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         NW_edge_type,gradient_xLR0_yLR0_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_NW_edge_xyLR0 = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_NW_edge_xLR0_yLR0: '',L1)', test_NW_edge_xyLR0
c$$$          end if
c$$$
c$$$          
c$$$          !test of the SW edge : I gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$
c$$$          call reflect_y(bf_nodes0)
c$$$          call reflect_y(bf_nodes1)
c$$$
c$$$          i1 = 2
c$$$          j1 = 2
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-4.992695567d0,3.40475565d0,-12.39524435d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         SW_edge_type,gradient_I_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_SW_edge_I = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_SW_edge_I: '',L1)', test_SW_edge_I
c$$$          end if
c$$$
c$$$
c$$$          !test of the SW edge : xLR0 gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$          !tested data
c$$$          newgrdpt_data = [-5.065722713d0,3.426851262d0,-12.37314874d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         SW_edge_type,gradient_xLR0_yI_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_SW_edge_xLR0 = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_SW_edge_xLR0: '',L1)', test_SW_edge_xLR0
c$$$          end if
c$$$
c$$$
c$$$          !test of the SW edge : yLR0 gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$          !tested data
c$$$          newgrdpt_data = [-5.624210148d0, 3.423421794d0,-12.37657821d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         SW_edge_type,gradient_xI_yLR0_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_SW_edge_yLR0 = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_SW_edge_yLR0: '',L1)', test_SW_edge_yLR0
c$$$          end if
c$$$
c$$$
c$$$          !test of the SW edge : xLR0_yLR0 gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$          !tested data
c$$$          newgrdpt_data = [-5.697237294d0, 3.445517406d0,-12.35448259d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         SW_edge_type,gradient_xLR0_yLR0_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_SW_edge_xyLR0 = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_SW_edge_xLR0_yLR0: '',L1)', test_SW_edge_xyLR0
c$$$          end if
c$$$
c$$$
c$$$          !test of the SE edge: I gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$          call reflect_x(bf_nodes0)
c$$$          call reflect_x(bf_nodes1)
c$$$
c$$$          i1 = 3
c$$$          j1 = 2
c$$$
c$$$          !tested data
c$$$          newgrdpt_data = [-4.992695567d0,-3.40475565d0,-12.39524435d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         SE_edge_type,gradient_I_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_SE_edge_I = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_SE_edge_I: '',L1)', test_SE_edge_I
c$$$          end if
c$$$
c$$$          !test of the SE edge: xLR0 gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$          !tested data
c$$$          newgrdpt_data = [-5.065722713d0,-3.426851262d0,-12.37314874d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         SE_edge_type,gradient_xLR0_yI_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_SE_edge_xLR0 = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_SE_edge_xLR0: '',L1)', test_SE_edge_xLR0
c$$$          end if
c$$$
c$$$          !test of the SE edge: yLR0 gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$          !tested data
c$$$          newgrdpt_data = [-5.624210148d0,-3.423421794d0,-12.37657821d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         SE_edge_type,gradient_xI_yLR0_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_SE_edge_yLR0 = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_SE_edge_yLR0: '',L1)', test_SE_edge_yLR0
c$$$          end if
c$$$
c$$$          !test of the SE edge: xLR0_yLR0 gradient
c$$$          !------------------------------------------------
c$$$          test_validated = .true.
c$$$          !tested data
c$$$          newgrdpt_data = [-5.697237294d0,-3.445517406d0,-12.35448259d0]
c$$$
c$$$          !test
c$$$          newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
c$$$     $         p_model, t, dt,
c$$$     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
c$$$     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
c$$$     $         i1,j1,
c$$$     $         SE_edge_type,gradient_xLR0_yLR0_type)
c$$$
c$$$          !comparison
c$$$          do k=1,ne
c$$$
c$$$             test_loc = is_test_validated(
c$$$     $            newgrdpt_data(k),
c$$$     $            newgrdpt(k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$
c$$$          end do
c$$$
c$$$          test_SE_edge_xyLR0 = test_validated
c$$$          if(detailled) then
c$$$             print '(''test_SE_edge_xyLR0: '',L1)', test_SE_edge_xyLR0
c$$$          end if          
c$$$
c$$$          test_validated =
c$$$     $         test_NE_corner.and.
c$$$     $         test_NW_corner.and.
c$$$     $         test_SW_corner.and.
c$$$     $         test_SE_corner.and.
c$$$     $         test_E_edge.and.
c$$$     $         test_W_edge.and.
c$$$     $         test_N_edge.and.
c$$$     $         test_S_edge.and.
c$$$     $         test_NE_edge_I.and.
c$$$     $         test_NE_edge_xLR0.and.
c$$$     $         test_NE_edge_yLR0.and.
c$$$     $         test_NE_edge_xyLR0.and.
c$$$     $         test_NW_edge_I.and.
c$$$     $         test_NW_edge_xLR0.and.
c$$$     $         test_NW_edge_yLR0.and.
c$$$     $         test_NW_edge_xyLR0.and.
c$$$     $         test_SW_edge_I.and.
c$$$     $         test_SW_edge_xLR0.and.
c$$$     $         test_SW_edge_yLR0.and.
c$$$     $         test_SW_edge_xyLR0.and.
c$$$     $         test_SE_edge_I.and.
c$$$     $         test_SE_edge_xLR0.and.
c$$$     $         test_SE_edge_yLR0.and.
c$$$     $         test_SE_edge_xyLR0
c$$$
c$$$        end function test_compute_newgrdpt
c$$$
c$$$
c$$$        function test_bf_compute_compute_newgrdpt(detailled)
c$$$     $     result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          logical, intent(in) :: detailled
c$$$          logical             :: test_validated
c$$$
c$$$
c$$$          integer    , dimension(:,:)  , allocatable :: align0
c$$$          integer    , dimension(:,:)  , allocatable :: grdpts_id0
c$$$          real(rkind), dimension(:)    , allocatable :: x_map0
c$$$          real(rkind), dimension(:)    , allocatable :: y_map0
c$$$          real(rkind), dimension(:,:,:), allocatable :: nodes0
c$$$
c$$$          integer    , dimension(:,:)  , allocatable :: align1
c$$$          integer    , dimension(:,:)  , allocatable :: grdpts_id1
c$$$          real(rkind), dimension(:)    , allocatable :: x_map1
c$$$          real(rkind), dimension(:)    , allocatable :: y_map1
c$$$          real(rkind), dimension(:,:,:), allocatable :: nodes1
c$$$
c$$$          integer(ikind)                             :: i1,j1
c$$$          integer                                    :: k,l
c$$$          logical                                    :: test_loc
c$$$
c$$$          integer    , dimension(:,:)  , allocatable :: compute_index
c$$$          real(rkind), dimension(:,:)  , allocatable :: newgrdpt_data
c$$$          real(rkind), dimension(ne)                 :: newgrdpt
c$$$
c$$$          type(bf_compute) :: bf_compute_used
c$$$          type(pmodel_eq)  :: p_model
c$$$
c$$$          real(rkind) :: t
c$$$          real(rkind) :: dt
c$$$
c$$$          !t and dt initialization
c$$$          t  = 0.0d0
c$$$          dt = 0.25d0
c$$$
c$$$
c$$$          !allocations
c$$$          allocate(align0(2,2))
c$$$          allocate(grdpts_id0(14,12))
c$$$          allocate(x_map0(14))
c$$$          allocate(y_map0(12))
c$$$          allocate(nodes0(14,12,ne))
c$$$
c$$$          allocate(align1(2,2))
c$$$          allocate(grdpts_id1(16,14))
c$$$          allocate(x_map1(16))
c$$$          allocate(y_map1(14))
c$$$          allocate(nodes1(16,14,ne))
c$$$
c$$$
c$$$          !initializations
c$$$
c$$$          !initialization of the alignments
c$$$          !--------------------------------------
c$$$          align0(1,1) = 4
c$$$          align0(1,2) = 13
c$$$          align0(2,1) = 3
c$$$          align0(2,2) = 10
c$$$
c$$$          align1(1,1) = 3
c$$$          align1(1,2) = 14
c$$$          align1(2,1) = 2
c$$$          align1(2,2) = 11
c$$$          
c$$$
c$$$          !initialization of the grdpts
c$$$          !--------------------------------------
c$$$          grdpts_id0 = reshape((/
c$$$     $         0,0,0,0,3,3,3,3,3,0,0,0,0,0,
c$$$     $         3,3,3,3,3,2,2,2,3,3,3,3,3,3,
c$$$     $         3,2,2,2,2,2,1,2,2,2,2,2,2,3,
c$$$     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
c$$$     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
c$$$     $         3,2,2,2,1,1,1,1,1,1,2,2,2,3,
c$$$     $         3,3,3,2,1,1,1,1,1,1,2,3,3,3,
c$$$     $         0,0,3,2,1,1,1,1,1,1,2,3,0,0,
c$$$     $         0,0,3,2,2,2,1,2,2,2,2,3,0,0,
c$$$     $         0,0,3,3,3,2,1,2,3,3,3,3,0,0,
c$$$     $         0,0,0,0,3,2,2,2,3,0,0,0,0,0,
c$$$     $         0,0,0,0,3,3,3,3,3,0,0,0,0,0/),
c$$$     $         (/14,12/))
c$$$
c$$$          grdpts_id1 = reshape((/
c$$$     $         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
c$$$     $         0,0,0,0,0,3,3,3,3,3,0,0,0,0,0,0,
c$$$     $         0,3,3,3,3,3,2,2,2,3,3,3,3,3,3,0,
c$$$     $         0,3,2,2,2,2,2,1,2,2,2,2,2,2,3,0,
c$$$     $         0,3,2,1,1,1,1,1,1,1,1,1,1,2,3,0,
c$$$     $         0,3,2,1,1,1,1,1,1,1,1,1,1,2,3,0,
c$$$     $         0,3,2,2,2,1,1,1,1,1,1,2,2,2,3,0,
c$$$     $         0,3,3,3,2,1,1,1,1,1,1,2,3,3,3,0,
c$$$     $         0,0,0,3,2,1,1,1,1,1,1,2,3,0,0,0,
c$$$     $         0,0,0,3,2,2,2,1,2,2,2,2,3,0,0,0,
c$$$     $         0,0,0,3,3,3,2,1,2,3,3,3,3,0,0,0,
c$$$     $         0,0,0,0,0,3,2,2,2,3,0,0,0,0,0,0,
c$$$     $         0,0,0,0,0,3,3,3,3,3,0,0,0,0,0,0,
c$$$     $         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/),
c$$$     $         (/16,14/))
c$$$
c$$$
c$$$          !initialization of the maps
c$$$          !--------------------------------------
c$$$          do k=1, size(x_map0,1)
c$$$             x_map0(k) = 1.0d0 + 0.25d0*(k-12)
c$$$          end do
c$$$
c$$$          do k=1, size(y_map0,1)
c$$$             y_map0(k) = 0.25d0 + 0.25d0*(k-9) !1.0
c$$$          end do
c$$$
c$$$          do k=1, size(x_map1,1)
c$$$             x_map1(k) = 0.75d0 + 0.25d0*(k-12)
c$$$          end do
c$$$
c$$$          do k=1, size(y_map1,1)
c$$$             y_map1(k) = 0.0d0 + 0.25d0*(k-9) !1.0
c$$$          end do
c$$$
c$$$
c$$$          !initialization of the nodes
c$$$          !--------------------------------------
c$$$          !S
c$$$          nodes0(8:9,1:2,:) = reshape((/
c$$$     $         3.0d0,  1.25d0,
c$$$     $         2.0d0,  -0.5d0,
c$$$     $         3.26d0, 9.23d0,
c$$$     $        -8.25d0, 7.85d0,
c$$$     $        -3.26d0,-6.15d0,
c$$$     $         0.75d0, 0.45d0/),
c$$$     $         (/2,2,3/))
c$$$
c$$$          nodes1(9:10,2:3,:) = reshape((/
c$$$     $         2.45d0, -0.26d0,
c$$$     $          1.0d0,   0.5d0,
c$$$     $        -2.15d0,  7.85d0,
c$$$     $         2.05d0,  9.26d0,
c$$$     $         0.75d0,  8.52d0,
c$$$     $        -0.25d0,  -0.1d0/),
c$$$     $         (/2,2,3/))
c$$$
c$$$          !SW_1
c$$$          nodes0(1:3,2:4,:) = reshape((/
c$$$     $          1.25d0,  -0.5d0,   0.5d0,
c$$$     $           3.0d0,   2.0d0,   1.0d0,
c$$$     $           0.8d0,   0.2d0,   0.6d0,
c$$$     $         -6.15d0,  0.45d0,  -0.1d0,
c$$$     $         -3.26d0,  0.75d0, -0.25d0,
c$$$     $         -7.15d0, -6.12d0,  3.25d0,
c$$$     $         -9.23d0, -7.85d0, -9.26d0,
c$$$     $         -3.26d0,  8.25d0, -2.05d0,
c$$$     $         -4.15d0, -3.25d0, -9.26d0/),
c$$$     $         (/3,3,3/))
c$$$
c$$$          nodes1(2:4,3:5,:) = reshape((/
c$$$     $        -0.26d0,   0.5d0,   0.2d0,
c$$$     $         2.45d0,   1.0d0,   8.9d0,
c$$$     $          1.2d0,   7.8d0,   2.3d0,
c$$$     $         8.52d0,  -0.1d0,  -2.3d0,
c$$$     $         0.75d0, -0.25d0, -9.26d0,
c$$$     $        -7.15d0, -1.23d0,   5.2d0,
c$$$     $        -7.85d0, -9.26d0, -7.26d0,
c$$$     $         2.15d0, -2.05d0, -1.25d0,
c$$$     $        -7.32d0,  -1.2d0, -9.63d0/),
c$$$     $         (/3,3,3/))
c$$$
c$$$          !NE_1
c$$$          nodes0(10:12,8:10,:) = reshape((/
c$$$     $          0.6d0,  0.2d0, 0.8d0,
c$$$     $          1.0d0,  2.0d0, 3.0d0,
c$$$     $          0.5d0, -0.5d0,1.25d0,
c$$$     $        -3.25d0, 6.12d0,7.15d0,
c$$$     $         0.25d0,-0.75d0,3.26d0,
c$$$     $          0.1d0,-0.45d0,6.15d0,
c$$$     $         9.26d0, 3.25d0,4.15d0,
c$$$     $         2.05d0,-8.25d0,3.26d0,
c$$$     $         9.26d0, 7.85d0,9.23d0/),
c$$$     $         (/3,3,3/))
c$$$
c$$$          nodes1(11:13,9:11,:) = reshape((/
c$$$     $         2.3d0,  7.8d0,  1.2d0,
c$$$     $         8.9d0,  1.0d0, 2.45d0,
c$$$     $         0.2d0,  0.5d0,-0.26d0,
c$$$     $        -5.2d0, 1.23d0, 7.15d0,
c$$$     $        9.26d0, 0.25d0,-0.75d0,
c$$$     $         2.3d0,  0.1d0,-8.52d0,
c$$$     $        9.63d0,  1.2d0, 7.32d0,
c$$$     $        1.25d0, 2.05d0,-2.15d0,
c$$$     $        7.26d0, 9.26d0, 7.85d0/),
c$$$     $         (/3,3,3/))
c$$$
c$$$
c$$$          !set the nodes in bf_compute object
c$$$          call bf_compute_used%set_alignment(align0)
c$$$          call bf_compute_used%set_grdpts_id(grdpts_id0)
c$$$          call bf_compute_used%set_x_map(x_map0)
c$$$          call bf_compute_used%set_y_map(y_map0)
c$$$          call bf_compute_used%set_nodes(nodes0)
c$$$
c$$$          
c$$$          !test the computation of the new grid point
c$$$          allocate(compute_index(2,3))
c$$$          allocate(newgrdpt_data(3,3))
c$$$
c$$$          compute_index(:,1) = [10,1]
c$$$          newgrdpt_data(:,1) = [-2.77171875d0,9.87d0,-4.37046875d0]
c$$$
c$$$          compute_index(:,2) = [1,2]
c$$$          newgrdpt_data(:,2) = [-11.75157856d0,-6.859983049d0,-7.799983049d0]
c$$$          
c$$$
c$$$          compute_index(:,3) = [14,12]
c$$$          newgrdpt_data(:,3) = [-11.75157856d0, 6.859983049d0, 7.799983049d0]
c$$$
c$$$
c$$$          !test the computation of the new grid points
c$$$          do k=1, size(compute_index,2)
c$$$             
c$$$             !indices of the new grid point
c$$$             i1 = compute_index(1,k)
c$$$             j1 = compute_index(2,k)
c$$$
c$$$             !get the new grid point
c$$$             newgrdpt = bf_compute_used%compute_newgrdpt(
c$$$     $            p_model,t,dt,
c$$$     $            align1, x_map1, y_map1, nodes1,
c$$$     $            i1, j1)
c$$$
c$$$             !compare with the data
c$$$             test_validated = .true.
c$$$             do l=1,ne
c$$$                test_loc = is_test_validated(
c$$$     $               newgrdpt(l),
c$$$     $               newgrdpt_data(l,k),
c$$$     $               detailled)
c$$$                test_validated = test_validated.and.test_loc
c$$$             end do
c$$$
c$$$             !print the comparison results
c$$$             if(detailled) then
c$$$                if(test_loc) then
c$$$                   print '(''test '',I2, '' validated'')', k
c$$$                else
c$$$                   print '(''**test '',I2, '' failed**'')', k
c$$$                end if
c$$$             end if
c$$$
c$$$          end do
c$$$
c$$$        end function test_bf_compute_compute_newgrdpt
c$$$
c$$$
c$$$        function test_bf_layer_compute_newgrdpt(detailled)
c$$$     $     result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          logical, intent(in) :: detailled
c$$$          logical             :: test_validated
c$$$
c$$$
c$$$          integer    , dimension(:,:)  , allocatable :: align0
c$$$          integer    , dimension(:,:)  , allocatable :: grdpts_id0
c$$$          real(rkind), dimension(:)    , allocatable :: x_map0
c$$$          real(rkind), dimension(:)    , allocatable :: y_map0
c$$$          real(rkind), dimension(:,:,:), allocatable :: nodes0
c$$$
c$$$          integer    , dimension(:,:)  , allocatable :: align1
c$$$          integer    , dimension(:,:)  , allocatable :: grdpts_id1
c$$$          real(rkind), dimension(:)    , allocatable :: x_map1
c$$$          real(rkind), dimension(:)    , allocatable :: y_map1
c$$$          real(rkind), dimension(:,:,:), allocatable :: nodes1
c$$$
c$$$          integer(ikind)                             :: i1,j1
c$$$          integer                                    :: k,l
c$$$          logical                                    :: test_loc
c$$$
c$$$          integer    , dimension(:,:)  , allocatable :: compute_index
c$$$          real(rkind), dimension(:,:)  , allocatable :: newgrdpt_data
c$$$          real(rkind), dimension(ne)                 :: newgrdpt
c$$$
c$$$          type(bf_layer)  :: bf_layer_used
c$$$          type(pmodel_eq) :: p_model
c$$$
c$$$          real(rkind) :: t
c$$$          real(rkind) :: dt
c$$$
c$$$          !t and dt initialization
c$$$          t  = 0.0d0
c$$$          dt = 0.25d0
c$$$
c$$$
c$$$          !allocations
c$$$          allocate(align0(2,2))
c$$$          allocate(grdpts_id0(14,12))
c$$$          allocate(x_map0(14))
c$$$          allocate(y_map0(12))
c$$$          allocate(nodes0(14,12,ne))
c$$$
c$$$          allocate(align1(2,2))
c$$$          allocate(grdpts_id1(16,14))
c$$$          allocate(x_map1(16))
c$$$          allocate(y_map1(14))
c$$$          allocate(nodes1(16,14,ne))
c$$$
c$$$
c$$$          !initializations
c$$$
c$$$          !initialization of the alignments
c$$$          !--------------------------------------
c$$$          align0(1,1) = 4
c$$$          align0(1,2) = 13
c$$$          align0(2,1) = 3
c$$$          align0(2,2) = 10
c$$$
c$$$          align1(1,1) = 3
c$$$          align1(1,2) = 14
c$$$          align1(2,1) = 2
c$$$          align1(2,2) = 11
c$$$          
c$$$
c$$$          !initialization of the grdpts
c$$$          !--------------------------------------
c$$$          grdpts_id0 = reshape((/
c$$$     $         0,0,0,0,3,3,3,3,3,0,0,0,0,0,
c$$$     $         3,3,3,3,3,2,2,2,3,3,3,3,3,3,
c$$$     $         3,2,2,2,2,2,1,2,2,2,2,2,2,3,
c$$$     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
c$$$     $         3,2,1,1,1,1,1,1,1,1,1,1,2,3,
c$$$     $         3,2,2,2,1,1,1,1,1,1,2,2,2,3,
c$$$     $         3,3,3,2,1,1,1,1,1,1,2,3,3,3,
c$$$     $         0,0,3,2,1,1,1,1,1,1,2,3,0,0,
c$$$     $         0,0,3,2,2,2,1,2,2,2,2,3,0,0,
c$$$     $         0,0,3,3,3,2,1,2,3,3,3,3,0,0,
c$$$     $         0,0,0,0,3,2,2,2,3,0,0,0,0,0,
c$$$     $         0,0,0,0,3,3,3,3,3,0,0,0,0,0/),
c$$$     $         (/14,12/))
c$$$
c$$$          grdpts_id1 = reshape((/
c$$$     $         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
c$$$     $         0,0,0,0,0,3,3,3,3,3,0,0,0,0,0,0,
c$$$     $         0,3,3,3,3,3,2,2,2,3,3,3,3,3,3,0,
c$$$     $         0,3,2,2,2,2,2,1,2,2,2,2,2,2,3,0,
c$$$     $         0,3,2,1,1,1,1,1,1,1,1,1,1,2,3,0,
c$$$     $         0,3,2,1,1,1,1,1,1,1,1,1,1,2,3,0,
c$$$     $         0,3,2,2,2,1,1,1,1,1,1,2,2,2,3,0,
c$$$     $         0,3,3,3,2,1,1,1,1,1,1,2,3,3,3,0,
c$$$     $         0,0,0,3,2,1,1,1,1,1,1,2,3,0,0,0,
c$$$     $         0,0,0,3,2,2,2,1,2,2,2,2,3,0,0,0,
c$$$     $         0,0,0,3,3,3,2,1,2,3,3,3,3,0,0,0,
c$$$     $         0,0,0,0,0,3,2,2,2,3,0,0,0,0,0,0,
c$$$     $         0,0,0,0,0,3,3,3,3,3,0,0,0,0,0,0,
c$$$     $         0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/),
c$$$     $         (/16,14/))
c$$$
c$$$
c$$$          !initialization of the maps
c$$$          !--------------------------------------
c$$$          do k=1, size(x_map0,1)
c$$$             x_map0(k) = 1.0d0 + 0.25d0*(k-12)
c$$$          end do
c$$$
c$$$          do k=1, size(y_map0,1)
c$$$             y_map0(k) = 0.25d0 + 0.25d0*(k-9) !1.0
c$$$          end do
c$$$
c$$$          do k=1, size(x_map1,1)
c$$$             x_map1(k) = 0.75d0 + 0.25d0*(k-12)
c$$$          end do
c$$$
c$$$          do k=1, size(y_map1,1)
c$$$             y_map1(k) = 0.0d0 + 0.25d0*(k-9) !1.0
c$$$          end do
c$$$
c$$$
c$$$          !initialization of the nodes
c$$$          !--------------------------------------
c$$$          !S
c$$$          nodes0(8:9,1:2,:) = reshape((/
c$$$     $         3.0d0,  1.25d0,
c$$$     $         2.0d0,  -0.5d0,
c$$$     $         3.26d0, 9.23d0,
c$$$     $        -8.25d0, 7.85d0,
c$$$     $        -3.26d0,-6.15d0,
c$$$     $         0.75d0, 0.45d0/),
c$$$     $         (/2,2,3/))
c$$$
c$$$          nodes1(9:10,2:3,:) = reshape((/
c$$$     $         2.45d0, -0.26d0,
c$$$     $          1.0d0,   0.5d0,
c$$$     $        -2.15d0,  7.85d0,
c$$$     $         2.05d0,  9.26d0,
c$$$     $         0.75d0,  8.52d0,
c$$$     $        -0.25d0,  -0.1d0/),
c$$$     $         (/2,2,3/))
c$$$
c$$$          !SW_1
c$$$          nodes0(1:3,2:4,:) = reshape((/
c$$$     $          1.25d0,  -0.5d0,   0.5d0,
c$$$     $           3.0d0,   2.0d0,   1.0d0,
c$$$     $           0.8d0,   0.2d0,   0.6d0,
c$$$     $         -6.15d0,  0.45d0,  -0.1d0,
c$$$     $         -3.26d0,  0.75d0, -0.25d0,
c$$$     $         -7.15d0, -6.12d0,  3.25d0,
c$$$     $         -9.23d0, -7.85d0, -9.26d0,
c$$$     $         -3.26d0,  8.25d0, -2.05d0,
c$$$     $         -4.15d0, -3.25d0, -9.26d0/),
c$$$     $         (/3,3,3/))
c$$$
c$$$          nodes1(2:4,3:5,:) = reshape((/
c$$$     $        -0.26d0,   0.5d0,   0.2d0,
c$$$     $         2.45d0,   1.0d0,   8.9d0,
c$$$     $          1.2d0,   7.8d0,   2.3d0,
c$$$     $         8.52d0,  -0.1d0,  -2.3d0,
c$$$     $         0.75d0, -0.25d0, -9.26d0,
c$$$     $        -7.15d0, -1.23d0,   5.2d0,
c$$$     $        -7.85d0, -9.26d0, -7.26d0,
c$$$     $         2.15d0, -2.05d0, -1.25d0,
c$$$     $        -7.32d0,  -1.2d0, -9.63d0/),
c$$$     $         (/3,3,3/))
c$$$
c$$$          !NE_1
c$$$          nodes0(10:12,8:10,:) = reshape((/
c$$$     $          0.6d0,  0.2d0, 0.8d0,
c$$$     $          1.0d0,  2.0d0, 3.0d0,
c$$$     $          0.5d0, -0.5d0,1.25d0,
c$$$     $        -3.25d0, 6.12d0,7.15d0,
c$$$     $         0.25d0,-0.75d0,3.26d0,
c$$$     $          0.1d0,-0.45d0,6.15d0,
c$$$     $         9.26d0, 3.25d0,4.15d0,
c$$$     $         2.05d0,-8.25d0,3.26d0,
c$$$     $         9.26d0, 7.85d0,9.23d0/),
c$$$     $         (/3,3,3/))
c$$$
c$$$          nodes1(11:13,9:11,:) = reshape((/
c$$$     $         2.3d0,  7.8d0,  1.2d0,
c$$$     $         8.9d0,  1.0d0, 2.45d0,
c$$$     $         0.2d0,  0.5d0,-0.26d0,
c$$$     $        -5.2d0, 1.23d0, 7.15d0,
c$$$     $        9.26d0, 0.25d0,-0.75d0,
c$$$     $         2.3d0,  0.1d0,-8.52d0,
c$$$     $        9.63d0,  1.2d0, 7.32d0,
c$$$     $        1.25d0, 2.05d0,-2.15d0,
c$$$     $        7.26d0, 9.26d0, 7.85d0/),
c$$$     $         (/3,3,3/))
c$$$
c$$$
c$$$          !set the nodes in bf_compute object
c$$$          call bf_layer_used%bf_compute_used%set_alignment(align0)
c$$$          call bf_layer_used%bf_compute_used%set_grdpts_id(grdpts_id0)
c$$$          call bf_layer_used%bf_compute_used%set_x_map(x_map0)
c$$$          call bf_layer_used%bf_compute_used%set_y_map(y_map0)
c$$$          call bf_layer_used%bf_compute_used%set_nodes(nodes0)
c$$$
c$$$          !set the nodes in the bf_layer_object
c$$$          call bf_layer_used%set_alignment_tab(align1)
c$$$          call bf_layer_used%set_x_map(x_map1)
c$$$          call bf_layer_used%set_y_map(y_map1)
c$$$          call bf_layer_used%set_nodes(nodes1)
c$$$
c$$$          
c$$$          !test the computation of the new grid point
c$$$          allocate(compute_index(2,3))
c$$$          allocate(newgrdpt_data(3,3))
c$$$
c$$$          compute_index(:,1) = [10,1]
c$$$          newgrdpt_data(:,1) = [-2.77171875d0,9.87d0,-4.37046875d0]
c$$$
c$$$          compute_index(:,2) = [1,2]
c$$$          newgrdpt_data(:,2) = [-11.75157856d0,-6.859983049d0,-7.799983049d0]
c$$$          
c$$$          compute_index(:,3) = [14,12]
c$$$          newgrdpt_data(:,3) = [-11.75157856d0, 6.859983049d0, 7.799983049d0]
c$$$
c$$$
c$$$          !test the computation of the new grid points
c$$$          do k=1, size(compute_index,2)
c$$$             
c$$$             !indices of the new grid point
c$$$             i1 = compute_index(1,k)
c$$$             j1 = compute_index(2,k)
c$$$
c$$$             !get the new grid point
c$$$             call bf_layer_used%compute_newgrdpt(
c$$$     $            p_model,t,dt,
c$$$     $            i1, j1)
c$$$
c$$$             !compare with the data
c$$$             test_validated = .true.
c$$$             newgrdpt = bf_layer_used%get_nodes([i1,j1])
c$$$             do l=1,ne
c$$$                test_loc = is_test_validated(
c$$$     $               newgrdpt(l),
c$$$     $               newgrdpt_data(l,k),
c$$$     $               detailled)
c$$$                test_validated = test_validated.and.test_loc
c$$$             end do
c$$$
c$$$             !print the comparison results
c$$$             if(detailled) then
c$$$                if(test_loc) then
c$$$                   print '(''test '',I2, '' validated'')', k
c$$$                else
c$$$                   print '(''**test '',I2, '' failed**'')', k
c$$$                end if
c$$$             end if
c$$$
c$$$          end do
c$$$
c$$$        end function test_bf_layer_compute_newgrdpt
c$$$
c$$$
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
          logical                      :: test_loc


          x_map = [1.0d0,1.0d0,2.0d0]
          y_map = [1.0d0,2.0d0,2.0d0]

          nodes = reshape( (/
     $         1.0d0,1.0d0,1.0d0,
     $         1.0d0,1.0d0,0.0d0,
     $         0.0d0,1.0d0,1.0d0,
     $         0.0d0,1.0d0,1.0d0
     $         /),
     $         (/3,ne/))

          nodes_inter_data = reshape( (/
     $           0.0d0, 0.0d0, 1.0d0,
     $          -1.0d0, 0.0d0, 2.0d0,
     $           0.0d0, 1.0d0,-1.0d0,
     $           0.0d0, 1.0d0,-1.0d0
     $         /),
     $         (/3,ne/))

          nodes_inter = bf_newgrdpt_used%get_interpolation_coeff_2D(
     $         x_map, y_map,nodes)

          test_validated = is_real_matrix_validated(
     $         nodes_inter,
     $         nodes_inter_data,
     $         detailled)

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


          x=0.25d0
          y=0.5d0
          inter_coeff = reshape( (/
     $          0.0d0, 0.0d0, 1.0d0,
     $         -1.0d0, 0.0d0, 1.0d0,
     $          0.0d0,-1.0d0, 1.0d0,
     $          0.0d0,-1.0d0, 1.0d0
     $         /),
     $         (/3,ne/))
          
          nodes_inter_data = [1.0d0,0.75d0,0.5d0,0.5d0]

          nodes_inter = bf_newgrdpt_used%interpolate_2D(
     $         x,y,
     $         inter_coeff)

          test_validated = is_real_vector_validated(
     $         nodes_inter,
     $         nodes_inter_data,
     $         detailled)

        end function test_interpolate_2D
c$$$
c$$$        
c$$$        subroutine reflect_x(nodes)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind), dimension(:,:,:), intent(inout) :: nodes
c$$$
c$$$          real(rkind), dimension(:,:,:), allocatable :: nodes_copy
c$$$          integer :: nx_nodes, ny_nodes
c$$$          integer :: i,j,k
c$$$
c$$$          nx_nodes = size(nodes,1)
c$$$          ny_nodes = size(nodes,2)
c$$$          allocate(nodes_copy(nx_nodes,ny_nodes,ne))
c$$$          nodes_copy(:,:,:) = nodes(:,:,:)
c$$$
c$$$          do k=1,ne
c$$$             do j=1,ny_nodes
c$$$                do i=1,nx_nodes
c$$$
c$$$                   if(k.ne.2) then
c$$$                      nodes(i,j,k) = nodes_copy(nx_nodes-i+1,j,k)
c$$$                   else
c$$$                      nodes(i,j,k) =-nodes_copy(nx_nodes-i+1,j,k)
c$$$                   end if
c$$$
c$$$                end do
c$$$             end do
c$$$          end do
c$$$
c$$$          deallocate(nodes_copy)
c$$$
c$$$        end subroutine reflect_x
c$$$
c$$$
c$$$        subroutine reflect_y(nodes)
c$$$
c$$$          implicit none
c$$$
c$$$          real(rkind), dimension(:,:,:), intent(inout) :: nodes
c$$$
c$$$          real(rkind), dimension(:,:,:), allocatable :: nodes_copy
c$$$          integer :: nx_nodes, ny_nodes
c$$$          integer :: i,j,k
c$$$
c$$$          nx_nodes = size(nodes,1)
c$$$          ny_nodes = size(nodes,2)
c$$$          allocate(nodes_copy(nx_nodes,ny_nodes,ne))
c$$$          nodes_copy(:,:,:) = nodes(:,:,:)
c$$$
c$$$          do k=1,ne
c$$$             do j=1,ny_nodes
c$$$                do i=1,nx_nodes
c$$$
c$$$                   if(k.ne.3) then
c$$$                      nodes(i,j,k) = nodes_copy(i,ny_nodes-j+1,k)
c$$$                   else
c$$$                      nodes(i,j,k) =-nodes_copy(i,ny_nodes-j+1,k)
c$$$                   end if
c$$$
c$$$                end do
c$$$             end do
c$$$          end do
c$$$
c$$$          deallocate(nodes_copy)
c$$$
c$$$        end subroutine reflect_y

      end program test_bf_newgrdpt_prim
