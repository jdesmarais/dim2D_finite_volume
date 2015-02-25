      program test_bf_newgrdpt_prim

        use bf_compute_class, only :
     $       bf_compute

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
     $       n1_direction,
     $       n2_direction,
     $       obc_eigenqties_bc,
     $       obc_eigenqties_lin,
     $       vector_x,
     $       vector_y

        use parameters_input, only :
     $       ne,
     $       obc_eigenqties_strategy,
     $       ic_choice

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
        logical           :: test_loc
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

        call p_model%initial_conditions%ini_far_field()

        far_field = p_model%get_far_field(0.0d0,1.0d0,1.0d0)

        if(.not.is_real_vector_validated(
     $       far_field,
     $       [1.46510213931996d0,0.146510214d0,0.0d0,2.84673289046992d0],
     $       .true.)) then
           print '(''the test requires p_model%get_far_field(t,x,y)='')'
           print '(''[1.465102139d0,0.14651021d0,0.0d0,2.84673289d0]'')'
           print '()'
           print '(''T0 should be 0.95'')'
           print '(''flow_direction should be x-direction'')'
           print '(''ic_choice should be newgrdpt_test'')'
           print '()'
           stop ''
           
        end if

        !test config
        detailled = .false.
        test_validated = .true.


        !test of get_interpolation_coeff_1D
        test_loc = test_get_interpolation_coeff_1D(bf_newgrdpt_used,detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_interpolation_coeff_1D: '', L1)', test_loc


        !test of interpolate_1D
        test_loc = test_interpolate_1D(bf_newgrdpt_used,detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_interpolate_1D: '', L1)', test_loc


        !test of compute_NewtonCotes_integration
        test_loc = test_compute_NewtonCotes_integration(bf_newgrdpt_used,detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_NewtonCotes_integration: '', L1)', test_loc

        !test of get_interpolation_coeff_2D
        test_loc = test_get_interpolation_coeff_2D(bf_newgrdpt_used,detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_interpolation_coeff_2D: '', L1)', test_loc


        !test of interpolate_2D
        test_loc = test_interpolate_2D(bf_newgrdpt_used,detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_interpolate_2D: '', L1)', test_loc
        print '()'

        detailled = .true.

        !test of compute_newgrdpt_x + symmetry check
        test_loc = test_sym_compute_newgrdpt_x(bf_newgrdpt_used,p_model,detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_sym_compute_newgrdpt_x: '', L1)', test_loc
        print '()'

        !to test this use, flow_direction=y_direction
        !test of compute_newgrdpt_y + symmetry check
        test_loc = test_sym_compute_newgrdpt_y(bf_newgrdpt_used,p_model,detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_sym_compute_newgrdpt_y: '', L1)', test_loc
        print '()'

        !test of compute_newgrdpt_xy + symmetry check
        test_loc = test_sym_compute_newgrdpt_xy(bf_newgrdpt_used,p_model,detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_sym_compute_newgrdpt_xy: '', L1)', test_loc
        print '()'

        !test of compute_newgrdpt
        test_loc = test_compute_newgrdpt(bf_newgrdpt_used,p_model,detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_newgrdpt: '', L1)', test_loc

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


        function test_sym_compute_newgrdpt_x(
     $     bf_newgrdpt_used,
     $     p_model,
     $     detailled)
     $     result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in)    :: bf_newgrdpt_used
          type(pmodel_eq)   , intent(inout) :: p_model
          logical           , intent(in)    :: detailled
          logical                           :: test_validated

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
          real(rkind), dimension(ne)     :: newgrdpt_sym
          integer(ikind)                 :: i1
          integer(ikind)                 :: j1
          logical                        :: side_x

          real(rkind), dimension(ne)     :: far_field

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
          bf_align1(2,1) = 0
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
          j1 = 2

          side_x = right

          newgrdpt = bf_newgrdpt_used%compute_newgrdpt_x(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,side_x,gradient_y_interior)

          select case(obc_eigenqties_strategy)

            case(obc_eigenqties_bc)
               newgrdpt_data = [
     $               1.22383078395524d0,
     $               0.39531842478603d0,
     $              -0.21050217290879d0,
     $               4.19684181018688d0]

            case(obc_eigenqties_lin)
               newgrdpt_data = [
     $              1.21167555521982d0,
     $              0.35901827468671d0,
     $             -0.20475732388332d0,
     $              4.20002914561351d0]
               
            case default
               print '(''test_bf_newgrdpt_prim'')'
               print '(''test_sym_compute_newgrdpt_x'')'
               print '(''obc_eigenqties_strategy not recognized'')'
               stop ''
               
          end select


          test_loc = is_real_vector_validated(
     $         newgrdpt,
     $         newgrdpt_data,
     $         detailled)
          test_validated = test_loc
          print '(''compute_newgrdpt_x: '',L1)', test_loc
          

          !compute the symmetrized newgrdpt
          bf_align0(1,1) = 1
          bf_align0(2,1) = 0
          bf_x_map0 = [-2.5d0, -1.5d0, -0.5d0]
          bf_y_map0 = [ 0.0d0, 0.25d0,  0.5d0]
          bf_nodes0 = reshape((/
     $         1.35d0, 1.30d0, 1.48d0, 
     $         1.4d0,  1.45d0, 1.26d0, 
     $         1.47d0, 1.27d0, 1.46d0, 
     $         
     $        -0.142d0,-0.127d0,-0.128d0,
     $        -0.132d0,-0.148d0,-1.138d0,
     $        -0.145d0,-0.143d0,-0.146d0,
     $         
     $         0.060d0, 0.020d0, 0.0050d0,
     $         0.015d0, 0.001d0, 0.0025d0,
     $         0.050d0, 0.002d0, 0.0100d0,
     $         
     $         4.855d0, 4.870d0, 4.88d0, 
     $         4.845d0, 4.865d0, 4.85d0, 
     $         4.860d0, 4.870d0, 4.89d0/),
     $         (/3,3,ne/))

          bf_align1(1,1) = 0
          bf_align1(2,1) = 0
          bf_x_map1 = [ -3.5d0,-2.5d0,-1.5d0,-0.5d0]
          bf_y_map1 = [  0.0d0,0.25d0, 0.5d0,0.75d0]
          bf_nodes1 = reshape((/
     $         0.0d0, 1.48d0, 1.455d0, 1.50d0,
     $         0.0d0, 1.25d0, 1.350d0, 1.20d0,
     $         0.0d0, 1.40d0, 1.250d0, 1.49d0,
     $         0.0d0, 0.00d0, 0.000d0, 0.00d0,
     $         
     $         0.0d0, -0.135d0, -0.450d0, -0.128d0,
     $         0.0d0, -0.122d0, -0.150d0, -0.148d0,
     $         0.0d0, -0.236d0, -1.152d0, -0.142d0,
     $         0.0d0,  0.000d0,  0.000d0,  0.000d0,
     $         
     $         0.0d0,  0.020d0, 0.0600d0, 0.006d0,
     $         0.0d0,  0.035d0, 0.0028d0, 0.000d0,
     $         0.0d0,  0.040d0, 0.0030d0, 0.020d0,
     $         0.0d0,  0.000d0, 0.0000d0, 0.000d0,
     $         
     $         0.0d0,  4.862d0,  4.825d0, 4.876d0,
     $         0.0d0,  4.892d0,  4.871d0, 4.890d0,
     $         0.0d0,  4.895d0,  4.757d0, 4.865d0,
     $         0.0d0,  0.000d0, 0.0000d0, 0.000d0
     $         /),
     $         (/4,4,ne/))

          i1 = 1
          j1 = 2

          side_x = left

          far_field = p_model%initial_conditions%get_far_field(0.0,1.0,0.0)
          far_field(2) = - far_field(2)
          call p_model%initial_conditions%set_far_field(far_field)

          newgrdpt_sym = bf_newgrdpt_used%compute_newgrdpt_x(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,side_x,gradient_y_interior)

          newgrdpt_sym(2) = -newgrdpt_sym(2)

          test_loc = is_real_vector_validated(
     $         newgrdpt_sym,
     $         newgrdpt_data,
     $         detailled)

          test_validated = test_validated.and.test_loc
          print '(''symmetry_newgrdpt_x: '',L1)', test_loc

          call p_model%initial_conditions%ini_far_field()

        end function test_sym_compute_newgrdpt_x


        function test_sym_compute_newgrdpt_y(
     $     bf_newgrdpt_used,
     $     p_model,
     $     detailled)
     $     result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in)    :: bf_newgrdpt_used
          type(pmodel_eq)   , intent(inout) :: p_model
          logical           , intent(in)    :: detailled
          logical                           :: test_validated

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
          real(rkind), dimension(ne)     :: newgrdpt_sym
          integer(ikind)                 :: i1
          integer(ikind)                 :: j1
          logical                        :: side_y

          logical                        :: test_loc

          real(rkind)                    :: temp
          real(rkind), dimension(ne)     :: far_field
          

          !initialization of the inputs
          t=0.0d0
          dt=0.25d0


          !computation of the new grdpt          
          bf_align0(1,1) = 0
          bf_align0(2,1) = 0
          bf_x_map0 = [0.0d0, 0.25d0, 0.5d0]
          bf_y_map0 = [0.5d0, 1.5d0 , 2.5d0]
          bf_nodes0 = reshape((/
     $         1.48d0, 1.26d0, 1.46d0,
     $         1.30d0, 1.45d0, 1.27d0,
     $         1.35d0,  1.4d0, 1.47d0,
     $         
     $         0.0050d0, 0.0025d0, 0.01d0,
     $         0.02d0 , 0.001d0 , 0.002d0,
     $         0.06d0  , 0.015d0 , 0.050d0,
     $         
     $         0.128d0, 1.138d0, 0.146d0,
     $         0.127d0, 0.148d0, 0.143d0,
     $         0.142d0, 0.132d0, 0.145d0,
     $         
     $         4.88d0, 4.850d0,	4.89d0,
     $         4.87d0, 4.865d0, 4.87d0,
     $         4.855d0,4.845d0, 4.86d0/),
     $         (/3,3,ne/))

          bf_align1(1,1) = 0
          bf_align1(2,1) = 0
          bf_x_map1 = [0.0d0,0.25d0,0.5d0,0.75d0]
          bf_y_map1 = [0.5d0, 1.5d0,2.5d0, 3.5d0]
          bf_nodes1 = reshape((/
     $         1.50d0 , 1.20d0 , 1.49d0, 0.0d0,
     $         1.455d0, 1.350d0, 1.25d0, 0.0d0,
     $         1.48d0 , 1.250d0, 1.40d0, 0.0d0,
     $         0.00d0 , 0.000d0, 0.00d0, 0.0d0,
     $         
     $         0.006d0,	0.0000d0, 0.020d0, 0.0d0,
     $         0.060d0,	0.0028d0, 0.003d0, 0.0d0,
     $         0.020d0,	0.0350d0, 0.040d0, 0.0d0,
     $         0.000d0, 0.0000d0, 0.000d0, 0.0d0,
     $         
     $         0.128d0, 0.148d0, 0.142d0, 0.0d0,
     $         0.450d0, 0.150d0, 1.152d0, 0.0d0,
     $         0.135d0, 0.122d0, 0.236d0, 0.0d0,
     $         0.000d0, 0.000d0, 0.000d0, 0.0d0,
     $         
     $         4.876d0, 4.890d0, 4.865d0, 0.0d0,
     $         4.825d0, 4.871d0, 4.757d0, 0.0d0,
     $         4.862d0, 4.892d0, 4.895d0, 0.0d0,
     $         0.000d0, 0.0000d0, 0.000d0, 0.0d0
     $         /),
     $         (/4,4,ne/))

          i1 = 2
          j1 = 4

          side_y = right

          far_field    = p_model%initial_conditions%get_far_field(0.0,1.0,0.0)
          temp         = far_field(2)
          far_field(2) = far_field(3)
          far_field(3) = temp
          call p_model%initial_conditions%set_far_field(far_field)

          newgrdpt = bf_newgrdpt_used%compute_newgrdpt_y(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,side_y,gradient_x_interior)

          select case(obc_eigenqties_strategy)

            case(obc_eigenqties_bc)
               newgrdpt_data = [
     $               1.22383078395524d0,
     $              -0.21050217290879d0,
     $               0.39531842478603d0,
     $               4.19684181018688d0]

            case(obc_eigenqties_lin)
               newgrdpt_data = [
     $              1.21167555521982d0,
     $             -0.20475732388332d0,
     $              0.35901827468671d0,
     $              4.20002914561351d0]
               
            case default
               print '(''test_bf_newgrdpt_prim'')'
               print '(''test_sym_compute_newgrdpt_y'')'
               print '(''obc_eigenqties_strategy not recognized'')'
               stop ''
               
          end select


          test_loc = is_real_vector_validated(
     $         newgrdpt,
     $         newgrdpt_data,
     $         detailled)
          test_validated = test_loc
          print '(''compute_newgrdpt_y: '',L1)', test_loc
          

          !compute the symmetrized newgrdpt
          bf_align0(1,1) = 0
          bf_align0(2,1) = 1
          bf_x_map0 = [ 0.0d0, 0.25d0,  0.5d0]
          bf_y_map0 = [-2.5d0, -1.5d0, -0.5d0]
          bf_nodes0 = reshape((/
     $         1.35d0,  1.4d0, 1.47d0,
     $         1.30d0, 1.45d0, 1.27d0,
     $         1.48d0, 1.26d0, 1.46d0,
     $         
     $         0.06d0  , 0.015d0 , 0.050d0,
     $         0.02d0  , 0.001d0 , 0.002d0,
     $         0.0050d0, 0.0025d0, 0.01d0,
     $         
     $        -0.142d0, -0.132d0, -0.145d0,
     $        -0.127d0, -0.148d0, -0.143d0,
     $        -0.128d0, -1.138d0, -0.146d0,
     $         
     $         4.855d0,4.845d0, 4.86d0,
     $         4.87d0, 4.865d0, 4.87d0,
     $         4.88d0, 4.850d0,	4.89d0/),
     $         (/3,3,ne/))


          bf_align1(1,1) = 0
          bf_align1(2,1) = 0
          bf_x_map1 = [  0.0d0,0.25d0, 0.5d0, 0.75d0]
          bf_y_map1 = [ -3.5d0,-2.5d0,-1.5d0,-0.500 ]
          bf_nodes1 = reshape((/
     $         0.00d0 , 0.000d0, 0.00d0, 0.0d0,
     $         1.48d0 , 1.250d0, 1.40d0, 0.0d0,
     $         1.455d0, 1.350d0, 1.25d0, 0.0d0,
     $         1.50d0 , 1.20d0 , 1.49d0, 0.0d0,
     $         
     $         0.000d0, 0.0000d0, 0.000d0, 0.0d0,
     $         0.020d0,	0.0350d0, 0.040d0, 0.0d0,
     $         0.060d0,	0.0028d0, 0.003d0, 0.0d0,
     $         0.006d0,	0.0000d0, 0.020d0, 0.0d0,
     $         
     $        -0.000d0, -0.000d0, -0.000d0, -0.0d0,
     $        -0.135d0, -0.122d0, -0.236d0, -0.0d0,
     $        -0.450d0, -0.150d0, -1.152d0, -0.0d0,
     $        -0.128d0, -0.148d0, -0.142d0, -0.0d0,
     $         
     $         0.000d0, 0.0000d0, 0.000d0, 0.0d0,
     $         4.862d0, 4.892d0, 4.895d0, 0.0d0,
     $         4.825d0, 4.871d0, 4.757d0, 0.0d0,
     $         4.876d0, 4.890d0, 4.865d0, 0.0d0
     $         /),
     $         (/4,4,ne/))

          i1 = 2
          j1 = 1

          side_y = left

          far_field    = p_model%initial_conditions%get_far_field(0.0,1.0,0.0)
          far_field(3) = - far_field(3)
          call p_model%initial_conditions%set_far_field(far_field)

          newgrdpt_sym = bf_newgrdpt_used%compute_newgrdpt_y(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,side_y,gradient_x_interior)

          newgrdpt_sym(3) = -newgrdpt_sym(3)

          test_loc = is_real_vector_validated(
     $              newgrdpt_sym,
     $              newgrdpt_data,
     $              detailled)
          test_validated = test_validated.and.test_loc
          print '(''symmetry_newgrdpt_x: '',L1)', test_loc

          call p_model%initial_conditions%ini_far_field()

        end function test_sym_compute_newgrdpt_y


        function test_sym_compute_newgrdpt_xy(bf_newgrdpt_used, p_model, detailled)
     $     result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in) :: bf_newgrdpt_used
          type(pmodel_eq)   , intent(in) :: p_model
          logical           , intent(in) :: detailled
          logical                        :: test_validated

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
          real(rkind), dimension(ne)     :: newgrdpt_sym
          integer(ikind)                 :: i1
          integer(ikind)                 :: j1
          integer                        :: n_direction
          logical                        :: side_n
          integer, dimension(2)          :: eigen_indices
          integer, dimension(2,3)        :: inter_indices1

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
          bf_align1(2,1) = 0
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
     $         gradient_x_x_oneside_R0,
     $         gradient_y_interior,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_R0,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_R0,
     $         eigen_indices,
     $         inter_indices1)

          select case(obc_eigenqties_strategy)

            case(obc_eigenqties_bc)
               newgrdpt_data = [
     $              1.55742091189301d0,
     $              0.56224898725725d0,
     $             -0.14504765609975d0,
     $              4.63647450016455d0]

            case(obc_eigenqties_lin)
               newgrdpt_data = [
     $              1.42742921947925d0,
     $              0.74488540650052d0,
     $             -0.03274712450770d0,
     $              4.27121569817512d0]
               
            case default
               print '(''test_bf_newgrdpt_prim'')'
               print '(''test_sym_compute_newgrdpt_xy'')'
               print '(''obc_eigenqties_strategy not recognized'')'
               stop ''
               
          end select


          test_loc = is_real_vector_validated(
     $         newgrdpt,
     $         newgrdpt_data,
     $         detailled)
          test_validated = test_loc
          print '(''compute_newgrdpt_xy: '',L1)', test_loc
          

          !compute the symmetrized newgrdpt
          bf_align0(1,1) = 0
          bf_align0(2,1) = 1
          bf_x_map0 = [ 0.5d0,   1.5d0, 2.5d0]
          bf_y_map0 = [-0.5d0, -0.25d0, 0.0d0]
          bf_nodes0 = reshape((/
     $         1.46d0, 1.27d0, 1.47d0,
     $         1.26d0, 1.45d0, 1.4d0,
     $         1.48d0, 1.30d0, 1.35d0,
     $         
     $         0.146d0, 0.143d0, 0.145d0,
     $         1.138d0, 0.148d0, 0.132d0,
     $         0.128d0, 0.127d0, 0.142d0,
     $         
     $        -0.0100d0,-0.002d0,-0.050d0,
     $        -0.0025d0,-0.001d0,-0.015d0,
     $        -0.0050d0,-0.020d0,-0.060d0,
     $         
     $         4.89d0, 4.870d0, 4.860d0,
     $         4.85d0, 4.865d0, 4.845d0,
     $         4.88d0, 4.870d0,	4.855d0/),
     $         (/3,3,ne/))

          bf_align1(1,1) = 0
          bf_align1(2,1) = 0
          bf_x_map1 = [  0.5d0, 1.5d0, 2.5d0, 3.5d0]
          bf_y_map1 = [-0.75d0,-0.5d0,-0.25d0,0.0d0]
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
     $         0.000d0, 0.0000d0, 0.000d0, 0.0d0,
     $         4.865d0, 4.757d0, 4.895d0, 0.0d0,
     $         4.890d0, 4.871d0, 4.892d0, 0.0d0,
     $         4.876d0, 4.825d0, 4.862d0, 0.0d0
     $         /),
     $         (/4,4,ne/))

          i1 = 4
          j1 = 1

          n_direction         = n1_direction
          side_n              = right
          eigen_indices       = [3,2]
          inter_indices1(:,1) = [2,2]
          inter_indices1(:,2) = [3,2]
          inter_indices1(:,3) = [3,3]

          newgrdpt_sym = bf_newgrdpt_used%compute_newgrdpt_xy(
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         n_direction,
     $         side_n,
     $         gradient_x_interior,
     $         gradient_y_y_oneside_L0,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_y_oneside_L0,
     $         gradient_x_x_oneside_R0,
     $         gradient_y_interior,
     $         eigen_indices,
     $         inter_indices1)


          newgrdpt_sym(3) = -newgrdpt_sym(3)

          test_loc = is_real_vector_validated(
     $         newgrdpt,
     $         newgrdpt_sym,
     $         detailled)

          test_validated = test_validated.and.test_loc
          print '(''symmetry_newgrdpt_xy: '',L1)', test_loc

        end function test_sym_compute_newgrdpt_xy


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


        function test_compute_newgrdpt(bf_newgrdpt_used, p_model, detailled)
     $     result(test_validated)

          implicit none

          class(bf_newgrdpt), intent(in)    :: bf_newgrdpt_used
          type(pmodel_eq)   , intent(inout) :: p_model
          logical           , intent(in)    :: detailled
          logical                           :: test_validated

          real(rkind)                        :: t
          real(rkind)                        :: dt
                                             
          integer(ikind), dimension(2,2)     :: bf_align0
          real(rkind), dimension(5)          :: bf_x_map0
          real(rkind), dimension(5)          :: bf_y_map0
          real(rkind), dimension(5,5,ne)     :: bf_nodes0
                                             
          integer(ikind), dimension(2,2)     :: bf_align1
          real(rkind), dimension(5)          :: bf_x_map1
          real(rkind), dimension(5)          :: bf_y_map1
          real(rkind), dimension(5,5,ne)     :: bf_nodes1
                                             
          integer(ikind)                     :: i1
          integer(ikind)                     :: j1
          
          integer, parameter                 :: nb_tests = 4
          integer, dimension(nb_tests)       :: procedure_type_test
          integer, dimension(nb_tests)       :: gradient_type_test
          real(rkind), dimension(4,nb_tests) :: newgrdpt_data_test
          character(12), dimension(nb_tests)  :: procedure_type_char_test
          character(12), dimension(nb_tests)  :: gradient_type_char_test

          test_validated = .true.


          !initialization of the inputs
          t=0.0d0
          dt=0.25d0


          !computation of the new grdpt
          bf_align0(1,1) = 0
          bf_align0(2,1) = 0

          bf_align1(1,1) = 0
          bf_align1(2,1) = 0

          !edges
          !------------------------------
          call initialize_nodes_and_maps(
     $         bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_x_map1, bf_y_map1, bf_nodes1)

          !E_edge
          i1 = 3
          j1 = 3

          !E_edge, gradient_I_type
          procedure_type_test(1)   = E_edge_type
          gradient_type_test(1)    = gradient_I_type
          newgrdpt_data_test(:,1)  = [
     $         1.04485111659085d0,
     $         0.30761070084187d0,
     $         0.46442897694658d0,
     $         4.27410998299761d0]

          procedure_type_char_test(1) = 'E_edge'
          gradient_type_char_test(1)  = 'I_grad'

          !E_edge, gradient_L0_type
          procedure_type_test(2)   = E_edge_type
          gradient_type_test(2)    = gradient_L0_type
          newgrdpt_data_test(:,2)  = [
     $         1.00164269892366d0,
     $         0.29896444891981d0,
     $        -0.37674225289292d0,
     $         4.01398521393387d0]

          procedure_type_char_test(2) = 'E_edge'
          gradient_type_char_test(2)  = 'L0_grad'

          !E_edge, gradient_R0_type
          procedure_type_test(3)   = E_edge_type
          gradient_type_test(3)    = gradient_R0_type
          newgrdpt_data_test(:,3)  = [
     $         1.08805953425803d0,
     $         0.31590541823443d0,
     $         1.37651536245884d0,
     $         5.24756704620429d0]

          procedure_type_char_test(3) = 'E_edge'
          gradient_type_char_test(3)  = 'R0_grad'

          !test procedure+gradient_type: E_type
          call test_procedure_and_gradient_type(
     $         bf_newgrdpt_used,
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         procedure_type_test(1:3),
     $         gradient_type_test(1:3),
     $         newgrdpt_data_test(:,1:3),
     $         procedure_type_char_test(1:3),
     $         gradient_type_char_test(1:3),
     $         test_validated,
     $         detailled)                  

          !W_type
          call reflect_x(
     $         p_model,
     $         bf_x_map0, bf_nodes0, 
     $         bf_x_map1, bf_nodes1,
     $         i1,
     $         newgrdpt_data_test(:,1:3))

          !W_edge, gradient_I_type
          procedure_type_test(1)      = W_edge_type
          gradient_type_test(1)       = gradient_I_type
          procedure_type_char_test(1) = 'W_edge'
          gradient_type_char_test(1)  = 'I_grad'

          !W_edge, gradient_L0_type
          procedure_type_test(2)      = W_edge_type
          gradient_type_test(2)       = gradient_L0_type
          procedure_type_char_test(2) = 'W_edge'
          gradient_type_char_test(2)  = 'L0_grad'

          !W_edge, gradient_R0_type
          procedure_type_test(3)      = W_edge_type
          gradient_type_test(3)       = gradient_R0_type
          procedure_type_char_test(3) = 'W_edge'
          gradient_type_char_test(3)  = 'R0_grad'

          !test procedure+gradient_type: W_type
          call test_procedure_and_gradient_type(
     $         bf_newgrdpt_used,
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         procedure_type_test(1:3),
     $         gradient_type_test(1:3),
     $         newgrdpt_data_test(:,1:3),
     $         procedure_type_char_test(1:3),
     $         gradient_type_char_test(1:3),
     $         test_validated,
     $         detailled)

          !S_type
          call transpose_data(
     $         p_model,
     $         bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         newgrdpt_data_test(:,1:3))

          !S_edge, gradient_I_type
          procedure_type_test(1)      = S_edge_type
          gradient_type_test(1)       = gradient_I_type
          procedure_type_char_test(1) = 'S_edge'
          gradient_type_char_test(1)  = 'I_grad'

          !S_edge, gradient_L0_type
          procedure_type_test(2)      = S_edge_type
          gradient_type_test(2)       = gradient_L0_type
          procedure_type_char_test(2) = 'S_edge'
          gradient_type_char_test(2)  = 'L0_grad'

          !S_edge, gradient_R0_type
          procedure_type_test(3)      = S_edge_type
          gradient_type_test(3)       = gradient_R0_type
          procedure_type_char_test(3) = 'S_edge'
          gradient_type_char_test(3)  = 'R0_grad'

          !test procedure+gradient_type: S_edge
          call test_procedure_and_gradient_type(
     $         bf_newgrdpt_used,
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         procedure_type_test(1:3),
     $         gradient_type_test(1:3),
     $         newgrdpt_data_test(:,1:3),
     $         procedure_type_char_test(1:3),
     $         gradient_type_char_test(1:3),
     $         test_validated,
     $         detailled)


          !N_type
          call reflect_y(
     $         p_model,
     $         bf_y_map0, bf_nodes0,
     $         bf_y_map1, bf_nodes1,
     $         j1,
     $         newgrdpt_data_test(:,1:3))

          !N_edge, gradient_I_type
          procedure_type_test(1)      = N_edge_type
          gradient_type_test(1)       = gradient_I_type
          procedure_type_char_test(1) = 'N_edge'
          gradient_type_char_test(1)  = 'I_grad'

          !N_edge, gradient_L0_type
          procedure_type_test(2)      = N_edge_type
          gradient_type_test(2)       = gradient_L0_type
          procedure_type_char_test(2) = 'N_edge'
          gradient_type_char_test(2)  = 'L0_grad'

          !N_edge, gradient_R0_type
          procedure_type_test(3)      = N_edge_type
          gradient_type_test(3)       = gradient_R0_type
          procedure_type_char_test(3) = 'N_edge'
          gradient_type_char_test(3)  = 'R0_grad'

          !test procedure+gradient_type: N_edge
          call test_procedure_and_gradient_type(
     $         bf_newgrdpt_used,
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         procedure_type_test(1:3),
     $         gradient_type_test(1:3),
     $         newgrdpt_data_test(:,1:3),
     $         procedure_type_char_test(1:3),
     $         gradient_type_char_test(1:3),
     $         test_validated,
     $         detailled)


          !corners
          !------------------------------
          call initialize_nodes_and_maps(
     $         bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_x_map1, bf_y_map1, bf_nodes1)

          !NE_corner
          i1 = 4
          j1 = 4

          procedure_type_test(1)   = NE_corner_type
          gradient_type_test(1)    = gradient_I_type
          newgrdpt_data_test(:,1)  = [
     $         1.55742091189301d0,
     $         0.56224898725725d0,
     $        -0.14504765609975d0,
     $         4.63647450016455d0]

          procedure_type_char_test(1) = 'NE_corner'
          gradient_type_char_test(1)  = 'I_grad'

          !test procedure+gradient_type: NE_corner
          call test_procedure_and_gradient_type(
     $         bf_newgrdpt_used,
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         procedure_type_test(1:1),
     $         gradient_type_test(1:1),
     $         newgrdpt_data_test(:,1:1),
     $         procedure_type_char_test(1:1),
     $         gradient_type_char_test(1:1),
     $         test_validated,
     $         detailled)

          !NW corner
          call reflect_x(
     $         p_model,
     $         bf_x_map0, bf_nodes0,
     $         bf_x_map1, bf_nodes1,
     $         i1,
     $         newgrdpt_data_test(:,1:1))

          procedure_type_test(1)      = NW_corner_type
          gradient_type_test(1)       = gradient_I_type
          procedure_type_char_test(1) = 'NW_corner'
          gradient_type_char_test(1)  = 'I_grad'

          !test procedure+gradient_type: NW_corner
          call test_procedure_and_gradient_type(
     $         bf_newgrdpt_used,
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         procedure_type_test(1:1),
     $         gradient_type_test(1:1),
     $         newgrdpt_data_test(:,1:1),
     $         procedure_type_char_test(1:1),
     $         gradient_type_char_test(1:1),
     $         test_validated,
     $         detailled)

          !SW corner
          call reflect_y(
     $         p_model,
     $         bf_y_map0, bf_nodes0,
     $         bf_y_map1, bf_nodes1,
     $         j1,
     $         newgrdpt_data_test(:,1:1))

          procedure_type_test(1)      = SW_corner_type
          gradient_type_test(1)       = gradient_I_type
          procedure_type_char_test(1) = 'SW_corner'
          gradient_type_char_test(1)  = 'I_grad'

          !test procedure+gradient_type: SW_corner
          call test_procedure_and_gradient_type(
     $         bf_newgrdpt_used,
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         procedure_type_test(1:1),
     $         gradient_type_test(1:1),
     $         newgrdpt_data_test(:,1:1),
     $         procedure_type_char_test(1:1),
     $         gradient_type_char_test(1:1),
     $         test_validated,
     $         detailled)

          !SE corner
          call reflect_x(
     $         p_model,
     $         bf_x_map0, bf_nodes0,
     $         bf_x_map1, bf_nodes1,
     $         i1,
     $         newgrdpt_data_test(:,1:1))

          procedure_type_test(1)      = SE_corner_type
          gradient_type_test(1)       = gradient_I_type
          procedure_type_char_test(1) = 'SE_corner'
          gradient_type_char_test(1)  = 'I_grad'

          !test procedure+gradient_type: SE_corner
          call test_procedure_and_gradient_type(
     $         bf_newgrdpt_used,
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         procedure_type_test(1:1),
     $         gradient_type_test(1:1),
     $         newgrdpt_data_test(:,1:1),
     $         procedure_type_char_test(1:1),
     $         gradient_type_char_test(1:1),
     $         test_validated,
     $         detailled)

          !anti_corners
          !------------------------------
          call initialize_nodes_and_maps(
     $         bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_x_map1, bf_y_map1, bf_nodes1)

          call p_model%initial_conditions%ini_far_field()

          !NE_anti_corner
          i1 = 3
          j1 = 3

          procedure_type_test(1)   = NE_edge_type
          gradient_type_test(1)    = gradient_I_type
          newgrdpt_data_test(:,1)  = [
     $         1.10410912269715d0,
     $         0.35097007239812d0,
     $         0.23207807274618d0,
     $         3.86145332394228d0]
          procedure_type_char_test(1) = 'NE_edge'
          gradient_type_char_test(1)  = 'I_grad'

          procedure_type_test(2)   = NE_edge_type
          gradient_type_test(2)    = gradient_xLR0_yI_type
          newgrdpt_data_test(:,2)  = [
     $         1.09198841545674d0,
     $         0.34791358035954d0,
     $         0.22873420436383d0,
     $         3.76108958832707d0]
          procedure_type_char_test(2) = 'NE_edge'
          gradient_type_char_test(2)  = 'xLR0_yI'

          procedure_type_test(3)   = NE_edge_type
          gradient_type_test(3)    = gradient_xI_yLR0_type
          newgrdpt_data_test(:,3)  = [
     $         1.34538489199497d0,
     $         0.47081951094764d0,
     $         0.47092928155778d0,
     $         5.00331418010598d0]
          procedure_type_char_test(3) = 'NE_edge'
          gradient_type_char_test(3)  = 'xI_yLR0'

          procedure_type_test(4)   = NE_edge_type
          gradient_type_test(4)    = gradient_xLR0_yLR0_type
          newgrdpt_data_test(:,4)  = [
     $         1.33326418475456d0,
     $         0.46755020750444d0,
     $         0.46571456281873d0,
     $         4.95180514591908d0]
          procedure_type_char_test(4) = 'NE_edge'
          gradient_type_char_test(4)  = 'xLR0_yLR0'

          !test procedure+gradient_type: NE_anti_corner
          call test_procedure_and_gradient_type(
     $         bf_newgrdpt_used,
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         procedure_type_test(1:4),
     $         gradient_type_test(1:4),
     $         newgrdpt_data_test(:,1:4),
     $         procedure_type_char_test(1:4),
     $         gradient_type_char_test(1:4),
     $         test_validated,
     $         detailled)

          !NW_anti_corner
          call reflect_x(
     $         p_model,
     $         bf_x_map0, bf_nodes0,
     $         bf_x_map1, bf_nodes1,
     $         i1,
     $         newgrdpt_data_test(:,1:4))

          procedure_type_test(1)   = NW_edge_type
          gradient_type_test(1)    = gradient_I_type
          procedure_type_char_test(1) = 'NW_edge'
          gradient_type_char_test(1)  = 'I_grad'

          procedure_type_test(2)   = NW_edge_type
          gradient_type_test(2)    = gradient_xLR0_yI_type
          procedure_type_char_test(2) = 'NW_edge'
          gradient_type_char_test(2)  = 'xLR0_yI'

          procedure_type_test(3)   = NW_edge_type
          gradient_type_test(3)    = gradient_xI_yLR0_type
          procedure_type_char_test(3) = 'NW_edge'
          gradient_type_char_test(3)  = 'xI_yLR0'

          procedure_type_test(4)   = NW_edge_type
          gradient_type_test(4)    = gradient_xLR0_yLR0_type
          procedure_type_char_test(4) = 'NW_edge'
          gradient_type_char_test(4)  = 'xLR0_yLR0'

          !test procedure+gradient_type: NW_anti_corner
          call test_procedure_and_gradient_type(
     $         bf_newgrdpt_used,
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         procedure_type_test(1:4),
     $         gradient_type_test(1:4),
     $         newgrdpt_data_test(:,1:4),
     $         procedure_type_char_test(1:4),
     $         gradient_type_char_test(1:4),
     $         test_validated,
     $         detailled)

          !SW_anti_corner
          call reflect_y(
     $         p_model,
     $         bf_y_map0, bf_nodes0,
     $         bf_y_map1, bf_nodes1,
     $         j1,
     $         newgrdpt_data_test(:,1:4))

          procedure_type_test(1)   = SW_edge_type
          gradient_type_test(1)    = gradient_I_type
          procedure_type_char_test(1) = 'SW_edge'
          gradient_type_char_test(1)  = 'I_grad'

          procedure_type_test(2)   = SW_edge_type
          gradient_type_test(2)    = gradient_xLR0_yI_type
          procedure_type_char_test(2) = 'SW_edge'
          gradient_type_char_test(2)  = 'xLR0_yI'

          procedure_type_test(3)   = SW_edge_type
          gradient_type_test(3)    = gradient_xI_yLR0_type
          procedure_type_char_test(3) = 'SW_edge'
          gradient_type_char_test(3)  = 'xI_yLR0'

          procedure_type_test(4)   = SW_edge_type
          gradient_type_test(4)    = gradient_xLR0_yLR0_type
          procedure_type_char_test(4) = 'SW_edge'
          gradient_type_char_test(4)  = 'xLR0_yLR0'

          !test procedure+gradient_type: SW_anti_corner
          call test_procedure_and_gradient_type(
     $         bf_newgrdpt_used,
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         procedure_type_test(1:4),
     $         gradient_type_test(1:4),
     $         newgrdpt_data_test(:,1:4),
     $         procedure_type_char_test(1:4),
     $         gradient_type_char_test(1:4),
     $         test_validated,
     $         detailled)

          !SE_anti_corner
          call reflect_X(
     $         p_model,
     $         bf_x_map0, bf_nodes0,
     $         bf_x_map1, bf_nodes1,
     $         i1,
     $         newgrdpt_data_test(:,1:4))

          procedure_type_test(1)   = SE_edge_type
          gradient_type_test(1)    = gradient_I_type
          procedure_type_char_test(1) = 'SE_edge'
          gradient_type_char_test(1)  = 'I_grad'

          procedure_type_test(2)   = SE_edge_type
          gradient_type_test(2)    = gradient_xLR0_yI_type
          procedure_type_char_test(2) = 'SE_edge'
          gradient_type_char_test(2)  = 'xLR0_yI'

          procedure_type_test(3)   = SE_edge_type
          gradient_type_test(3)    = gradient_xI_yLR0_type
          procedure_type_char_test(3) = 'SE_edge'
          gradient_type_char_test(3)  = 'xI_yLR0'

          procedure_type_test(4)   = SE_edge_type
          gradient_type_test(4)    = gradient_xLR0_yLR0_type
          procedure_type_char_test(4) = 'SE_edge'
          gradient_type_char_test(4)  = 'xLR0_yLR0'

          !test procedure+gradient_type: SE_anti_corner
          call test_procedure_and_gradient_type(
     $         bf_newgrdpt_used,
     $         p_model, t,dt,
     $         bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $         bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $         i1,j1,
     $         procedure_type_test(1:4),
     $         gradient_type_test(1:4),
     $         newgrdpt_data_test(:,1:4),
     $         procedure_type_char_test(1:4),
     $         gradient_type_char_test(1:4),
     $         test_validated,
     $         detailled)

        end function test_compute_newgrdpt


        subroutine test_procedure_and_gradient_type(
     $     bf_newgrdpt_used,
     $     p_model, t,dt,
     $     bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $     bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $     i1,j1,
     $     procedure_type_test,
     $     gradient_type_test,
     $     newgrdpt_data_test,
     $     procedure_type_char_test,
     $     gradient_type_char_test,
     $     test_validated,
     $     detailled)

          implicit none

          class(bf_newgrdpt)              , intent(in)    :: bf_newgrdpt_used
          type(pmodel_eq)                 , intent(in)    :: p_model
          real(rkind)                     , intent(in)    :: t
          real(rkind)                     , intent(in)    :: dt
          integer(ikind), dimension(2,2)  , intent(in)    :: bf_align0
          real(rkind)   , dimension(:)    , intent(in)    :: bf_x_map0
          real(rkind)   , dimension(:)    , intent(in)    :: bf_y_map0
          real(rkind)   , dimension(:,:,:), intent(in)    :: bf_nodes0
          integer(ikind), dimension(2,2)  , intent(in)    :: bf_align1
          real(rkind)   , dimension(:)    , intent(in)    :: bf_x_map1
          real(rkind)   , dimension(:)    , intent(in)    :: bf_y_map1
          real(rkind)   , dimension(:,:,:), intent(in)    :: bf_nodes1
          integer(ikind)                  , intent(in)    :: i1
          integer(ikind)                  , intent(in)    :: j1
          integer       , dimension(:)    , intent(in)    :: procedure_type_test
          integer       , dimension(:)    , intent(in)    :: gradient_type_test
          real(rkind)   , dimension(:,:)  , intent(in)    :: newgrdpt_data_test
          character*(*) , dimension(:)    , intent(in)    :: procedure_type_char_test
          character*(*) , dimension(:)    , intent(in)    :: gradient_type_char_test
          logical                         , intent(inout) :: test_validated
          logical                         , intent(in)    :: detailled

          
          integer                    :: k
          real(rkind), dimension(ne) :: newgrdpt
          real(rkind), dimension(ne) :: newgrdpt_data
          logical                    :: test_loc


          do k=1,size(procedure_type_test,1)

             newgrdpt = bf_newgrdpt_used%compute_newgrdpt(
     $            p_model, t,dt,
     $            bf_align0, bf_x_map0, bf_y_map0, bf_nodes0,
     $            bf_align1, bf_x_map1, bf_y_map1, bf_nodes1,
     $            i1,j1,
     $            procedure_type_test(k),
     $            gradient_type_test(k))

             newgrdpt_data = newgrdpt_data_test(:,k)

             select case(obc_eigenqties_strategy)

               case(obc_eigenqties_bc)

                  test_loc = is_real_vector_validated(
     $                 newgrdpt,
     $                 newgrdpt_data,
     $                 detailled)

                  test_validated = test_validated.and.test_loc

                  print '(''test_compute_newgrdpt ('',A12,'',''A12''): '',L1)',
     $                 procedure_type_char_test(k),
     $                 gradient_type_char_test(k),
     $                 test_loc

               case(obc_eigenqties_lin)
                  stop 'not implemented'
               
               case default
                  print '(''test_bf_newgrdpt_prim'')'
                  print '(''test_sym_compute_newgrdpt_xy'')'
                  print '(''obc_eigenqties_strategy not recognized'')'
                  stop ''
               
             end select

          end do

        end subroutine test_procedure_and_gradient_type


        subroutine initialize_nodes_and_maps(
     $     bf_x_map0, bf_y_map0, bf_nodes0,
     $     bf_x_map1, bf_y_map1, bf_nodes1)

          implicit none

          real(rkind), dimension(5)     , intent(out) :: bf_x_map0
          real(rkind), dimension(5)     , intent(out) :: bf_y_map0
          real(rkind), dimension(5,5,ne), intent(out) :: bf_nodes0
          real(rkind), dimension(5)     , intent(out) :: bf_x_map1
          real(rkind), dimension(5)     , intent(out) :: bf_y_map1
          real(rkind), dimension(5,5,ne), intent(out) :: bf_nodes1
          

          bf_x_map0 = [0.5d0, 1.5d0 , 2.5d0,  3.5d0, 4.5d0]
          bf_y_map0 = [0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0]
          bf_nodes0 = reshape((/
     $         1.48d0, 1.30d0, 1.35d0, 1.31d0, 1.43d0,
     $         1.26d0, 1.45d0,  1.4d0, 1.29d0, 1.37d0,
     $         1.46d0, 1.27d0, 1.47d0, 1.28d0, 1.25d0,
     $         1.48d0, 1.26d0, 1.41d0, 1.34d0, 1.31d0,
     $         1.42d0, 1.46d0, 1.38d0, 1.26d0, 1.37d0,
     $         
     $         0.128d0, 0.127d0, 0.142d0, 0.129d0, 0.136d0,
     $         1.138d0, 0.148d0, 0.132d0, 0.125d0, 0.175d0,
     $         0.146d0, 0.143d0, 0.145d0, 0.182d0, 0.135d0,
     $         0.123d0, 0.129d0, 0.124d0, 0.162d0, 0.152d0,
     $         0.168d0, 0.198d0, 0.186d0, 0.163d0, 0.126d0,
     $         
     $         0.0050d0, 0.020d0, 0.060d0, 0.056d0, 0.062d0,
     $         0.0025d0, 0.001d0, 0.015d0,  0.07d0, 0.085d0,
     $         0.0100d0, 0.002d0, 0.050d0,  0.08d0, 0.015d0,
     $           0.08d0, 0.015d0,  0.09d0, 0.065d0, 0.042d0,
     $          0.026d0,  0.03d0, 0.045d0, 0.052d0, 0.023d0,
     $         
     $         4.88d0, 4.870d0,	4.855d0, 4.834d0, 4.592d0,
     $         4.85d0, 4.865d0, 4.845d0, 4.875d0, 4.815d0,
     $         4.89d0, 4.870d0, 4.860d0, 4.826d0, 4.723d0,
     $         4.83d0,  4.95d0,  4.62d0, 4.952d0, 4.852d0,
     $         4.81d0, 4.758d0, 4.762d0,  4.95d0, 4.703d0
     $         /),
     $         (/5,5,ne/))

          bf_x_map1 = [0.5d0, 1.5d0,2.5d0, 3.5d0,4.5d0]
          bf_y_map1 = [0.0d0,0.25d0,0.5d0,0.75d0, 1.0d0]
          bf_nodes1 = reshape((/
     $         1.50d0, 1.455d0, 1.48d0, 1.29d0, 1.45d0,
     $         1.20d0, 1.350d0, 1.25d0, 1.25d0, 1.36d0,
     $         1.49d0, 1.250d0, 1.40d0, 1.54d0, 1.256d0,
     $         1.52d0, 1.28d0,  1.45d0, 1.52d0, 1.29d0,
     $         1.59d0, 1.26d0,  1.26d0, 1.37d0, 1.48d0,
     $         
     $         0.128d0, 0.450d0, 0.135d0, 0.14d0,  0.18d0,
     $         0.148d0, 0.150d0, 0.122d0, 0.136d0, 0.1742d0,
     $         0.142d0, 1.152d0, 0.236d0, 0.12d0,  0.1592d0,
     $         0.136d0, 0.185d0, 0.296d0, 0.192d0, 0.1263d0,
     $         0.158d0, 0.1265d0, 0.13d0, 0.182d0, 0.1426d0,
     $         
     $         0.006d0,	0.0600d0, 0.020d0,  0.052d0, 0.562d0,
     $         0.000d0,	0.0028d0, 0.035d0, 0.4825d0, 0.365d0,
     $         0.020d0,	0.0030d0, 0.040d0,  0.085d0, 0.0269d0,
     $          0.06d0,  0.085d0, 0.015d0, 0.0256d0, 0.0365d0,
     $         0.059d0, 0.0425d0,  0.05d0, 0.0259d0, 0.0152d0,
     $         
     $         4.876d0, 4.825d0, 4.862d0, 4.821d0, 4.932d0,
     $         4.890d0, 4.871d0, 4.892d0, 4.875d0, 4.926d0,
     $         4.865d0, 4.757d0, 4.895d0, 4.856d0, 4.518d0,
     $         4.625d0, 4.785d0, 4.825d0, 4.826d0, 4.165d0,
     $         4.852d0, 4.926d0, 4.815d0, 4.925d0, 4.291d0
     $         
     $         /),
     $         (/5,5,ne/))

        end subroutine initialize_nodes_and_maps


        subroutine reflect_x(
     $     p_model,
     $     bf_x_map0, bf_nodes0, 
     $     bf_x_map1, bf_nodes1,
     $     i1,
     $     newgrdpt_data_test)

          implicit none

          type(pmodel_eq)               , intent(inout) :: p_model
          real(rkind), dimension(5)     , intent(inout) :: bf_x_map0
          real(rkind), dimension(5,5,ne), intent(inout) :: bf_nodes0
          real(rkind), dimension(5)     , intent(inout) :: bf_x_map1
          real(rkind), dimension(5,5,ne), intent(inout) :: bf_nodes1
          integer(ikind)                , intent(inout) :: i1
          real(rkind), dimension(:,:)   , intent(inout) :: newgrdpt_data_test

          integer :: i
          integer :: i_r
          integer :: j
          integer :: k

          integer    , dimension(ne)     :: var_type
          real(rkind), dimension(ne)     :: far_field
          real(rkind), dimension(5)      :: bf_x_map0_r
          real(rkind), dimension(5,5,ne) :: bf_nodes0_r
          real(rkind), dimension(5)      :: bf_x_map1_r
          real(rkind), dimension(5,5,ne) :: bf_nodes1_r


          !get the type of variables
          var_type = p_model%get_var_type()


          !reflect the far field
          far_field = p_model%initial_conditions%get_far_field(0.0,1.0,0.0)
          do k=1,ne
             if(var_type(k).eq.vector_x) then
                far_field(k) = -far_field(k)
             end if
          end do
          call p_model%initial_conditions%set_far_field(far_field)


          !reflect the x_maps
          do i=1, size(bf_x_map0,1)

             i_r = size(bf_x_map0)-(i-1)

             bf_x_map0_r(i) = -bf_x_map0(i_r)
             bf_x_map1_r(i) = -bf_x_map1(i_r)

          end do
          bf_x_map0 = bf_x_map0_r
          bf_x_map1 = bf_x_map1_r


          !reflect the nodes
          do k=1, ne

             if(var_type(k).eq.vector_x) then

                do j=1, size(bf_nodes0,2)
                   do i=1, size(bf_nodes0,1)
                      i_r = size(bf_nodes0,1)-(i-1)
                      bf_nodes0_r(i,j,k) = - bf_nodes0(i_r,j,k)
                      bf_nodes1_r(i,j,k) = - bf_nodes1(i_r,j,k)
                   end do
                end do

             else

                do j=1, size(bf_nodes0,2)
                   do i=1, size(bf_nodes0,1)
                      i_r = size(bf_nodes0,1)-(i-1)
                      bf_nodes0_r(i,j,k) = bf_nodes0(i_r,j,k)
                      bf_nodes1_r(i,j,k) = bf_nodes1(i_r,j,k)
                   end do
                end do

             end if             

          end do
          bf_nodes0 = bf_nodes0_r
          bf_nodes1 = bf_nodes1_r


          !reflect the i1,j1
          i1 = size(bf_x_map0,1) - (i1-1)


          !reflect the newgrdpt_data_test
          do i=1, size(newgrdpt_data_test,2)
             do k=1, ne
                if(var_type(k).eq.vector_x) then
                   newgrdpt_data_test(k,i) = - newgrdpt_data_test(k,i)
                end if
             end do
          end do

        end subroutine reflect_x


        subroutine reflect_y(
     $     p_model,
     $     bf_y_map0, bf_nodes0, 
     $     bf_y_map1, bf_nodes1,
     $     j1,
     $     newgrdpt_data_test)

          implicit none

          type(pmodel_eq)               , intent(inout) :: p_model
          real(rkind), dimension(5)     , intent(inout) :: bf_y_map0
          real(rkind), dimension(5,5,ne), intent(inout) :: bf_nodes0
          real(rkind), dimension(5)     , intent(inout) :: bf_y_map1
          real(rkind), dimension(5,5,ne), intent(inout) :: bf_nodes1
          integer(ikind)                , intent(inout) :: j1
          real(rkind), dimension(:,:)   , intent(inout) :: newgrdpt_data_test

          integer :: i
          integer :: j
          integer :: j_r
          integer :: k

          integer    , dimension(ne)     :: var_type
          real(rkind), dimension(ne)     :: far_field
          real(rkind), dimension(5)      :: bf_y_map0_r
          real(rkind), dimension(5,5,ne) :: bf_nodes0_r
          real(rkind), dimension(5)      :: bf_y_map1_r
          real(rkind), dimension(5,5,ne) :: bf_nodes1_r


          !get the type of variables
          var_type = p_model%get_var_type()


          !reflect the far field
          far_field = p_model%initial_conditions%get_far_field(0.0,1.0,0.0)
          do k=1,ne
             if(var_type(k).eq.vector_y) then
                far_field(k) = -far_field(k)
             end if
          end do
          call p_model%initial_conditions%set_far_field(far_field)


          !reflect the y_maps
          do j=1, size(bf_y_map0,1)

             j_r = size(bf_y_map0,1)-(j-1)

             bf_y_map0_r(j) = -bf_y_map0(j_r)
             bf_y_map1_r(j) = -bf_y_map1(j_r)

          end do
          bf_y_map0 = bf_y_map0_r
          bf_y_map1 = bf_y_map1_r


          !reflect the nodes
          do k=1, ne

             if(var_type(k).eq.vector_y) then

                do j=1, size(bf_nodes0,2)

                   j_r = size(bf_nodes0,2)-(j-1)

                   do i=1, size(bf_nodes0,1)
                      bf_nodes0_r(i,j,k) = - bf_nodes0(i,j_r,k)
                      bf_nodes1_r(i,j,k) = - bf_nodes1(i,j_r,k)
                   end do

                end do

             else

                do j=1, size(bf_nodes0,2)

                   j_r = size(bf_nodes0,2)-(j-1)

                   do i=1, size(bf_nodes0,1)
                      bf_nodes0_r(i,j,k) = bf_nodes0(i,j_r,k)
                      bf_nodes1_r(i,j,k) = bf_nodes1(i,j_r,k)
                   end do

                end do

             end if             

          end do
          bf_nodes0 = bf_nodes0_r
          bf_nodes1 = bf_nodes1_r


          !reflect the j1
          j1 = size(bf_y_map0,1) - (j1-1)


          !reflect the newgrdpt_data_test
          do i=1, size(newgrdpt_data_test,2)
             do k=1, ne
                if(var_type(k).eq.vector_y) then
                   newgrdpt_data_test(k,i) = - newgrdpt_data_test(k,i)
                end if
             end do
          end do

        end subroutine reflect_y


        subroutine transpose_data(
     $     p_model,
     $     bf_x_map0, bf_y_map0, bf_nodes0, 
     $     bf_x_map1, bf_y_map1, bf_nodes1,
     $     i1,j1,
     $     newgrdpt_data_test)

          implicit none

          type(pmodel_eq)               , intent(inout) :: p_model
          real(rkind), dimension(5)     , intent(inout) :: bf_x_map0
          real(rkind), dimension(5)     , intent(inout) :: bf_y_map0
          real(rkind), dimension(5,5,ne), intent(inout) :: bf_nodes0
          real(rkind), dimension(5)     , intent(inout) :: bf_x_map1
          real(rkind), dimension(5)     , intent(inout) :: bf_y_map1
          real(rkind), dimension(5,5,ne), intent(inout) :: bf_nodes1
          integer(ikind)                , intent(inout) :: i1
          integer(ikind)                , intent(inout) :: j1
          real(rkind), dimension(:,:)   , intent(inout) :: newgrdpt_data_test

          integer :: i
          integer :: k

          real(rkind), dimension(ne)     :: far_field
          real(rkind)                    :: temp

          real(rkind), dimension(5)      :: bf_x_map0_t
          real(rkind), dimension(5)      :: bf_y_map0_t
          real(rkind), dimension(5,5,ne) :: bf_nodes0_t
          real(rkind), dimension(5)      :: bf_x_map1_t
          real(rkind), dimension(5)      :: bf_y_map1_t
          real(rkind), dimension(5,5,ne) :: bf_nodes1_t


          !transpose the far field
          far_field    = p_model%initial_conditions%get_far_field(0.0,1.0,0.0)
          temp         = far_field(2)
          far_field(2) = far_field(3)
          far_field(3) = temp
          call p_model%initial_conditions%set_far_field(far_field)


          !transpose the x_map and y_map
          bf_x_map0_t = bf_y_map0
          bf_y_map0_t = bf_x_map0
          bf_x_map1_t = bf_y_map1
          bf_y_map1_t = bf_x_map1

          bf_x_map0 = bf_x_map0_t
          bf_y_map0 = bf_y_map0_t
          bf_x_map1 = bf_x_map1_t
          bf_y_map1 = bf_y_map1_t


          !transpose the nodes
          do k=1,ne
             
             bf_nodes0_t(:,:,k) = transpose(bf_nodes0(:,:,k))
             bf_nodes1_t(:,:,k) = transpose(bf_nodes1(:,:,k))

          end do

          bf_nodes0(:,:,1) = bf_nodes0_t(:,:,1)
          bf_nodes0(:,:,2) = bf_nodes0_t(:,:,3)
          bf_nodes0(:,:,3) = bf_nodes0_t(:,:,2)
          bf_nodes0(:,:,4) = bf_nodes0_t(:,:,4)

          bf_nodes1(:,:,1) = bf_nodes1_t(:,:,1)
          bf_nodes1(:,:,2) = bf_nodes1_t(:,:,3)
          bf_nodes1(:,:,3) = bf_nodes1_t(:,:,2)
          bf_nodes1(:,:,4) = bf_nodes1_t(:,:,4)


          !transpose the indices
          k  = i1
          i1 = j1
          j1 = k


          !transpose the newgrdpt_data_test
          do i=1, size(newgrdpt_data_test,2)

             temp                    = newgrdpt_data_test(2,i)
             newgrdpt_data_test(2,i) = newgrdpt_data_test(3,i)
             newgrdpt_data_test(3,i) = temp

          end do

        end subroutine transpose_data

      end program test_bf_newgrdpt_prim
