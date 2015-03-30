      program test_bc_operators_openbc_normal

        use bc_operators_class, only :
     $       bc_operators

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated,
     $       is_real_matrix3D_validated
        
        use hedstrom_xy_module, only :
     $       compute_timedev_y_edge_local

        use openbc_operators_module, only :
     $       incoming_left,
     $       incoming_right,
     $       compute_fluxes_at_the_edges_2ndorder

        use parameters_bf_layer, only :
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       no_overlap,
     $       align_N, align_S,
     $       align_E, align_W

        use parameters_constant, only :
     $       vector_x,
     $       vector_y,
     $       left,
     $       right

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_fd_module, only :
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_R0

        use sd_operators_x_oneside_L0_class, only :
     $       sd_operators_x_oneside_L0

        use sd_operators_x_oneside_L1_class, only :
     $       sd_operators_x_oneside_L1

        use sd_operators_x_oneside_R1_class, only :
     $       sd_operators_x_oneside_R1

        use sd_operators_x_oneside_R0_class, only :
     $       sd_operators_x_oneside_R0

        use sd_operators_y_oneside_L0_class, only :
     $       sd_operators_y_oneside_L0

        use sd_operators_y_oneside_L1_class, only :
     $       sd_operators_y_oneside_L1

        use sd_operators_y_oneside_R1_class, only :
     $       sd_operators_y_oneside_R1

        use sd_operators_y_oneside_R0_class, only :
     $       sd_operators_y_oneside_R0

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_compute_fluxes_x_for_bc_y_edge(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_fluxes_x_for_bc_y_edge: '',L1)', test_loc
        print '()'


        test_loc = test_compute_fluxes_y_for_bc_x_edge(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_fluxes_y_for_bc_x_edge: '',L1)', test_loc
        print '()'


        test_loc = test_apply_bc_on_timedev_edge(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_apply_bc_on_timedev_edge: '',L1)', test_loc
        print '()'

        print '(''test_validated: '',L1)', test_validated


        contains


        function test_compute_fluxes_x_for_bc_y_edge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bc_operators)                 :: bc_operators_openbc_normal_used
          type(pmodel_eq)                    :: p_model
          real(rkind)                        :: t          
          type(sd_operators_x_oneside_L0)    :: s_x_L0
          type(sd_operators_x_oneside_L1)    :: s_x_L1
          type(sd_operators_x_oneside_R1)    :: s_x_R1
          type(sd_operators_x_oneside_R0)    :: s_x_R0
          type(sd_operators_y_oneside_L0)    :: s_y_L0
          type(sd_operators_y_oneside_L1)    :: s_y_L1
          type(sd_operators_y_oneside_R1)    :: s_y_R1
          type(sd_operators_y_oneside_R0)    :: s_y_R0  
          real(rkind)                        :: dx
          real(rkind)                        :: dy
          integer(ikind)                     :: i_min
          integer(ikind)                     :: i_max
          integer(ikind)                     :: j_min
          
          real(rkind)   , dimension(nx)         :: interior_x_map
          real(rkind)   , dimension(ny)         :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne)   :: interior_nodes
          real(rkind)   , dimension(nx+1,ny,ne) :: flux_x_ref
          real(rkind)   , dimension(nx,ny+1,ne) :: flux_y_ref
          
          integer(ikind), dimension(2,2)    :: bf_alignment
          integer       , dimension(6,5)    :: bf_grdpts_id
          real(rkind)   , dimension(6)      :: bf_x_map
          real(rkind)   , dimension(5)      :: bf_y_map
          real(rkind)   , dimension(6,5,ne) :: bf_nodes
          real(rkind)   , dimension(7,5,ne) :: bf_flux_x


          integer :: i,j,k

          test_validated = .true.


          !N_edge: left missing grid-points
          !============================================================
          ! in this test, we want to compare the computation of the
          ! fluxes by the buffer layer where only a limited number of
          ! grid-points are available (and so require to extract nodes
          ! from the interior domain) with the reference function
          ! computing the fluxes directly from the interior nodes
          !           _________________
          !          |     |     |     | ____ north buffer layer
          !     _____|_____|_____|_____|
          !    |/////////////////|__________ interior_domain
          !    |//|__|__|__|__|//|
          !    |//|           |//|
          !    |//|           |//|
          !
          !initialize the nodes
          call p_model%initial_conditions%ini_far_field()

          t=0.0d0
          
          interior_x_map = [0.5d0, 1.5d0 , 2.5d0,  3.5d0, 4.5d0, 5.5d0 ]
          interior_y_map = [0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.25d0]

          interior_nodes = reshape((/
     $         1.48d0, 1.30d0, 1.35d0, 1.31d0, 1.43d0, 1.31d0,
     $         1.26d0, 1.45d0, 1.40d0, 1.29d0, 1.37d0, 1.41d0,
     $         1.46d0, 1.27d0, 1.47d0, 1.28d0, 1.25d0, 1.43d0,
     $         1.48d0, 1.26d0, 1.41d0, 1.34d0, 1.31d0, 1.39d0,
     $         1.42d0, 1.46d0, 1.38d0, 1.26d0, 1.37d0, 1.33d0,
     $         1.41d0, 1.22d0, 1.42d0, 1.23d0, 1.21d0, 1.40d0,
     $         
     $         0.128d0, 0.127d0, 0.142d0, 0.129d0, 0.136d0, 0.124d0,
     $         1.138d0, 0.148d0, 0.132d0, 0.125d0, 0.175d0, 0.123d0,
     $         0.146d0, 0.143d0, 0.145d0, 0.182d0, 0.135d0, 0.154d0,
     $         0.123d0, 0.129d0, 0.124d0, 0.162d0, 0.152d0, 0.142d0,
     $         0.168d0, 0.198d0, 0.186d0, 0.163d0, 0.126d0, 0.168d0,
     $         0.164d0, 0.134d0, 0.154d0, 0.128d0, 0.153d0, 0.145d0,
     $         
     $         0.0050d0, 0.020d0, 0.060d0, 0.056d0, 0.062d0, 0.062d0,
     $         0.0025d0, 0.001d0, 0.015d0, 0.070d0, 0.085d0, 0.011d0,
     $         0.0100d0, 0.002d0, 0.050d0, 0.080d0, 0.015d0, 0.057d0,
     $         0.0800d0, 0.015d0, 0.090d0, 0.065d0, 0.042d0, 0.067d0,
     $         0.0260d0, 0.030d0, 0.045d0, 0.052d0, 0.023d0, 0.051d0,
     $         0.0200d0, 0.012d0, 0.098d0, 0.056d0, 0.024d0, 0.090d0,
     $         
     $         4.88d0, 4.870d0,	4.855d0, 4.834d0, 4.592d0, 4.834d0,
     $         4.85d0, 4.865d0, 4.845d0, 4.875d0, 4.815d0, 4.875d0,
     $         4.89d0, 4.870d0, 4.860d0, 4.826d0, 4.723d0, 4.826d0,
     $         4.83d0, 4.950d0, 4.620d0, 4.952d0, 4.852d0, 4.952d0,
     $         4.81d0, 4.758d0, 4.762d0, 4.950d0, 4.703d0, 4.950d0,
     $         4.98d0, 4.780d0, 4.608d0, 4.628d0, 4.237d0, 4.862d0
     $         /),
     $         (/6,6,ne/))

          bf_alignment = reshape((/
     $         align_E,align_N,align_E+1,align_N/),
     $         (/2,2/))

          bf_grdpts_id = reshape((/
     $         1,1,2,3,0,0,
     $         1,1,2,3,0,0,
     $         2,2,2,3,0,0,
     $         3,3,3,3,0,0,
     $         0,0,0,0,0,0/),
     $         (/6,5/))

          bf_x_map = [2.5d0,  3.5d0, 4.5d0, 5.5d0, 6.5d0, 7.5d0 ]
          bf_y_map = [0.5d0, 0.75d0, 1.0d0, 1.25d0, 1.50d0]

          bf_nodes = reshape((/
     $         1.47d0, 1.28d0, 1.25d0, 1.43d0, -99.0d0, -99.0d0,
     $         1.41d0, 1.34d0, 1.31d0, 1.39d0, -99.0d0, -99.0d0,
     $         1.38d0, 1.26d0, 1.37d0, 1.33d0, -99.0d0, -99.0d0,
     $         1.42d0, 1.23d0, 1.21d0, 1.40d0, -99.0d0, -99.0d0,
     $        -99.0d0,-99.0d0,-99.0d0,-99.0d0, -99.0d0, -99.0d0,
     $         
     $         0.145d0, 0.182d0, 0.135d0, 0.154d0, -99.0d0, -99.0d0,
     $         0.124d0, 0.162d0, 0.152d0, 0.142d0, -99.0d0, -99.0d0,
     $         0.186d0, 0.163d0, 0.126d0, 0.168d0, -99.0d0, -99.0d0,
     $         0.154d0, 0.128d0, 0.153d0, 0.145d0, -99.0d0, -99.0d0,
     $         -99.0d0, -99.0d0, -99.0d0, -99.0d0, -99.0d0, -99.0d0,
     $         
     $         0.050d0,  0.08d0, 0.015d0, 0.057d0, -99.0d0, -99.0d0,
     $         0.090d0, 0.065d0, 0.042d0, 0.067d0, -99.0d0, -99.0d0,
     $         0.045d0, 0.052d0, 0.023d0, 0.051d0, -99.0d0, -99.0d0,
     $         0.098d0, 0.056d0, 0.024d0, 0.090d0, -99.0d0, -99.0d0,
     $         -99.0d0, -99.0d0, -99.0d0, -99.0d0, -99.0d0, -99.0d0,
     $         
     $         4.860d0, 4.826d0, 4.723d0, 4.826d0, -99.0d0, -99.0d0,
     $         4.620d0, 4.952d0, 4.852d0, 4.952d0, -99.0d0, -99.0d0,
     $         4.762d0, 4.950d0, 4.703d0, 4.950d0, -99.0d0, -99.0d0,
     $         4.608d0, 4.628d0, 4.237d0, 4.862d0, -99.0d0, -99.0d0,
     $         -99.0d0, -99.0d0, -99.0d0, -99.0d0, -99.0d0, -99.0d0
     $         /),
     $         (/6,5,ne/))

          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)


          ! verify that the nodes are correctly initialized
          test_loc = is_real_matrix3D_validated(
     $         interior_nodes(3:6,3:6,:),
     $         bf_nodes(1:4,1:4,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''bf_nodes and interior_nodes not matching'')'
          end if


          ! compute the reference fluxes
          call compute_fluxes_at_the_edges_2ndorder(
     $         interior_nodes,dx,dy,
     $         s_x_L0,s_x_L1,s_x_R1,s_x_R0,
     $         s_y_L0,s_y_L1,s_y_R1,s_y_R0,
     $         p_model,
     $         flux_x_ref,flux_y_ref)


          ! computation of the fluxes using the buffer layer
          i_min = 1
          i_max = 3
          j_min = 3

          call bc_operators_openbc_normal_used%compute_fluxes_x_for_bc_y_edge(
     $         bf_alignment,
     $         bf_grdpts_id,
     $         bf_nodes,
     $         interior_nodes,
     $         dx,dy,
     $         s_y_R1,s_y_R0,
     $         p_model,
     $         i_min, i_max, j_min,
     $         [.true.,.true.],
     $         bf_flux_x)

          ! validation
          test_loc = is_real_matrix3D_validated(
     $         bf_flux_x(1:3,3:4,:),
     $         flux_x_ref(3:5,5:6,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''flux N left failed'')'
          end if


          !N_edge: right missing grid-points
          !============================================================
          ! in this test, we want to compare the computation of the
          ! fluxes by the buffer layer where only a limited number of
          ! grid-points are available (and so require to extract nodes
          ! from the interior domain) with the reference function
          ! computing the fluxes directly from the interior nodes
          !        _________________
          !       |     |     |     | ____ north buffer layer
          !       |_____|_____|_____|______
          !             |/////////////////|__________ interior_domain
          !             |//|__|__|__|__|//|
          !             |//|           |//|
          !             |//|           |//|
          !
          !initialize the nodes
          bf_alignment = reshape((/
     $         align_W-1,align_N,align_W,align_N/),
     $         (/2,2/))

          bf_grdpts_id = reshape((/
     $         0,0,3,2,1,1,
     $         0,0,3,2,1,1,
     $         0,0,3,2,2,2,
     $         0,0,3,3,3,3,
     $         0,0,0,0,0,0/),
     $         (/6,5/))

          bf_x_map = [-1.5d0, -0.5d0, 0.5d0, 1.5d0, 2.5d0, 3.5d0]
          bf_y_map = [0.5d0, 0.75d0, 1.0d0, 1.25d0, 1.50d0]

          bf_nodes = reshape((/
     $         (((-99.0d0,i=1,6),j=1,5),k=1,ne)/),
     $         (/6,5,ne/))
          
          bf_nodes(3:6,1:4,:) = interior_nodes(1:4,3:6,:)

          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)


          ! computation of the fluxes using the buffer layer
          i_min = 5
          i_max = 7
          j_min = 3

          call bc_operators_openbc_normal_used%compute_fluxes_x_for_bc_y_edge(
     $         bf_alignment,
     $         bf_grdpts_id,
     $         bf_nodes,
     $         interior_nodes,
     $         dx,dy,
     $         s_y_R1,s_y_R0,
     $         p_model,
     $         i_min, i_max, j_min,
     $         [.true.,.true.],
     $         bf_flux_x)

          ! validation
          test_loc = is_real_matrix3D_validated(
     $         bf_flux_x(5:7,3:4,:),
     $         flux_x_ref(3:5,5:6,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''flux N right failed'')'
          end if


          !S_edge: left missing grid-points
          !============================================================
          ! in this test, we want to compare the computation of the
          ! fluxes by the buffer layer where only a limited number of
          ! grid-points are available (and so require to extract nodes
          ! from the interior domain) with the reference function
          ! computing the fluxes directly from the interior nodes
          !                             
          !    |//|           |//|__________ interior_domain
          !    |//|__ __ __ __|//|
          !    |//|__|__|__|__|//|
          !    |/////////////////|________
          !             |     |     |     |__south buffer layer
          !             |_____|_____|_____|
          !
          !initialize the nodes
          bf_alignment = reshape((/
     $         align_E,align_S,align_E+1,align_S/),
     $         (/2,2/))

          bf_grdpts_id = reshape((/
     $         0,0,0,0,0,0,
     $         3,3,3,3,0,0,
     $         2,2,2,3,0,0,
     $         1,1,2,3,0,0,
     $         1,1,2,3,0,0/),
     $         (/6,5/))

          bf_x_map = [2.5d0,  3.5d0, 4.5d0, 5.5d0, 6.5d0, 7.5d0 ]
          bf_y_map = [0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.25d0]

          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)

          bf_nodes(1:4,2:5,:) = interior_nodes(3:6,1:4,:)

          ! computation of the fluxes using the buffer layer
          i_min = 1
          i_max = 3
          j_min = 2

          call bc_operators_openbc_normal_used%compute_fluxes_x_for_bc_y_edge(
     $         bf_alignment,
     $         bf_grdpts_id,
     $         bf_nodes,
     $         interior_nodes,
     $         dx,dy,
     $         s_y_L0,s_y_L1,
     $         p_model,
     $         i_min, i_max, j_min,
     $         [.true.,.true.],
     $         bf_flux_x)

          ! validation
          test_loc = is_real_matrix3D_validated(
     $         bf_flux_x(1:3,2:3,:),
     $         flux_x_ref(3:5,1:2,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''flux S left failed'')'
          end if


          !S_edge: right missing grid-points
          !============================================================
          ! in this test, we want to compare the computation of the
          ! fluxes by the buffer layer where only a limited number of
          ! grid-points are available (and so require to extract nodes
          ! from the interior domain) with the reference function
          ! computing the fluxes directly from the interior nodes
          !                                   
          !          |//|           |//|__________ interior_domain
          !          |//|__ __ __ __|//|
          !          |//|__|__|__|__|//|
          !     _____|/////////////////|
          !    |     |     |     |________________south buffer layer
          !    |_____|_____|_____|
          !          
          !initialize the nodes
          bf_alignment = reshape((/
     $         align_W-1,align_S,align_W,align_S/),
     $         (/2,2/))

          bf_grdpts_id = reshape((/
     $         0,0,0,0,0,0,
     $         0,0,3,3,3,3,
     $         0,0,3,2,2,2,
     $         0,0,3,2,1,1,
     $         0,0,3,2,1,1/),
     $         (/6,5/))

          bf_x_map = [-1.5d0, -0.5d0, 0.5d0, 1.5d0, 2.5d0,  3.5d0]
          bf_y_map = [0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.25d0]

          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)

          bf_nodes(3:6,2:5,:) = interior_nodes(1:4,1:4,:)

          ! computation of the fluxes using the buffer layer
          i_min = 5
          i_max = 7
          j_min = 2

          call bc_operators_openbc_normal_used%compute_fluxes_x_for_bc_y_edge(
     $         bf_alignment,
     $         bf_grdpts_id,
     $         bf_nodes,
     $         interior_nodes,
     $         dx,dy,
     $         s_y_L0,s_y_L1,
     $         p_model,
     $         i_min, i_max, j_min,
     $         [.true.,.true.],
     $         bf_flux_x)

          ! validation
          test_loc = is_real_matrix3D_validated(
     $         bf_flux_x(5:7,2:3,:),
     $         flux_x_ref(3:5,1:2,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''flux S right failed'')'
          end if

        end function test_compute_fluxes_x_for_bc_y_edge


        function test_compute_fluxes_y_for_bc_x_edge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bc_operators)                 :: bc_operators_openbc_normal_used
          type(pmodel_eq)                    :: p_model
          real(rkind)                        :: t          
          type(sd_operators_x_oneside_L0)    :: s_x_L0
          type(sd_operators_x_oneside_L1)    :: s_x_L1
          type(sd_operators_x_oneside_R1)    :: s_x_R1
          type(sd_operators_x_oneside_R0)    :: s_x_R0
          type(sd_operators_y_oneside_L0)    :: s_y_L0
          type(sd_operators_y_oneside_L1)    :: s_y_L1
          type(sd_operators_y_oneside_R1)    :: s_y_R1
          type(sd_operators_y_oneside_R0)    :: s_y_R0  
          real(rkind)                        :: dx
          real(rkind)                        :: dy
          integer(ikind)                     :: i_min
          integer(ikind)                     :: j_min
          integer(ikind)                     :: j_max

          real(rkind)   , dimension(nx)         :: interior_x_map
          real(rkind)   , dimension(ny)         :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne)   :: interior_nodes
          real(rkind)   , dimension(nx+1,ny,ne) :: flux_x_ref
          real(rkind)   , dimension(nx,ny+1,ne) :: flux_y_ref
          
          integer(ikind), dimension(2,2)    :: bf_alignment
          integer       , dimension(5,6)    :: bf_grdpts_id
          real(rkind)   , dimension(5)      :: bf_x_map
          real(rkind)   , dimension(6)      :: bf_y_map
          real(rkind)   , dimension(5,6,ne) :: bf_nodes
          real(rkind)   , dimension(5,7,ne) :: bf_flux_y


          test_validated = .true.


          !E_edge: left missing grid-points
          !============================================================
          ! in this test, we want to compare the computation of the
          ! fluxes by the buffer layer where only a limited number of
          ! grid-points are available (and so require to extract nodes
          ! from the interior domain) with the reference function
          ! computing the fluxes directly from the interior nodes
          !                 _____      
          !                |     | 
          !     ___________|_____|
          !    ////////////|     |____ east buffer layer
          !    |__|__|__|//|_____|
          !             |//|     |
          !             |//|_____|
          !             |//|
          !             |//|__________ interior_domain
          !    
          !
          !initialize the nodes
          call p_model%initial_conditions%ini_far_field()

          t=0.0d0
          
          interior_x_map = [0.5d0, 1.5d0 , 2.5d0,  3.5d0, 4.5d0, 5.5d0 ]
          interior_y_map = [0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.25d0]

          interior_nodes = reshape((/
     $         1.48d0, 1.30d0, 1.35d0, 1.31d0, 1.43d0, 1.31d0,
     $         1.26d0, 1.45d0, 1.40d0, 1.29d0, 1.37d0, 1.41d0,
     $         1.46d0, 1.27d0, 1.47d0, 1.28d0, 1.25d0, 1.43d0,
     $         1.48d0, 1.26d0, 1.41d0, 1.34d0, 1.31d0, 1.39d0,
     $         1.42d0, 1.46d0, 1.38d0, 1.26d0, 1.37d0, 1.33d0,
     $         1.41d0, 1.22d0, 1.42d0, 1.23d0, 1.21d0, 1.40d0,
     $         
     $         0.128d0, 0.127d0, 0.142d0, 0.129d0, 0.136d0, 0.124d0,
     $         1.138d0, 0.148d0, 0.132d0, 0.125d0, 0.175d0, 0.123d0,
     $         0.146d0, 0.143d0, 0.145d0, 0.182d0, 0.135d0, 0.154d0,
     $         0.123d0, 0.129d0, 0.124d0, 0.162d0, 0.152d0, 0.142d0,
     $         0.168d0, 0.198d0, 0.186d0, 0.163d0, 0.126d0, 0.168d0,
     $         0.164d0, 0.134d0, 0.154d0, 0.128d0, 0.153d0, 0.145d0,
     $         
     $         0.0050d0, 0.020d0, 0.060d0, 0.056d0, 0.062d0, 0.062d0,
     $         0.0025d0, 0.001d0, 0.015d0, 0.070d0, 0.085d0, 0.011d0,
     $         0.0100d0, 0.002d0, 0.050d0, 0.080d0, 0.015d0, 0.057d0,
     $         0.0800d0, 0.015d0, 0.090d0, 0.065d0, 0.042d0, 0.067d0,
     $         0.0260d0, 0.030d0, 0.045d0, 0.052d0, 0.023d0, 0.051d0,
     $         0.0200d0, 0.012d0, 0.098d0, 0.056d0, 0.024d0, 0.090d0,
     $         
     $         4.88d0, 4.870d0,	4.855d0, 4.834d0, 4.592d0, 4.834d0,
     $         4.85d0, 4.865d0, 4.845d0, 4.875d0, 4.815d0, 4.875d0,
     $         4.89d0, 4.870d0, 4.860d0, 4.826d0, 4.723d0, 4.826d0,
     $         4.83d0, 4.950d0, 4.620d0, 4.952d0, 4.852d0, 4.952d0,
     $         4.81d0, 4.758d0, 4.762d0, 4.950d0, 4.703d0, 4.950d0,
     $         4.98d0, 4.780d0, 4.608d0, 4.628d0, 4.237d0, 4.862d0
     $         /),
     $         (/6,6,ne/))

          bf_alignment = reshape((/
     $         align_E,align_N,align_E,align_N+1/),
     $         (/2,2/))

          bf_grdpts_id = reshape((/
     $         1,1,2,3,0,
     $         1,1,2,3,0,
     $         2,2,2,3,0,
     $         3,3,3,3,0,
     $         0,0,0,0,0,
     $         0,0,0,0,0/),
     $         (/5,6/))

          bf_x_map = [3.5d0, 4.5d0, 5.5d0, 6.5d0, 7.5d0 ]
          bf_y_map = [0.5d0, 0.75d0, 1.0d0, 1.25d0, 1.50d0, 1.75d0]

          bf_nodes(1:4,1:4,:) = interior_nodes(3:6,3:6,:)

          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)

          ! compute the reference fluxes
          call compute_fluxes_at_the_edges_2ndorder(
     $         interior_nodes,dx,dy,
     $         s_x_L0,s_x_L1,s_x_R1,s_x_R0,
     $         s_y_L0,s_y_L1,s_y_R1,s_y_R0,
     $         p_model,
     $         flux_x_ref,flux_y_ref)


          ! computation of the fluxes using the buffer layer
          i_min = 3
          j_min = 1
          j_max = 3


          call bc_operators_openbc_normal_used%compute_fluxes_y_for_bc_x_edge(
     $         bf_alignment,
     $         bf_grdpts_id,
     $         bf_nodes,
     $         interior_nodes,
     $         dx,dy,
     $         s_x_R1,s_x_R0,
     $         p_model,
     $         i_min, j_min, j_max,
     $         [.true.,.true.],
     $         bf_flux_y)

          ! validation
          test_loc = is_real_matrix3D_validated(
     $         bf_flux_y(3:4,1:3,:),
     $         flux_y_ref(5:6,3:5,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''flux E left failed'')'
          end if


          !E_edge: right missing grid-points
          !============================================================
          ! in this test, we want to compare the computation of the
          ! fluxes by the buffer layer where only a limited number of
          ! grid-points are available (and so require to extract nodes
          ! from the interior domain) with the reference function
          ! computing the fluxes directly from the interior nodes
          !                 
          !             |//|
          !             |//|_____ 
          !             |//|     |
          !     ________|//|_____|
          !    |  |  |  |//|     |__________ interior_domain
          !    |///////////|_____|
          !                |     |
          !                |_____|
          !
          !initialize the nodes
          call p_model%initial_conditions%ini_far_field()

          t=0.0d0

          bf_alignment = reshape((/
     $         align_E,align_S-1,align_E,align_S/),
     $         (/2,2/))

          bf_grdpts_id = reshape((/
     $         0,0,0,0,0,
     $         0,0,0,0,0,
     $         3,3,3,3,0,
     $         2,2,2,3,0,
     $         1,1,2,3,0,
     $         1,1,2,3,0/),
     $         (/5,6/))

          bf_x_map = [3.5d0, 4.50d0, 5.5d0, 6.50d0, 7.5d0]
          bf_y_map = [0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.25d0]

          bf_nodes(1:4,3:6,:) = interior_nodes(3:6,1:4,:)

          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)

          ! compute the reference fluxes
          call compute_fluxes_at_the_edges_2ndorder(
     $         interior_nodes,dx,dy,
     $         s_x_L0,s_x_L1,s_x_R1,s_x_R0,
     $         s_y_L0,s_y_L1,s_y_R1,s_y_R0,
     $         p_model,
     $         flux_x_ref,flux_y_ref)


          ! computation of the fluxes using the buffer layer
          i_min = 3
          j_min = 5
          j_max = 7

          call bc_operators_openbc_normal_used%compute_fluxes_y_for_bc_x_edge(
     $         bf_alignment,
     $         bf_grdpts_id,
     $         bf_nodes,
     $         interior_nodes,
     $         dx,dy,
     $         s_x_R1,s_x_R0,
     $         p_model,
     $         i_min, j_min, j_max,
     $         [.true.,.true.],
     $         bf_flux_y)

          ! validation
          test_loc = is_real_matrix3D_validated(
     $         bf_flux_y(3:4,5:7,:),
     $         flux_y_ref(5:6,3:5,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''flux E right failed'')'
          end if


          !W_edge: left missing grid-points
          !============================================================
          ! in this test, we want to compare the computation of the
          ! fluxes by the buffer layer where only a limited number of
          ! grid-points are available (and so require to extract nodes
          ! from the interior domain) with the reference function
          ! computing the fluxes directly from the interior nodes
          !                 _____      
          !                |     | 
          !                |_____|_____________
          !                |     |/////////////
          !    W_edge______|_____|//|
          !                |     |//|  
          !                |_____|//|___ interior_domain
          !                      |//|
          !                      |//|
          !    
          !
          !initialize the nodes
          call p_model%initial_conditions%ini_far_field()

          bf_alignment = reshape((/
     $         align_W,align_N,align_W,align_N+1/),
     $         (/2,2/))

          bf_grdpts_id = reshape((/
     $         0,3,2,1,1,
     $         0,3,2,1,1,
     $         0,3,2,2,2,
     $         0,3,3,3,3,
     $         0,0,0,0,0,
     $         0,0,0,0,0/),
     $         (/5,6/))

          bf_x_map = [0.5d0, 1.5d0, 2.5d0, 3.5d0, 4.5d0]
          bf_y_map = [0.5d0, 0.75d0, 1.0d0, 1.25d0, 1.50d0, 1.75d0]

          bf_nodes(2:5,1:4,:) = interior_nodes(1:4,3:6,:)

          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)

          ! compute the reference fluxes
          call compute_fluxes_at_the_edges_2ndorder(
     $         interior_nodes,dx,dy,
     $         s_x_L0,s_x_L1,s_x_R1,s_x_R0,
     $         s_y_L0,s_y_L1,s_y_R1,s_y_R0,
     $         p_model,
     $         flux_x_ref,flux_y_ref)


          ! computation of the fluxes using the buffer layer
          i_min = 2
          j_min = 1
          j_max = 3

          call bc_operators_openbc_normal_used%compute_fluxes_y_for_bc_x_edge(
     $         bf_alignment,
     $         bf_grdpts_id,
     $         bf_nodes,
     $         interior_nodes,
     $         dx,dy,
     $         s_x_L0,s_x_L1,
     $         p_model,
     $         i_min, j_min, j_max,
     $         [.true.,.true.],
     $         bf_flux_y)

          ! validation
          test_loc = is_real_matrix3D_validated(
     $         bf_flux_y(2:3,1:3,:),
     $         flux_y_ref(1:2,3:5,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''flux W left failed'')'
          end if


          !W_edge: left missing grid-points
          !============================================================
          ! in this test, we want to compare the computation of the
          ! fluxes by the buffer layer where only a limited number of
          ! grid-points are available (and so require to extract nodes
          ! from the interior domain) with the reference function
          ! computing the fluxes directly from the interior nodes
          !                 _____      
          !                |     |//| 
          !                |_____|//|__________
          !                |     |/////////////
          !    W_edge______|_____|/////////////
          !                |     |
          !                |_____|
          !    
          !
          !initialize the nodes
          bf_alignment = reshape((/
     $         align_W,align_S-1,align_W,align_S/),
     $         (/2,2/))

          bf_grdpts_id = reshape((/
     $         0,0,0,0,0,
     $         0,0,0,0,0,
     $         0,3,3,3,3,
     $         0,3,2,2,2,
     $         0,3,2,1,1,
     $         0,3,2,1,1/),
     $         (/5,6/))

          bf_x_map = [0.5d0, 1.5d0, 2.5d0, 3.5d0, 4.5d0]
          bf_y_map = [-0.5d0, -0.25d0, 0.0d0, 0.25d0, 0.5d0, 0.75d0]

          bf_nodes(2:5,3:6,:) = interior_nodes(1:4,1:4,:)

          dx = bf_x_map(2) - bf_x_map(1)
          dy = bf_y_map(2) - bf_y_map(1)

          ! compute the reference fluxes
          call compute_fluxes_at_the_edges_2ndorder(
     $         interior_nodes,dx,dy,
     $         s_x_L0,s_x_L1,s_x_R1,s_x_R0,
     $         s_y_L0,s_y_L1,s_y_R1,s_y_R0,
     $         p_model,
     $         flux_x_ref,flux_y_ref)


          ! computation of the fluxes using the buffer layer
          i_min = 2
          j_min = 5
          j_max = 7

          call bc_operators_openbc_normal_used%compute_fluxes_y_for_bc_x_edge(
     $         bf_alignment,
     $         bf_grdpts_id,
     $         bf_nodes,
     $         interior_nodes,
     $         dx,dy,
     $         s_x_L0,s_x_L1,
     $         p_model,
     $         i_min, j_min, j_max,
     $         [.true.,.true.],
     $         bf_flux_y)

          ! validation
          test_loc = is_real_matrix3D_validated(
     $         bf_flux_y(2:3,5:7,:),
     $         flux_y_ref(1:2,3:5,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''flux W right failed'')'
          end if

        end function test_compute_fluxes_y_for_bc_x_edge


        function test_apply_bc_on_timedev_edge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          type(bc_operators)                 :: bc_operators_openbc_normal_used
          type(pmodel_eq)                    :: p_model
          real(rkind)                        :: t
          real(rkind), dimension(6,6,ne)     :: nodes
          real(rkind), dimension(6)          :: x_map
          real(rkind), dimension(6)          :: y_map
          real(rkind), dimension(7,6,ne)     :: flux_x
          real(rkind), dimension(6,7,ne)     :: flux_y
          type(sd_operators_x_oneside_L0)    :: s_x_L0
          type(sd_operators_x_oneside_L1)    :: s_x_L1
          type(sd_operators_x_oneside_R1)    :: s_x_R1
          type(sd_operators_x_oneside_R0)    :: s_x_R0
          type(sd_operators_y_oneside_L0)    :: s_y_L0
          type(sd_operators_y_oneside_L1)    :: s_y_L1
          type(sd_operators_y_oneside_R1)    :: s_y_R1
          type(sd_operators_y_oneside_R0)    :: s_y_R0          
          real(rkind)                        :: dx
          real(rkind)                        :: dy
          integer(ikind)                     :: i_min
          integer(ikind)                     :: i_max
          integer(ikind)                     :: j_min
          integer(ikind)                     :: j_max
          real(rkind), dimension(6,6,ne)     :: timedev

          real(rkind), dimension(nx+1,ny,ne) :: flux_x_ref
          real(rkind), dimension(nx,ny+1,ne) :: flux_y_ref
          real(rkind), dimension(2,2,ne)     :: timedev_ref

          integer(ikind) :: i,j
          integer        :: k

          integer(ikind), dimension(2,2)      :: bf_alignment
          integer       , dimension(1,1)      :: bf_grdpts_id
          real(rkind)   , dimension(nx,ny,ne) :: interior_nodes


          test_validated = .true.


          !N_edge
          !============================================================

          !initialize the nodes
          call p_model%initial_conditions%ini_far_field()

          t=0.0d0
          
          x_map = [0.5d0, 1.5d0 , 2.5d0,  3.5d0, 4.5d0, 5.5d0 ]
          y_map = [0.0d0, 0.25d0, 0.5d0, 0.75d0, 1.0d0, 1.25d0]
          nodes = reshape((/
     $         1.48d0, 1.30d0, 1.35d0, 1.31d0, 1.43d0, 1.31d0,
     $         1.26d0, 1.45d0,  1.4d0, 1.29d0, 1.37d0, 1.41d0,
     $         1.46d0, 1.27d0, 1.47d0, 1.28d0, 1.25d0, 1.43d0,
     $         1.48d0, 1.26d0, 1.41d0, 1.34d0, 1.31d0, 1.39d0,
     $         1.42d0, 1.46d0, 1.38d0, 1.26d0, 1.37d0, 1.33d0,
     $         1.41d0, 1.22d0, 1.42d0, 1.23d0, 1.21d0, 1.40d0,
     $         
     $         0.128d0, 0.127d0, 0.142d0, 0.129d0, 0.136d0, 0.124d0,
     $         1.138d0, 0.148d0, 0.132d0, 0.125d0, 0.175d0, 0.123d0,
     $         0.146d0, 0.143d0, 0.145d0, 0.182d0, 0.135d0, 0.154d0,
     $         0.123d0, 0.129d0, 0.124d0, 0.162d0, 0.152d0, 0.142d0,
     $         0.168d0, 0.198d0, 0.186d0, 0.163d0, 0.126d0, 0.168d0,
     $         0.164d0, 0.134d0, 0.154d0, 0.128d0, 0.153d0, 0.145d0,
     $         
     $         0.0050d0, 0.020d0, 0.060d0, 0.056d0, 0.062d0, 0.062d0,
     $         0.0025d0, 0.001d0, 0.015d0,  0.07d0, 0.085d0, 0.011d0,
     $         0.0100d0, 0.002d0, 0.050d0,  0.08d0, 0.015d0, 0.057d0,
     $           0.08d0, 0.015d0,  0.09d0, 0.065d0, 0.042d0, 0.067d0,
     $          0.026d0,  0.03d0, 0.045d0, 0.052d0, 0.023d0, 0.051d0,
     $           0.02d0, 0.012d0, 0.098d0, 0.056d0, 0.024d0, 0.090d0,
     $         
     $         4.88d0, 4.870d0,	4.855d0, 4.834d0, 4.592d0, 4.834d0,
     $         4.85d0, 4.865d0, 4.845d0, 4.875d0, 4.815d0, 4.875d0,
     $         4.89d0, 4.870d0, 4.860d0, 4.826d0, 4.723d0, 4.826d0,
     $         4.83d0,  4.95d0,  4.62d0, 4.952d0, 4.852d0, 4.952d0,
     $         4.81d0, 4.758d0, 4.762d0,  4.95d0, 4.703d0, 4.95d0,
     $         4.98d0, 4.780d0, 4.608d0, 4.628d0, 4.237d0, 4.862d0
     $         /),
     $         (/6,6,ne/))

          nodes(1:6,1:2,:) = reshape((/
     $         (((-99.0d0,i=1,6),j=1,2),k=1,ne)/),
     $         (/6,2,ne/))

          flux_x = reshape(
     $         (/ (((-99.0d0,i=1,7), j=1,6), k=1,ne) /),
     $         (/ 7,6,ne /))

          timedev = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,6), k=1,ne) /),
     $         (/ 6,6,ne /))

          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)

          !output
          i_min = 3
          i_max = 4
          j_min = 5
          call bc_operators_openbc_normal_used%apply_bc_on_timedev_N_edge(
     $         t,
     $         bf_alignment,
     $         bf_grdpts_id,
     $         x_map,y_map,nodes,
     $         interior_nodes,
     $         s_y_R1, s_y_R0,
     $         p_model,
     $         i_min, i_max, j_min,
     $         no_overlap,
     $         flux_x,
     $         timedev)

          call compute_fluxes_at_the_edges_2ndorder(
     $         nodes,dx,dy,
     $         s_x_L0,s_x_L1,s_x_R1,s_x_R0,
     $         s_y_L0,s_y_L1,s_y_R1,s_y_R0,
     $         p_model,
     $         flux_x_ref,flux_y_ref)

          do j=5,6
             do i=3,4
                timedev_ref(i-2,j-4,:) = compute_timedev_y_edge_local(
     $               t,
     $               x_map,
     $               y_map,
     $               nodes,
     $               p_model,
     $               gradient_y_y_oneside_R0,dy,
     $               i,j,
     $               flux_x_ref,dx,
     $               right)

             end do
          end do

          !validation
          test_loc = check_around_timedev(
     $         timedev,
     $         reshape((/3,5,4,6/),(/2,2/)),
     $         -99.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''check_around_timedev(N): failed'')'
          end if

          test_loc = is_real_matrix3D_validated(timedev_ref,timedev(3:4,5:6,:),detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''timedev(N) failed'')'
          end if


          !S_edge
          !============================================================
          !input
          call reflect_y(
     $         p_model,
     $         y_map,
     $         nodes,
     $         timedev_ref)

          flux_x = reshape(
     $         (/ (((-99.0d0,i=1,7), j=1,6), k=1,ne) /),
     $         (/ 7,6,ne /))

          timedev = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,6), k=1,ne) /),
     $         (/ 6,6,ne /))

          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)

          !output
          i_min = 3
          i_max = 4
          j_min = 1
          call bc_operators_openbc_normal_used%apply_bc_on_timedev_S_edge(
     $         t,
     $         bf_alignment,
     $         bf_grdpts_id,
     $         x_map, y_map, nodes,
     $         interior_nodes,
     $         s_y_L0, s_y_L1,
     $         p_model,
     $         i_min, i_max, j_min,
     $         no_overlap,
     $         flux_x,
     $         timedev)
          
          !validation
          test_loc = check_around_timedev(
     $         timedev,
     $         reshape((/3,1,4,2/),(/2,2/)),
     $         -99.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''check_around_timedev(S): failed'')'
          end if

          test_loc = is_real_matrix3D_validated(timedev_ref,timedev(3:4,1:2,:),detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''timedev(S) failed'')'
          end if


          !W_edge
          !============================================================
          !input
          call transpose_data(
     $         p_model,
     $         x_map,
     $         y_map,
     $         nodes,
     $         timedev_ref)

          flux_x = reshape(
     $         (/ (((-99.0d0,i=1,7), j=1,6), k=1,ne) /),
     $         (/ 7,6,ne /))

          flux_y = reshape(
     $         (/ (((-99.0d0,i=1,7), j=1,6), k=1,ne) /),
     $         (/ 6,7,ne /))

          timedev = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,6), k=1,ne) /),
     $         (/ 6,6,ne /))

          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)

          !output
          j_min = 3
          j_max = 4
          i_min = 1
          call bc_operators_openbc_normal_used%apply_bc_on_timedev_W_edge(
     $         t,
     $         bf_alignment,
     $         bf_grdpts_id,
     $         x_map, y_map, nodes,
     $         interior_nodes,
     $         s_x_L0, s_x_L1,
     $         p_model,
     $         i_min, j_min, j_max,
     $         no_overlap,
     $         flux_y,
     $         timedev)
          
          !validation
          test_loc = check_around_timedev(
     $         timedev,
     $         reshape((/1,3,2,4/),(/2,2/)),
     $         -99.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''check_around_timedev(W): failed'')'
          end if

          test_loc = is_real_matrix3D_validated(timedev_ref,timedev(1:2,3:4,:),detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''timedev(W) failed'')'
          end if


          !E_edge
          !============================================================
          !input
          call reflect_x(
     $         p_model,
     $         x_map,
     $         nodes,
     $         timedev_ref)

          flux_x = reshape(
     $         (/ (((-99.0d0,i=1,7), j=1,6), k=1,ne) /),
     $         (/ 7,6,ne /))

          flux_y = reshape(
     $         (/ (((-99.0d0,i=1,7), j=1,6), k=1,ne) /),
     $         (/ 6,7,ne /))

          timedev = reshape(
     $         (/ (((-99.0d0,i=1,6), j=1,6), k=1,ne) /),
     $         (/ 6,6,ne /))

          dx = x_map(2) - x_map(1)
          dy = y_map(2) - y_map(1)

          !output
          j_min = 3
          j_max = 4
          i_min = 5
          call bc_operators_openbc_normal_used%apply_bc_on_timedev_E_edge(
     $         t,
     $         bf_alignment,
     $         bf_grdpts_id,
     $         x_map, y_map, nodes,
     $         interior_nodes,
     $         s_x_R1, s_x_R0,
     $         p_model,
     $         i_min, j_min, j_max,
     $         no_overlap,
     $         flux_y,
     $         timedev)
          
          !validation
          test_loc = check_around_timedev(
     $         timedev,
     $         reshape((/5,3,6,4/),(/2,2/)),
     $         -99.0d0,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''check_around_timedev(E): failed'')'
          end if

          test_loc = is_real_matrix3D_validated(timedev_ref,timedev(5:6,3:4,:),detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''timedev(E) failed'')'
          end if

        end function test_apply_bc_on_timedev_edge


        subroutine check_inputs()

          implicit none

          type(pmodel_eq)            :: p_model
          real(rkind), dimension(ne) :: far_field

          call p_model%initial_conditions%ini_far_field()
          far_field = p_model%initial_conditions%get_far_field(0.0d0,1.0d0,1.0d0)

          if(.not.(
     $         is_real_vector_validated(
     $         far_field,
     $         [1.46510213931996d0,0.146510214d0,0.0d0,2.84673289046992d0],
     $         .true.).and.
     $         (nx.eq.6).and.
     $         (ny.eq.6))) then
             print '(''the test requires p_model%get_far_field(t,x,y)='')'
             print '(''[1.465102139d0,0.14651021d0,0.0d0,2.84673289d0]'')'
             print '()'
             print '(''T0 should be 0.95'')'
             print '(''flow_direction should be x-direction'')'
             print '(''ic_choice should be newgrdpt_test'')'
             print '()'
             stop ''
             
          end if

        end subroutine check_inputs


        function check_around_timedev(
     $     timedev,
     $     gen_coords,
     $     cst,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind)   , dimension(:,:,:), intent(in) :: timedev
          integer(ikind), dimension(2,2)  , intent(in) :: gen_coords
          real(rkind)                     , intent(in) :: cst
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated


          integer(ikind) :: i,j
          integer        :: k

          
          test_validated = .true.


          do k=1, size(timedev,3)
             do j=1, size(timedev,2)
                do i=1, size(timedev,1)

                   test_loc = is_real_validated(timedev(i,j,k),cst,.false.)

                   if(.not.test_loc) then
                      test_loc =
     $                     (i.ge.gen_coords(1,1)).and.
     $                     (i.le.gen_coords(1,2)).and.
     $                     (j.ge.gen_coords(2,1)).and.
     $                     (j.le.gen_coords(2,2))
                      
                   end if

                   if(detailled.and.(.not.test_loc)) then
                      print *, timedev(i,j,k)
                      print *, cst
                      print '(''['',3I4,'']'')', i,j,k
                   end if

                   test_validated = test_validated.and.test_loc

                end do
             end do
          end do

        end function check_around_timedev


        subroutine reflect_x(
     $     p_model,
     $     bf_x_map,
     $     bf_nodes,
     $     timedev_ref)

          implicit none

          type(pmodel_eq)               , intent(inout) :: p_model
          real(rkind), dimension(6)     , intent(inout) :: bf_x_map
          real(rkind), dimension(6,6,ne), intent(inout) :: bf_nodes
          real(rkind), dimension(2,2,ne), intent(inout) :: timedev_ref

          integer :: i
          integer :: i_r
          integer :: j
          integer :: k

          integer    , dimension(ne)     :: var_type
          real(rkind), dimension(ne)     :: far_field
          real(rkind), dimension(6)      :: bf_x_map_r
          real(rkind), dimension(6,6,ne) :: bf_nodes_r
          real(rkind), dimension(2,2,ne) :: timedev_ref_r


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
          do i=1, size(bf_x_map,1)

             i_r = size(bf_x_map)-(i-1)
             bf_x_map_r(i) = -bf_x_map(i_r)

          end do
          bf_x_map = bf_x_map_r


          !reflect the nodes
          do k=1, ne

             if(var_type(k).eq.vector_x) then

                do j=1, size(bf_nodes,2)
                   do i=1, size(bf_nodes,1)
                      i_r = size(bf_nodes,1)-(i-1)
                      bf_nodes_r(i,j,k) = - bf_nodes(i_r,j,k)
                   end do
                end do

             else

                do j=1, size(bf_nodes,2)
                   do i=1, size(bf_nodes,1)
                      i_r = size(bf_nodes,1)-(i-1)
                      bf_nodes_r(i,j,k) = bf_nodes(i_r,j,k)
                   end do
                end do

             end if             

          end do
          bf_nodes = bf_nodes_r


          !reflect the timedev_ref
          do k=1, ne
             if(var_type(k).eq.vector_x) then
                do j=1, size(timedev_ref,2)
                   do i=1, size(timedev_ref,1)
                      i_r = size(timedev_ref,1)-(i-1)
                      timedev_ref_r(i,j,k) = - timedev_ref(i_r,j,k)
                   end do
                end do
             else
                do j=1, size(timedev_ref,2)
                   do i=1, size(timedev_ref,1)
                      i_r = size(timedev_ref,1)-(i-1)
                      timedev_ref_r(i,j,k) = timedev_ref(i_r,j,k)
                   end do
                end do
             end if
          end do
          timedev_ref = timedev_ref_r

        end subroutine reflect_x


        subroutine reflect_y(
     $     p_model,
     $     bf_y_map,
     $     bf_nodes,
     $     timedev_ref)

          implicit none

          type(pmodel_eq)               , intent(inout) :: p_model
          real(rkind), dimension(6)     , intent(inout) :: bf_y_map
          real(rkind), dimension(6,6,ne), intent(inout) :: bf_nodes
          real(rkind), dimension(2,2,ne), intent(inout) :: timedev_ref

          integer :: i
          integer :: j
          integer :: j_r
          integer :: k

          integer    , dimension(ne)     :: var_type
          real(rkind), dimension(ne)     :: far_field
          real(rkind), dimension(6)      :: bf_y_map_r
          real(rkind), dimension(6,6,ne) :: bf_nodes_r
          real(rkind), dimension(2,2,ne) :: timedev_ref_r


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
          do j=1, size(bf_y_map,1)

             j_r = size(bf_y_map,1)-(j-1)
             bf_y_map_r(j) = -bf_y_map(j_r)

          end do
          bf_y_map = bf_y_map_r


          !reflect the nodes
          do k=1, ne

             if(var_type(k).eq.vector_y) then

                do j=1, size(bf_nodes,2)

                   j_r = size(bf_nodes,2)-(j-1)

                   do i=1, size(bf_nodes,1)
                      bf_nodes_r(i,j,k) = - bf_nodes(i,j_r,k)
                   end do

                end do

             else

                do j=1, size(bf_nodes,2)

                   j_r = size(bf_nodes,2)-(j-1)

                   do i=1, size(bf_nodes,1)
                      bf_nodes_r(i,j,k) = bf_nodes(i,j_r,k)
                   end do

                end do

             end if             

          end do
          bf_nodes = bf_nodes_r


          !reflect the timedev_ref
          do k=1, ne

             if(var_type(k).eq.vector_y) then 

                do j=1, size(timedev_ref,2)

                   j_r = size(timedev_ref,2)-(j-1)

                   do i=1, size(timedev_ref,1)
                      timedev_ref_r(i,j,k) = - timedev_ref(i,j_r,k)
                   end do
                end do

             else
                
                do j=1, size(timedev_ref,2)

                   j_r = size(timedev_ref,2)-(j-1)

                   do i=1, size(timedev_ref,1)
                      timedev_ref_r(i,j,k) = timedev_ref(i,j_r,k)
                   end do
                end do

             end if
          end do
          timedev_ref = timedev_ref_r

        end subroutine reflect_y


        subroutine transpose_data(
     $     p_model,
     $     bf_x_map, bf_y_map, bf_nodes,
     $     timedev_ref)

          implicit none

          type(pmodel_eq)               , intent(inout) :: p_model
          real(rkind), dimension(6)     , intent(inout) :: bf_x_map
          real(rkind), dimension(6)     , intent(inout) :: bf_y_map
          real(rkind), dimension(6,6,ne), intent(inout) :: bf_nodes
          real(rkind), dimension(2,2,ne), intent(inout) :: timedev_ref

          integer :: k

          real(rkind), dimension(ne)     :: far_field
          real(rkind)                    :: temp

          real(rkind), dimension(6)      :: bf_x_map_t
          real(rkind), dimension(6)      :: bf_y_map_t
          real(rkind), dimension(6,6,ne) :: bf_nodes_t
          real(rkind), dimension(2,2,ne) :: timedev_ref_t


          !transpose the far field
          far_field    = p_model%initial_conditions%get_far_field(0.0,1.0,0.0)
          temp         = far_field(2)
          far_field(2) = far_field(3)
          far_field(3) = temp
          call p_model%initial_conditions%set_far_field(far_field)


          !transpose the x_map and y_map
          bf_x_map_t = bf_y_map
          bf_y_map_t = bf_x_map
          bf_x_map   = bf_x_map_t
          bf_y_map   = bf_y_map_t


          !transpose the nodes
          do k=1,ne             
             bf_nodes_t(:,:,k) = transpose(bf_nodes(:,:,k))
          end do

          bf_nodes(:,:,1) = bf_nodes_t(:,:,1)
          bf_nodes(:,:,2) = bf_nodes_t(:,:,3)
          bf_nodes(:,:,3) = bf_nodes_t(:,:,2)
          bf_nodes(:,:,4) = bf_nodes_t(:,:,4)


          !transpose the newgrdpt_data_test
          do k=1,ne
             timedev_ref_t(:,:,k) = transpose(timedev_ref(:,:,k))
          end do

          timedev_ref(:,:,1) = timedev_ref_t(:,:,1)
          timedev_ref(:,:,2) = timedev_ref_t(:,:,3)
          timedev_ref(:,:,3) = timedev_ref_t(:,:,2)
          timedev_ref(:,:,4) = timedev_ref_t(:,:,4)          

        end subroutine transpose_data

      end program test_bc_operators_openbc_normal
