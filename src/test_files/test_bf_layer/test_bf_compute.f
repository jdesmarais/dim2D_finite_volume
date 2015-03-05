      !test the function td_operators%compute_time_dev_nopt()
      program test_bf_compute

        use bc_operators_class, only :
     $       bc_operators

        use bf_compute_class, only :
     $       bf_compute

        use check_data_module, only :
     $       is_real_vector_validated,
     $       is_real_matrix3D_validated,
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       no_pt,
     $       bc_pt,
     $       bc_interior_pt,
     $       interior_pt,
     $       
     $       align_S,align_E,
     $       
     $       SW_corner_type,
     $       S_edge_type,
     $       SE_corner_type,
     $       W_edge_type,
     $       E_edge_type,
     $       NW_corner_type,
     $       N_edge_type,
     $       NE_corner_type,
     $       
     $       no_overlap

        use parameters_constant, only :
     $       bc_timedev_choice,
     $       hedstrom_xy_choice
        
        use parameters_input, only :
     $       x_min, x_max,
     $       y_min, y_max,
     $       nx,ny,ne,
     $       bc_size,
     $       bc_choice,
     $       bc_N_type_choice,
     $       bc_S_type_choice,
     $       bc_E_type_choice,
     $       bc_W_type_choice
        
        use parameters_kind, only :
     $       ikind,
     $       rkind
        
        use pmodel_eq_class, only :
     $       pmodel_eq

        use rk3tvd_steps_module, only :
     $       compute_1st_step,
     $       compute_1st_step_nopt
        
        use sd_operators_class, only :
     $       sd_operators
        
        use td_operators_class, only :
     $       td_operators
        
        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        test_loc = test_does_previous_timestep_exist(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_does_previous_timestep_exist: '',L1)', test_loc
        print '()'

        
        test_loc = test_allocate_tables(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_allocate_tables: '',L1)', test_loc
        print '()'


        test_loc = test_deallocate_tables(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_deallocate_tables: '',L1)', test_loc
        print '()'


        test_loc = test_compute_time_dev(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_time_dev: '',L1)', test_loc
        print '()'

        
        test_loc = test_compute_integration_step(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_integration_step: '',L1)', test_loc
        print '()'


        test_loc = test_get_grdpts_id_part(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_grdpts_id_part: '',L1)', test_loc
        print '()'


        contains


        function test_does_previous_timestep_exist(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_compute) :: bf_compute_used          

          test_validated = .not.(bf_compute_used%does_previous_timestep_exist())

          if(detailled.and.(.not.test_validated)) then
             print '(''test failed: '',L1,'' -> '',L1)', .true., .false.
          end if

        end function test_does_previous_timestep_exist


        function test_allocate_tables(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_compute)                      :: bf_compute_used
          integer    , parameter                :: size_x=6
          integer    , parameter                :: size_y=10
          integer    , dimension(2,2)           :: alignment
          real(rkind)                           :: dx
          real(rkind)                           :: dy
          real(rkind), dimension(size_x)        :: x_map
          real(rkind), dimension(size_y)        :: y_map
          integer    , dimension(size_x,size_y) :: grdpts_id

          integer :: i,j

          test_validated = .true.

          !input
          alignment = reshape((/
     $         align_E,
     $         align_S+1,
     $         align_E+size_x-2*bc_size,
     $         align_S+1+size_y-2*bc_size/),
     $         (/2,2/))

          dx = (x_max-x_min)/size_x
          dy = (y_max-y_min)/size_y

          x_map = (/ (x_min+(i-1)*dx,i=1,size_x) /)
          y_map = (/ (y_min+(j-1)*dy,j=1,size_y) /)
          grdpts_id = reshape(
     $         (/ ((interior_pt, i=1,size_x),j=1,size_y) /),
     $         (/size_x,size_y/))

          !output
          call bf_compute_used%allocate_tables(
     $         size_x,size_y,
     $         alignment,
     $         x_map,
     $         y_map,
     $         grdpts_id)

          !validation
          !alignment
          test_loc = is_int_matrix_validated(
     $         bf_compute_used%alignment_tmp,alignment,detailled)
          test_validated = test_validated.and.test_loc

          !x_map
          test_loc = is_real_vector_validated(
     $         bf_compute_used%x_map_tmp,x_map,detailled)
          test_validated = test_validated.and.test_loc

          !y_map
          test_loc = is_real_vector_validated(
     $         bf_compute_used%y_map_tmp,y_map,detailled)
          test_validated = test_validated.and.test_loc

          !grdpts_id
          test_loc = is_int_matrix_validated(
     $         bf_compute_used%grdpts_id_tmp,grdpts_id,detailled)
          test_validated = test_validated.and.test_loc

          !nodes
          test_loc = size_x.eq.size(bf_compute_used%nodes_tmp,1)
          test_validated = test_validated.and.test_loc

          test_loc = size_y.eq.size(bf_compute_used%nodes_tmp,2)
          test_validated = test_validated.and.test_loc

          test_loc = ne.eq.size(bf_compute_used%nodes_tmp,3)
          test_validated = test_validated.and.test_loc
          
          !timedev
          test_loc = size_x.eq.size(bf_compute_used%time_dev,1)
          test_validated = test_validated.and.test_loc

          test_loc = size_y.eq.size(bf_compute_used%time_dev,2)
          test_validated = test_validated.and.test_loc

          test_loc = ne.eq.size(bf_compute_used%time_dev,3)
          test_validated = test_validated.and.test_loc

        end function test_allocate_tables

      
        function test_deallocate_tables(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_compute) :: bf_compute_used
          integer    , parameter                :: size_x=6
          integer    , parameter                :: size_y=10
          integer    , dimension(2,2)           :: alignment
          real(rkind), dimension(size_x)        :: x_map
          real(rkind), dimension(size_y)        :: y_map
          integer    , dimension(size_x,size_y) :: grdpts_id

          logical :: test_loc

          test_validated = .true.

          call bf_compute_used%deallocate_tables()

          test_loc = .not.allocated(bf_compute_used%alignment_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%grdpts_id_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%x_map_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%y_map_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%nodes_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%time_dev)
          test_validated = test_validated.and.test_loc
          
          if(detailled.and.(.not.test_validated)) then
             print '(''test_deallocate_tables(1) failed'')'
          end if

          call bf_compute_used%allocate_tables(
     $         size_x,size_y,
     $         alignment,
     $         x_map,
     $         y_map,
     $         grdpts_id)

          call bf_compute_used%deallocate_tables()

          test_loc = .not.allocated(bf_compute_used%alignment_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%grdpts_id_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%x_map_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%y_map_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%nodes_tmp)
          test_validated = test_validated.and.test_loc
          test_loc = .not.allocated(bf_compute_used%time_dev)
          test_validated = test_validated.and.test_loc

        end function test_deallocate_tables


        function test_compute_time_dev(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_compute)                                 :: bf_compute_used
          type(td_operators)                               :: td_operators_used
          real(rkind)                                      :: t
          real(rkind)   , dimension(nx,ny,ne)              :: nodes
          real(rkind)                                      :: dx
          real(rkind)                                      :: dy
          real(rkind)   , dimension(nx)                    :: x_map
          real(rkind)   , dimension(ny)                    :: y_map
          type(sd_operators)                               :: s
          type(pmodel_eq)                                  :: p_model
          type(bc_operators)                               :: bc_used
          integer(ikind), dimension(2,2)                   :: bf_alignment
          integer       , dimension(nx,ny)                 :: grdpts_id
          real(rkind)   , dimension(nx,ny,ne)              :: interior_nodes
          integer       , dimension(:,:)     , allocatable :: bc_sections
          integer(ikind), dimension(2)                     :: x_borders
          integer(ikind), dimension(2)                     :: y_borders

          real(rkind)   , dimension(nx,ny,ne) :: timedev_test

          integer(ikind) :: i,j
          integer        :: k

          !input
          t = 0.0d0

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

          dx = 0.1d0
          dy = 0.2d0

          x_map = (/ (x_min+(i-1)*dx,i=1,nx) /)
          y_map = (/ (y_min+(j-1)*dy,j=1,ny) /)

          bf_alignment = reshape((/
     $         align_E,
     $         align_S+1,
     $         align_E+nx-2*bc_size,
     $         align_S+1+ny-2*bc_size/),
     $         (/2,2/))

          grdpts_id = reshape(
     $         (/ ((interior_pt, i=1,nx),j=1,ny) /),
     $         (/nx,ny/))

          interior_nodes = reshape((/
     $         (((100*(k-1) + 10*(j-1) + (i-1),i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          allocate(bc_sections(5,8))

          bc_sections = reshape((/
     $         SW_corner_type, 1           , 1           , no_overlap, no_overlap,
     $         S_edge_type   , bc_size+1   , 1           , nx-bc_size, no_overlap,
     $         SE_corner_type, nx-bc_size+1, 1           , no_overlap, no_overlap,
     $         W_edge_type   , 1           , bc_size+1   , ny-bc_size, no_overlap,
     $         E_edge_type   , nx-bc_size+1, bc_size+1   , ny-bc_size, no_overlap,
     $         NW_corner_type, 1           , ny-bc_size+1, no_overlap, no_overlap,
     $         N_edge_type   , bc_size+1   , ny-bc_size+1, nx-bc_size, no_overlap,
     $         NE_corner_type, nx-bc_size+1, ny-bc_size+1, no_overlap, no_overlap/),
     $         (/5,8/))

          x_borders = [bc_size+1,nx-bc_size]
          y_borders = [bc_size+1,ny-bc_size]

          call bf_compute_used%allocate_tables(
     $         nx,ny,
     $         bf_alignment,
     $         x_map,
     $         y_map,
     $         grdpts_id)

          !output
          call bf_compute_used%compute_time_dev(
     $         td_operators_used,
     $         t, nodes, x_map, y_map,
     $         s,
     $         p_model,bc_used,
     $         bf_alignment,
     $         grdpts_id,
     $         interior_nodes,
     $         bc_sections,
     $         x_borders, y_borders)

          !validation
          timedev_test = td_operators_used%compute_time_dev(
     $         t, nodes, x_map, y_map,
     $         s,p_model,bc_used)

          test_validated = is_real_matrix3D_validated(
     $         bf_compute_used%time_dev,
     $         timedev_test,
     $         detailled)

        end function test_compute_time_dev


      function test_compute_integration_step(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_compute)                                 :: bf_compute_used
          type(td_operators)                               :: td_operators_used
          real(rkind)                                      :: t
          real(rkind)   , dimension(nx,ny,ne)              :: nodes
          real(rkind)                                      :: dx
          real(rkind)                                      :: dy
          real(rkind)   , dimension(nx)                    :: x_map
          real(rkind)   , dimension(ny)                    :: y_map
          type(sd_operators)                               :: s
          type(pmodel_eq)                                  :: p_model
          type(bc_operators)                               :: bc_used
          integer(ikind), dimension(2,2)                   :: bf_alignment
          integer       , dimension(nx,ny)                 :: grdpts_id
          real(rkind)   , dimension(nx,ny,ne)              :: interior_nodes
          real(rkind)   , dimension(nx,ny,ne)              :: interior_nodes_tmp
          integer       , dimension(:,:)     , allocatable :: bc_sections
          integer(ikind), dimension(2)                     :: x_borders
          integer(ikind), dimension(2)                     :: y_borders

          real(rkind)                         :: dt
          real(rkind)   , dimension(nx,ny,ne) :: interior_timedev

          integer(ikind) :: i,j
          logical        :: test_loc


          test_validated = .true.


          !input
          t  = 0.0d0
          dt = 0.01d0

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

          dx = 0.1d0
          dy = 0.2d0

          x_map = (/ (x_min+(i-1)*dx,i=1,nx) /)
          y_map = (/ (y_min+(j-1)*dy,j=1,ny) /)

          bf_alignment = reshape((/
     $         align_E,
     $         align_S+1,
     $         align_E+nx-2*bc_size,
     $         align_S+1+ny-2*bc_size/),
     $         (/2,2/))

          grdpts_id = reshape(
     $         (/ ((interior_pt, i=1,nx),j=1,ny) /),
     $         (/nx,ny/))

          interior_nodes = nodes

          allocate(bc_sections(5,8))

          bc_sections = reshape((/
     $         SW_corner_type, 1           , 1           , no_overlap, no_overlap,
     $         S_edge_type   , bc_size+1   , 1           , nx-bc_size, no_overlap,
     $         SE_corner_type, nx-bc_size+1, 1           , no_overlap, no_overlap,
     $         W_edge_type   , 1           , bc_size+1   , ny-bc_size, no_overlap,
     $         E_edge_type   , nx-bc_size+1, bc_size+1   , ny-bc_size, no_overlap,
     $         NW_corner_type, 1           , ny-bc_size+1, no_overlap, no_overlap,
     $         N_edge_type   , bc_size+1   , ny-bc_size+1, nx-bc_size, no_overlap,
     $         NE_corner_type, nx-bc_size+1, ny-bc_size+1, no_overlap, no_overlap/),
     $         (/5,8/))

          x_borders = [bc_size+1,nx-bc_size]
          y_borders = [bc_size+1,ny-bc_size]

          call bf_compute_used%allocate_tables(
     $         nx,ny,
     $         bf_alignment,
     $         x_map,
     $         y_map,
     $         grdpts_id)

          call bf_compute_used%compute_time_dev(
     $         td_operators_used,
     $         t, nodes, x_map, y_map,
     $         s,
     $         p_model,bc_used,
     $         bf_alignment,
     $         grdpts_id,
     $         interior_nodes,
     $         bc_sections,
     $         x_borders, y_borders)

          !output
          call bf_compute_used%compute_integration_step(
     $         grdpts_id, nodes, dt,
     $         x_borders, y_borders,
     $         compute_1st_step_nopt)

          !validation
          interior_timedev = td_operators_used%compute_time_dev(
     $         t, interior_nodes, x_map, y_map,
     $         s,p_model,bc_used)

          call compute_1st_step(
     $         interior_nodes,dt,
     $         interior_nodes_tmp,
     $         interior_timedev,
     $         [1,nx],
     $         [1,ny])

          test_loc = is_real_matrix3D_validated(
     $         bf_compute_used%nodes_tmp,
     $         interior_nodes_tmp,
     $         detailled)
          test_validated = test_validated.and.test_loc

          test_loc = is_real_matrix3D_validated(
     $         nodes,
     $         interior_nodes,
     $         detailled)
          test_validated = test_validated.and.test_loc
          

        end function test_compute_integration_step

        
        function test_get_grdpts_id_part(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_compute)               :: bf_compute_used
          integer(ikind), dimension(2,2) :: bf_alignment
          integer(ikind), dimension(2,2) :: gen_coords
          integer       , dimension(6,3) :: tmp_grdpts_id

          integer(ikind) :: i,j


          test_validated = .true.


          bf_alignment = reshape(
     $         (/align_E,align_S+5,align_E+2,align_S+6/),
     $         (/2,2/))

          gen_coords = reshape(
     $         (/align_E-2,align_S+6,align_E+3,align_S+8/),
     $         (/2,2/))

          tmp_grdpts_id = reshape(
     $         (/ ((no_pt,i=1,6),j=1,3) /),
     $         (/6,3/))


          !test w/o allocation
          !------------------------------------------------------------
          call bf_compute_used%get_grdpts_id_part(
     $         tmp_grdpts_id,
     $         gen_coords)

          test_loc = is_int_matrix_validated(
     $         tmp_grdpts_id,
     $         reshape(
     $         (/ ((no_pt,i=1,6),j=1,3) /),
     $         (/6,3/)),
     $         detailled)
          test_validated = test_validated.and.test_loc


          !test w/ allocation
          !------------------------------------------------------------
          !input
          call bf_compute_used%allocate_tables(
     $         7,6,
     $         bf_alignment,
     $         (/((i-1)*0.1d0,i=1,7)/),
     $         (/((j-1)*0.1d0,j=1,6)/),
     $         reshape((/
     $             1,1,2,3,3,3,3,
     $             1,1,2,2,2,2,3,
     $             1,1,1,1,1,2,3,
     $             1,1,1,1,1,2,3,
     $             1,1,2,2,2,2,3,
     $             1,1,2,3,3,3,3/),
     $         (/7,6/)))
          
          !output
          call bf_compute_used%get_grdpts_id_part(
     $         tmp_grdpts_id,
     $         gen_coords)

          !validation
          test_loc = is_int_matrix_validated(
     $         tmp_grdpts_id,
     $         reshape((/
     $            1,1,1,1,1,2,
     $            1,1,2,2,2,2,
     $            1,1,2,3,3,3/),
     $            (/6,3/)),
     $         detailled)
          test_validated = test_validated.and.test_loc

        end function test_get_grdpts_id_part


        subroutine check_inputs()

          implicit none

          logical :: test_parameters

          test_parameters = bc_choice.eq.hedstrom_xy_choice

          test_parameters = test_parameters.and.(bc_N_type_choice.eq.bc_timedev_choice)
          test_parameters = test_parameters.and.(bc_S_type_choice.eq.bc_timedev_choice)
          test_parameters = test_parameters.and.(bc_E_type_choice.eq.bc_timedev_choice)
          test_parameters = test_parameters.and.(bc_W_type_choice.eq.bc_timedev_choice)
          test_parameters = test_parameters.and.(nx.eq.6)
          test_parameters = test_parameters.and.(ny.eq.6)
          
        end subroutine check_inputs

      end program test_bf_compute
