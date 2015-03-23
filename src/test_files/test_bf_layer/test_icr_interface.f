      program test_icr_interface

        use bf_interface_coords_class, only :
     $       bf_interface_coords

        use bf_sublayer_class, only :
     $       bf_sublayer

        use dim2d_parameters, only :
     $       cv_r

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated,
     $       is_int_vector_validated,
     $       is_int_matrix_validated

        use icr_interface_class, only :
     $       icr_interface

        use icr_path_class, only :
     $       icr_path

        use icr_path_chain_class, only :
     $       icr_path_chain

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W

        use parameters_constant, only :
     $       N,S,E,W,
     $       obc_eigenqties_bc,
     $       obc_eigenqties_lin

        use parameters_input, only :
     $       nx,ny,ne,
     $       obc_eigenqties_strategy

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_ini(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_ini: '',L1)', test_loc
        print '()'


        test_loc = test_stage(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_stage: '',L1)', test_loc
        print '()'


        test_loc = test_commit(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_commit: '',L1)', test_loc
        print '()'


        test_loc = test_finalize_domain_increase(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_finalize_domain_increase: '',L1)', test_loc
        print '()'


        contains


        function test_ini(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(icr_interface) :: icr_interface_used

          
          test_validated = .true.


          call icr_interface_used%ini()

          test_loc = icr_interface_used%paths%get_nb_paths().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''paths%get_nb_paths failed'')'
          end if


          test_loc = .not.associated(icr_interface_used%paths%get_head_path())
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''paths%get_head_path failed'')'
          end if


          test_loc = .not.associated(icr_interface_used%paths%get_tail_path())
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''paths%get_tail_path failed'')'
          end if


          test_loc = .not.associated(icr_interface_used%current_path)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''current_path failed'')'
          end if

        end function test_ini


        function test_stage(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(icr_interface)        :: icr_interface_used
          type(bf_interface_coords)  :: bf_interface_used
          real(rkind), dimension(nx) :: interior_x_map
          real(rkind), dimension(ny) :: interior_y_map

          
          type(icr_path)                :: icr_path_test
          type(icr_path_chain), pointer :: icr_path_tested
          logical :: test_loc


          test_validated = .true.


          !input
          !============================================================
          call bf_interface_used%ini(interior_x_map,interior_y_map)
          call icr_interface_used%ini()


          !output
          !============================================================
          !this first grid point should be saved in paths(1)
          call icr_interface_used%stage(
     $         [align_E,align_S+8],
     $         bf_interface_used)

          !this second grid-point should be saved in paths(2)
          !as it is too far from the first grid-point
          call icr_interface_used%stage(
     $         [align_E,align_S+15],
     $         bf_interface_used)

          !this third grid-point should be saved in paths(1)
          call icr_interface_used%stage(
     $         [align_E,align_S+7],
     $         bf_interface_used)

          !this fourth grid-point should be saved in paths(2)
          call icr_interface_used%stage(
     $         [align_E,align_S+16],
     $         bf_interface_used)


          !validation
          !============================================================
          !verify the number of activated paths
          !------------------------------------------------------------
          test_loc = icr_interface_used%paths%get_nb_paths().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test paths%nb_paths failed'')'
          end if


          !verify the content of the paths(1)
          !------------------------------------------------------------
          icr_path_tested => icr_interface_used%paths%get_head_path()
          call icr_path_test%ini()
          call icr_path_test%stage([align_E,align_S+8],bf_interface_used)
          call icr_path_test%stage([align_E,align_S+7],bf_interface_used)
          test_loc = is_path_validated(
     $         icr_path_test,
     $         icr_path_tested,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test paths(1) failed'')'
          end if


          !verify the content of the paths(2)
          !------------------------------------------------------------
          icr_path_tested => icr_interface_used%paths%get_head_path()
          icr_path_tested => icr_path_tested%get_next()
          call icr_path_test%ini()
          call icr_path_test%stage([align_E,align_S+15],bf_interface_used)
          call icr_path_test%stage([align_E,align_S+16],bf_interface_used)
          test_loc = is_path_validated(
     $         icr_path_test,
     $         icr_path_tested,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test paths(2) failed'')'
          end if

        end function test_stage


        function test_commit(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(icr_interface)              :: icr_interface_used
          type(bf_interface_coords)        :: bf_interface_used
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes0
          real(rkind), dimension(nx,ny,ne) :: interior_nodes1
          real(rkind)                      :: t
          real(rkind)                      :: dt
          type(pmodel_eq)                  :: p_model



          test_validated = .true.


          test_loc = test_commit_wo_merge(
     $         detailled,
     $         icr_interface_used,
     $         bf_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1,
     $         t,dt,
     $         p_model)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test_commit_wo_merge failed'')'
          end if


          test_loc = test_commit_w_merge(
     $         detailled,
     $         icr_interface_used,
     $         bf_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1,
     $         t,dt,
     $         p_model)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test_commit_w_merge failed'')'
          end if    

        end function test_commit


        function test_commit_wo_merge(
     $     detailled,
     $     icr_interface_used,
     $     bf_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1,
     $     t,dt,
     $     p_model)
     $     result(test_validated)

          implicit none

          logical                         , intent(in)    :: detailled
          type(icr_interface)             , intent(inout) :: icr_interface_used
          type(bf_interface_coords)       , intent(inout) :: bf_interface_used
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes1
          real(rkind)                     , intent(in)    :: t
          real(rkind)                     , intent(in)    :: dt
          type(pmodel_eq)                 , intent(in)    :: p_model
          logical                                         :: test_validated

          type(bf_sublayer), pointer :: bf_sublayer_ptr
          logical :: test_loc


          test_validated = .true.


          !input
          !============================================================
          call bf_interface_used%ini(interior_x_map,interior_y_map)
          call icr_interface_used%ini()


          !this first grid point should be saved in paths(1)
          call icr_interface_used%stage(
     $         [align_E,align_S+8],
     $         bf_interface_used)

          !this second grid-point should be saved in paths(2)
          !as it is too far from the first grid-point
          call icr_interface_used%stage(
     $         [align_E,align_S+15],
     $         bf_interface_used)

          !this third grid-point should be saved in paths(1)
          call icr_interface_used%stage(
     $         [align_E,align_S+7],
     $         bf_interface_used)

          !this fourth grid-point should be saved in paths(2)
          call icr_interface_used%stage(
     $         [align_E,align_S+16],
     $         bf_interface_used)


          !output
          !============================================================
          call icr_interface_used%commit(
     $         N,
     $         bf_interface_used,
     $         p_model,
     $         t,
     $         dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)

          call icr_interface_used%commit(
     $         E,
     $         bf_interface_used,
     $         p_model,
     $         t,
     $         dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)          


          !validation
          !============================================================
          !verify the number of sublayers in N,S,E,W
          !------------------------------------------------------------
          test_loc = bf_interface_used%mainlayer_pointers(N)%get_nb_sublayers().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test mainlayer_pointers(N)%get_nb_sublayers failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(S)%get_nb_sublayers().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test mainlayer_pointers(S)%get_nb_sublayers failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(E)%get_nb_sublayers().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test mainlayer_pointers(E)%get_nb_sublayers failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(W)%get_nb_sublayers().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test mainlayer_pointers(W)%get_nb_sublayers failed'')'
          end if


          !verify the configuration of the buffer layer E
          !------------------------------------------------------------
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
          
          !alignment (E,1)
          !............................................................
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%alignment,
     $         reshape((/align_E, align_S+7, align_E, align_S+8/),(/2,2/)))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment(E,1) failed'')'
          end if

          !grdpts_id (E,1)
          !............................................................
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%grdpts_id,
     $         reshape((/
     $           1,1,2,3,3,
     $           1,1,2,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,2,2,3,
     $           1,1,2,3,3/),
     $         (/5,6/)))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test grdpts_id(E,1) failed'')'
          end if


          bf_sublayer_ptr => bf_sublayer_ptr%get_next()
          
          !alignment (E,2)
          !............................................................
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%alignment,
     $         reshape((/align_E, align_S+15, align_E, align_S+16/),(/2,2/)))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment(E,2) failed'')'
          end if

          !grdpts_id (E,2)
          !............................................................
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%grdpts_id,
     $         reshape((/
     $           1,1,2,3,3,
     $           1,1,2,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,2,2,3,
     $           1,1,2,3,3/),
     $         (/5,6/)))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test grdpts_id(E,2) failed'')'
          end if
          
        end function test_commit_wo_merge


        function test_commit_w_merge(
     $     detailled,
     $     icr_interface_used,
     $     bf_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1,
     $     t,dt,
     $     p_model)
     $     result(test_validated)

          implicit none

          logical                         , intent(in)    :: detailled
          type(icr_interface)             , intent(inout) :: icr_interface_used
          type(bf_interface_coords)       , intent(inout) :: bf_interface_used
          real(rkind), dimension(nx)      , intent(in)    :: interior_x_map
          real(rkind), dimension(ny)      , intent(in)    :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes0
          real(rkind), dimension(nx,ny,ne), intent(in)    :: interior_nodes1
          real(rkind)                     , intent(in)    :: t
          real(rkind)                     , intent(in)    :: dt
          type(pmodel_eq)                 , intent(in)    :: p_model
          logical                                         :: test_validated

          type(bf_sublayer), pointer :: bf_sublayer_ptr
          logical :: test_loc


          test_validated = .true.


          !input
          !============================================================

          !this first grid point should be saved in paths(1)
          call icr_interface_used%stage(
     $         [align_E,align_S+9],
     $         bf_interface_used)

          !this second grid-point should be saved in paths(2)
          !as it is too far from the first grid-point
          call icr_interface_used%stage(
     $         [align_E,align_S+17],
     $         bf_interface_used)

          !this third grid-point should be saved in paths(1)
          call icr_interface_used%stage(
     $         [align_E,align_S+10],
     $         bf_interface_used)

          !this fourth grid-point should be saved in paths(3)
          call icr_interface_used%stage(
     $         [align_E,align_S+24],
     $         bf_interface_used)


          !output
          !============================================================
          call icr_interface_used%commit(
     $         N,
     $         bf_interface_used,
     $         p_model,
     $         t,
     $         dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)

          call icr_interface_used%commit(
     $         E,
     $         bf_interface_used,
     $         p_model,
     $         t,
     $         dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)          


          !validation
          !============================================================
          !verify the number of sublayers in N,S,E,W
          !------------------------------------------------------------
          test_loc = bf_interface_used%mainlayer_pointers(N)%get_nb_sublayers().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test mainlayer_pointers(N)%get_nb_sublayers failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(S)%get_nb_sublayers().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test mainlayer_pointers(S)%get_nb_sublayers failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(E)%get_nb_sublayers().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test mainlayer_pointers(E)%get_nb_sublayers failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(W)%get_nb_sublayers().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test mainlayer_pointers(W)%get_nb_sublayers failed'')'
          end if


          !verify the configuration of the buffer layer E
          !------------------------------------------------------------
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
          
          !alignment (E,1)
          !............................................................
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%alignment,
     $         reshape((/align_E, align_S+7, align_E, align_S+17/),(/2,2/)))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment(E,1) failed'')'
          end if

          !grdpts_id (E,1)
          !............................................................
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%grdpts_id,
     $         reshape((/
     $           1,1,2,3,3,
     $           1,1,2,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,2,2,3,
     $           1,1,2,3,3/),
     $         (/5,15/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test grdpts_id(E,1) failed'')'
          end if


          bf_sublayer_ptr => bf_sublayer_ptr%get_next()
          
          !alignment (E,2)
          !............................................................
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%alignment,
     $         reshape((/align_E, align_S+24, align_E, align_S+24/),(/2,2/)))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment(E,2) failed'')'
          end if

          !grdpts_id (E,2)
          !............................................................
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%grdpts_id,
     $         reshape((/
     $           1,1,2,3,3,
     $           1,1,2,2,3,
     $           1,1,1,2,3,
     $           1,1,2,2,3,
     $           1,1,2,3,3/),
     $         (/5,5/)))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test grdpts_id(E,2) failed'')'
          end if
          
        end function test_commit_w_merge



        function test_finalize_domain_increase(
     $     detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in)    :: detailled
          logical                :: test_validated

          type(icr_interface)              :: icr_interface_used
          type(bf_interface_coords)        :: bf_interface_used
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes0
          real(rkind), dimension(nx,ny,ne) :: interior_nodes1
          real(rkind)                      :: t
          real(rkind)                      :: dt
          type(pmodel_eq)                  :: p_model


          type(bf_sublayer), pointer :: bf_sublayer_ptr
          logical :: test_loc


          test_validated = .true.


          !input
          !============================================================
          call bf_interface_used%ini(interior_x_map,interior_y_map)
          call icr_interface_used%ini()


          !this first grid point should be saved in paths(1)
          call icr_interface_used%stage(
     $         [align_E,align_S+8],
     $         bf_interface_used)

          !this second grid-point should be saved in paths(2)
          !as it is too far from the first grid-point
          call icr_interface_used%stage(
     $         [align_E,align_S+15],
     $         bf_interface_used)

          !this third grid-point should be saved in paths(1)
          call icr_interface_used%stage(
     $         [align_E,align_S+7],
     $         bf_interface_used)

          !this fourth grid-point should be saved in paths(2)
          call icr_interface_used%stage(
     $         [align_E,align_S+16],
     $         bf_interface_used)


          !output
          !============================================================
          call icr_interface_used%finalize_domain_increase(
     $         bf_interface_used,
     $         p_model,
     $         t,
     $         dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)


          !validation
          !============================================================
          !verify the number of sublayers in N,S,E,W
          !------------------------------------------------------------
          test_loc = bf_interface_used%mainlayer_pointers(N)%get_nb_sublayers().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test mainlayer_pointers(N)%get_nb_sublayers failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(S)%get_nb_sublayers().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test mainlayer_pointers(S)%get_nb_sublayers failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(E)%get_nb_sublayers().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test mainlayer_pointers(E)%get_nb_sublayers failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(W)%get_nb_sublayers().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test mainlayer_pointers(W)%get_nb_sublayers failed'')'
          end if


          !verify the configuration of the buffer layer E
          !------------------------------------------------------------
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
          
          !alignment (E,1)
          !............................................................
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%alignment,
     $         reshape((/align_E, align_S+7, align_E, align_S+8/),(/2,2/)))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment(E,1) failed'')'
          end if

          !grdpts_id (E,1)
          !............................................................
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%grdpts_id,
     $         reshape((/
     $           1,1,2,3,3,
     $           1,1,2,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,2,2,3,
     $           1,1,2,3,3/),
     $         (/5,6/)))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test grdpts_id(E,1) failed'')'
          end if


          bf_sublayer_ptr => bf_sublayer_ptr%get_next()
          
          !alignment (E,2)
          !............................................................
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%alignment,
     $         reshape((/align_E, align_S+15, align_E, align_S+16/),(/2,2/)))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment(E,2) failed'')'
          end if

          !grdpts_id (E,2)
          !............................................................
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%grdpts_id,
     $         reshape((/
     $           1,1,2,3,3,
     $           1,1,2,2,3,
     $           1,1,1,2,3,
     $           1,1,1,2,3,
     $           1,1,2,2,3,
     $           1,1,2,3,3/),
     $         (/5,6/)))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test grdpts_id(E,2) failed'')'
          end if
          
        end function test_finalize_domain_increase


        function is_path_validated(path1,path2,detailled)
     $     result(test_validated)

          implicit none

          class(icr_path), intent(in) :: path1
          class(icr_path), intent(in) :: path2
          logical        , intent(in) :: detailled
          logical                     :: test_validated

          logical :: test_loc


          test_validated = .true.


          !nb_pts
          test_loc = path1%nb_pts.eq.path2%nb_pts
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_pts: '',I2,''->'',I2)', 
     $            path1%nb_pts,
     $            path2%nb_pts
          end if

          if(path1%nb_pts.gt.0) then

             !mainlayer_id
             test_loc = path1%mainlayer_id.eq.path2%mainlayer_id
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''mainlayer_id: '',I2,''->'',I2)', 
     $               path1%mainlayer_id,
     $               path2%mainlayer_id
             end if
   
             !alignment
             test_loc = is_int_matrix_validated(
     $            path1%alignment,
     $            path2%alignment,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''alignment failed'')'
             end if
   
             !matching sublayer
             if(
     $            associated(path1%matching_sublayer).and.
     $            associated(path2%matching_sublayer)) then
             
                test_loc = associated(
     $               path1%matching_sublayer,
     $               path2%matching_sublayer)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''matching_sublayer failed'')'
                end if
   
             end if

             !pts
             if(path1%nb_pts.eq.path2%nb_pts) then
                test_loc = is_int_matrix_validated(
     $               path1%pts(:,1:path1%nb_pts),
     $               path2%pts(:,1:path2%nb_pts),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''pts failed'')'
                end if
             end if

          end if

        end function is_path_validated
        

        subroutine check_inputs()

          implicit none

          type(pmodel_eq)            :: p_model
          real(rkind), dimension(ne) :: far_field

          if(ne.ne.4) then
             stop 'the test requires ne=4: DIM2D model'
          end if

          if(.not.is_real_validated(cv_r,2.5d0,.false.)) then
             stop 'the test requires c_v_r=2.5'
          end if

          call p_model%initial_conditions%ini_far_field()

          far_field = p_model%get_far_field(0.0d0,1.0d0,1.0d0)

          if(.not.is_real_vector_validated(
     $         far_field,
     $         [1.46510213931996d0,0.146510214d0,0.0d0,2.84673289046992d0],
     $         .true.)) then
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

      end program test_icr_interface
