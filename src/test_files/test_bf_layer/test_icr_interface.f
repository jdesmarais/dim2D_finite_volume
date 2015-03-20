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

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W

        use parameters_constant, only :
     $       N,S,E,
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


        test_loc = test_finalize_domain_increase(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_finalize_domain_extension: '',L1)', test_loc
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

          test_loc = icr_interface_used%current_path_is_head_path
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''current_path_is_head_path failed'')'
          end if


          test_loc = icr_interface_used%head_path%nb_pts.eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''head_path%nb_pts failed'')'
          end if


          test_loc = icr_interface_used%current_path%nb_pts.eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''current_path%nb_pts failed'')'
          end if

        end function test_ini


        function test_stage(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(icr_interface)              :: icr_interface_used
          type(bf_interface_coords)        :: bf_interface_used
          type(pmodel_eq)                  :: p_model
          real(rkind)                      :: t
          real(rkind)                      :: dt
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes0
          real(rkind), dimension(nx,ny,ne) :: interior_nodes1

          
          type(icr_path)             :: icr_path_test
          type(bf_sublayer), pointer :: new_sublayer
          logical :: test_loc


          test_validated = .true.


          !input
          !============================================================
          call bf_interface_used%ini(interior_x_map,interior_y_map)
          call icr_interface_used%ini()


          !output
          !============================================================
          !this first grid point should be saved in the head_path
          call icr_interface_used%stage(
     $         [align_E,align_S+8],
     $         bf_interface_used,
     $         p_model,
     $         t,
     $         dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)

          !this second grid-point should be saved in the current_path
          !as it is too far from the first grid-point
          call icr_interface_used%stage(
     $         [align_E,align_S+15],
     $         bf_interface_used,
     $         p_model,
     $         t,
     $         dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)

          !this third grid-point should force the current_path to
          !commit its changes to the domain extension. It should be
          !saved in a reinitialized current_path, once its update
          !operations have been commit
          call icr_interface_used%stage(
     $         [align_E,align_S+7],
     $         bf_interface_used,
     $         p_model,
     $         t,
     $         dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)

          !this fourth grid-point should force the current path to
          !commit its changes to the domain extension. As the third
          !grid point is closed to the first grid-point in the
          !head_path, the head_path and the current_path should be
          !merged for the update
          call icr_interface_used%stage(
     $         [align_E,align_S+16],
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
          !verify the content of the head path
          !------------------------------------------------------------
          call icr_path_test%ini()
          call icr_path_test%stage([align_E,align_S+16],bf_interface_used)
          test_loc = is_path_validated(
     $         icr_path_test,
     $         icr_interface_used%head_path,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test head_path failed'')'
          end if


          !verify the content of the current_path
          !------------------------------------------------------------
          call icr_path_test%ini()
          test_loc = is_path_validated(
     $         icr_path_test,
     $         icr_interface_used%current_path,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test head_path failed'')'
          end if


          !verify the content of the mainlayer_pointers(E)
          !------------------------------------------------------------
          !nb_sublayers
          !............................................................
          test_loc = bf_interface_used%mainlayer_pointers(E)%get_nb_sublayers().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_sublayers failed'')'
          end if

          !first sublayer: alignment
          !............................................................
          new_sublayer => bf_interface_used%mainlayer_pointers(E)%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         new_sublayer%alignment,
     $         reshape((/
     $            align_E,align_S+7,align_E,align_S+8/),
     $            (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment(1) failed'')'
          end if

          !first sublayer: grdpts_id
          !............................................................
          test_loc = is_int_matrix_validated(
     $         new_sublayer%grdpts_id,
     $         reshape((/
     $            1,1,2,3,3,
     $            1,1,2,2,3,
     $            1,1,1,2,3,
     $            1,1,1,2,3,
     $            1,1,2,2,3,
     $            1,1,2,3,3/),
     $            (/5,6/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment(1) failed'')'
          end if

          !second sublayer: alignment
          !............................................................
          new_sublayer => new_sublayer%get_next()
          test_loc = is_int_matrix_validated(
     $         new_sublayer%alignment,
     $         reshape((/
     $            align_E,align_S+15,align_E,align_S+15/),
     $            (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment(2) failed'')'
          end if

          !second sublayer: grdpts_id
          !............................................................
          test_loc = is_int_matrix_validated(
     $         new_sublayer%grdpts_id,
     $         reshape((/
     $            1,1,2,3,3,
     $            1,1,2,2,3,
     $            1,1,1,2,3,
     $            1,1,2,2,3,
     $            1,1,2,3,3/),
     $            (/5,5/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment(2) failed'')'
          end if

        end function test_stage


        function test_finalize_domain_increase(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          logical :: test_loc


          test_validated = .true.


          test_loc = test_finalize_domain_increase_test1(detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test_finalize_domain_increase_test1 failed'')'
          end if

          test_loc = test_finalize_domain_increase_test2(detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test_finalize_domain_increase_test2 failed'')'
          end if

        end function test_finalize_domain_increase


        function test_finalize_domain_increase_test1(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(icr_interface)              :: icr_interface_used
          type(bf_interface_coords)        :: bf_interface_used
          type(pmodel_eq)                  :: p_model
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes0
          real(rkind), dimension(nx,ny,ne) :: interior_nodes1

          real(rkind) :: t
          real(rkind) :: dt

          type(bf_sublayer), pointer :: new_sublayer

          real(rkind), dimension(ne) :: test_newgrdpt
          

          test_validated = .true.

          
          !input
          !============================================================
          call ini_for_test_finalize(
     $         1,
     $         icr_interface_used,
     $         bf_interface_used,
     $         p_model,
     $         t,dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1,
     $         test_newgrdpt)

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
          new_sublayer => bf_interface_used%mainlayer_pointers(E)%get_head_sublayer()
          
          !test the number of sublayers
          !------------------------------------------------------------
          test_loc = bf_interface_used%mainlayer_pointers(E)%get_nb_sublayers().eq.1
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_sublayers failed'')'
          end if

          !test the alignment
          !------------------------------------------------------------
          test_loc = is_int_matrix_validated(
     $         new_sublayer%alignment,
     $         reshape((/
     $            align_E, align_S+8, align_E, align_S+14/),
     $            (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment failed'')'
          end if

          !test the configuration of the grid points
          !------------------------------------------------------------
          test_loc = is_int_matrix_validated(
     $         new_sublayer%grdpts_id,
     $         reshape((/
     $            1,1,2,3,3,
     $            1,1,2,2,3,
     $            1,1,1,2,3,
     $            1,1,2,2,3,
     $            1,1,2,3,3,
     $         
     $            1,1,2,3,0,
     $         
     $            1,1,2,3,3,
     $            1,1,2,2,3,
     $            1,1,1,2,3,
     $            1,1,2,2,3,
     $            1,1,2,3,3/),
     $            (/5,11/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment failed'')'
          end if


          !test some new grid points
          !------------------------------------------------------------
          test_loc = is_real_vector_validated(
     $         new_sublayer%nodes(5,3,:),
     $         test_newgrdpt,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test newgrdpt(1) failed'')'
          end if
          
           test_loc = is_real_vector_validated(
     $         new_sublayer%nodes(5,9,:),
     $         test_newgrdpt,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test newgrdpt(2) failed'')'
          end if


       end function test_finalize_domain_increase_test1


       function test_finalize_domain_increase_test2(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(icr_interface)              :: icr_interface_used
          type(bf_interface_coords)        :: bf_interface_used
          type(pmodel_eq)                  :: p_model
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes0
          real(rkind), dimension(nx,ny,ne) :: interior_nodes1

          real(rkind) :: t
          real(rkind) :: dt

          type(bf_sublayer), pointer :: new_sublayer

          real(rkind), dimension(ne) :: test_newgrdpt
          

          test_validated = .true.

          
          !input
          !============================================================
          call ini_for_test_finalize(
     $         2,
     $         icr_interface_used,
     $         bf_interface_used,
     $         p_model,
     $         t,dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1,
     $         test_newgrdpt)

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
          !test the number of sublayers
          !------------------------------------------------------------
          test_loc = bf_interface_used%mainlayer_pointers(E)%get_nb_sublayers().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_sublayers failed'')'
          end if

          !test of the first sublayer
          !------------------------------------------------------------
          new_sublayer => bf_interface_used%mainlayer_pointers(E)%get_head_sublayer()

          !test the alignment
          !............................................................
          test_loc = is_int_matrix_validated(
     $         new_sublayer%alignment,
     $         reshape((/
     $            align_E, align_S+8, align_E, align_S+8/),
     $            (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment failed'')'
          end if

          !test the configuration of the grid points
          !............................................................
          test_loc = is_int_matrix_validated(
     $         new_sublayer%grdpts_id,
     $         reshape((/
     $            1,1,2,3,3,
     $            1,1,2,2,3,
     $            1,1,1,2,3,
     $            1,1,2,2,3,
     $            1,1,2,3,3/),
     $            (/5,5/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment failed'')'
          end if


          !test some new grid points
          !............................................................
          test_loc = is_real_vector_validated(
     $         new_sublayer%nodes(5,3,:),
     $         test_newgrdpt,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test newgrdpt failed'')'
          end if


          !test of the second sublayer
          !------------------------------------------------------------
          new_sublayer => new_sublayer%get_next()

          !test the alignment
          !............................................................
          test_loc = is_int_matrix_validated(
     $         new_sublayer%alignment,
     $         reshape((/
     $            align_E, align_S+15, align_E, align_S+15/),
     $            (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment failed'')'
          end if

          !test the configuration of the grid points
          !............................................................
          test_loc = is_int_matrix_validated(
     $         new_sublayer%grdpts_id,
     $         reshape((/
     $            1,1,2,3,3,
     $            1,1,2,2,3,
     $            1,1,1,2,3,
     $            1,1,2,2,3,
     $            1,1,2,3,3/),
     $            (/5,5/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment failed'')'
          end if


          !test some new grid points
          !............................................................
          test_loc = is_real_vector_validated(
     $         new_sublayer%nodes(5,3,:),
     $         test_newgrdpt,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test newgrdpt failed'')'
          end if

       end function test_finalize_domain_increase_test2


       subroutine ini_for_test_finalize(
     $     test_id,
     $     icr_interface_used,
     $     bf_interface_used,
     $     p_model,
     $     t,dt,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes0,
     $     interior_nodes1,
     $     test_newgrdpt)

         implicit none

         integer                         , intent(in)  :: test_id
         type(icr_interface)             , intent(out) :: icr_interface_used
         type(bf_interface_coords)       , intent(out) :: bf_interface_used
         type(pmodel_eq)                 , intent(out) :: p_model
         real(rkind)                     , intent(out) :: t
         real(rkind)                     , intent(out) :: dt
         real(rkind), dimension(nx)      , intent(out) :: interior_x_map
         real(rkind), dimension(ny)      , intent(out) :: interior_y_map
         real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes0
         real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes1
         real(rkind), dimension(ne)      , intent(out) :: test_newgrdpt

         integer(ikind) :: i,j


         call icr_interface_used%ini()
         
         call bf_interface_used%ini(interior_x_map,interior_y_map)
         
         call p_model%initial_conditions%ini_far_field()

         t  = 0.0d0
         dt = 0.25d0

         interior_x_map = (/ ((i-1)*1.00d0,i=1,nx) /)
         interior_y_map = (/ ((j-1)*0.25d0,j=1,ny) /)

         interior_nodes0(align_E-1:align_E+1,
     $        align_S+7:align_S+9,:)
     $        = reshape((/
     $        1.48d0, 1.30d0, 1.35d0,
     $        1.26d0, 1.45d0, 1.40d0,
     $        1.46d0, 1.27d0, 1.47d0,
     $        
     $        0.128d0, 0.127d0, 0.142d0,
     $        1.138d0, 0.148d0, 0.132d0,
     $        0.146d0, 0.143d0, 0.145d0,
     $        
     $        0.0050d0, 0.020d0, 0.060d0,
     $        0.0025d0, 0.001d0, 0.015d0,
     $        0.0100d0, 0.002d0, 0.050d0,
     $        
     $        4.88d0, 4.870d0, 4.855d0,
     $        4.85d0, 4.865d0, 4.845d0,
     $        4.89d0, 4.870d0, 4.860d0/),
     $        (/3,3,ne/))               
               
         interior_nodes1(align_E-1:align_E+1,
     $        align_S+7:align_S+9,:)
     $        = reshape((/
     $        1.50d0, 1.455d0, 1.48d0,
     $        1.20d0, 1.350d0, 1.25d0,
     $        1.49d0, 1.250d0, 1.40d0,
     $        
     $        0.128d0, 0.450d0, 0.135d0,
     $        0.148d0, 0.150d0, 0.122d0,
     $        0.142d0, 1.152d0, 0.236d0,
     $        
     $        0.006d0, 0.0600d0, 0.020d0,
     $        0.000d0, 0.0028d0, 0.035d0,
     $        0.020d0, 0.0030d0, 0.040d0,
     $        
     $        4.876d0, 4.825d0, 4.862d0,
     $        4.890d0, 4.871d0, 4.892d0,
     $        4.865d0, 4.757d0, 4.895d0/),
     $        (/3,3,ne/))         

         select case(obc_eigenqties_strategy)

            case(obc_eigenqties_bc)
               test_newgrdpt = [
     $              1.22383078395524d0,
     $              0.39531842478603d0,
     $              -0.21050217290879d0,
     $              4.19684181018688d0]
               
            case(obc_eigenqties_lin)
               test_newgrdpt = [
     $              1.21167555521982d0,
     $              0.35901827468671d0,
     $              -0.20475732388332d0,
     $              4.20002914561351d0]
               
            case default
               print '(''test_bf_newgrdpt_prim'')'
               print '(''test_sym_compute_newgrdpt_x'')'
               print '(''obc_eigenqties_strategy not recognized'')'
               stop ''
               
          end select


          select case(test_id)
            case(1)

               !interior nodes
               !----------------------------------------
               interior_nodes0(align_E-1:align_E+1,
     $              align_S+13:align_S+15,:)
     $         = interior_nodes0(align_E-1:align_E+1,
     $              align_S+7:align_S+9,:)
               
               interior_nodes1(align_E-1:align_E+1,
     $              align_S+13:align_S+15,:)
     $         = interior_nodes1(align_E-1:align_E+1,
     $              align_S+7:align_S+9,:)

               !head path
               !----------------------------------------
               call icr_interface_used%head_path%ini()
               call icr_interface_used%head_path%stage(
     $              [align_E,align_S+8],
     $              bf_interface_used)
               
               !currrent path
               !----------------------------------------
               call icr_interface_used%current_path%ini()
               call icr_interface_used%current_path%stage(
     $              [align_E,align_S+14],
     $              bf_interface_used)

            case(2)

               !interior nodes
               !----------------------------------------
               interior_nodes0(align_E-1:align_E+1,
     $              align_S+14:align_S+16,:)
     $         = interior_nodes0(align_E-1:align_E+1,
     $              align_S+7:align_S+9,:)
               
               interior_nodes1(align_E-1:align_E+1,
     $              align_S+14:align_S+16,:)
     $         = interior_nodes1(align_E-1:align_E+1,
     $              align_S+7:align_S+9,:)

               !head path
               !----------------------------------------
               call icr_interface_used%head_path%ini()
               call icr_interface_used%head_path%stage(
     $              [align_E,align_S+8],
     $              bf_interface_used)
               
               !currrent path
               !----------------------------------------
               call icr_interface_used%current_path%ini()
               call icr_interface_used%current_path%stage(
     $              [align_E,align_S+15],
     $              bf_interface_used)
               
         end select

        end subroutine ini_for_test_finalize


        function is_path_validated(path1,path2,detailled)
     $     result(test_validated)

          implicit none

          type(icr_path), intent(in) :: path1
          type(icr_path), intent(in) :: path2
          logical       , intent(in) :: detailled
          logical                    :: test_validated

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
