      program test_mainlayer_interface_grdpts_id_update

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated,
     $       is_int_matrix_validated

        use dim2d_parameters, only :
     $       cv_r

        use mainlayer_interface_grdpts_id_update_class, only :
     $       mainlayer_interface_grdpts_id_update

        use parameters_bf_layer, only :
     $       interior_pt,
     $       bc_interior_pt,
     $       bc_pt,
     $       no_pt,
     $       
     $       align_N,
     $       align_S,
     $       align_E,
     $       
     $       NE_interface_type,
     $       SE_interface_type

        use parameters_constant, only :
     $       obc_eigenqties_bc,
     $       obc_eigenqties_lin,
     $       N,S,E

        use parameters_input, only :
     $       nx,ny,ne,
     $       obc_eigenqties_strategy,
     $       debug_geometry_update

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


        test_loc = test_update_grdpts_id_around(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_update_grdpts_id_around_new_interior_pt: '',L1)', test_loc
        print '()'


        test_loc = test_detect_bc_interior_pt_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_detect_bc_interior_pt_crenel: '',L1)', test_loc
        print '()'


        test_loc = test_curb_bc_interior_pt_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_can_interior_crenel_be_curbed: '',L1)', test_loc
        print '()'


        test_loc = test_finalize_for_bc_interior_pt_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_finalize_for_bc_interior_pt_crenel: '',L1)', test_loc
        print '()'


        test_loc = test_detect_and_curb_bc_pt_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_detect_and_curb_bc_pt_crenel: '',L1)', test_loc
        print '()'


        test_loc = test_finalize_for_bc_pt_crenel_local(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_finalize_for_bc_pt_crenel_local: '',L1)', test_loc
        print '()'


        test_loc = test_finalize_for_bc_pt_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_finalize_for_bc_pt_crenel: '',L1)', test_loc
        print '()'


        test_loc = test_update_grdpts_id_in_bf_layer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_update_grdpts_id_in_bf_layer: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_update_grdpts_id_in_bf_layer(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(mainlayer_interface_grdpts_id_update) :: mainlayer_interface_used
          type(pmodel_eq)                            :: p_model
          real(rkind)                                :: t
          real(rkind)                                :: dt
          real(rkind), dimension(nx)                 :: interior_x_map
          real(rkind), dimension(ny)                 :: interior_y_map
          real(rkind), dimension(nx,ny,ne)           :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)           :: interior_nodes1
          type(bf_sublayer)                          :: bf_sublayer_ptr
          integer       , dimension(6,5)             :: test_grdpts_id
          real(rkind)   , dimension(ne)              :: test_newgrdpt

          logical :: test_loc
          integer :: k
          

          test_validated = .true.

          
          !input
          !------------------------------------------------------------
          call ini_for_test_update_grdpts_id_around_edge(
     $         t,dt,
     $         p_model,
     $         bf_sublayer_ptr,
     $         test_newgrdpt,
     $         test_grdpts_id)


          !output
          !------------------------------------------------------------
          call mainlayer_interface_used%update_grdpts_id_in_bf_layer(
     $         p_model,
     $         t,dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1,
     $         bf_sublayer_ptr,
     $         reshape((/align_E+1,align_S+10/),(/2,1/)))


          !validation
          !------------------------------------------------------------
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%grdpts_id,
     $         test_grdpts_id,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test grdpts_id failed'')'
             print '(''grdpts_id:'')'
             do k=1,5
                print '(6I2)', bf_sublayer_ptr%grdpts_id(:,5-(k-1))
             end do
             print '()'

             print '(''test_grdpts_id:'')'
             do k=1,5
                print '(6I2)', test_grdpts_id(:,5-(k-1))
             end do
             print '()'
          end if

          if(.not.debug_geometry_update) then
             test_loc = is_real_vector_validated(
     $            bf_sublayer_ptr%nodes(6,3,:),
     $            test_newgrdpt,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test newgrdpt failed'')'
             end if
          end if

        end function test_update_grdpts_id_in_bf_layer

        
        function test_finalize_for_bc_pt_crenel(detailled)
     $     result(test_validated)
          
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(mainlayer_interface_grdpts_id_update) :: mainlayer_interface_used
          type(bf_sublayer), pointer                 :: N_sublayer_ptr
          type(bf_sublayer), pointer                 :: S_sublayer_ptr
          type(bf_sublayer), pointer                 :: E_sublayer_ptr
          integer(ikind), dimension(2)               :: gen_coords
          integer(ikind), dimension(2)               :: match_table
          integer       , dimension(5,ny)            :: test_grdpts_id

          logical :: test_loc
          integer :: k


          test_validated = .true.


          call ini_mainlayer_interface_for_bc_pt_crenel(
     $         mainlayer_interface_used,
     $         N_sublayer_ptr,
     $         S_sublayer_ptr,
     $         E_sublayer_ptr)


          do k=1,4

             !input
             call get_test_param_bc_pt_crenel_nl(
     $            k,
     $            N_sublayer_ptr,
     $            S_sublayer_ptr,
     $            E_sublayer_ptr,
     $            gen_coords,
     $            match_table,
     $            test_grdpts_id)

             !output
             call mainlayer_interface_used%finalize_for_bc_pt_crenel(
     $            E_sublayer_ptr,
     $            gen_coords,
     $            match_table)

             !validation
             test_loc = is_int_matrix_validated(
     $            E_sublayer_ptr%grdpts_id,
     $            test_grdpts_id,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')'
             end if

          end do

        end function test_finalize_for_bc_pt_crenel


        function test_finalize_for_bc_pt_crenel_local(detailled)
     $     result(test_validated)
          
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(mainlayer_interface_grdpts_id_update) :: mainlayer_interface_used
          type(bf_sublayer), pointer                 :: N_sublayer_ptr
          type(bf_sublayer), pointer                 :: S_sublayer_ptr
          type(bf_sublayer), pointer                 :: E_sublayer_ptr
          integer(ikind), dimension(2)               :: gen_coords
          integer(ikind), dimension(2)               :: match_table
          integer       , dimension(5,ny)            :: test_grdpts_id

          logical :: test_loc
          integer :: k


          test_validated = .true.


          call ini_mainlayer_interface_for_bc_pt_crenel(
     $         mainlayer_interface_used,
     $         N_sublayer_ptr,
     $         S_sublayer_ptr,
     $         E_sublayer_ptr)


          do k=1,4

             !input
             call get_test_param_bc_pt_crenel(
     $            k,
     $            N_sublayer_ptr,
     $            S_sublayer_ptr,
     $            E_sublayer_ptr,
     $            gen_coords,
     $            match_table,
     $            test_grdpts_id)

             !output
             call mainlayer_interface_used%finalize_for_bc_pt_crenel_local(
     $            E_sublayer_ptr,
     $            gen_coords,
     $            match_table)

             !validation
             test_loc = is_int_matrix_validated(
     $            E_sublayer_ptr%grdpts_id,
     $            test_grdpts_id,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')'
             end if

          end do

        end function test_finalize_for_bc_pt_crenel_local


        function test_detect_and_curb_bc_pt_crenel(detailled)
     $     result(test_validated)
          
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(mainlayer_interface_grdpts_id_update) :: mainlayer_interface_used
          type(bf_sublayer), pointer                 :: N_sublayer_ptr
          type(bf_sublayer), pointer                 :: S_sublayer_ptr
          type(bf_sublayer), pointer                 :: E_sublayer_ptr
          integer(ikind), dimension(2)               :: gen_coords
          integer(ikind), dimension(2)               :: match_table
          integer       , dimension(5,ny)            :: test_grdpts_id

          logical :: test_loc
          integer :: k


          test_validated = .true.


          call ini_mainlayer_interface_for_bc_pt_crenel(
     $         mainlayer_interface_used,
     $         N_sublayer_ptr,
     $         S_sublayer_ptr,
     $         E_sublayer_ptr)


          do k=1,3

             !input
             call get_test_param_bc_pt_crenel(
     $            k,
     $            N_sublayer_ptr,
     $            S_sublayer_ptr,
     $            E_sublayer_ptr,
     $            gen_coords,
     $            match_table,
     $            test_grdpts_id)

             !output
             call mainlayer_interface_used%detect_and_curb_bc_pt_crenel(
     $            E_sublayer_ptr,
     $            gen_coords,
     $            match_table)

             !validation
             test_loc = is_int_matrix_validated(
     $            E_sublayer_ptr%grdpts_id,
     $            test_grdpts_id,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')'
             end if

          end do

        end function test_detect_and_curb_bc_pt_crenel


        subroutine ini_mainlayer_interface_for_bc_pt_crenel(
     $     mainlayer_interface_used,
     $     N_sublayer_ptr,
     $     S_sublayer_ptr,
     $     E_sublayer_ptr)

          implicit none

          type(mainlayer_interface_grdpts_id_update), intent(out) :: mainlayer_interface_used
          type(bf_sublayer), pointer                , intent(out) :: N_sublayer_ptr
          type(bf_sublayer), pointer                , intent(out) :: S_sublayer_ptr
          type(bf_sublayer), pointer                , intent(out) :: E_sublayer_ptr


          integer :: j


          !North buffer layer
          !------------------------------------------------------------
          allocate(N_sublayer_ptr)
          N_sublayer_ptr%localization = N
          N_sublayer_ptr%alignment = reshape((/
     $         align_E-5, align_N, align_E, align_N+3/),
     $         (/2,2/))
          allocate(N_sublayer_ptr%grdpts_id(10,8))
          N_sublayer_ptr%grdpts_id = reshape((/
     $         1,1,1,1,1,1,1,1,2,3,
     $         1,1,1,1,1,1,1,1,2,3,
     $         2,2,1,1,1,1,1,1,2,3,          
     $         3,2,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,2,3,
     $         3,2,2,2,2,2,2,2,2,3,
     $         3,3,3,3,3,3,3,3,3,3/),
     $         (/10,8/))

          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         NE_interface_type,
     $         N_sublayer_ptr)


          !South buffer layer
          !------------------------------------------------------------
          allocate(S_sublayer_ptr)
          S_sublayer_ptr%localization = S
          S_sublayer_ptr%alignment = reshape((/
     $         align_E-5, align_S-3, align_E, align_S/),
     $         (/2,2/))
          allocate(S_sublayer_ptr%grdpts_id(10,8))
          S_sublayer_ptr%grdpts_id = reshape((/
     $         3,3,3,3,3,3,3,3,3,3,
     $         3,2,2,2,2,2,2,2,2,3,
     $         3,2,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,2,3,
     $         2,2,1,1,1,1,1,1,2,3,          
     $         1,1,1,1,1,1,1,1,2,3,
     $         1,1,1,1,1,1,1,1,2,3/),
     $         (/10,8/))
          
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         S_sublayer_ptr)

          
          !East buffer layer
          !------------------------------------------------------------
          allocate(E_sublayer_ptr)
          E_sublayer_ptr%localization = E
          E_sublayer_ptr%alignment = reshape((/
     $         align_E,align_S+1,align_E,align_N-1/),
     $         (/2,2/))
          allocate(E_sublayer_ptr%grdpts_id(5,ny))
          
          do j=1,ny
             E_sublayer_ptr%grdpts_id(:,j) = [1,1,1,2,3]
          end do

          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         E_sublayer_ptr)

          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         NE_interface_type,
     $         E_sublayer_ptr)


        end subroutine ini_mainlayer_interface_for_bc_pt_crenel


        subroutine get_test_param_bc_pt_crenel(
     $     test_id,
     $     N_sublayer_ptr,
     $     S_sublayer_ptr,
     $     E_sublayer_ptr,
     $     gen_coords,
     $     match_table,
     $     test_grdpts_id)

          implicit none

          integer                        , intent(in)    :: test_id
          type(bf_sublayer)              , intent(inout) :: N_sublayer_ptr
          type(bf_sublayer)              , intent(inout) :: S_sublayer_ptr
          type(bf_sublayer)              , intent(inout) :: E_sublayer_ptr
          integer(ikind), dimension(2)   , intent(out)   :: gen_coords
          integer(ikind), dimension(2)   , intent(out)   :: match_table
          integer       , dimension(5,ny), intent(out)   :: test_grdpts_id

          integer :: j

          match_table = [align_E-3,align_S+1-3]
          
          do j=1,ny
             test_grdpts_id(:,j) = [1,1,1,2,3]
          end do

          select case(test_id)

            case(1)
               E_sublayer_ptr%grdpts_id(1:5,align_S+4:align_S+8) =
     $              reshape((/
     $              1,1,1,2,3,
     $              1,1,2,2,3,
     $              1,1,2,3,3,
     $              1,1,2,2,3,
     $              1,1,1,2,3/),
     $              (/5,5/))

               gen_coords = [align_E+1,align_S+6]

            case(2)
               E_sublayer_ptr%grdpts_id(1:5,1:3) =
     $              reshape((/
     $              1,1,2,3,3,
     $              1,1,2,2,3,
     $              1,1,1,2,3/),
     $              (/5,3/))

               gen_coords = [align_E+1,align_S-1]

               S_sublayer_ptr%grdpts_id = reshape((/
     $              3,3,3,3,3,3,3,3,3,3,
     $              3,2,2,2,2,2,2,2,2,3,
     $              3,2,1,1,1,1,1,1,2,3,
     $              3,2,1,1,1,1,1,2,2,3,
     $              3,2,1,1,1,1,1,2,3,3,
     $              2,2,1,1,1,1,1,2,2,3,          
     $              1,1,1,1,1,1,1,1,2,3,
     $              1,1,1,1,1,1,1,1,2,3/),
     $              (/10,8/))

            case(3)
               E_sublayer_ptr%grdpts_id(1:5,ny-2:ny) =
     $              reshape((/
     $              1,1,1,2,3,
     $              1,1,2,2,3,
     $              1,1,2,3,3/),
     $              (/5,3/))

               gen_coords = [align_E+1,align_N+1]

               N_sublayer_ptr%grdpts_id = reshape((/
     $              1,1,1,1,1,1,1,1,2,3,
     $              1,1,1,1,1,1,1,1,2,3,
     $              2,2,1,1,1,1,1,2,2,3,          
     $              3,2,1,1,1,1,1,2,3,3,
     $              3,2,1,1,1,1,1,2,2,3,
     $              3,2,1,1,1,1,1,1,2,3,
     $              3,2,2,2,2,2,2,2,2,3,
     $              3,3,3,3,3,3,3,3,3,3/),
     $              (/10,8/))

            case(4)
               gen_coords = [align_E+2,align_N-1]

            end select

        end subroutine get_test_param_bc_pt_crenel


        subroutine get_test_param_bc_pt_crenel_nl(
     $     test_id,
     $     N_sublayer_ptr,
     $     S_sublayer_ptr,
     $     E_sublayer_ptr,
     $     gen_coords,
     $     match_table,
     $     test_grdpts_id)

          implicit none

          integer                        , intent(in)    :: test_id
          type(bf_sublayer)              , intent(inout) :: N_sublayer_ptr
          type(bf_sublayer)              , intent(inout) :: S_sublayer_ptr
          type(bf_sublayer)              , intent(inout) :: E_sublayer_ptr
          integer(ikind), dimension(2)   , intent(out)   :: gen_coords
          integer(ikind), dimension(2)   , intent(out)   :: match_table
          integer       , dimension(5,ny), intent(out)   :: test_grdpts_id

          integer :: j

          match_table = [align_E-3,align_S+1-3]
          
          do j=1,ny
             test_grdpts_id(:,j) = [1,1,1,2,3]
          end do

          select case(test_id)

            case(1)
               E_sublayer_ptr%grdpts_id(1:5,align_S+4:align_S+8) =
     $              reshape((/
     $              1,1,1,2,3,
     $              1,1,2,2,3,
     $              1,1,2,3,3,
     $              1,1,2,2,3,
     $              1,1,1,2,3/),
     $              (/5,5/))

               gen_coords = [align_E,align_S+4]

            case(2)
               E_sublayer_ptr%grdpts_id(1:5,1:3) =
     $              reshape((/
     $              1,1,2,3,3,
     $              1,1,2,2,3,
     $              1,1,1,2,3/),
     $              (/5,3/))

               gen_coords = [align_E,align_S+1]

               S_sublayer_ptr%grdpts_id = reshape((/
     $              3,3,3,3,3,3,3,3,3,3,
     $              3,2,2,2,2,2,2,2,2,3,
     $              3,2,1,1,1,1,1,1,2,3,
     $              3,2,1,1,1,1,1,2,2,3,
     $              3,2,1,1,1,1,1,2,3,3,
     $              2,2,1,1,1,1,1,2,2,3,          
     $              1,1,1,1,1,1,1,1,2,3,
     $              1,1,1,1,1,1,1,1,2,3/),
     $              (/10,8/))

            case(3)
               E_sublayer_ptr%grdpts_id(1:5,ny-2:ny) =
     $              reshape((/
     $              1,1,1,2,3,
     $              1,1,2,2,3,
     $              1,1,2,3,3/),
     $              (/5,3/))

               gen_coords = [align_E,align_N-1]

               N_sublayer_ptr%grdpts_id = reshape((/
     $              1,1,1,1,1,1,1,1,2,3,
     $              1,1,1,1,1,1,1,1,2,3,
     $              2,2,1,1,1,1,1,2,2,3,          
     $              3,2,1,1,1,1,1,2,3,3,
     $              3,2,1,1,1,1,1,2,2,3,
     $              3,2,1,1,1,1,1,1,2,3,
     $              3,2,2,2,2,2,2,2,2,3,
     $              3,3,3,3,3,3,3,3,3,3/),
     $              (/10,8/))

            case(4)
               gen_coords = [align_E,align_N-1]

            end select

        end subroutine get_test_param_bc_pt_crenel_nl


        function test_finalize_for_bc_interior_pt_crenel(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mainlayer_interface_grdpts_id_update) :: mainlayer_interface_used
          type(bf_sublayer), pointer                 :: bf_sublayer_ptr
          integer(ikind)                             :: gen_i
          integer(ikind)                             :: gen_j
          integer(ikind), dimension(2)               :: match_table
          integer       , dimension(5,6)             :: test_grdpts_id


          !input
          call get_param_finalize_for_interior_crenel(
     $         mainlayer_interface_used,
     $         bf_sublayer_ptr,
     $         gen_i,
     $         gen_j,
     $         match_table,
     $         test_grdpts_id)


          !output
          call mainlayer_interface_used%finalize_for_bc_interior_pt_crenel(
     $         bf_sublayer_ptr,
     $         [gen_i,gen_j],
     $         match_table)


          !validation
          test_validated = is_int_matrix_validated(
     $         bf_sublayer_ptr%grdpts_id,
     $         test_grdpts_id,
     $         detailled)


        end function test_finalize_for_bc_interior_pt_crenel

      
        subroutine get_param_finalize_for_interior_crenel(
     $     mainlayer_interface_used,
     $     bf_sublayer_ptr,
     $     gen_i,
     $     gen_j,
     $     match_table,
     $     test_grdpts_id)

          implicit none

          type(mainlayer_interface_grdpts_id_update), intent(inout) :: mainlayer_interface_used
          type(bf_sublayer), pointer                , intent(inout) :: bf_sublayer_ptr
          integer(ikind)                            , intent(out)   :: gen_i
          integer(ikind)                            , intent(out)   :: gen_j
          integer(ikind), dimension(2)              , intent(out)   :: match_table
          integer       , dimension(5,6)            , intent(out)   :: test_grdpts_id
        
          type(bf_sublayer), pointer :: nbf_sublayer_ptr


          gen_i = align_E
          gen_j = align_S+1
          match_table = [align_E-3,align_S+1-3]          

          test_grdpts_id = reshape((/
     $         1,1,1,2,3,
     $         1,1,1,2,3,
     $         1,1,1,2,3,
     $         1,1,1,2,3,
     $         1,1,1,2,3,
     $         1,1,1,2,3/),
     $         (/5,6/))


          !East buffer layer
          allocate(bf_sublayer_ptr)

          bf_sublayer_ptr%localization = E

          bf_sublayer_ptr%alignment = reshape((/
     $         align_E,align_S+1,align_E,align_S+2/),
     $         (/2,2/))

          allocate(bf_sublayer_ptr%grdpts_id(5,6))
          bf_sublayer_ptr%grdpts_id = reshape((/
     $         1,1,1,2,3,
     $         1,1,2,2,3,
     $         1,1,1,2,3,
     $         1,1,2,2,3,
     $         1,1,1,2,3,
     $         1,1,1,2,3/),
     $         (/5,6/))

          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         bf_sublayer_ptr)


          !South buffer layer
          allocate(nbf_sublayer_ptr)

          nbf_sublayer_ptr%localization = S

          nbf_sublayer_ptr%alignment = reshape((/
     $         align_E-4,align_S-2,align_E,align_S/),
     $         (/2,2/))

          allocate(nbf_sublayer_ptr%grdpts_id(9,7))
          nbf_sublayer_ptr%grdpts_id = reshape((/
     $         3,3,3,3,3,3,3,3,3,
     $         3,2,2,2,2,2,2,2,3,
     $         3,2,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,2,3,
     $         2,2,1,1,1,1,1,2,3/),
     $         (/9,7/))

          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         nbf_sublayer_ptr)

        end subroutine get_param_finalize_for_interior_crenel


        function test_curb_bc_interior_pt_crenel(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mainlayer_interface_grdpts_id_update) :: mainlayer_interface_used
          type(bf_sublayer), pointer                 :: bf_sublayer_ptr
          integer(ikind)                             :: gen_i
          integer(ikind)                             :: gen_j
          integer(ikind), dimension(2)               :: match_table
          logical                                    :: test_can_be_curbed
          logical                                    :: can_be_curbed

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,2

             !input
             call get_test_param_curb_interior_crenel(
     $            k,
     $            mainlayer_interface_used,
     $            bf_sublayer_ptr,
     $            gen_i,
     $            gen_j,
     $            match_table,
     $            test_can_be_curbed)

             !output
             can_be_curbed = mainlayer_interface_used%can_bc_interior_pt_crenel_be_curbed(
     $            bf_sublayer_ptr,
     $            gen_i,
     $            gen_j,
     $            match_table)

             !validation
             test_loc = can_be_curbed.eqv.test_can_be_curbed
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_curb_bc_interior_pt_crenel


        subroutine get_test_param_curb_interior_crenel(
     $     test_id,
     $     mainlayer_interface_used,
     $     bf_sublayer_ptr,
     $     gen_i,
     $     gen_j,
     $     match_table,
     $     test_can_be_curbed)

          implicit none

          integer                                   , intent(in)    :: test_id
          type(mainlayer_interface_grdpts_id_update), intent(inout) :: mainlayer_interface_used
          type(bf_sublayer), pointer                , intent(inout) :: bf_sublayer_ptr
          integer(ikind)                            , intent(out)   :: gen_i
          integer(ikind)                            , intent(out)   :: gen_j
          integer(ikind), dimension(2)              , intent(out)   :: match_table
          logical                                   , intent(out)   :: test_can_be_curbed
        
          type(bf_sublayer), pointer :: nbf_sublayer_ptr


          !East buffer layer
          allocate(bf_sublayer_ptr)

          bf_sublayer_ptr%localization = E

          bf_sublayer_ptr%alignment = reshape((/
     $         align_E,align_S+1,align_E,align_S+1/),
     $         (/2,2/))

          allocate(bf_sublayer_ptr%grdpts_id(5,5))
          bf_sublayer_ptr%grdpts_id = reshape((/
     $         1,1,1,2,3,
     $         1,1,2,2,3,
     $         1,1,2,3,3,
     $         1,1,2,3,0,
     $         1,1,2,3,0/),
     $         (/5,5/))

          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         bf_sublayer_ptr)

          match_table(1) = align_E-3
          match_table(2) = align_S+1-3


          !South buffer layer
          allocate(nbf_sublayer_ptr)

          nbf_sublayer_ptr%localization = S

          nbf_sublayer_ptr%alignment = reshape((/
     $         align_E-4,align_S-2,align_E,align_S/),
     $         (/2,2/))

          allocate(nbf_sublayer_ptr%grdpts_id(9,7))
          nbf_sublayer_ptr%grdpts_id = reshape((/
     $         3,3,3,3,3,3,3,3,3,
     $         3,2,2,2,2,2,2,2,3,
     $         3,2,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,2,2,3,
     $         3,2,1,1,1,1,2,3,3,
     $         2,2,1,1,1,1,2,3,0/),
     $         (/9,7/))

          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         nbf_sublayer_ptr)


          select case(test_id)
            case(1)
               gen_i = align_E
               gen_j = align_S+1
               test_can_be_curbed = .false.
            case(2)
               gen_i = align_E
               gen_j = align_S-1
               test_can_be_curbed = .true.
          end select

        end subroutine get_test_param_curb_interior_crenel


        function test_detect_bc_interior_pt_crenel(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mainlayer_interface_grdpts_id_update) :: mainlayer_interface_used
          type(bf_sublayer), pointer                 :: bf_sublayer_ptr
          integer(ikind)                             :: gen_i
          integer(ikind)                             :: gen_j
          integer(ikind), dimension(2)               :: match_table
          logical                                    :: test_crenel
          logical                                    :: crenel          

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,2

             !input
             call get_test_param_detect_interior_crenel(
     $            k,
     $            mainlayer_interface_used,
     $            bf_sublayer_ptr,
     $            gen_i,
     $            gen_j,
     $            match_table,
     $            test_crenel)

             !output
             crenel = mainlayer_interface_used%detect_bc_interior_pt_crenel(
     $            bf_sublayer_ptr,
     $            gen_i,
     $            gen_j,
     $            match_table)

             !validation
             test_loc = crenel.eqv.test_crenel
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_detect_bc_interior_pt_crenel


        subroutine get_test_param_detect_interior_crenel(
     $     test_id,
     $     mainlayer_interface_used,
     $     bf_sublayer_ptr,
     $     gen_i,
     $     gen_j,
     $     match_table,
     $     test_crenel)

          implicit none

          integer                                   , intent(in)    :: test_id
          type(mainlayer_interface_grdpts_id_update), intent(inout) :: mainlayer_interface_used
          type(bf_sublayer), pointer                , intent(inout) :: bf_sublayer_ptr
          integer(ikind)                            , intent(out)   :: gen_i
          integer(ikind)                            , intent(out)   :: gen_j
          integer(ikind), dimension(2)              , intent(out)   :: match_table
          logical                                   , intent(out)   :: test_crenel
        
          type(bf_sublayer), pointer :: nbf_sublayer_ptr


          !East buffer layer
          allocate(bf_sublayer_ptr)

          bf_sublayer_ptr%localization = E

          bf_sublayer_ptr%alignment = reshape((/
     $         align_E,align_S+1,align_E,align_S+1/),
     $         (/2,2/))

          allocate(bf_sublayer_ptr%grdpts_id(5,5))
          bf_sublayer_ptr%grdpts_id = reshape((/
     $         1,1,1,2,3,
     $         1,1,2,2,3,
     $         1,1,2,3,3,
     $         1,1,2,3,0,
     $         1,1,2,3,0/),
     $         (/5,5/))

          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         bf_sublayer_ptr)

          match_table(1) = align_E-3
          match_table(2) = align_S+1-3


          !South buffer layer
          allocate(nbf_sublayer_ptr)

          nbf_sublayer_ptr%localization = S

          nbf_sublayer_ptr%alignment = reshape((/
     $         align_E-4,align_S-2,align_E,align_S/),
     $         (/2,2/))

          allocate(nbf_sublayer_ptr%grdpts_id(9,7))
          nbf_sublayer_ptr%grdpts_id = reshape((/
     $         3,3,3,3,3,3,3,3,3,
     $         3,2,2,2,2,2,2,2,3,
     $         3,2,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,2,2,3,
     $         3,2,1,1,1,1,2,3,3,
     $         2,2,1,1,1,1,2,3,0/),
     $         (/9,7/))

          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         nbf_sublayer_ptr)


          select case(test_id)
            case(1)
               gen_i = align_E
               gen_j = align_S+1
               test_crenel = .false.
            case(2)
               gen_i = align_E
               gen_j = align_S-1
               test_crenel = .true.
          end select

        end subroutine get_test_param_detect_interior_crenel


        function test_update_grdpts_id_around(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          if(debug_geometry_update) then
             test_validated = test_update_grdpts_id_around_geometry(detailled)
          else
             test_validated = test_update_grdpts_id_around_edge(detailled)
          end if

        end function test_update_grdpts_id_around


        function test_update_grdpts_id_around_geometry(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(mainlayer_interface_grdpts_id_update) :: mainlayer_interface_used
          type(pmodel_eq)                            :: p_model
          real(rkind)                                :: t
          real(rkind)                                :: dt
          real(rkind), dimension(nx)                 :: interior_x_map
          real(rkind), dimension(ny)                 :: interior_y_map
          real(rkind), dimension(nx,ny,ne)           :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)           :: interior_nodes1
          type(bf_sublayer)                          :: bf_sublayer_ptr
          integer(ikind)                             :: i_prev
          integer(ikind)                             :: j_prev
          integer(ikind)                             :: i
          integer(ikind)                             :: j
          integer(ikind), dimension(2)               :: match_table
          integer       , dimension(17,17)           :: test_grdpts_id
          
          
          test_validated = .true.


          allocate(bf_sublayer_ptr%grdpts_id(17,17))
          i_prev = 9
          j_prev = 9
          match_table = [0,0]
          

          do j=j_prev-6,j_prev+6
             do i=i_prev-6, i_prev+6
                
                !input
                call get_param_test_update_grdpts_id(
     $               i_prev,j_prev,
     $               i,j,
     $               bf_sublayer_ptr%grdpts_id,
     $               test_grdpts_id)

                !output
                call mainlayer_interface_used%update_grdpts_id_around_new_interior_pt(
     $               p_model,
     $               t,dt,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes0,
     $               interior_nodes1,
     $               bf_sublayer_ptr,
     $               i_prev,
     $               j_prev,
     $               i,j,
     $               match_table)

                !validation
                test_loc = is_int_matrix_validated(
     $               bf_sublayer_ptr%grdpts_id,
     $               test_grdpts_id,
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test ('',3I3,'') failed'')', i,j
                end if

             end do
          end do

        end function test_update_grdpts_id_around_geometry


        function test_update_grdpts_id_around_edge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(mainlayer_interface_grdpts_id_update) :: mainlayer_interface_used
          type(pmodel_eq)                            :: p_model
          real(rkind)                                :: t
          real(rkind)                                :: dt
          real(rkind), dimension(nx)                 :: interior_x_map
          real(rkind), dimension(ny)                 :: interior_y_map
          real(rkind), dimension(nx,ny,ne)           :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)           :: interior_nodes1
          type(bf_sublayer)                          :: bf_sublayer_ptr
          integer(ikind), dimension(2)               :: match_table
          integer       , dimension(6,5)             :: test_grdpts_id
          real(rkind)   , dimension(ne)              :: test_newgrdpt

          logical :: test_loc
          integer :: k
          

          test_validated = .true.

          
          !input
          !------------------------------------------------------------
          call ini_for_test_update_grdpts_id_around_edge(
     $         t,dt,
     $         p_model,
     $         bf_sublayer_ptr,
     $         test_newgrdpt,
     $         test_grdpts_id)

          match_table(1) = align_E+3-6
          match_table(2) = align_S+10-3


          !output
          !------------------------------------------------------------
          call mainlayer_interface_used%update_grdpts_id_around_new_interior_pt(
     $         p_model,
     $         t,dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1,
     $         bf_sublayer_ptr,
     $         nx/2,
     $         ny/2,
     $         align_E+1,align_S+10,
     $         match_table)


          !validation
          !------------------------------------------------------------
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%grdpts_id,
     $         test_grdpts_id,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test grdpts_id failed'')'
             print '(''grdpts_id:'')'
             do k=1,5
                print '(6I2)', bf_sublayer_ptr%grdpts_id(:,5-(k-1))
             end do
             print '()'

             print '(''test_grdpts_id:'')'
             do k=1,5
                print '(6I2)', test_grdpts_id(:,5-(k-1))
             end do
             print '()'
          end if


          test_loc = is_real_vector_validated(
     $         bf_sublayer_ptr%nodes(6,3,:),
     $         test_newgrdpt,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test newgrdpt failed'')'
          end if

        end function test_update_grdpts_id_around_edge


        subroutine ini_for_test_update_grdpts_id_around_edge(
     $     t,dt,
     $     p_model,
     $     bf_sublayer_ptr,
     $     test_newgrdpt,
     $     test_grdpts_id)

          implicit none

          real(rkind)                , intent(out) :: t
          real(rkind)                , intent(out) :: dt
          type(pmodel_eq)            , intent(out) :: p_model
          type(bf_sublayer)          , intent(out) :: bf_sublayer_ptr
          real(rkind), dimension(ne) , intent(out) :: test_newgrdpt
          integer    , dimension(6,5), intent(out) :: test_grdpts_id

          
          !initialization
          !------------------------------------------------------------
          t  = 0.0d0
          dt = 0.25d0
          call p_model%initial_conditions%ini_far_field()


          !initialization of the data at t
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%grdpts_id(6,5))
          allocate(bf_sublayer_ptr%x_map(6))
          allocate(bf_sublayer_ptr%y_map(5))
          allocate(bf_sublayer_ptr%nodes(6,5,ne))

          bf_sublayer_ptr%localization = E

          bf_sublayer_ptr%alignment = reshape((/
     $         align_E, align_S+10, align_E+1, align_S+10/),(/2,2/))

          bf_sublayer_ptr%grdpts_id = reshape((/
     $         1,1,1,2,3,0,
     $         1,1,1,2,3,0,
     $         1,1,1,2,3,0,
     $         1,1,1,2,3,0,
     $         1,1,1,2,3,0/),
     $         (/6,5/))

          bf_sublayer_ptr%x_map(1:6) = [-1.5d0, -0.50d0, 0.5d0, 1.50d0, 2.5d0, 3.50d0]
          bf_sublayer_ptr%y_map(1:5) = [-0.25d0, 0.0d0, 0.25d0, 0.5d0, 0.75d0]

          bf_sublayer_ptr%nodes(3:6,2:5,:) = reshape((/
     $           1.50d0, 1.455d0, 1.48d0, 0.0d0,
     $           1.20d0, 1.350d0, 1.25d0, 0.0d0,
     $           1.49d0, 1.250d0, 1.40d0, 0.0d0,
     $           0.00d0, 0.000d0, 0.00d0, 0.0d0,
     $           
     $           0.128d0, 0.450d0, 0.135d0, 0.0d0,
     $           0.148d0, 0.150d0, 0.122d0, 0.0d0,
     $           0.142d0, 1.152d0, 0.236d0, 0.0d0,
     $           0.000d0, 0.000d0, 0.000d0, 0.0d0,
     $           
     $           0.006d0, 0.0600d0, 0.020d0, 0.0d0,
     $           0.000d0, 0.0028d0, 0.035d0, 0.0d0,
     $           0.020d0, 0.0030d0, 0.040d0, 0.0d0,
     $           0.000d0, 0.0000d0, 0.000d0, 0.0d0,
     $           
     $           4.876d0, 4.825d0, 4.862d0, 0.0d0,
     $           4.890d0, 4.871d0, 4.892d0, 0.0d0,
     $           4.865d0, 4.757d0, 4.895d0, 0.0d0,
     $           0.000d0, 0.000d0, 0.000d0, 0.0d0
     $           /),
     $           (/4,4,ne/))


          !initialization of the data at t-dt
          !------------------------------------------------------------
          allocate(bf_sublayer_ptr%bf_compute_used%alignment_tmp(2,2))
          allocate(bf_sublayer_ptr%bf_compute_used%grdpts_id_tmp(5,5))
          allocate(bf_sublayer_ptr%bf_compute_used%x_map_tmp(5))
          allocate(bf_sublayer_ptr%bf_compute_used%y_map_tmp(5))
          allocate(bf_sublayer_ptr%bf_compute_used%nodes_tmp(5,5,ne))

          bf_sublayer_ptr%bf_compute_used%alignment_tmp = reshape((/
     $         align_E, align_S+10, align_E, align_S+10/),(/2,2/))
          
          bf_sublayer_ptr%bf_compute_used%grdpts_id_tmp = reshape((/
     $         1,1,1,2,3,
     $         1,1,1,2,3,
     $         1,1,1,2,3,
     $         1,1,1,2,3,
     $         1,1,1,2,3/),
     $         (/5,5/))

          bf_sublayer_ptr%bf_compute_used%x_map_tmp(1:5) = [-1.50d0, -0.50d0, 0.50d0, 1.50d0, 2.50d0]
          bf_sublayer_ptr%bf_compute_used%y_map_tmp(1:5) = [-0.25d0,  0.0d0 , 0.25d0, 0.50d0, 0.75d0]

          bf_sublayer_ptr%bf_compute_used%nodes_tmp(3:5,2:4,:) = reshape((/
     $           1.48d0, 1.30d0, 1.35d0,
     $           1.26d0, 1.45d0, 1.4d0,
     $           1.46d0, 1.27d0, 1.47d0,
     $           
     $           0.128d0, 0.127d0, 0.142d0,
     $           1.138d0, 0.148d0, 0.132d0,
     $           0.146d0, 0.143d0, 0.145d0,
     $           
     $           0.0050d0, 0.020d0, 0.060d0,
     $           0.0025d0, 0.001d0, 0.015d0,
     $           0.0100d0, 0.002d0, 0.050d0,
     $           
     $           4.88d0, 4.870d0, 4.855d0,
     $           4.85d0, 4.865d0, 4.845d0,
     $           4.89d0, 4.870d0, 4.860d0/),
     $           (/3,3,ne/))

          !initialization of the new grid-point
          !------------------------------------------------------------
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
               print '(''test_mainlayer_interface_grdpts_id_update'')'
               print '(''test_update_grdpts_id_around'')'
               print '(''obc_eigenqties_strategy not recognized'')'
               stop ''
               
          end select
          

          test_grdpts_id = reshape((/
     $         1,1,1,2,3,3,
     $         1,1,1,2,2,3,
     $         1,1,1,1,2,3,
     $         1,1,1,2,2,3,
     $         1,1,1,2,3,3/),
     $         (/6,5/))

        end subroutine ini_for_test_update_grdpts_id_around_edge


        subroutine get_param_test_update_grdpts_id(
     $     i_prev,j_prev,
     $     i_center,j_center,
     $     grdpts_id,
     $     test_grdpts_id)

          implicit none

          integer(ikind)                  , intent(in)  :: i_prev
          integer(ikind)                  , intent(in)  :: j_prev
          integer(ikind)                  , intent(in)  :: i_center
          integer(ikind)                  , intent(in)  :: j_center
          integer       , dimension(17,17), intent(out) :: grdpts_id
          integer       , dimension(17,17), intent(out) :: test_grdpts_id

          integer(ikind) :: i,j


          grdpts_id = reshape((/
     $         ((no_pt, i=1,17),j=1,17)/),
     $         (/17,17/))
          
          do j=j_prev-2,j_prev+2
             do i=i_prev-2,i_prev+2
                grdpts_id(i,j) = bc_pt
             end do
          end do

          do j=j_prev-1,j_prev+1
             do i=i_prev-1,i_prev+1
                grdpts_id(i,j) = bc_interior_pt
             end do
          end do

          grdpts_id(i_prev,j_prev) = interior_pt


          test_grdpts_id = reshape((/
     $         ((no_pt, i=1,17),j=1,17)/),
     $         (/17,17/))

          do j=j_prev-2,j_prev+2
             do i=i_prev-2,i_prev+2
                test_grdpts_id(i,j) = bc_pt
             end do
          end do

          do j=j_center-2,j_center+2
             do i=i_center-2,i_center+2
                test_grdpts_id(i,j) = bc_pt
             end do
          end do

          do j=j_prev-1,j_prev+1
             do i=i_prev-1,i_prev+1
                test_grdpts_id(i,j) = bc_interior_pt
             end do
          end do

          do j=j_center-1,j_center+1
             do i=i_center-1,i_center+1
                test_grdpts_id(i,j) = bc_interior_pt
             end do
          end do

          test_grdpts_id(i_prev,j_prev)     = interior_pt
          test_grdpts_id(i_center,j_center) = interior_pt

        end subroutine get_param_test_update_grdpts_id


        subroutine check_inputs()

          implicit none

          real(rkind), dimension(4) :: far_field
          type(pmodel_eq)           :: p_model

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

      end program test_mainlayer_interface_grdpts_id_update
