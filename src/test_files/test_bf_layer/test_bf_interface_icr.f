      program test_bf_interface_icr

        use bf_interface_icr_class, only :
     $       bf_interface_icr

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_int_matrix_validated,
     $       is_int_matrix3D_validated

        use icr_interface_class, only :
     $       icr_interface

        use icr_path_chain_class, only :
     $       icr_path_chain

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W,
     $       
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $       
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       
     $       dct_icr_distance,
     $       no_overlap,
     $       N_overlap,
     $       S_overlap,
     $       E_overlap,
     $       W_overlap,
     $       NS_overlap,
     $       EW_overlap,
     $       
     $       cpt2not_and_cpt3normal,
     $       cpt2not_and_cpt3not

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne

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


        test_loc = test_is_interior_node_activated(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_is_interior_node_activated: '',L1)', test_loc
        print ''


        test_loc = test_is_bf_layer_node_activated(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_is_bf_layer_node_activated: '',L1)', test_loc
        print ''


        test_loc = test_is_node_activated(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_is_node_activated: '',L1)', test_loc
        print ''


        test_loc = test_stage(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_stage: '',L1)', test_loc
        print ''


        test_loc = test_analyze_bc_section_edge_y(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_bc_section_edge_y: '',L1)', test_loc
        print ''


        test_loc = test_analyze_bc_section_edge_x(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_bc_section_edge_x: '',L1)', test_loc
        print ''


        test_loc = test_analyze_bc_section_square_bounds(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_bc_section_square_bounds: '',L1)', test_loc
        print ''


        test_loc = test_analyze_bc_section_square_xy(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_bc_section_square_xy: '',L1)', test_loc
        print '()'


        test_loc = test_analyze_bc_section(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_bc_section: '',L1)', test_loc
        print '()'


        test_loc = test_get_interior_bcs_mainlayer_id(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_interior_bcs_mainlayer_id: '',L1)', test_loc
        print '()'


        test_loc = test_analyze_interior_bc_sections(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_interior_bc_section: '',L1)', test_loc
        print '()'


        test_loc = test_analyze_bf_layer_bc_sections(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_bf_layer_bc_section: '',L1)', test_loc
        print '()'


        test_loc = test_analyze_and_update_boundary(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_analyze_and_update_boundary: '',L1)', test_loc
        print '()'


        test_loc = test_adapt_domain_extension(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_adapt_domain_extension: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_adapt_domain_extension(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_icr)                      :: bf_interface_used
          type(icr_interface)                         :: icr_interface_used
          real(rkind), dimension(nx)                  :: interior_x_map
          real(rkind), dimension(ny)                  :: interior_y_map
          real(rkind), dimension(nx,ny,ne)            :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)            :: interior_nodes1
          integer                                     :: mainlayer_id
          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections
          real(rkind)                                 :: t
          real(rkind)                                 :: dt
          type(pmodel_eq)                             :: p_model

          logical                        :: test_loc
          integer(ikind), dimension(2,2) :: test_alignment
          type(bf_sublayer), pointer     :: bf_sublayer_ptr


          test_validated = .true.


          !input
          call get_test_param_analyze_and_update_boundary(
     $         bf_interface_used,
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes1,
     $         mainlayer_id,
     $         interior_bc_sections,
     $         t,
     $         dt,
     $         test_alignment)

          !output
          call bf_interface_used%adapt_domain_extension(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1,
     $         p_model,
     $         t,
     $         dt,
     $         interior_bc_sections)

          !validation
          test_loc = bf_interface_used%mainlayer_pointers(N)%get_nb_sublayers().eq.1
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_sublayers(N) failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(S)%get_nb_sublayers().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_sublayers(S) failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(E)%get_nb_sublayers().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_sublayers(E) failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(W)%get_nb_sublayers().eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_sublayers(W) failed'')'
          end if
          
          if(test_loc) then
             bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%get_head_sublayer()
             test_loc = is_int_matrix_validated(
     $            bf_sublayer_ptr%get_alignment_tab(),
     $            test_alignment,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''alignment failed'')'
             end if
          end if

        end function test_adapt_domain_extension


        function test_analyze_and_update_boundary(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_icr)                      :: bf_interface_used
          type(icr_interface)                         :: icr_interface_used
          real(rkind), dimension(nx)                  :: interior_x_map
          real(rkind), dimension(ny)                  :: interior_y_map
          real(rkind), dimension(nx,ny,ne)            :: interior_nodes0
          real(rkind), dimension(nx,ny,ne)            :: interior_nodes1
          integer                                     :: mainlayer_id
          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections
          real(rkind)                                 :: t
          real(rkind)                                 :: dt
          type(pmodel_eq)                             :: p_model

          logical                        :: test_loc
          integer(ikind), dimension(2,2) :: test_alignment
          type(bf_sublayer), pointer     :: bf_sublayer_ptr


          test_validated = .true.


          !input
          call get_test_param_analyze_and_update_boundary(
     $         bf_interface_used,
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes1,
     $         mainlayer_id,
     $         interior_bc_sections,
     $         t,
     $         dt,
     $         test_alignment)

          !output
          call bf_interface_used%analyze_and_update_boundary(
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1,
     $         p_model,
     $         t,
     $         dt,
     $         interior_bc_sections,
     $         mainlayer_id)

          !validation
          test_loc = bf_interface_used%mainlayer_pointers(N)%get_nb_sublayers().eq.1
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_sublayers failed'')'
          end if
          
          if(test_loc) then
             bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%get_head_sublayer()
             test_loc = is_int_matrix_validated(
     $            bf_sublayer_ptr%get_alignment_tab(),
     $            test_alignment,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''alignment failed'')'
             end if
          end if

        end function test_analyze_and_update_boundary


        subroutine get_test_param_analyze_and_update_boundary(
     $     bf_interface_used,
     $     icr_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     mainlayer_id,
     $     interior_bc_sections,
     $     t,
     $     dt,
     $     test_alignment)

          implicit none


          type(bf_interface_icr)                     , intent(inout) :: bf_interface_used
          type(icr_interface)                        , intent(inout) :: icr_interface_used
          real(rkind), dimension(nx)                 , intent(inout) :: interior_x_map
          real(rkind), dimension(ny)                 , intent(inout) :: interior_y_map
          real(rkind), dimension(nx,ny,ne)           , intent(inout) :: interior_nodes
          integer                                    , intent(inout) :: mainlayer_id
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: interior_bc_sections
          real(rkind)                                , intent(out)   :: t
          real(rkind)                                , intent(out)   :: dt
          integer(ikind), dimension(2,2)             , intent(out)   :: test_alignment
 

          integer(ikind), dimension(2,2) :: bf_alignment_tmp
          type(bf_sublayer), pointer     :: added_sublayer
          real(rkind), parameter :: A = -1.0d0
          real(rkind), parameter :: D =  1.0d0
          

          t=0.0d0
          dt=0.1d0

          mainlayer_id = N


          call bf_interface_used%ini(interior_x_map,interior_y_map)
          call icr_interface_used%ini()


          !first buffer layer
          bf_alignment_tmp = reshape((/
     $         align_W+5,align_N,align_W+6,align_N/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          allocate(added_sublayer%bc_sections(5,1))

          added_sublayer%bc_sections(:,1) = [N_edge_type,1,3,6,no_overlap]

          interior_nodes(align_W+3:align_W+8,align_N-dct_icr_distance,1) = 
     $         [A,A,A,A,A,A]


          !second buffer layer
          bf_alignment_tmp = reshape((/
     $         align_W+17,align_N,align_W+19,align_N/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          allocate(added_sublayer%bc_sections(5,1))

          added_sublayer%bc_sections(:,1) = [N_edge_type,1,3,7,no_overlap]

          interior_nodes(align_W+15:align_W+21,align_N-dct_icr_distance,1) = 
     $         [A,A,A,A,A,A,A]


          !interior_bc_sections
          allocate(interior_bc_sections(5,1))
          interior_bc_sections(:,1) = [N_edge_type,align_W+9,align_N,align_W+14,no_overlap]
          
          interior_nodes(align_W+9:align_W+14,align_N-dct_icr_distance,1) = 
     $         [A,D,A,D,A,D]


          !test_alignment
          test_alignment = reshape((/
     $         align_W+1,align_N,align_W+21,align_N/),
     $         (/2,2/))

        end subroutine get_test_param_analyze_and_update_boundary


        function test_analyze_bf_layer_bc_sections(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_icr)            :: bf_interface_used
          type(icr_interface)               :: icr_interface_used
          real(rkind), dimension(nx)        :: interior_x_map
          real(rkind), dimension(ny)        :: interior_y_map
          real(rkind), dimension(nx,ny,ne)  :: interior_nodes
          integer                           :: mainlayer_id
          type(pmodel_eq)                   :: p_model

          logical                           :: test_loc
          integer       , dimension(2)      :: test_nb_pts
          integer(ikind), dimension(2,18)   :: test_pts
          type(icr_path_chain), pointer     :: icr_path_used


          test_validated = .true.


          !input
          call get_test_param_analyze_bf_layer_bc_sections(
     $         bf_interface_used,
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         mainlayer_id,
     $         test_nb_pts,
     $         test_pts)

          !output
          call bf_interface_used%analyze_bf_layer_bc_sections(
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         N)

          !validation
          test_loc = icr_interface_used%paths%get_nb_paths().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test_nb_pts failed '')'
          end if

          if(test_loc) then

             !test the first path
             icr_path_used => icr_interface_used%paths%get_head_path()
             test_loc = is_int_matrix_validated(
     $            icr_path_used%pts(:,1:test_nb_pts(1)),
     $            test_pts(:,1:test_nb_pts(1)),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test_pts(1) failed'')'
             end if

             !test the second path
             icr_path_used => icr_path_used%get_next()
             test_loc = is_int_matrix_validated(
     $            icr_path_used%pts(:,1:test_nb_pts(2)),
     $            test_pts(:,test_nb_pts(1)+1:test_nb_pts(1)+test_nb_pts(2)),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test_pts(2) failed'')'
             end if

          end if

        end function test_analyze_bf_layer_bc_sections


        subroutine get_test_param_analyze_bf_layer_bc_sections(
     $     bf_interface_used,
     $     icr_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     mainlayer_id,
     $     test_nb_pts,
     $     test_pts)

          implicit none

          
          type(bf_interface_icr)          , intent(inout) :: bf_interface_used
          type(icr_interface)             , intent(inout) :: icr_interface_used
          real(rkind), dimension(nx)      , intent(inout) :: interior_x_map
          real(rkind), dimension(ny)      , intent(inout) :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(inout) :: interior_nodes
          integer                         , intent(out)   :: mainlayer_id
          integer       , dimension(2)    , intent(out)   :: test_nb_pts
          integer(ikind), dimension(2,18) , intent(out)   :: test_pts

          
          integer(ikind), dimension(2,2) :: bf_alignment_tmp
          real(rkind)   , parameter      :: A = -1.0d0
          real(rkind)   , parameter      :: D =  1.0d0
          type(bf_sublayer), pointer     :: added_sublayer


          mainlayer_id = N


          call bf_interface_used%ini(interior_x_map,interior_y_map)
          call icr_interface_used%ini()


          !first buffer layer
          bf_alignment_tmp = reshape((/
     $         align_W+5,align_N,align_W+6,align_N/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          allocate(added_sublayer%bc_sections(5,1))

          added_sublayer%bc_sections(:,1) = [N_edge_type,1,3,6,no_overlap]

          test_nb_pts(1) = 6
          test_pts(:,1:6) = reshape((/
     $         align_W+3,align_N,
     $         align_W+4,align_N,
     $         align_W+5,align_N,
     $         align_W+6,align_N,
     $         align_W+7,align_N,
     $         align_W+8,align_N/),
     $         (/2,6/))

          interior_nodes(align_W+3:align_W+8,align_N-dct_icr_distance,1) = 
     $         [A,A,A,A,A,A]


          !second buffer layer
          bf_alignment_tmp = reshape((/
     $         align_W+17,align_N,align_W+19,align_N/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          allocate(added_sublayer%bc_sections(5,1))

          added_sublayer%bc_sections(:,1) = [N_edge_type,1,3,7,no_overlap]

          test_nb_pts(2) = 7
          test_pts(:,7:13) = reshape((/
     $         align_W+15,align_N,
     $         align_W+16,align_N,
     $         align_W+17,align_N,
     $         align_W+18,align_N,
     $         align_W+19,align_N,
     $         align_W+20,align_N,
     $         align_W+21,align_N/),
     $         (/2,7/))


          interior_nodes(align_W+15:align_W+21,align_N-dct_icr_distance,1) = 
     $         [A,A,A,A,A,A,A]
          

        end subroutine get_test_param_analyze_bf_layer_bc_sections


        function test_analyze_interior_bc_sections(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_icr)                        :: bf_interface_used
          type(icr_interface)                           :: icr_interface_used
          real(rkind), dimension(nx)                    :: interior_x_map
          real(rkind), dimension(ny)                    :: interior_y_map
          real(rkind), dimension(nx,ny,ne)              :: interior_nodes
          type(pmodel_eq)                               :: p_model
          integer(ikind), dimension(:,:)  , allocatable :: interior_bc_sections

          logical                         :: test_loc
          integer       , dimension(2)    :: test_nb_pts
          integer(ikind), dimension(2,18) :: test_pts
          type(icr_path_chain), pointer   :: icr_path_used


          allocate(interior_bc_sections(5,2))


          test_validated = .true.


          !input
          call get_test_param_analyze_bc_section_edge_x(
     $         1,1,
     $         bf_interface_used,
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         interior_bc_sections(:,1),
     $         test_nb_pts(1),
     $         test_pts(:,1:9))

          call get_test_param_analyze_bc_section_edge_y(
     $         1,1,
     $         bf_interface_used,
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         interior_bc_sections(:,2),
     $         test_nb_pts(2),
     $         test_pts(:,10:18))

          !output
          call bf_interface_used%analyze_interior_bc_sections(
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         N,
     $         interior_bc_sections)

          !validation
          test_loc = icr_interface_used%paths%get_nb_paths().eq.1
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test_nb_pts failed '')'
          end if

          if(test_loc) then
             icr_path_used => icr_interface_used%paths%get_head_path()
             test_loc = is_int_matrix_validated(
     $            icr_path_used%pts(:,1:test_nb_pts(2)),
     $            test_pts(:,test_nb_pts(1)+1:test_nb_pts(1)+test_nb_pts(2)),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test_pts failed'')'
             end if
          end if

          icr_path_used => icr_interface_used%paths%get_head_path()
          if(associated(icr_path_used)) then
             call icr_interface_used%paths%remove_path(icr_path_used)
          end if

        end function test_analyze_interior_bc_sections


        function test_get_interior_bcs_mainlayer_id(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated 


          type(bf_interface_icr)       :: bf_interface_used
          integer(ikind), dimension(5) :: bc_section
          integer                      :: mainlayer_id
          integer                      :: test_mainlayer_id

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,4
             
             !input
             call get_param_test_interior_bcs(
     $            k,
     $            bc_section,
     $            test_mainlayer_id)

             !output
             mainlayer_id = bf_interface_used%get_interior_bcs_mainlayer_id(
     $            bc_section)

             !validation
             test_loc = mainlayer_id.eq.test_mainlayer_id
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test ('',I2,'') failed'')',k
             end if             

          end do          

        end function test_get_interior_bcs_mainlayer_id


        subroutine get_param_test_interior_bcs(
     $     test_id,
     $     bc_section,
     $     test_mainlayer_id)

          implicit none

          integer                     , intent(in)  :: test_id
          integer(ikind), dimension(5), intent(out) :: bc_section
          integer                     , intent(out) :: test_mainlayer_id


          select case(test_id)
            case(1)
               bc_section = [N_edge_type,5,ny-1,no_overlap,no_overlap]
               test_mainlayer_id = N

            case(2)
               bc_section = [S_edge_type,4,1,no_overlap,no_overlap]
               test_mainlayer_id = S

            case(3)
               bc_section = [E_edge_type,align_E,5,no_overlap,no_overlap]
               test_mainlayer_id = E

            case(4)
               bc_section = [W_edge_type,1,6,no_overlap,no_overlap]
               test_mainlayer_id = W
               
          end select

        end subroutine get_param_test_interior_bcs


        function test_analyze_bc_section(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated 


          type(bf_interface_icr)           :: bf_interface_used
          type(icr_interface)              :: icr_interface_used
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          type(pmodel_eq)                  :: p_model
          integer(ikind), dimension(5)     :: bc_section

          logical                        :: test_loc
          integer                        :: test_nb_pts
          integer(ikind), dimension(2,9) :: test_pts
          type(icr_path_chain), pointer  :: icr_path_used

          type(bf_sublayer), pointer   :: bf_sublayer_ptr
          integer(ikind), dimension(2) :: match_table
          logical                      :: interior_domain


          test_validated = .true.


          interior_domain = .true.


          !test edge_x bc_section
          !============================================================
          !input
          call get_test_param_analyze_bc_section_edge_x(
     $         1,1,
     $         bf_interface_used,
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bc_section,
     $         test_nb_pts,
     $         test_pts)

          !output
          call bf_interface_used%analyze_bc_section(
     $         icr_interface_used,
     $         bf_sublayer_ptr,
     $         match_table,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         interior_domain,
     $         bc_section)

          !validation
          if(test_nb_pts.gt.0) then
             test_loc = icr_interface_used%paths%get_nb_paths().eq.1
          else
             test_loc = icr_interface_used%paths%get_nb_paths().eq.0
          end if
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test_nb_pts(1) failed '')'
          end if

          if(test_loc.and.(test_nb_pts.gt.0)) then
             icr_path_used => icr_interface_used%paths%get_head_path()
             test_loc = is_int_matrix_validated(
     $            icr_path_used%pts(:,1:test_nb_pts),
     $            test_pts(:,1:test_nb_pts),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test_pts(1)  failed'')'
             end if
          end if

          icr_path_used => icr_interface_used%paths%get_head_path()
          if(associated(icr_path_used)) then
             call icr_interface_used%paths%remove_path(icr_path_used)
          end if


          !test edge_y bc_section
          !============================================================
          !input
          call get_test_param_analyze_bc_section_edge_y(
     $         1,1,
     $         bf_interface_used,
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bc_section,
     $         test_nb_pts,
     $         test_pts)

          !output
          call bf_interface_used%analyze_bc_section(
     $         icr_interface_used,
     $         bf_sublayer_ptr,
     $         match_table,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         interior_domain,
     $         bc_section)

          !validation
          if(test_nb_pts.gt.0) then
             test_loc = icr_interface_used%paths%get_nb_paths().eq.1
          else
             test_loc = icr_interface_used%paths%get_nb_paths().eq.0
          end if
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test_nb_pts(2) failed '')'
          end if

          if(test_loc.and.(test_nb_pts.gt.0)) then
             icr_path_used => icr_interface_used%paths%get_head_path()
             test_loc = is_int_matrix_validated(
     $            icr_path_used%pts(:,1:test_nb_pts),
     $            test_pts(:,1:test_nb_pts),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test_pts(2) failed'')'
             end if
          end if

          icr_path_used => icr_interface_used%paths%get_head_path()
          if(associated(icr_path_used)) then
             call icr_interface_used%paths%remove_path(icr_path_used)
          end if


          !test square bc_section
          !============================================================
          !input
          call get_test_param_analyze_square_xy(
     $         5,
     $         bf_interface_used,
     $         icr_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bc_section,
     $         test_nb_pts,
     $         test_pts)

          !output
          call bf_interface_used%analyze_bc_section(
     $         icr_interface_used,
     $         bf_sublayer_ptr,
     $         match_table,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         interior_domain,
     $         bc_section)

         !validation
         if(test_nb_pts.gt.0) then
            test_loc = icr_interface_used%paths%get_nb_paths().eq.1
         else
            test_loc = icr_interface_used%paths%get_nb_paths().eq.0
         end if
         test_validated = test_validated.and.test_loc
         if(detailled.and.(.not.test_loc)) then
            print '(''test_nb_pts(3) failed '')'
            print '(''nb_paths: '',I2)', icr_interface_used%paths%get_nb_paths()
         end if

         if(test_loc.and.(test_nb_pts.gt.0)) then
            icr_path_used => icr_interface_used%paths%get_head_path()
            test_loc = is_int_matrix_validated(
     $           icr_path_used%pts(:,1:test_nb_pts),
     $           test_pts(:,1:test_nb_pts),
     $           detailled)
            test_validated = test_validated.and.test_loc
            if(detailled.and.(.not.test_loc)) then
               print '(''test_pts(3) failed'')'
            end if
         end if

         icr_path_used => icr_interface_used%paths%get_head_path()
         if(associated(icr_path_used)) then
            call icr_interface_used%paths%remove_path(icr_path_used)
         end if         

        end function test_analyze_bc_section


        function test_analyze_bc_section_square_xy(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated 


          type(bf_interface_icr)           :: bf_interface_used
          type(icr_interface)              :: icr_interface_used
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          type(pmodel_eq)                  :: p_model
          integer(ikind), dimension(5)     :: bc_section

          integer                        :: k
          logical                        :: test_loc
          integer                        :: test_nb_pts
          integer(ikind), dimension(2,9) :: test_pts
          type(icr_path_chain), pointer  :: icr_path_used

          type(bf_sublayer), pointer   :: bf_sublayer_ptr
          integer(ikind), dimension(2) :: match_table
          logical                      :: interior_domain


          test_validated = .true.


          interior_domain = .true.

          do k=1,11

             !input
             call get_test_param_analyze_square_xy(
     $            k,
     $            bf_interface_used,
     $            icr_interface_used,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            bc_section,
     $            test_nb_pts,
     $            test_pts)

             !output
             call bf_interface_used%analyze_bc_section_square_xy(
     $            icr_interface_used,
     $            bf_sublayer_ptr,
     $            match_table,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            p_model,
     $            interior_domain,
     $            bc_section)

             !validation
             if(test_nb_pts.gt.0) then
                test_loc = icr_interface_used%paths%get_nb_paths().eq.1
             else
                test_loc = icr_interface_used%paths%get_nb_paths().eq.0
             end if
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test_nb_pts('',I2,'') failed '')',k
                print '(''nb_paths: '',I2)', icr_interface_used%paths%get_nb_paths()
             end if

             if(test_loc.and.(test_nb_pts.gt.0)) then
                icr_path_used => icr_interface_used%paths%get_head_path()
                test_loc = is_int_matrix_validated(
     $               icr_path_used%pts(:,1:test_nb_pts),
     $               test_pts(:,1:test_nb_pts),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test_pts('',I2,'') failed'')',k
                end if
             end if

             icr_path_used => icr_interface_used%paths%get_head_path()
             if(associated(icr_path_used)) then
                call icr_interface_used%paths%remove_path(icr_path_used)
             end if

           end do

        end function test_analyze_bc_section_square_xy


        subroutine get_test_param_analyze_square_xy(
     $     test_id,
     $     bf_interface_used,
     $     icr_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     bc_section,
     $     test_nb_pts,
     $     test_pts)

          implicit none

          integer                         , intent(in)  :: test_id
          type(bf_interface_icr)          , intent(out) :: bf_interface_used
          type(icr_interface)             , intent(out) :: icr_interface_used
          real(rkind), dimension(nx)      , intent(out) :: interior_x_map
          real(rkind), dimension(ny)      , intent(out) :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes
          integer(ikind), dimension(5)    , intent(out) :: bc_section
          integer                         , intent(out) :: test_nb_pts
          integer(ikind), dimension(2,9)  , intent(out) :: test_pts

          integer(ikind)                 :: i,j,k
          integer(ikind), dimension(2,3) :: gen_coords
          real(rkind)   , dimension(3)   :: nodes_set
          
          real(rkind), parameter :: A = -1.0d0
          real(rkind), parameter :: D =  1.0d0


          call bf_interface_used%ini(interior_x_map,interior_y_map)
          call icr_interface_used%ini()


          interior_x_map = (/ ((i-1)*1.0d0, i=1,nx) /)
          interior_y_map = (/ ((j-1)*1.0d0, j=1,ny) /)

          interior_nodes = reshape((/
     $         (((D,i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          bc_section = [NW_corner_type,1,ny-1,no_overlap,no_overlap]

          gen_coords(:,1) = [ 2+dct_icr_distance, ny-1-dct_icr_distance-1]
          gen_coords(:,2) = [ 2+dct_icr_distance, ny-1-dct_icr_distance  ]
          gen_coords(:,3) = [ 3+dct_icr_distance, ny-1-dct_icr_distance  ]
          

          select case(test_id)

            !all pts activated
            case(1)
               nodes_set = [A,A,A]

               test_nb_pts = 1
               test_pts(:,1) = [2,ny-1]

            !all pts desactivated
            case(2)
               nodes_set = [D,D,D]
               
               test_nb_pts = 0

            !only one point activated
            case(3)
               nodes_set = [A,D,D]
               
               test_nb_pts = 1
               test_pts(:,1) = [2,ny-1]

            !only one point activated
            case(4)
               nodes_set = [D,A,D]
               
               test_nb_pts = 1
               test_pts(:,1) = [2,ny-1]

            !only one point activated
            case(5)
               nodes_set = [D,D,A]
               
               test_nb_pts = 1
               test_pts(:,1) = [2,ny-1]

            !all activated in spite of the overlap
            case(6)
               nodes_set = [A,A,A]
               bc_section(5) = N_overlap

               test_nb_pts = 1
               test_pts(:,1) = [2,ny-1]

            !all activated in spite of the overlap
            case(7)
               nodes_set = [A,A,A]
               bc_section(5) = W_overlap

               test_nb_pts = 1
               test_pts(:,1) = [2,ny-1]

            !no point activated because of the overlap
            case(8)
               nodes_set = [A,A,A]
               bc_section(5) = E_overlap

               test_nb_pts = 0

            !no point activated because of the overlap
            case(9)
               nodes_set = [A,A,A]
               bc_section(5) = S_overlap

               test_nb_pts = 0

            !no point activated because of the overlap
            case(10)
               nodes_set = [A,A,A]
               bc_section(4) = cpt2not_and_cpt3normal

               test_nb_pts = 0

            !no point activated because of the overlap
            case(11)
               nodes_set = [A,A,A]
               bc_section(4) = cpt2not_and_cpt3not

               test_nb_pts = 0

          end select

          do k=1,3
             interior_nodes(gen_coords(1,k),gen_coords(2,k),1) = nodes_set(k)
          end do


        end subroutine get_test_param_analyze_square_xy

        
        function test_analyze_bc_section_square_bounds(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_interface_icr) :: bf_interface_used
          real(rkind), dimension(nx) :: interior_x_map
          real(rkind), dimension(ny) :: interior_y_map          

          integer(ikind), dimension(5)     :: bc_section

          integer(ikind), dimension(2,2,2) :: analyzed_grdpts_bounds
          integer(ikind), dimension(2,3)   :: activated_grdpts
          integer                          :: nb_activated_grdpts
          logical                          :: no_activation

          integer(ikind), dimension(2,2,2) :: test_analyzed_grdpts_bounds
          integer(ikind), dimension(2,3)   :: test_activated_grdpts
          integer                          :: test_nb_activated_grdpts
          logical                          :: test_no_activation

          integer :: k
          logical :: test_loc


          test_validated = .true.


          call bf_interface_used%ini(interior_x_map,interior_y_map)


          do k=1,8

             !input
             call get_test_param_square_bounds(
     $            k,
     $            bc_section,
     $            test_analyzed_grdpts_bounds,
     $            test_activated_grdpts,
     $            test_nb_activated_grdpts,
     $            test_no_activation)

             !output
             call bf_interface_used%analyze_bc_section_square_bounds(
     $            bc_section,
     $            analyzed_grdpts_bounds,
     $            activated_grdpts,
     $            nb_activated_grdpts,
     $            no_activation)

             !validation
             !no_activation
             test_loc = no_activation.eqv.test_no_activation
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test no_activation failed'')'
             end if

             if(test_loc.and.(.not.no_activation)) then

                !analyze_grdpts_bounds
                test_loc = is_int_matrix3D_validated(
     $               analyzed_grdpts_bounds,
     $               test_analyzed_grdpts_bounds,
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test analyze_grdpts_bounds failed'')'
                end if

                !nb_activated_grdpts
                test_loc = nb_activated_grdpts.eq.test_nb_activated_grdpts
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test nb_activated_grdpts failed'')'
                end if

                !activated_grdpts
                test_loc = is_int_matrix_validated(
     $               activated_grdpts(:,1:test_nb_activated_grdpts),
     $               test_activated_grdpts(:,1:test_nb_activated_grdpts),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test activated_grdpts failed'')'
                end if

             end if

          end do

        end function test_analyze_bc_section_square_bounds


        subroutine get_test_param_square_bounds(
     $     test_id,
     $     bc_section,
     $     test_analyzed_grdpts_bounds,
     $     test_activated_grdpts,
     $     test_nb_activated_grdpts,
     $     test_no_activation)

          implicit none

          integer                         , intent(in)  :: test_id
          integer(ikind), dimension(5)    , intent(out) :: bc_section
          integer(ikind), dimension(2,2,2), intent(out) :: test_analyzed_grdpts_bounds
          integer(ikind), dimension(2,3)  , intent(out) :: test_activated_grdpts
          integer                         , intent(out) :: test_nb_activated_grdpts
          logical                         , intent(out) :: test_no_activation

          
          bc_section(4) = no_overlap
          bc_section(5) = no_overlap

          test_no_activation = .false.


          select case(test_id)
            case(1)
               bc_section(1) = NW_corner_type
               bc_section(2) = 1
               bc_section(3) = ny-1

               test_analyzed_grdpts_bounds = reshape((/
     $              2+dct_icr_distance,
     $              ny-1-dct_icr_distance-1,
     $              2+dct_icr_distance,
     $              ny-1-dct_icr_distance-1,
     $              2+dct_icr_distance,
     $              ny-1-dct_icr_distance,
     $              2+dct_icr_distance+1,
     $              ny-1-dct_icr_distance/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 1

               test_activated_grdpts(:,1) = [2,ny-1]


            case(2)
               bc_section(1) = NE_corner_type
               bc_section(2) = nx-1
               bc_section(3) = ny-1

               test_analyzed_grdpts_bounds = reshape((/
     $              nx-1-dct_icr_distance,
     $              ny-1-dct_icr_distance-1,
     $              nx-1-dct_icr_distance,
     $              ny-1-dct_icr_distance-1,
     $              nx-1-dct_icr_distance-1,
     $              ny-1-dct_icr_distance,
     $              nx-1-dct_icr_distance,
     $              ny-1-dct_icr_distance/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 1

               test_activated_grdpts(:,1) = [nx-1,ny-1]

            case(3)
               bc_section(1) = SW_corner_type
               bc_section(2) = 1
               bc_section(3) = 1

               test_analyzed_grdpts_bounds = reshape((/
     $              2+dct_icr_distance,
     $              2+dct_icr_distance,
     $              2+dct_icr_distance+1,
     $              2+dct_icr_distance,
     $              
     $              2+dct_icr_distance,
     $              2+dct_icr_distance+1,
     $              2+dct_icr_distance,
     $              2+dct_icr_distance+1/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 1

               test_activated_grdpts(:,1) = [2,2]

            case(4)
               bc_section(1) = SE_corner_type
               bc_section(2) = nx-1
               bc_section(3) = 1

               test_analyzed_grdpts_bounds = reshape((/
     $              nx-1-dct_icr_distance-1,
     $              2+dct_icr_distance,
     $              nx-1-dct_icr_distance,
     $              2+dct_icr_distance,
     $              
     $              nx-1-dct_icr_distance,
     $              2+dct_icr_distance+1,
     $              nx-1-dct_icr_distance,
     $              2+dct_icr_distance+1/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 1

               test_activated_grdpts(:,1) = [nx-1,2]

            case(5)
               bc_section(1) = NW_edge_type
               bc_section(2) = 1
               bc_section(3) = ny-1

               test_analyzed_grdpts_bounds = reshape((/
     $              0,
     $              ny-1-dct_icr_distance,
     $              2+dct_icr_distance,
     $              ny-1-dct_icr_distance,
     $              
     $              2+dct_icr_distance,
     $              ny-1-dct_icr_distance+1,
     $              2+dct_icr_distance,
     $              ny+1/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 3

               test_activated_grdpts(:,1) = [1,ny-1]
               test_activated_grdpts(:,2) = [2,ny-1]
               test_activated_grdpts(:,3) = [2,ny  ]
               

            case(6)
               bc_section(1) = NE_edge_type
               bc_section(2) = nx-1
               bc_section(3) = ny-1

               test_analyzed_grdpts_bounds = reshape((/
     $              nx-1-dct_icr_distance,
     $              ny-1-dct_icr_distance,
     $              nx+1,
     $              ny-1-dct_icr_distance,
     $              
     $              nx-1-dct_icr_distance,
     $              ny-1-dct_icr_distance+1,
     $              nx-1-dct_icr_distance,
     $              ny+1/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 3

               test_activated_grdpts(:,1) = [nx-1,ny-1]
               test_activated_grdpts(:,2) = [nx  ,ny-1]
               test_activated_grdpts(:,3) = [nx-1,ny  ]

            case(7)
               bc_section(1) = SW_edge_type
               bc_section(2) = 1
               bc_section(3) = 1

               test_analyzed_grdpts_bounds = reshape((/
     $              2+dct_icr_distance,
     $              0,
     $              2+dct_icr_distance,
     $              2+dct_icr_distance-1,
     $              
     $              0,
     $              2+dct_icr_distance,
     $              2+dct_icr_distance,
     $              2+dct_icr_distance/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 3

               test_activated_grdpts(:,1) = [2,1]
               test_activated_grdpts(:,2) = [1,2]
               test_activated_grdpts(:,3) = [2,2]

            case(8)
               bc_section(1) = SE_edge_type
               bc_section(2) = nx-1
               bc_section(3) = 1

               test_analyzed_grdpts_bounds = reshape((/
     $              nx-1-dct_icr_distance,
     $              0,
     $              nx-1-dct_icr_distance,
     $              2+dct_icr_distance-1,
     $              
     $              nx-1-dct_icr_distance,
     $              2+dct_icr_distance,
     $              nx+1,
     $              2+dct_icr_distance/),
     $              (/2,2,2/))

               test_nb_activated_grdpts = 3

               test_activated_grdpts(:,1) = [nx-1,1]
               test_activated_grdpts(:,2) = [nx-1,2]
               test_activated_grdpts(:,3) = [nx  ,2]

          end select

        end subroutine get_test_param_square_bounds
        

        function test_analyze_bc_section_edge_x(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_icr)           :: bf_interface_used
          type(icr_interface)              :: icr_interface_used
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          type(pmodel_eq)                  :: p_model
          integer(ikind), dimension(5)     :: bc_section

          type(bf_sublayer), pointer   :: bf_sublayer_ptr
          integer(ikind), dimension(2) :: match_table
          logical                      :: interior_domain

          integer                        :: k,l
          logical                        :: test_loc
          integer                        :: test_nb_pts
          integer(ikind), dimension(2,9) :: test_pts
          type(icr_path_chain), pointer  :: icr_path_used


          test_validated = .true.


          interior_domain = .true.

          do l=1,2
             do k=1,8

                !input
                call get_test_param_analyze_bc_section_edge_x(
     $               k,l,
     $               bf_interface_used,
     $               icr_interface_used,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               bc_section,
     $               test_nb_pts,
     $               test_pts)

                !output
                call bf_interface_used%analyze_bc_section_edge_x(
     $               icr_interface_used,
     $               bf_sublayer_ptr,
     $               match_table,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               p_model,
     $               interior_domain,
     $               bc_section)

                !validation
                if(test_nb_pts.gt.0) then
                   test_loc = icr_interface_used%paths%get_nb_paths().eq.1
                else
                   test_loc = icr_interface_used%paths%get_nb_paths().eq.0
                end if
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test_nb_pts('',2I2,'') failed '')',k,l
                end if

                if(test_loc.and.(test_nb_pts.gt.0)) then
                   icr_path_used => icr_interface_used%paths%get_head_path()
                   test_loc = is_int_matrix_validated(
     $                  icr_path_used%pts(:,1:test_nb_pts),
     $                  test_pts(:,1:test_nb_pts),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print '(''test_pts('',2I2,'') failed'')',k,l
                   end if
                end if

                icr_path_used => icr_interface_used%paths%get_head_path()
                if(associated(icr_path_used)) then
                   call icr_interface_used%paths%remove_path(icr_path_used)
                end if

             end do
          end do

        end function test_analyze_bc_section_edge_x


        subroutine get_test_param_analyze_bc_section_edge_x(
     $     test_id,
     $     edge_type,
     $     bf_interface_used,
     $     icr_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     bc_section,
     $     test_nb_pts,
     $     test_pts)

          implicit none

          integer                         , intent(in)  :: test_id
          integer                         , intent(in)  :: edge_type
          type(bf_interface_icr)          , intent(out) :: bf_interface_used
          type(icr_interface)             , intent(out) :: icr_interface_used
          real(rkind), dimension(nx)      , intent(out) :: interior_x_map
          real(rkind), dimension(ny)      , intent(out) :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes
          integer(ikind), dimension(5)    , intent(out) :: bc_section
          integer                         , intent(out) :: test_nb_pts
          integer(ikind), dimension(2,9)  , intent(out) :: test_pts

          real(rkind), dimension(ne) :: unactivated_node
          real(rkind), dimension(ne) :: activated_node

          integer(ikind) :: i,j,k
          integer(ikind) :: j_s
          integer :: i_pt
          integer :: i_activated
          
          real(rkind), parameter :: A = -1.0d0
          real(rkind), parameter :: D =  1.0d0


          call bf_interface_used%ini(interior_x_map,interior_y_map)
          call icr_interface_used%ini()


          interior_x_map = (/ ((i-1)*1.0d0, i=1,nx) /)
          interior_y_map = (/ ((j-1)*1.0d0, j=1,ny) /)

          interior_nodes = reshape((/
     $         (((0.0d0,i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          j_s = ny/2

          select case(edge_type)
            case(1)
               i_pt = nx-1
               bc_section = [E_edge_type,i_pt  ,j_s,j_s+8,no_overlap]
               i_activated = i_pt-dct_icr_distance
            case(2)
               i_pt = 2
               bc_section = [W_edge_type,i_pt-1,j_s,j_s+8,no_overlap]
               i_activated = i_pt+dct_icr_distance
          end select

          unactivated_node = [ D,0.0d0,0.0d0]
          activated_node   = [ A,0.0d0,0.0d0]
          

          select case(test_id)

            !all pts activated
            case(1)
               interior_nodes(i_activated,j_s:j_s+8,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 9

               test_pts(1,:) = i_pt
               test_pts(2,:) = [j_s,j_s+1,j_s+2,j_s+3,j_s+4,
     $                          j_s+5,j_s+6,j_s+7,j_s+8]

            !all pts desactivated
            case(2)
               interior_nodes(i_activated,j_s:j_s+8,1) = [D,D,D,D,D,D,D,D,D]
               
               test_nb_pts = 0

            !one out of two points activated
            case(3)
               interior_nodes(i_activated,j_s:j_s+8,1) = [A,D,A,D,A,D,A,D,A]
               
               test_nb_pts = 9

               test_pts(1,:) = i_pt
               test_pts(2,:) = [j_s,j_s+1,j_s+2,j_s+3,j_s+4,
     $                          j_s+5,j_s+6,j_s+7,j_s+8]

            !one out of three points activated
            case(4)
               interior_nodes(i_activated,j_s:j_s+8,1) = [A,D,D,A,D,D,A,D,D]
               
               test_nb_pts = 8
               test_pts(1,:)   = i_pt
               test_pts(2,1:8) = [j_s,j_s+1,j_s+2,j_s+3,j_s+4,
     $                            j_s+5,j_s+6,j_s+7]
               

            !one out of four points activated
            case(5)
               interior_nodes(i_activated,j_s:j_s+8,1) = [A,D,D,D,A,D,D,D,A]
               
               test_nb_pts = 7

               test_pts(1,:)   = i_pt
               test_pts(2,1:7) = [j_s,j_s+1,j_s+3,j_s+4,
     $                            j_s+5,j_s+7,j_s+8]

            !all activated in spite of the overlap
            case(6)
               select case(edge_type)
                 case(1)
                    bc_section(5) = E_overlap
                 case(2)
                    bc_section(5) = W_overlap
               end select

               interior_nodes(i_activated,j_s:j_s+8,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 9

               test_pts(1,:) = i_pt
               test_pts(2,:) = [j_s,j_s+1,j_s+2,j_s+3,j_s+4,
     $                          j_s+5,j_s+6,j_s+7,j_s+8]

            !no point activated because of the overlap
            case(7)
               select case(edge_type)
                 case(1)
                    bc_section(5) = W_overlap
                 case(2)
                    bc_section(5) = E_overlap
               end select

               interior_nodes(i_activated,j_s:j_s+8,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 0

            !no point activated because of the overlap
            case(8)
               bc_section(5) = EW_overlap

               interior_nodes(i_activated,j_s:j_s+8,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 0

          end select               

        end subroutine get_test_param_analyze_bc_section_edge_x


        function test_analyze_bc_section_edge_y(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_icr)           :: bf_interface_used
          type(icr_interface)              :: icr_interface_used
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          type(pmodel_eq)                  :: p_model
          integer(ikind), dimension(5)     :: bc_section

          type(bf_sublayer), pointer   :: bf_sublayer_ptr
          integer(ikind), dimension(2) :: match_table
          logical                      :: interior_domain

          integer                        :: k,l
          logical                        :: test_loc
          integer                        :: test_nb_pts
          integer(ikind), dimension(2,9) :: test_pts
          type(icr_path_chain), pointer  :: icr_path_used


          test_validated = .true.


          interior_domain = .true.

          do l=1,2
             do k=1,8

                !input
                call get_test_param_analyze_bc_section_edge_y(
     $               k,l,
     $               bf_interface_used,
     $               icr_interface_used,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               bc_section,
     $               test_nb_pts,
     $               test_pts)

                !output
                call bf_interface_used%analyze_bc_section_edge_y(
     $               icr_interface_used,
     $               bf_sublayer_ptr,
     $               match_table,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               p_model,
     $               interior_domain,
     $               bc_section)

                !validation
                if(test_nb_pts.gt.0) then
                   test_loc = icr_interface_used%paths%get_nb_paths().eq.1
                else
                   test_loc = icr_interface_used%paths%get_nb_paths().eq.0
                end if
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test_nb_pts('',2I2,'') failed '')',k,l
                end if

                if(test_loc.and.(test_nb_pts.gt.0)) then
                   icr_path_used => icr_interface_used%paths%get_head_path()
                   test_loc = is_int_matrix_validated(
     $                  icr_path_used%pts(:,1:test_nb_pts),
     $                  test_pts(:,1:test_nb_pts),
     $                  detailled)
                   test_validated = test_validated.and.test_loc
                   if(detailled.and.(.not.test_loc)) then
                      print '(''test_pts('',2I2,'') failed'')',k,l
                   end if
                end if

                icr_path_used => icr_interface_used%paths%get_head_path()
                if(associated(icr_path_used)) then
                   call icr_interface_used%paths%remove_path(icr_path_used)
                end if

             end do
          end do

        end function test_analyze_bc_section_edge_y


        subroutine get_test_param_analyze_bc_section_edge_y(
     $     test_id,
     $     edge_type,
     $     bf_interface_used,
     $     icr_interface_used,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     bc_section,
     $     test_nb_pts,
     $     test_pts)

          implicit none

          integer                         , intent(in)  :: test_id
          integer                         , intent(in)  :: edge_type
          type(bf_interface_icr)          , intent(out) :: bf_interface_used
          type(icr_interface)             , intent(out) :: icr_interface_used
          real(rkind), dimension(nx)      , intent(out) :: interior_x_map
          real(rkind), dimension(ny)      , intent(out) :: interior_y_map
          real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes
          integer(ikind), dimension(5)    , intent(out) :: bc_section
          integer                         , intent(out) :: test_nb_pts
          integer(ikind), dimension(2,9)  , intent(out) :: test_pts

          real(rkind), dimension(ne) :: unactivated_node
          real(rkind), dimension(ne) :: activated_node

          integer(ikind) :: i,j,k
          integer(ikind) :: i_s
          integer :: j_pt
          integer :: j_activated
          
          real(rkind), parameter :: A = -1.0d0
          real(rkind), parameter :: D =  1.0d0


          call bf_interface_used%ini(interior_x_map,interior_y_map)
          call icr_interface_used%ini()


          interior_x_map = (/ ((i-1)*1.0d0, i=1,nx) /)
          interior_y_map = (/ ((j-1)*1.0d0, j=1,ny) /)

          interior_nodes = reshape((/
     $         (((0.0d0,i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          i_s = nx/2

          select case(edge_type)
            case(1)
               j_pt = ny-1
               bc_section = [N_edge_type,i_s,j_pt,i_s+8,no_overlap]
               j_activated = j_pt-dct_icr_distance
            case(2)
               j_pt = 2
               bc_section = [S_edge_type,i_s,j_pt-1,i_s+8,no_overlap]
               j_activated = j_pt+dct_icr_distance
          end select

          unactivated_node = [ D,0.0d0,0.0d0]
          activated_node   = [ A,0.0d0,0.0d0]
          

          select case(test_id)

            !all pts activated
            case(1)
               interior_nodes(i_s:i_s+8,j_activated,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 9
               test_pts(1,:) = [i_s,i_s+1,i_s+2,i_s+3,i_s+4,
     $                          i_s+5,i_s+6,i_s+7,i_s+8]
               test_pts(2,:) = j_pt

            !all pts desactivated
            case(2)
               interior_nodes(i_s:i_s+8,j_activated,1) = [D,D,D,D,D,D,D,D,D]
               
               test_nb_pts = 0

            !one out of two points activated
            case(3)
               interior_nodes(i_s:i_s+8,j_activated,1) = [A,D,A,D,A,D,A,D,A]
               
               test_nb_pts = 9
               test_pts(1,:) = [i_s,i_s+1,i_s+2,i_s+3,i_s+4,
     $                          i_s+5,i_s+6,i_s+7,i_s+8]
               test_pts(2,:) = j_pt

            !one out of three points activated
            case(4)
               interior_nodes(i_s:i_s+8,j_activated,1) = [A,D,D,A,D,D,A,D,D]
               
               test_nb_pts = 8
               test_pts(1,1:8) = [i_s,i_s+1,i_s+2,i_s+3,i_s+4,
     $                          i_s+5,i_s+6,i_s+7]
               test_pts(2,:) = j_pt

            !one out of four points activated
            case(5)
               interior_nodes(i_s:i_s+8,j_activated,1) = [A,D,D,D,A,D,D,D,A]
               
               test_nb_pts = 7
               test_pts(1,1:7) = [i_s,i_s+1,i_s+3,i_s+4,
     $                            i_s+5,i_s+7,i_s+8]
               test_pts(2,:) = j_pt

            !all activated in spite of the overlap
            case(6)
               select case(edge_type)
                 case(1)
                    bc_section(5) = N_overlap
                 case(2)
                    bc_section(5) = S_overlap
               end select

               interior_nodes(i_s:i_s+8,j_activated,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 9

               test_pts(1,:) = [i_s,i_s+1,i_s+2,i_s+3,i_s+4,
     $                          i_s+5,i_s+6,i_s+7,i_s+8]
               test_pts(2,:) = j_pt

            !no point activated because of the overlap
            case(7)
               select case(edge_type)
                 case(1)
                    bc_section(5) = S_overlap
                 case(2)
                    bc_section(5) = N_overlap
               end select

               interior_nodes(i_s:i_s+8,j_activated,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 0

           !no point activated because of the overlap
            case(8)
               bc_section(5) = NS_overlap

               interior_nodes(i_s:i_s+8,j_activated,1) = [A,A,A,A,A,A,A,A,A]

               test_nb_pts = 0

          end select               

        end subroutine get_test_param_analyze_bc_section_edge_y


        function test_stage(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_interface_icr) :: bf_interface_used
          type(icr_interface)    :: icr_interface_used

          real(rkind)   , dimension(nx)       :: interior_x_map
          real(rkind)   , dimension(ny)       :: interior_y_map
          integer(ikind), dimension(2)        :: match_table
          type(icr_path_chain), pointer       :: head_path


          !input
          call bf_interface_used%ini(interior_x_map, interior_y_map)
          call icr_interface_used%ini()
          match_table = [-1,2]


          !output
          call bf_interface_used%stage(
     $         icr_interface_used,
     $         [2,3],
     $         match_table,
     $         .true.)

          call bf_interface_used%stage(
     $         icr_interface_used,
     $         [2,3],
     $         match_table,
     $         .false.)

          
          !validation
          head_path => icr_interface_used%paths%get_head_path()
          test_validated = is_int_matrix_validated(
     $         head_path%pts(:,1:2),
     $         reshape((/
     $            2,3,
     $            1,5/),
     $         (/2,2/)),
     $         detailled)

        end function test_stage


        function test_is_node_activated(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_interface_icr)              :: bf_interface_used
          type(bf_sublayer), pointer          :: bf_sublayer_ptr
          integer(ikind), dimension(2)        :: match_table
          type(pmodel_eq)                     :: p_model
          real(rkind)   , dimension(nx)       :: interior_x_map
          real(rkind)   , dimension(ny)       :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
          integer(ikind), dimension(2)        :: loc_coords

          integer(ikind), dimension(2,2) :: bf_alignment_tmp

          logical :: test_loc
          
          test_validated = .true.


          !input
          call bf_interface_used%ini(interior_x_map,interior_y_map)
           
          bf_alignment_tmp = reshape((/
     $         align_E,align_S+5,align_E,align_S+5/),
     $         (/2,2/))
          
          bf_sublayer_ptr => bf_interface_used%allocate_sublayer(
     $         E,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)
          
          match_table = bf_sublayer_ptr%get_general_to_local_coord_tab()
           
          interior_nodes(3,3,1)        =  1.0d0
          bf_sublayer_ptr%nodes(3,3,1) = -1.0d0
          loc_coords = [3,3]
          

          !output
          test_loc = .not.(bf_interface_used%is_node_activated(
     $         loc_coords,
     $         bf_sublayer_ptr,
     $         match_table,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         .true.))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test interior failed'')'
          end if


          test_loc = bf_interface_used%is_node_activated(
     $         loc_coords,
     $         bf_sublayer_ptr,
     $         match_table,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         .false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test bf_layer failed'')'
          end if

        end function test_is_node_activated

        
        function test_is_bf_layer_node_activated(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_interface_icr)              :: bf_interface_used
          type(bf_sublayer), pointer          :: bf_sublayer_ptr
          integer(ikind), dimension(2)        :: match_table
          type(pmodel_eq)                     :: p_model
          real(rkind)   , dimension(nx)       :: interior_x_map
          real(rkind)   , dimension(ny)       :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
          integer(ikind), dimension(2)        :: loc_coords
          logical                             :: node_activated
          logical                             :: test_node_activated

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,4

             !input
             call get_param_test_bf_node_activated(
     $            k,
     $            bf_interface_used,
     $            bf_sublayer_ptr,
     $            match_table,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            loc_coords,
     $            test_node_activated)

             !output
             node_activated = bf_interface_used%is_bf_layer_node_activated(
     $            loc_coords,
     $            bf_sublayer_ptr,
     $            match_table,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            p_model)

             !validation
             test_loc = node_activated.eqv.test_node_activated
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test ('',I2,'') failed'')', k
             end if

          end do

        end function test_is_bf_layer_node_activated


        subroutine get_param_test_bf_node_activated(
     $     test_id,
     $     bf_interface_used,
     $     bf_sublayer_ptr,
     $     match_table,
     $     interior_x_map,
     $     interior_y_map,
     $     interior_nodes,
     $     loc_coords,
     $     test_node_activated)

          implicit none

          integer                            , intent(in)    :: test_id
          type(bf_interface_icr)             , intent(inout) :: bf_interface_used
          type(bf_sublayer), pointer         , intent(inout) :: bf_sublayer_ptr
          integer(ikind), dimension(2)       , intent(inout) :: match_table
          real(rkind)   , dimension(nx)      , intent(inout) :: interior_x_map
          real(rkind)   , dimension(ny)      , intent(inout) :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne), intent(inout) :: interior_nodes
          integer(ikind), dimension(2)       , intent(out)   :: loc_coords
          logical                            , intent(out)   :: test_node_activated


          integer(ikind), dimension(2,2) :: bf_alignment_tmp
          type(bf_sublayer), pointer     :: bf_sublayer2_ptr


          select case(test_id)

            !check whether the node is activated with the buffer layer
            case(1)
               call bf_interface_used%ini(interior_x_map,interior_y_map)

               bf_alignment_tmp = reshape((/
     $              align_E,align_S+5,align_E,align_S+5/),
     $              (/2,2/))

               bf_sublayer_ptr => bf_interface_used%allocate_sublayer(
     $              E,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               match_table = bf_sublayer_ptr%get_general_to_local_coord_tab()

               bf_sublayer_ptr%nodes(3,3,1) = -1.0d0

               loc_coords = [3,3]

               test_node_activated = .true.


            !check whether the node is activated by the interior nodes
            case(2)
               
               loc_coords = [1,1]

               interior_nodes(align_E-2,align_S+3,1) = -1.0d0

               test_node_activated = .true.


            !check whether the node is activated by another buffer layer
            case(3)

               bf_alignment_tmp = reshape((/
     $              align_W+5,align_N,align_W+5,align_N/),
     $              (/2,2/))

               bf_sublayer2_ptr => bf_interface_used%allocate_sublayer(
     $              N,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               bf_sublayer2_ptr%nodes(3,3,1) = -1.0d0
               
               loc_coords = [align_W+5-match_table(1),align_N-match_table(2)]
               
               test_node_activated = .true.

            !check whether the node is not activated by another buffer layer
            case(4)
               
               loc_coords = [align_W+7-match_table(1),align_N-match_table(2)]
               
               test_node_activated = .false.
               
          end select

        end subroutine get_param_test_bf_node_activated


        function test_is_interior_node_activated(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_interface_icr)              :: bf_interface_used
          type(pmodel_eq)                     :: p_model
          real(rkind)   , dimension(nx)       :: interior_x_map
          real(rkind)   , dimension(ny)       :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
          integer(ikind), dimension(2)        :: loc_coords

          logical :: test_loc


          test_validated = .true.


          loc_coords = [2,2]

          interior_nodes(2,2,1) = 1.0d0

          test_loc = .not.(bf_interface_used%is_interior_node_activated(
     $         loc_coords,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test desactivated failed'')'
          end if

          interior_nodes(2,2,1) = -1.0d0

          test_loc = (bf_interface_used%is_interior_node_activated(
     $         loc_coords,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test activated failed'')'
          end if

        end function test_is_interior_node_activated


        subroutine check_inputs()

          implicit none

          type(pmodel_eq)                :: p_model
          real(rkind), dimension(3)      :: x_map
          real(rkind), dimension(3)      :: y_map
          real(rkind), dimension(3,3,ne) :: nodes
          
          logical :: test_loc


          nodes(2,2,1) = -1.0d0
          test_loc = p_model%are_openbc_undermined(x_map,y_map,nodes)
          if(.not.test_loc) then
             print '(''the test requires: '')'
             print '(''openbc_undermined if nodes(2,2,1)<0'')'
             stop ''
          end if

          
          nodes(2,2,1) = 1.0d0
          test_loc = .not.p_model%are_openbc_undermined(x_map,y_map,nodes)
          if(.not.test_loc) then
             print '(''the test requires: '')'
             print '(''.not.openbc_undermined if nodes(2,2,1)>0'')'
             stop ''
          end if

        end subroutine check_inputs


      end program test_bf_interface_icr
