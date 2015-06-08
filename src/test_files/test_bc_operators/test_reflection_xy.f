      program test_reflection_xy

        use bc_operators_reflection_xy_class, only :
     $       bc_operators_reflection_xy

        use check_data_module, only :
     $       is_int_vector_validated,
     $       is_real_matrix3D_validated

        use parameters_bf_layer, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NW_corner_type,
     $       NE_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       no_overlap

        use parameters_constant, only :
     $       reflection_xy_choice,
     $       bc_timedev_choice,
     $       bc_nodes_choice

        use parameters_input, only :
     $       nx,ny,ne,bc_size,
     $       bc_N_type_choice,bc_S_type_choice,
     $       bc_E_type_choice,bc_W_type_choice,
     $       bc_choice

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

        test_loc = test_apply_bc_on_nodes(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_apply_bc_on_nodes: '',L1)', test_loc
        print '()'


c$$$        test_loc = test_apply_bc_on_nodes_nopt(detailled)
c$$$        test_validated = test_validated.and.test_loc
c$$$        print '(''test_apply_bc_on_nodes_nopt: '',L1)', test_loc
c$$$        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_apply_bc_on_nodes(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bc_operators_reflection_xy) :: bc_operators_used
          type(pmodel_eq)                  :: p_model
          real(rkind)                      :: t
          real(rkind), dimension(nx)       :: x_map
          real(rkind), dimension(ny)       :: y_map
          real(rkind), dimension(nx,ny,ne) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(nx,ny,ne) :: test_nodes

          integer(ikind) :: i,j
          integer        :: k

          
          test_validated = .true.


          !input

          nodes = reshape((/
     $         ((((i-1)+6*(j-1)+100*(k-1), i=1,6), j=1,8),k=1,ne)/),
     $         (/nx,ny,ne/))

          test_nodes = nodes

          !W_edge
          test_nodes(1:2,3:6,:) = reshape((/
     $         15 ,  14,   21, 20,   27,  26,   33,  32,
     $        -115,-114, -121,-120,-127,-126, -133,-132,
     $         215, 214,  221, 220, 227, 226,  233, 232
     $         /),(/2,4,ne/))

          !E edge
          test_nodes(5:6,3:6,:) = reshape((/
     $         15 ,  14,   21, 20,   27,  26,   33,  32,
     $        -115,-114, -121,-120,-127,-126, -133,-132,
     $         215, 214,  221, 220, 227, 226,  233, 232
     $         /),(/2,4,ne/))

          !S edge
          test_nodes(1:6,1:2,:) = reshape((/
     $          21,  20,   20,  21,  21,  20,
     $          15,  14,   14,  15,  15,  14,
     $        -121,-120,  120, 121,-121,-120,
     $        -115,-114,  114, 115,-115,-114,
     $        -221,-220, -220,-221,-221,-220,
     $        -215,-214, -214,-215,-215,-214
     $         /),(/6,2,ne/))

          !N edge
          test_nodes(1:6,7:8,:) = reshape((/
     $          33,  32,   32,  33,  33,  32,
     $          27,  26,   26,  27,  27,  26,
     $        -133,-132,  132, 133,-133,-132,
     $        -127,-126,  126, 127,-127,-126,
     $        -233,-232, -232,-233,-233,-232,
     $        -227,-226, -226,-227,-227,-226
     $         /),(/6,2,ne/))

          !output
          call bc_operators_used%apply_bc_on_nodes(
     $         [W_edge_type,1,bc_size+1,ny-bc_size+1],
     $         t,x_map,y_map,nodes_tmp,
     $         p_model,
     $         nodes)

          call bc_operators_used%apply_bc_on_nodes(
     $         [E_edge_type,nx-bc_size+1,bc_size+1,ny-bc_size+1],
     $         t,x_map,y_map,nodes_tmp,
     $         p_model,
     $         nodes)

          call bc_operators_used%apply_bc_on_nodes(
     $         [S_edge_type,1,1,nx],
     $         t,x_map,y_map,nodes_tmp,
     $         p_model,
     $         nodes)

          call bc_operators_used%apply_bc_on_nodes(
     $         [N_edge_type,1,ny-bc_size+1,nx],
     $         t,x_map,y_map,nodes_tmp,
     $         p_model,
     $         nodes)

          !validation
          test_loc = is_real_matrix3D_validated(
     $         nodes,
     $         test_nodes,
     $         detailled)
          test_validated = test_validated.and.test_loc

        end function test_apply_bc_on_nodes

c$$$
c$$$        function test_apply_bc_on_nodes_nopt(detailled)
c$$$     $       result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          logical, intent(in) :: detailled
c$$$          logical             :: test_validated
c$$$
c$$$
c$$$          type(bc_operators)               :: bc_operators_used
c$$$          type(pmodel_eq)                  :: p_model
c$$$          real(rkind), dimension(nx,ny,ne) :: nodes
c$$$          real(rkind), dimension(nx,ny,ne) :: test_nodes
c$$$
c$$$          integer(ikind) :: i,j
c$$$          integer        :: k
c$$$
c$$$          integer(ikind), dimension(:,:), allocatable :: bc_sections
c$$$
c$$$          
c$$$          !input
c$$$          call bc_operators_used%ini(p_model)
c$$$
c$$$          nodes = reshape((/
c$$$     $         ((((i-1)+6*(j-1)+100*(k-1), i=1,6), j=1,8),k=1,ne)/),
c$$$     $         (/nx,ny,ne/))
c$$$
c$$$          test_nodes = nodes
c$$$
c$$$          call bc_operators_used%apply_bc_on_nodes(test_nodes)
c$$$          
c$$$          allocate(bc_sections(5,8))
c$$$          bc_sections(:,1) = [W_edge_type   ,1,3,6         ,no_overlap]
c$$$          bc_sections(:,2) = [E_edge_type   ,5,3,6         ,no_overlap]
c$$$          bc_sections(:,3) = [SW_corner_type,1,1,no_overlap,no_overlap]
c$$$          bc_sections(:,4) = [S_edge_type   ,3,1,4         ,no_overlap]
c$$$          bc_sections(:,5) = [SE_corner_type,5,1,no_overlap,no_overlap]
c$$$          bc_sections(:,6) = [NW_corner_type,1,7,no_overlap,no_overlap]
c$$$          bc_sections(:,7) = [N_edge_type   ,3,7,4         ,no_overlap]
c$$$          bc_sections(:,8) = [NE_corner_type,5,7,no_overlap,no_overlap]
c$$$
c$$$
c$$$          !output
c$$$          call bc_operators_used%apply_bc_on_nodes_nopt(nodes,bc_sections)
c$$$
c$$$
c$$$          !validation
c$$$          test_validated = is_real_matrix3D_validated(
c$$$     $         nodes,
c$$$     $         test_nodes,
c$$$     $         detailled)
c$$$          
c$$$
c$$$        end function test_apply_bc_on_nodes_nopt


        subroutine check_inputs()

          implicit none


          if(.not.(
     $         (nx.eq.6).and.
     $         (ny.eq.8).and.
     $         (ne.eq.3).and.
     $         (bc_choice.eq.reflection_xy_choice).and.
     $         (bc_N_type_choice.eq.bc_nodes_choice).and.
     $         (bc_S_type_choice.eq.bc_nodes_choice).and.
     $         (bc_E_type_choice.eq.bc_nodes_choice).and.
     $         (bc_W_type_choice.eq.bc_nodes_choice))) then

             print '(''the test requires:'')'
             print '(''     - nx=6'')'
             print '(''     - ny=8'')'
             print '(''     - ne=3'')'
             print '(''     - bc_choice=reflection_xy_choice'')'
             print '(''     - bc_N_type=bc_nodes_choice'')'
             print '(''     - bc_S_type=bc_nodes_choice'')'
             print '(''     - bc_E_type=bc_nodes_choice'')'
             print '(''     - bc_W_type=bc_nodes_choice'')'
             stop ''

          end if

        end subroutine check_inputs

      end program test_reflection_xy
