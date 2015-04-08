      program test_periodic_xy

        use bc_operators_class, only :
     $       bc_operators

        use check_data_module, only :
     $       is_int_validated,
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
     $       periodic_xy_choice,
     $       bc_timedev_choice,
     $       bc_nodes_choice

        use parameters_input, only :
     $       nx,ny,ne,
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


        test_loc = test_ini(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_ini: '',L1)', test_loc
        print '()'


        test_loc = test_apply_bc_on_nodes(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_apply_bc_on_nodes: '',L1)', test_loc
        print '()'


        test_loc = test_apply_bc_on_nodes_nopt(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_apply_bc_on_nodes_nopt: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_ini(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bc_operators) :: bc_operators_used
          type(pmodel_eq)    :: p_model


          test_validated = .true.


          ! output
          call bc_operators_used%ini(p_model)


          ! validation
          test_loc = is_int_validated(
     $         bc_operators_used%period_x,
     $         4,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test period_x failed'')'
          end if

          test_loc = is_int_validated(
     $         bc_operators_used%period_y,
     $         6,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test period_y failed'')'
          end if

          test_loc = is_int_vector_validated(
     $         bc_operators_used%bc_type,
     $         [bc_nodes_choice,bc_nodes_choice,
     $          bc_nodes_choice,bc_nodes_choice],
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test bc_type failed'')'
          end if

        end function test_ini


        function test_apply_bc_on_nodes(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bc_operators)               :: bc_operators_used
          type(pmodel_eq)                  :: p_model
          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(nx,ny,ne) :: test_nodes

          integer(ikind) :: i,j
          integer        :: k

          
          test_validated = .true.


          !input
          call bc_operators_used%ini(p_model)

          nodes = reshape((/
     $         ((((i-1)+8*(j-1)+100*(k-1), i=1,8), j=1,10),k=1,ne)/),
     $         (/nx,ny,ne/))

          test_nodes = nodes

          !W_edge
          test_nodes(1:2,3:8,:) = reshape((/
     $           20,  21,  28,  29,  36,  37,  44,  45,  52,  53,  60,  61,
     $          120, 121, 128, 129, 136, 137, 144, 145, 152, 153, 160, 161,
     $          220, 221, 228, 229, 236, 237, 244, 245, 252, 253, 260, 261
     $         /),(/2,6,ne/))

          !E edge
          test_nodes(7:8,3:8,:) = reshape((/
     $           18,  19,  26,  27,  34,  35,  42,  43,  50,  51,  58,  59,
     $          118, 119, 126, 127, 134, 135, 142, 143, 150, 151, 158, 159,
     $          218, 219, 226, 227, 234, 235, 242, 243, 250, 251, 258, 259
     $         /),(/2,6,ne/))

          !S edge
          test_nodes(1:8,1:2,:) = reshape((/
     $         52,53,50,51,52,53,50,51,
     $         60,61,58,59,60,61,58,59,
     $         152,153,150,151,152,153,150,151,
     $         160,161,158,159,160,161,158,159,
     $         252,253,250,251,252,253,250,251,
     $         260,261,258,259,260,261,258,259
     $         /),(/8,2,ne/))

          !N edge
          test_nodes(1:8,9:10,:) = reshape((/
     $         20,21,18,19,20,21,18,19,
     $         28,29,26,27,28,29,26,27,
     $         120,121,118,119,120,121,118,119,
     $         128,129,126,127,128,129,126,127,
     $         220,221,218,219,220,221,218,219,
     $         228,229,226,227,228,229,226,227
     $         /),(/8,2,ne/))

          !output
          call bc_operators_used%apply_bc_on_nodes(nodes)

          !validation
          test_loc = is_real_matrix3D_validated(
     $         nodes,
     $         test_nodes,
     $         detailled)
          test_validated = test_validated.and.test_loc

        end function test_apply_bc_on_nodes


        function test_apply_bc_on_nodes_nopt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bc_operators)               :: bc_operators_used
          type(pmodel_eq)                  :: p_model
          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(nx,ny,ne) :: test_nodes

          integer(ikind) :: i,j
          integer        :: k

          integer(ikind), dimension(:,:), allocatable :: bc_sections

          
          !input
          call bc_operators_used%ini(p_model)

          nodes = reshape((/
     $         ((((i-1)+8*(j-1)+100*(k-1), i=1,8), j=1,10),k=1,ne)/),
     $         (/nx,ny,ne/))

          test_nodes = nodes

          call bc_operators_used%apply_bc_on_nodes(test_nodes)
          
          allocate(bc_sections(5,8))
          bc_sections(:,1) = [W_edge_type   ,1,3,8         ,no_overlap]
          bc_sections(:,2) = [E_edge_type   ,7,3,8         ,no_overlap]
          bc_sections(:,3) = [SW_corner_type,1,1,no_overlap,no_overlap]
          bc_sections(:,4) = [S_edge_type   ,3,1,6         ,no_overlap]
          bc_sections(:,5) = [SE_corner_type,7,1,no_overlap,no_overlap]
          bc_sections(:,6) = [NW_corner_type,1,9,no_overlap,no_overlap]
          bc_sections(:,7) = [N_edge_type   ,3,9,6         ,no_overlap]
          bc_sections(:,8) = [NE_corner_type,7,9,no_overlap,no_overlap]


          !output
          call bc_operators_used%apply_bc_on_nodes_nopt(nodes,bc_sections)


          !validation
          test_validated = is_real_matrix3D_validated(
     $         nodes,
     $         test_nodes,
     $         detailled)
          

        end function test_apply_bc_on_nodes_nopt


        subroutine check_inputs()

          implicit none


          if(.not.(
     $         (nx.eq.8).and.
     $         (ny.eq.10).and.
     $         (ne.eq.3).and.
     $         (bc_choice.eq.periodic_xy_choice).and.
     $         (bc_N_type_choice.eq.bc_nodes_choice).and.
     $         (bc_S_type_choice.eq.bc_nodes_choice).and.
     $         (bc_E_type_choice.eq.bc_nodes_choice).and.
     $         (bc_W_type_choice.eq.bc_nodes_choice))) then

             print '(''the test requires:'')'
             print '(''     - nx=8'')'
             print '(''     - ny=10'')'
             print '(''     - ne=3'')'
             print '(''     - bc_choice=periodic_xy_choice'')'
             print '(''     - bc_N_type=bc_nodes_choice'')'
             print '(''     - bc_S_type=bc_nodes_choice'')'
             print '(''     - bc_E_type=bc_nodes_choice'')'
             print '(''     - bc_W_type=bc_nodes_choice'')'
             stop ''

          end if

        end subroutine check_inputs

      end program test_periodic_xy
