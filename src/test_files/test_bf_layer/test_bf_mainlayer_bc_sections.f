       program test_bf_mainlayer_bc_sections

        use bf_mainlayer_bc_sections_module, only :
     $       update_interior_bc_sections_from_mainlayers

        use bf_mainlayer_pointer_class, only :
     $       bf_mainlayer_pointer

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       SW_corner_type,
     $       SE_corner_type,
     $       NW_corner_type,
     $       NE_corner_type,
     $       
     $       align_N,align_S,
     $       align_E,align_W,
     $       no_overlap

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only : 
     $       ikind,
     $       rkind


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .false.
        test_validated = .true.


        call check_inputs()


        test_loc = test_update_interior_bc_sections_from_mainlayers(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_update_interior_bc_sections_from_mainlayers: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_update_interior_bc_sections_from_mainlayers(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,2

             test_loc = perform_test_update_interior_bc_sections(k,detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test ('',I2,'') failed'')', k
             end if

          end do

        end function test_update_interior_bc_sections_from_mainlayers


        function perform_test_update_interior_bc_sections(test_id,detailled)
     $     result(test_validated)

          implicit none

          integer, intent(in) :: test_id
          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_mainlayer_pointer), dimension(4)    :: mainlayer_pointers
          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections
          integer(ikind), dimension(:,:), allocatable :: interior_bc_sections_test


          !input
          call get_param_test_update_interior_bc_sections(
     $         test_id,
     $         mainlayer_pointers,
     $         interior_bc_sections_test)

          
          !output
          call update_interior_bc_sections_from_mainlayers(
     $         mainlayer_pointers,
     $         interior_bc_sections)


          !validation
          test_validated = is_int_matrix_validated(
     $         interior_bc_sections,
     $         interior_bc_sections_test,
     $         detailled)

        end function perform_test_update_interior_bc_sections


        subroutine get_param_test_update_interior_bc_sections(
     $     test_id,
     $     mainlayer_pointers,
     $     interior_bc_sections_test)

          implicit none

          integer                                    , intent(in)    :: test_id
          type(bf_mainlayer_pointer), dimension(4)   , intent(inout) :: mainlayer_pointers
          integer(ikind), dimension(:,:), allocatable, intent(inout) :: interior_bc_sections_test

          type(bf_sublayer), pointer       :: added_sublayer
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          integer(ikind), dimension(2,2)   :: bf_alignment_tmp
          integer :: k


          !initialize the mainlayers
          do k=1,4
             call mainlayer_pointers(k)%ini()
             call mainlayer_pointers(k)%ini_mainlayer(k)
          end do

          
          select case(test_id)
            case(1)

               !initialization of the mainlayer_pointers
               !------------------------------------------------------------
               bf_alignment_tmp = reshape((/
     $              align_W+5,align_N,align_W+8,align_N/),
     $              (/2,2/))

               added_sublayer => mainlayer_pointers(N)%add_sublayer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               bf_alignment_tmp = reshape((/
     $              align_E,align_S+5,align_E+5,align_S+10/),
     $              (/2,2/))

               added_sublayer => mainlayer_pointers(E)%add_sublayer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               !interior_bc_sections
               !------------------------------------------------------------
               allocate(interior_bc_sections_test(5,10))

               interior_bc_sections_test(:,1) = [
     $              SW_corner_type, 1,1, no_overlap, no_overlap]

               interior_bc_sections_test(:,2) = [
     $              S_edge_type, align_W+1, 1, align_E-1, no_overlap]

               interior_bc_sections_test(:,3) = [
     $              SE_corner_type, align_E, 1, no_overlap, no_overlap]

               interior_bc_sections_test(:,4) = [
     $              W_edge_type, 1, align_S+1, align_N-1, no_overlap]

               interior_bc_sections_test(:,5) = [
     $              E_edge_type, align_E, align_S+1, align_S+2, no_overlap]

               interior_bc_sections_test(:,6) = [
     $              E_edge_type, align_E, align_S+13, align_N-1, no_overlap]

               interior_bc_sections_test(:,7) = [
     $              NW_corner_type, 1, align_N, no_overlap, no_overlap]

               interior_bc_sections_test(:,8) = [
     $              N_edge_type, bc_size+1, align_N, align_W+2, no_overlap]
               
               interior_bc_sections_test(:,9) = [
     $              N_edge_type, align_W+11, align_N, align_E-1, no_overlap]

               interior_bc_sections_test(:,10) = [
     $              NE_corner_type, align_E, align_N, no_overlap, no_overlap]

            case(2)

               !initialization of the mainlayer_pointers
               !------------------------------------------------------------
               !North buffer layers
               bf_alignment_tmp = reshape((/
     $              align_W-6,align_N,align_W+4,align_N+1/),
     $              (/2,2/))

               added_sublayer => mainlayer_pointers(N)%add_sublayer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               bf_alignment_tmp = reshape((/
     $              align_W+11,align_N,align_W+15,align_N+1/),
     $              (/2,2/))

               added_sublayer => mainlayer_pointers(N)%add_sublayer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               !South buffer layers
               bf_alignment_tmp = reshape((/
     $              align_E-4,align_S-1,align_E+6,align_S/),
     $              (/2,2/))

               added_sublayer => mainlayer_pointers(S)%add_sublayer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               bf_alignment_tmp = reshape((/
     $              align_E-15,align_S-1,align_E-11,align_S/),
     $              (/2,2/))

               added_sublayer => mainlayer_pointers(S)%add_sublayer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               !West buffer layers
               bf_alignment_tmp = reshape((/
     $              align_W-5,align_S+5,align_W,align_S+6/),
     $              (/2,2/))

               added_sublayer => mainlayer_pointers(W)%add_sublayer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               bf_alignment_tmp = reshape((/
     $              align_W-6,align_N-2,align_W,align_N-1/),
     $              (/2,2/))

               added_sublayer => mainlayer_pointers(W)%add_sublayer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               !East buffer layers
               bf_alignment_tmp = reshape((/
     $              align_E,align_N-6,align_E+5,align_N-5/),
     $              (/2,2/))

               added_sublayer => mainlayer_pointers(E)%add_sublayer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               bf_alignment_tmp = reshape((/
     $              align_E,align_S+1,align_E+6,align_S+2/),
     $              (/2,2/))

               added_sublayer => mainlayer_pointers(E)%add_sublayer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               allocate(interior_bc_sections_test(5,10))

               interior_bc_sections_test(:,1) = [
     $              SW_corner_type, 1,1, no_overlap, no_overlap]
               
               interior_bc_sections_test(:,2) = [
     $              S_edge_type, align_W+1, 1, align_E-18, no_overlap]
               
               interior_bc_sections_test(:,3) = [
     $              S_edge_type, align_E-8, 1, align_E-7, no_overlap]
               
               interior_bc_sections_test(:,4) = [
     $              W_edge_type, 1, align_S+1, align_S+2, no_overlap]

               interior_bc_sections_test(:,5) = [
     $              W_edge_type, 1, align_S+9, align_N-5, no_overlap]

               interior_bc_sections_test(:,6) = [
     $              E_edge_type, align_E, align_S+5, align_N-9, no_overlap]

               interior_bc_sections_test(:,7) = [
     $              E_edge_type, align_E, align_N-2, align_N-1, no_overlap]

               interior_bc_sections_test(:,8) = [
     $              N_edge_type, align_W+7, align_N, align_W+8, no_overlap]

               interior_bc_sections_test(:,9) = [
     $              N_edge_type, align_W+18, align_N, align_E-1, no_overlap]

               interior_bc_sections_test(:,10) = [
     $              NE_corner_type, align_E, align_N, no_overlap, no_overlap]

          end select               

        end subroutine get_param_test_update_interior_bc_sections


        subroutine check_inputs()

          implicit none

          if(.not.(
     $         (nx.eq.30).and.
     $         (ny.eq.35))) then

             print '(''the test is designed for'')'
             print '(''  - nx = 30'')'
             print '(''  - ny = 35'')'
             stop ''

          end if

        end subroutine check_inputs

      end program test_bf_mainlayer_bc_sections
