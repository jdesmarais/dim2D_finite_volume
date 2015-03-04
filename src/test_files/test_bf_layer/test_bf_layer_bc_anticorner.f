      program test_bf_layer_bc_anticorner

        use bf_layer_bc_anticorner_module, only :
     $       get_coords_from_pattern,
     $       are_grdpts_needed_for_flux_x,
     $       are_grdpts_needed_for_flux_y,
     $       extract_grdpts_to_compute_anticorner_fluxes

        use check_data_module, only :
     $       is_int_vector_validated,
     $       is_int_matrix_validated ,
     $       is_real_matrix3D_validated

        use parameters_bf_layer, only :
     $       align_E

        use parameters_constant, only :
     $       sd_interior_type,
     $       sd_R1_type

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.

        
        test_loc = test_get_coords_from_pattern(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_coords_from_pattern: '',L1)', test_loc
        print '()'

        test_loc = test_are_grdpts_needed_for_flux_x(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_are_grdpts_needed_for_flux_x: '',L1)', test_loc
        print '()'

        test_loc = test_are_grdpts_needed_for_flux_y(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_are_grdpts_needed_for_flux_y: '',L1)', test_loc
        print '()'

        test_loc = test_extract_grdpts_to_compute_anticorner_fluxes(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_extract_grdpts_to_compute_anticorner_fluxes: '',L1)', test_loc
        print '()'

        contains


        function test_get_coords_from_pattern(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer                   :: i,j
          integer, dimension(2,2,2) :: pattern
          integer                   :: size_x,size_y

          integer, dimension(2,2,2) :: border_coords_test
          integer, dimension(2,2)   :: cpt_coords_test
          logical, dimension(2)     :: grdpts_needed_test

          integer, dimension(2,2)   :: border_coords
          integer, dimension(2)     :: cpt_coords
          logical                   :: grdpts_needed

          integer :: k
          logical :: test_loc


          test_validated = .true.


          !input
          i=3
          j=4

          size_x = 10
          size_y = 6

          pattern = reshape((/
     $         -1,-2,2,1,
     $         -3,-3,3,2/),
     $         (/2,2,2/))

          border_coords_test = reshape((/
     $          2,2,5,5,
     $          0,1,6,6/),
     $         (/2,2,2/))

          cpt_coords_test = reshape((/
     $         2,3,
     $         4,4/),
     $         (/2,2/))

          grdpts_needed_test = [.false.,.true.]


          do k=1,size(pattern,3)

             !output
             call get_coords_from_pattern(
     $            pattern(:,:,k),
     $            i,j,
     $            size_x, size_y,
     $            border_coords,
     $            cpt_coords,
     $            grdpts_needed)

             !validation+detailled
             test_loc = is_int_matrix_validated(border_coords,border_coords_test(:,:,k),detailled)
             if(detailled.and.(.not.test_loc)) then
                print '(''borders_coords failed'')'
             end if
             test_validated = test_validated.and.test_loc

             test_loc = is_int_vector_validated(cpt_coords,cpt_coords_test(:,k),detailled)
             if(detailled.and.(.not.test_loc)) then
                print '(''cpt_coords failed'')'
             end if
             test_validated = test_validated.and.test_loc

             test_loc = grdpts_needed.eqv.grdpts_needed_test(k)
             if(detailled.and.(.not.test_loc)) then
                print '(''grdpts_needed failed'')'
             end if
             test_validated = test_validated.and.test_loc

          end do

        end function test_get_coords_from_pattern


        function test_are_grdpts_needed_for_flux_x(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(pmodel_eq)           :: p_model
          integer, dimension(2)     :: operator_type_test
          integer                   :: i
          integer                   :: j
          integer                   :: size_x
          integer                   :: size_y

          integer, dimension(2,2,2) :: border_coords_test
          integer, dimension(2,2)   :: cpt_coords_test
          logical, dimension(2)     :: grdpts_needed_test

          integer, dimension(2,2)   :: border_coords
          integer, dimension(2)     :: cpt_coords
          logical                   :: grdpts_needed
          
          integer                   :: k
          logical                   :: test_loc


          test_validated = .true.


          !input
          i=3
          j=2
          operator_type_test = [sd_interior_type,sd_R1_type]
          size_x = 5
          size_y = 3
          border_coords_test = reshape((/
     $         1,0,5,4,
     $         1,1,5,3/),
     $         (/2,2,2/))
          cpt_coords_test = reshape((/
     $         3,3,
     $         3,2/),
     $         (/2,2/))
          grdpts_needed_test = [.true.,.false.]
          

          do k=1, size(border_coords_test,3)

             !output
             grdpts_needed = are_grdpts_needed_for_flux_x(
     $            p_model,
     $            operator_type_test(k),
     $            i,j,
     $            size_x,size_y,
     $            border_coords,
     $            cpt_coords)

             !validation+detailled
             test_loc = is_int_matrix_validated(border_coords,border_coords_test(:,:,k),detailled)
             if(detailled.and.(.not.test_loc)) then
                print '(''borders_coords failed'')'
             end if
             test_validated = test_validated.and.test_loc

             test_loc = is_int_vector_validated(cpt_coords,cpt_coords_test(:,k),detailled)
             if(detailled.and.(.not.test_loc)) then
                print '(''cpt_coords failed'')'
             end if
             test_validated = test_validated.and.test_loc

             test_loc = grdpts_needed.eqv.grdpts_needed_test(k)
             if(detailled.and.(.not.test_loc)) then
                print '(''grdpts_needed failed'')'
             end if
             test_validated = test_validated.and.test_loc

          end do


        end function test_are_grdpts_needed_for_flux_x


        function test_are_grdpts_needed_for_flux_y(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(pmodel_eq)           :: p_model
          integer, dimension(2)     :: operator_type_test
          integer                   :: i
          integer                   :: j
          integer                   :: size_x
          integer                   :: size_y

          integer, dimension(2,2,2) :: border_coords_test
          integer, dimension(2,2)   :: cpt_coords_test
          logical, dimension(2)     :: grdpts_needed_test

          integer, dimension(2,2)   :: border_coords
          integer, dimension(2)     :: cpt_coords
          logical                   :: grdpts_needed
          
          integer                   :: k
          logical                   :: test_loc


          test_validated = .true.


          !input
          i=2
          j=3
          operator_type_test = [sd_interior_type,sd_R1_type]
          size_x = 3
          size_y = 5
          border_coords_test = reshape((/
     $         0,1,4,5,
     $         1,1,3,5/),
     $         (/2,2,2/))
          cpt_coords_test = reshape((/
     $         3,3,
     $         2,3/),
     $         (/2,2/))
          grdpts_needed_test = [.true.,.false.]
          

          do k=1, size(border_coords_test,3)

             !output
             grdpts_needed = are_grdpts_needed_for_flux_y(
     $            p_model,
     $            operator_type_test(k),
     $            i,j,
     $            size_x,size_y,
     $            border_coords,
     $            cpt_coords)

             !validation+detailled
             test_loc = is_int_matrix_validated(border_coords,border_coords_test(:,:,k),detailled)
             if(detailled.and.(.not.test_loc)) then
                print '(''borders_coords failed'')'
             end if
             test_validated = test_validated.and.test_loc

             test_loc = is_int_vector_validated(cpt_coords,cpt_coords_test(:,k),detailled)
             if(detailled.and.(.not.test_loc)) then
                print '(''cpt_coords failed'')'
             end if
             test_validated = test_validated.and.test_loc

             test_loc = grdpts_needed.eqv.grdpts_needed_test(k)
             if(detailled.and.(.not.test_loc)) then
                print '(''grdpts_needed failed'')'
             end if
             test_validated = test_validated.and.test_loc

          end do


        end function test_are_grdpts_needed_for_flux_y


        function test_extract_grdpts_to_compute_anticorner_fluxes(
     $     detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer(ikind), dimension(2,2)      :: bf_alignment
          integer       , dimension(6,6)      :: bf_grdpts_id
          real(rkind)   , dimension(6,6,ne)   :: bf_nodes
          real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
          integer(ikind), dimension(2,2)      :: gen_coords
          real(rkind)   , dimension(4,5,ne)   :: tmp_nodes_test
          real(rkind)   , dimension(4,5,ne)   :: tmp_nodes

          integer(ikind) :: i,j
          integer        :: k


          test_validated = .true.


          !input
          bf_alignment = reshape((/
     $         align_E,5,align_E+1,6/),
     $         (/2,2/))
          
          bf_grdpts_id = reshape((/
     $         1,1,2,3,3,3,
     $         1,1,2,2,2,3,
     $         1,1,1,1,2,3,
     $         1,1,1,1,2,3,
     $         1,1,2,2,2,3,
     $         1,1,2,3,3,3/),
     $         (/6,6/))

          do k=1, ne
             do j=1, size(bf_nodes,2)
                do i=1, size(bf_nodes,1)
                   bf_nodes(i,j,k) = (align_E-3+i-1) + 10*(2+j-1) + 100*k
                end do
             end do
          end do

          do k=1, ne
             do j=1, size(interior_nodes,2)
                do i=1, size(interior_nodes,1)
                   interior_nodes(i,j,k) = (i-1) + 10*(j-1) + 100*k
                end do
             end do
          end do

          do k=1, ne
             do j=1, size(tmp_nodes,2)
                do i=1, size(tmp_nodes,1)
                   tmp_nodes_test(i,j,k) = (align_E-3+i-1) + 10*(5+j-1) + 100*k
                end do
             end do
          end do

          gen_coords = reshape((/
     $         align_E-2, 6, align_E+1, 10/),
     $         (/2,2/))


          !output
          call extract_grdpts_to_compute_anticorner_fluxes(
     $         bf_alignment,
     $         bf_grdpts_id,
     $         bf_nodes,
     $         interior_nodes,
     $         gen_coords,
     $         tmp_nodes)


          !validation+detailled
          test_loc = is_real_matrix3D_validated(
     $         tmp_nodes, tmp_nodes_test, detailled)
          test_validated = test_validated.and.test_loc


        end function test_extract_grdpts_to_compute_anticorner_fluxes

      end program test_bf_layer_bc_anticorner
