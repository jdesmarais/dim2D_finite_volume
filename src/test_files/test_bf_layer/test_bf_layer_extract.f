      program test_bf_layer_extract

        use bf_layer_extract_module, only :
     $     get_indices_to_extract_interior_data,
     $     get_indices_to_extract_bf_layer_data,
     $     get_bf_layer_match_table,
     $     get_grdpts_id_from_interior

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_input, only :
     $       nx,ny

        use parameters_kind, only :
     $       ikind


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled      = .true.
        test_validated = .true.
        

        test_loc = test_get_indices_to_extract_interior_data(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_indices_to_extract_interior_data: '',L1)', test_loc
        print '()'

        test_loc = test_get_indices_to_extract_bf_layer_data(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_indices_to_extract_bf_layer_data: '',L1)', test_loc
        print '()'

        test_loc = test_get_bf_layer_match_table(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_bf_layer_match_table: '',L1)', test_loc
        print '()'

        test_loc = test_get_grdpts_id_from_interior(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_grdpts_id_from_interior: '',L1)', test_loc
        print '()'


        contains


        function test_get_indices_to_extract_interior_data(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer                        :: k
          integer(ikind), dimension(2,2) :: gen_coords
          integer(ikind), dimension(6)   :: output_test
          integer(ikind), dimension(6)   :: output
          character(6)  , dimension(6)   :: char_test

          !input data
          gen_coords(1,1) = 3
          gen_coords(1,2) = nx-4
          gen_coords(2,1) = -2
          gen_coords(2,2) = ny+5
          
          output_test(1) = nx-6
          output_test(2) = ny
          output_test(3) = 1
          output_test(4) = 4
          output_test(5) = 3
          output_test(6) = 1          

          char_test = [
     $         'size_x',
     $         'size_y',
     $         'i_recv',
     $         'j_recv',
     $         'i_send',
     $         'j_send']
          
          !output data
          call get_indices_to_extract_interior_data(
     $         gen_coords,
     $         output(1), output(2),
     $         output(3), output(4),
     $         output(5), output(6))
          
          !validation
          test_validated = .true.
          do k=1,6
             test_validated = test_validated.and.(
     $            output(k).eq.output_test(k))
          end do

          if((.not.test_validated).and.detailled) then
             do k=1,6
                print '(I2,'' -> '',I2)', output(k), output_test(k)
             end do
          end if

        end function test_get_indices_to_extract_interior_data


        function test_get_indices_to_extract_bf_layer_data(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer                        :: k
          integer(ikind), dimension(2,2) :: bf_align
          integer(ikind), dimension(2,2) :: gen_coords
          integer(ikind), dimension(6)   :: output_test
          integer(ikind), dimension(6)   :: output
          character(6)  , dimension(6)   :: char_test

          !input data
          bf_align(1,1) = 3
          bf_align(1,2) = nx-2
          bf_align(2,1) = 3
          bf_align(2,2) = ny-2

          gen_coords(1,1) = 3
          gen_coords(1,2) = nx-4
          gen_coords(2,1) = -2
          gen_coords(2,2) = ny+5
          
          output_test(1) = nx-6
          output_test(2) = ny
          output_test(3) = 1
          output_test(4) = 4
          output_test(5) = 3
          output_test(6) = 1          

          char_test = [
     $         'size_x',
     $         'size_y',
     $         'i_recv',
     $         'j_recv',
     $         'i_send',
     $         'j_send']
          
          !output data
          call get_indices_to_extract_bf_layer_data(
     $         bf_align,
     $         gen_coords,
     $         output(1), output(2),
     $         output(3), output(4),
     $         output(5), output(6))
          
          !validation
          test_validated = .true.
          do k=1,6
             test_validated = test_validated.and.(
     $            output(k).eq.output_test(k))
          end do

          if((.not.test_validated).and.detailled) then
             do k=1,6
                print '(I2,'' -> '',I2)', output(k), output_test(k)
             end do
          end if

        end function test_get_indices_to_extract_bf_layer_data


        function test_get_bf_layer_match_table(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer(ikind), dimension(2,2) :: alignment
          integer(ikind), dimension(2)   :: match_table_test
          integer(ikind), dimension(2)   :: match_table

          !input data
          alignment(1,1) = 3
          alignment(1,2) = 5
          alignment(2,1) = 4
          alignment(2,2) = 6
          
          match_table_test = [0,1]

          !output data
          match_table = get_bf_layer_match_table(alignment)

          !validation
          test_validated = 
     $         (match_table(1).eq.match_table_test(1)).and.
     $         (match_table(2).eq.match_table_test(2))

          if((.not.test_validated).and.detailled) then
             
             print '(I2,'' -> '',I2)', match_table(1), match_table_test(1)
             print '(I2,'' -> '',I2)', match_table(2), match_table_test(2)
             
          end if

        end function test_get_bf_layer_match_table


        function test_get_grdpts_id_from_interior(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(5,5,3) :: grdpts_id_test
          integer, dimension(2,2,3) :: gen_coords_test
          
          integer, dimension(5,5)   :: grdpts_id
          integer                   :: k
          logical                   :: test_loc
          
          test_validated = .true.

          !input
          grdpts_id_test = reshape((/
     $         0,0,0,3,2,
     $         0,0,0,3,3,
     $         0,0,0,0,0,
     $         0,0,0,0,0,
     $         0,0,0,0,0,
     $         
     $         1,1,2,3,0,
     $         1,1,2,3,0,
     $         2,2,2,3,0,
     $         3,3,3,3,0,
     $         0,0,0,0,0,
     $         
     $         0,0,0,0,0,
     $         3,3,3,3,3,
     $         2,2,2,2,2,
     $         1,1,1,1,1,
     $         1,1,1,1,1/),
     $         (/5,5,3/))

          gen_coords_test = reshape((/
     $         -2,ny-1,2,ny+3,
     $         nx-3,ny-3,nx+1,ny+1,
     $         3,0,7,4/),
     $         (/2,2,3/))
          

          do k=1, size(gen_coords_test,3)

             !output
             call get_grdpts_id_from_interior(
     $            grdpts_id,
     $            gen_coords_test(:,:,k))

             !validation
             test_loc = is_int_matrix_validated(
     $            grdpts_id,grdpts_id_test(:,:,k),detailled)
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')'
             end if

          end do

        end function test_get_grdpts_id_from_interior

      end program test_bf_layer_extract
