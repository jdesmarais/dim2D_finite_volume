      program test_bf_layer_extract

        use bf_layer_extract_module, only :
     $     get_indices_to_extract_interior_data,
     $     get_indices_to_extract_bf_layer_data,
     $     get_bf_layer_match_table,
     $     get_grdpts_id_from_interior,
     $     get_grdpts_id_from_bf_layer,
     $     get_nodes_from_interior,
     $     get_nodes_from_bf_layer,
     $     get_map_from_interior,
     $     set_grdpts_id_in_bf_layer

        use check_data_module, only :
     $       is_int_matrix_validated,
     $       is_real_vector_validated,
     $       is_real_matrix3D_validated

        use parameters_bf_layer, only :
     $       align_N,align_W

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind


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

        test_loc = test_get_grdpts_id_from_bf_layer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_grdpts_id_from_bf_layer: '',L1)', test_loc
        print '()'

        test_loc = test_get_nodes_from_interior(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_nodes_from_interior: '',L1)', test_loc
        print '()'

        test_loc = test_get_nodes_from_bf_layer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_nodes_from_bf_layer: '',L1)', test_loc
        print '()'

        test_loc = test_get_map_from_interior(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_map_from_interior: '',L1)', test_loc
        print '()'

        test_loc = test_set_grdpts_id_in_bf_layer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_set_grdpts_id_in_bf_layer: '',L1)', test_loc
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


        function test_get_grdpts_id_from_bf_layer(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer                         :: i,j
          integer       , dimension(4,7)  :: tmp_grdpts_id
          integer       , dimension(4,7)  :: tmp_grdpts_id_test
          integer       , dimension(2,2)  :: gen_coords
          integer(ikind), dimension(2,2)  :: bf_alignment
          integer       , dimension(12,9) :: bf_grdpts_id
          

          !input
          tmp_grdpts_id = reshape((/
     $         ((0, i=1,4),j=1,7)/),
     $         (/4,7/))

          tmp_grdpts_id_test = reshape((/
     $         1,1,1,1,
     $         1,2,2,2,
     $         1,2,3,3,
     $         2,2,3,0,
     $         3,3,3,0,
     $         0,0,0,0,
     $         0,0,0,0/),
     $         (/4,7/))

          gen_coords = reshape((/
     $         align_W+9,align_N+2,align_W+12,align_N+8/),
     $         (/2,2/))

          bf_alignment = reshape((/
     $         align_W+5,align_N,align_W+12,align_N+4/),
     $         (/2,2/))

          bf_grdpts_id = reshape((/
     $         1,1,1,1,1,1,1,1,1,1,1,1,
     $         1,1,1,1,1,1,1,1,1,1,1,1,
     $         2,2,1,1,1,1,1,1,1,1,2,2,
     $         3,2,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,2,2,2,2,3,
     $         3,2,1,1,1,1,1,2,3,3,3,3,
     $         3,2,2,2,2,2,2,2,3,0,0,0,
     $         3,3,3,3,3,3,3,3,3,0,0,0/),
     $         (/12,9/))

          !output
          call get_grdpts_id_from_bf_layer(
     $         tmp_grdpts_id,
     $         gen_coords,
     $         bf_alignment,
     $         bf_grdpts_id)

          !validation
          test_validated = is_int_matrix_validated(
     $         tmp_grdpts_id,
     $         tmp_grdpts_id_test,
     $         detailled)

        end function test_get_grdpts_id_from_bf_layer


        function test_get_nodes_from_interior(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind)   , dimension(7,6,ne)   :: tmp_nodes
          real(rkind)   , dimension(7,6,ne)   :: tmp_nodes_test
          integer(ikind), dimension(2,2)      :: gen_coords
          real(rkind)   , dimension(nx,ny,ne) :: interior_nodes

          integer(ikind), dimension(6)        :: extract_param
          
          logical        :: test_loc
          integer(ikind) :: i,j
          integer        :: k

          test_validated = .true.


          gen_coords = reshape(
     $         (/align_W,align_N-4,align_W+6,align_N+1/),(/2,2/))

          interior_nodes = reshape((/
     $         (((200*(k-1)+20*(j-1)+(i-1),i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          tmp_nodes_test = reshape((/
     $         (((200*(k-1)+20*(align_N-5+j-1)+(align_W-1+i-1),i=1,7),j=1,6),k=1,ne)/),
     $         (/7,6,ne/))


          !no optional arguments test
          !------------------------------------------------------------
          call get_nodes_from_interior(
     $         tmp_nodes,
     $         gen_coords,
     $         interior_nodes)

          test_loc = is_real_matrix3D_validated(
     $         tmp_nodes,
     $         tmp_nodes_test,
     $         detailled)
          test_loc = test_loc.and.test_validated
          if(detailled.and.(.not.test_loc)) then
             print '(''no optional arg failed'')'
          end if


          !extract_param_out optional argument
          !------------------------------------------------------------
          tmp_nodes = reshape((/
     $         (((-99.0d0,i=1,7),j=1,6),k=1,ne)/),
     $         (/7,6,ne/))

          call get_nodes_from_interior(
     $         tmp_nodes,
     $         gen_coords,
     $         interior_nodes,
     $         extract_param_out=extract_param)

          test_loc = is_real_matrix3D_validated(
     $         tmp_nodes,
     $         tmp_nodes_test,
     $         detailled)
          test_loc = test_loc.and.test_validated
          if(detailled.and.(.not.test_loc)) then
             print '(''extract_param_out optional arg failed'')'
          end if


          !extract_param_in optional argument
          !------------------------------------------------------------
          tmp_nodes = reshape((/
     $         (((-99.0d0,i=1,7),j=1,6),k=1,ne)/),
     $         (/7,6,ne/))

          call get_nodes_from_interior(
     $         tmp_nodes,
     $         gen_coords,
     $         interior_nodes,
     $         extract_param_in=extract_param)

          test_loc = is_real_matrix3D_validated(
     $         tmp_nodes,
     $         tmp_nodes_test,
     $         detailled)
          test_loc = test_loc.and.test_validated
          if(detailled.and.(.not.test_loc)) then
             print '(''extract_param_in optional arg failed'')'
          end if

        end function test_get_nodes_from_interior


        function test_get_nodes_from_bf_layer(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind)   , dimension(7,6,ne) :: tmp_nodes
          real(rkind)   , dimension(7,6,ne) :: tmp_nodes_test
          integer(ikind), dimension(2,2)    :: gen_coords
          integer(ikind), dimension(2,2)    :: bf_alignment
          real(rkind)   , dimension(9,8,ne) :: bf_nodes

          integer(ikind), dimension(6)      :: extract_param
          
          logical        :: test_loc
          integer(ikind) :: i,j
          integer        :: k

          test_validated = .true.


          gen_coords   = reshape((/align_W-4, align_N-1, align_W+2, align_N+4/),(/2,2/))

          bf_alignment = reshape((/align_W-3, align_N  , align_W+1, align_N+3/),(/2,2/))

          bf_nodes = reshape((/
     $         (((200*(k-1)+20*(align_N-3+j-1)+(align_W-6+i-1),i=1,9),j=1,8),k=1,ne)/),
     $         (/9,8,ne/))

          tmp_nodes_test = reshape((/
     $         (((200*(k-1)+20*(align_N-2+j-1)+(align_W-5+i-1),i=1,7),j=1,6),k=1,ne)/),
     $         (/7,6,ne/))


          !no optional arguments test
          !------------------------------------------------------------
          call get_nodes_from_bf_layer(
     $         tmp_nodes,
     $         gen_coords,
     $         bf_alignment,
     $         bf_nodes)

          test_loc = is_real_matrix3D_validated(
     $         tmp_nodes,
     $         tmp_nodes_test,
     $         detailled)
          test_loc = test_loc.and.test_validated
          if(detailled.and.(.not.test_loc)) then
             print '(''no optional arg failed'')'
          end if


          !extract_param_out optional argument
          !------------------------------------------------------------
          tmp_nodes = reshape((/
     $         (((-99.0d0,i=1,7),j=1,6),k=1,ne)/),
     $         (/7,6,ne/))

          call get_nodes_from_bf_layer(
     $         tmp_nodes,
     $         gen_coords,
     $         bf_alignment,
     $         bf_nodes,
     $         extract_param_out=extract_param)

          test_loc = is_real_matrix3D_validated(
     $         tmp_nodes,
     $         tmp_nodes_test,
     $         detailled)
          test_loc = test_loc.and.test_validated
          if(detailled.and.(.not.test_loc)) then
             print '(''extract_param_out optional arg failed'')'
          end if


          !extract_param_in optional argument
          !------------------------------------------------------------
          tmp_nodes = reshape((/
     $         (((-99.0d0,i=1,7),j=1,6),k=1,ne)/),
     $         (/7,6,ne/))

          call get_nodes_from_bf_layer(
     $         tmp_nodes,
     $         gen_coords,
     $         bf_alignment,
     $         bf_nodes,
     $         extract_param_in=extract_param)

          test_loc = is_real_matrix3D_validated(
     $         tmp_nodes,
     $         tmp_nodes_test,
     $         detailled)
          test_loc = test_loc.and.test_validated
          if(detailled.and.(.not.test_loc)) then
             print '(''extract_param_in optional arg failed'')'
          end if

        end function test_get_nodes_from_bf_layer


        function test_get_map_from_interior(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind)   , dimension(5)  :: tmp_map_test
          real(rkind)   , dimension(5)  :: tmp_map
          integer(ikind), dimension(2)  :: gen_coords
          real(rkind)   , dimension(20) :: interior_map

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,5

             !input
             call get_param_test_get_map_from_interior(
     $            k,
     $            tmp_map_test,
     $            gen_coords,
     $            interior_map)


             !output
             call get_map_from_interior(
     $            tmp_map,
     $            gen_coords,
     $            interior_map)


             !validation
             test_loc = is_real_vector_validated(
     $            tmp_map,
     $            tmp_map_test,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_validated)) then
                print '(''test('',I1,'') failed'')',k
             end if

          end do

        end function test_get_map_from_interior


        subroutine get_param_test_get_map_from_interior(
     $     test_id,
     $     tmp_map_test,
     $     gen_coords,
     $     interior_map)

          implicit none

          integer                      , intent(in)  :: test_id
          real(rkind)   , dimension(5) , intent(out) :: tmp_map_test
          integer(ikind), dimension(2) , intent(out) :: gen_coords
          real(rkind)   , dimension(20), intent(out) :: interior_map

          integer :: i


          interior_map = (/(i,i=1,20)/)

          select case(test_id)

            case(1)
               gen_coords = [-10,-6]

            case(2)
               gen_coords = [-2,2]

            case(3)
               gen_coords = [5,9]

            case(4)
               gen_coords = [8,12]

            case(5)
               gen_coords = [13,17]

          end select

          tmp_map_test = (/ (gen_coords(1)-1+i,i=1,5) /)

        end subroutine get_param_test_get_map_from_interior


        function test_set_grdpts_id_in_bf_layer(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer                         :: i,j
          integer       , dimension(4,7)  :: tmp_grdpts_id
          integer       , dimension(2,2)  :: gen_coords
          integer(ikind), dimension(2,2)  :: bf_alignment
          integer       , dimension(12,9) :: bf_grdpts_id
          integer       , dimension(12,9) :: bf_grdpts_id_test

          !input
          tmp_grdpts_id = reshape((/
     $         ((i+4*(j-1), i=1,4),j=1,7)/),
     $         (/4,7/))

          gen_coords = reshape((/
     $         align_W+9,align_N+2,align_W+12,align_N+8/),
     $         (/2,2/))

          bf_alignment = reshape((/
     $         align_W+5,align_N,align_W+12,align_N+4/),
     $         (/2,2/))

          bf_grdpts_id = reshape((/
     $         1,1,1,1,1,1,1,1,1,1,1,1,
     $         1,1,1,1,1,1,1,1,1,1,1,1,
     $         2,2,1,1,1,1,1,1,1,1,2,2,
     $         3,2,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,1,1,1,2,3,
     $         3,2,1,1,1,1,1,2,2,2,2,3,
     $         3,2,1,1,1,1,1,2,3,3,3,3,
     $         3,2,2,2,2,2,2,2,3,0,0,0,
     $         3,3,3,3,3,3,3,3,3,0,0,0/),
     $         (/12,9/))

          bf_grdpts_id_test = reshape((/
     $         1,1,1,1,1,1, 1, 1, 1, 1,1,1,
     $         1,1,1,1,1,1, 1, 1, 1, 1,1,1,
     $         2,2,1,1,1,1, 1, 1, 1, 1,2,2,
     $         3,2,1,1,1,1, 1, 1, 1, 1,2,3,
     $         3,2,1,1,1,1, 1, 2, 3, 4,2,3,
     $         3,2,1,1,1,1, 5, 6, 7, 8,2,3,
     $         3,2,1,1,1,1, 9,10,11,12,3,3,
     $         3,2,2,2,2,2,13,14,15,16,0,0,
     $         3,3,3,3,3,3,17,18,19,20,0,0/),
     $         (/12,9/))

          !output
          call set_grdpts_id_in_bf_layer(
     $         tmp_grdpts_id,
     $         gen_coords,
     $         bf_alignment,
     $         bf_grdpts_id)

          !validation
          test_validated = is_int_matrix_validated(
     $         bf_grdpts_id,
     $         bf_grdpts_id_test,
     $         detailled)

        end function test_set_grdpts_id_in_bf_layer

      end program test_bf_layer_extract
