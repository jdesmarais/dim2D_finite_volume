      program test_bf_layer_basic

        use bf_layer_basic_class, only :
     $     bf_layer_basic

        use parameters_input, only :
     $       nx,ny

        use parameters_kind, only :
     $       ikind
        

        logical              :: detailled
        logical              :: test_loc
        logical              :: test_validated

        detailled      = .true.
        test_validated = .true.

        test_loc = test_get_local_coords(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_local_coord: '',L1)', test_loc
        print '()'

        test_loc = test_get_general_to_local_coord_tab(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_general_to_local_coord_tab: '',L1)', test_loc
        print '()'

        
        contains


        function test_get_local_coords(detailled)
     $       result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_basic)           :: bf_layer_used
          integer(ikind), dimension(2,2) :: bf_alignment
          integer(ikind), dimension(2)   :: gen_coords
          integer(ikind), dimension(2)   :: local_coords
          integer(ikind), dimension(2)   :: local_coords_test
          
          integer :: k


          !input
          bf_alignment(1,1) = 3
          bf_alignment(1,2) = 6
          bf_alignment(2,1) = 4
          bf_alignment(2,2) = 8
          call bf_layer_used%set_alignment_tab(bf_alignment)
          gen_coords(1) = 4
          gen_coords(2) = 10
          local_coords_test = [4,9]


          !output
          local_coords = bf_layer_used%get_local_coord(gen_coords)


          !validation
          test_validated = 
     $         (local_coords(1).eq.local_coords_test(1)).and.
     $         (local_coords(2).eq.local_coords_test(2))

          if((.not.test_validated).and.detailled) then
             do k=1,2
                print '(I2,'' -> '',I2)', local_coords(k), local_coords_test(k)
             end do
          end if          

        end function test_get_local_coords


        function test_get_general_to_local_coord_tab(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_basic)           :: bf_layer_used
          integer(ikind), dimension(2,2) :: alignment
          integer(ikind), dimension(2)   :: match_table_test
          integer(ikind), dimension(2)   :: match_table

          !input data
          alignment(1,1) = 3
          alignment(1,2) = 5
          alignment(2,1) = 4
          alignment(2,2) = 6
          call bf_layer_used%set_alignment_tab(alignment)

          match_table_test = [0,1]

          !output data
          match_table = bf_layer_used%get_general_to_local_coord_tab()

          !validation
          test_validated = 
     $         (match_table(1).eq.match_table_test(1)).and.
     $         (match_table(2).eq.match_table_test(2))

          if((.not.test_validated).and.detailled) then
             
             print '(I2,'' -> '',I2)', match_table(1), match_table_test(1)
             print '(I2,'' -> '',I2)', match_table(2), match_table_test(2)
             
          end if

        end function test_get_general_to_local_coord_tab

      end program test_bf_layer_basic
