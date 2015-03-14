      program test_bf_layer_newgrdpt

        use bf_layer_newgrdpt_class, only :
     $       bf_layer_newgrdpt

        use parameters_kind, only :
     $       ikind

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled      = .true.
        test_validated = .true.


        test_loc = test_get_grdpts_id_part(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_grdpts_id_part: '',L1)', test_loc
        print '()'


        test_loc = test_set_grdpts_id_part(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_set_grdpts_id_part: '',L1)', test_loc
        print '()'        


        contains


        function test_get_grdpts_id_part(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_layer_newgrdpt)        :: bf_layer_used
          integer(ikind), dimension(2,2) :: bf_alignment0
          integer(ikind), dimension(2,2) :: bf_alignment1
          integer(ikind), dimension(2,2) :: gen_coords
          integer       , dimension(6,3) :: tmp_grdpts_id

          integer(ikind) :: i,j
          logical        :: test_loc


          test_validated = .true.


          !bf layer at t
          bf_alignment1 = reshape(
     $         (/align_E,align_S+5,align_E+2,align_S+7/),
     $         (/2,2/))

           bf_layer_used%alignment = bf_alignment1

           allocate(bf_layer_used%grdpts_id(7,7))
           
           bf_layer_used%grdpts_id = reshape((/
     $          1,1,2,3,3,3,0,
     $          1,1,2,2,2,3,3,
     $          1,1,1,1,2,2,3,
     $          1,1,1,1,1,2,3,
     $          1,1,1,1,2,2,3,
     $          1,1,2,2,2,3,3,
     $          1,1,2,3,3,3,0/),
     $          (/7,7/))

          !bf_layer at t-dt
          bf_alignment0 = reshape(
     $         (/align_E,align_S+5,align_E+2,align_S+6/),
     $         (/2,2/))

          call bf_layer_used%bf_compute_used%allocate_tables(
     $         7,6,
     $         bf_alignment0,
     $         (/((i-1)*0.1d0,i=1,7)/),
     $         (/((j-1)*0.1d0,j=1,6)/),
     $         reshape((/
     $             1,1,2,3,3,3,3,
     $             1,1,2,2,2,2,3,
     $             1,1,1,1,1,2,3,
     $             1,1,1,1,1,2,3,
     $             1,1,2,2,2,2,3,
     $             1,1,2,3,3,3,3/),
     $         (/7,6/)))

          !tmp_grdpts_id
          gen_coords = reshape(
     $         (/align_E-2,align_S+6,align_E+3,align_S+8/),
     $         (/2,2/))


          !test at t-dt
          !------------------------------------------------------------

          !extraction at t-dt
          tmp_grdpts_id = reshape(
     $         (/ ((no_pt,i=1,6),j=1,3) /),
     $         (/6,3/))

          call bf_layer_used%get_grdpts_id_part(
     $         tmp_grdpts_id,
     $         gen_coords,
     $         previous_step=.true.)

          test_loc = is_int_matrix_validated(
     $         tmp_grdpts_id,
     $         reshape((/
     $            1,1,1,1,1,2,
     $            1,1,2,2,2,2,
     $            1,1,2,3,3,3/),
     $            (/6,3/)),
     $         detailled)
          test_validated = test_validated.and.test_loc

          if(detailled.and.(.not.test_loc)) then
             print '(''extraction at t-dt failed'')'
          end if


          !test at t
          !------------------------------------------------------------

          !extraction at t
          tmp_grdpts_id = reshape(
     $         (/ ((no_pt,i=1,6),j=1,3) /),
     $         (/6,3/))

          call bf_layer_used%get_grdpts_id_part(
     $         tmp_grdpts_id,
     $         gen_coords)

          test_loc = is_int_matrix_validated(
     $         tmp_grdpts_id,
     $         reshape((/
     $            1,1,1,1,1,2,
     $            1,1,1,1,2,2,
     $            1,1,2,2,2,3/),
     $            (/6,3/)),
     $         detailled)
          test_validated = test_validated.and.test_loc

          if(detailled.and.(.not.test_loc)) then
             print '(''extraction at t failed'')'
          end if

        end function test_get_grdpts_id_part


        function test_set_grdpts_id_part(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_layer_newgrdpt)        :: bf_layer_used
          integer(ikind), dimension(2,2) :: bf_alignment1
          integer(ikind), dimension(2,2) :: gen_coords
          integer       , dimension(6,3) :: tmp_grdpts_id

          logical :: test_loc


          test_validated = .true.


          !input
          bf_alignment1 = reshape(
     $         (/align_E,align_S+5,align_E+2,align_S+7/),
     $         (/2,2/))

           bf_layer_used%alignment = bf_alignment1

           allocate(bf_layer_used%grdpts_id(7,7))
           
           bf_layer_used%grdpts_id = reshape((/
     $          1,1,2,3,3,3,0,
     $          1,1,2,2,2,3,3,
     $          1,1,1,1,2,2,3,
     $          1,1,1,1,1,2,3,
     $          1,1,1,1,2,2,3,
     $          1,1,2,2,2,3,3,
     $          1,1,2,3,3,3,0/),
     $          (/7,7/))          

          gen_coords = reshape(
     $         (/align_E-2,align_S+6,align_E+3,align_S+8/),
     $         (/2,2/))

          tmp_grdpts_id = reshape((/
     $          1,1,1,1,1,2,
     $          1,1,2,2,2,2,
     $          1,1,2,3,3,3/),
     $         (/6,3/))

          !output
          call bf_layer_used%set_grdpts_id_part(
     $         tmp_grdpts_id,
     $         gen_coords)

          !validation
          test_loc = is_int_matrix_validated(
     $         bf_layer_used%grdpts_id,
     $         reshape((/
     $            1,1,2,3,3,3,0,
     $            1,1,2,2,2,3,3,
     $            1,1,1,1,2,2,3,
     $            1,1,1,1,1,2,3,
     $            1,1,2,2,2,2,3,
     $            1,1,2,3,3,3,3,
     $            1,1,2,3,3,3,0/),
     $         (/7,7/)),
     $         detailled)
          test_validated = test_validated.and.test_loc

        end function test_set_grdpts_id_part


      end program test_bf_layer_newgrdpt
