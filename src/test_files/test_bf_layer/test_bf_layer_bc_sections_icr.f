      program test_bf_layer_bc_sections_icr

        use bf_layer_bc_sections_icr_module, only :
     $     get_edge_crenel_id_param

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       no_bc_procedure_type,
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $       no_overlap,
     $       N_overlap,
     $       S_overlap,
     $       E_overlap,
     $       W_overlap

        use parameters_constant, only :
     $       left, right

        use parameters_kind, only :
     $       ikind


        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled      = .true.
        test_validated = .true.


        test_loc = test_get_edge_crenel_id_param(detailled)
        test_validated=test_validated.and.test_loc
        print '(''test_get_edge_crenel_id_param: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated


        contains


        function test_get_edge_crenel_id_param(detailled)
     $       result(test_validated)
        
          implicit none
          
          logical, intent(in) :: detailled
          logical             :: test_validated


          logical       , dimension(2)   :: side
          integer(ikind), dimension(5,4) :: bc_sections

          integer(ikind), dimension(2,2,2,4) :: test_grdpts_ex_borders
          integer(ikind), dimension(2,2,4)   :: test1_loc_borders
          integer(ikind), dimension(3,3,2,4) :: test1_array
          integer(ikind), dimension(2,2,4)   :: test2_loc_borders
          integer(ikind), dimension(3,3,2,4) :: test2_array

          integer(ikind), dimension(2,2) :: grdpts_ex_borders
          integer(ikind), dimension(2,2) :: loc_borders1
          integer(ikind), dimension(3,3) :: array1
          integer(ikind), dimension(2,2) :: loc_borders2
          integer(ikind), dimension(3,3) :: array2

          integer :: k,l
          logical :: test_loc


          side = [left,right]
          
          bc_sections(:,1) = [N_edge_type,1,2,5,no_overlap]
          bc_sections(:,2) = [S_edge_type,1,2,5,no_overlap]
          bc_sections(:,3) = [E_edge_type,2,1,5,no_overlap]
          bc_sections(:,4) = [W_edge_type,2,1,5,no_overlap]


          test_grdpts_ex_borders(:,:,1,1) = reshape((/-1, 2,0,4/),(/2,2/))
          test_grdpts_ex_borders(:,:,2,1) = reshape((/ 6, 2,7,4/),(/2,2/))
          test_grdpts_ex_borders(:,:,1,2) = reshape((/-1, 1,0,3/),(/2,2/))
          test_grdpts_ex_borders(:,:,2,2) = reshape((/ 6, 1,7,3/),(/2,2/))
          test_grdpts_ex_borders(:,:,1,3) = reshape((/ 2,-1,4,0/),(/2,2/))
          test_grdpts_ex_borders(:,:,2,3) = reshape((/ 2, 6,4,7/),(/2,2/))
          test_grdpts_ex_borders(:,:,1,4) = reshape((/ 1,-1,3,0/),(/2,2/))
          test_grdpts_ex_borders(:,:,2,4) = reshape((/ 1, 6,3,7/),(/2,2/))

          test1_loc_borders(:,:,1) = reshape((/1,1,2,3/),(/2,2/))
          test1_loc_borders(:,:,2) = reshape((/1,1,2,3/),(/2,2/))
          test1_loc_borders(:,:,3) = reshape((/1,1,3,2/),(/2,2/))
          test1_loc_borders(:,:,4) = reshape((/1,1,3,2/),(/2,2/))

          test1_array(:,:,1,1) = reshape((/2,2,0,2,3,0,2,3,0/),(/3,3/))
          test1_array(:,:,2,1) = reshape((/2,2,0,3,2,0,3,2,0/),(/3,3/))
          test1_array(:,:,1,2) = reshape((/2,3,0,2,3,0,2,2,0/),(/3,3/))
          test1_array(:,:,2,2) = reshape((/3,2,0,3,2,0,2,2,0/),(/3,3/))
          test1_array(:,:,1,3) = reshape((/2,2,2,2,3,3,0,0,0/),(/3,3/))
          test1_array(:,:,2,3) = reshape((/2,3,3,2,2,2,0,0,0/),(/3,3/))
          test1_array(:,:,1,4) = reshape((/2,2,2,3,3,2,0,0,0/),(/3,3/))
          test1_array(:,:,2,4) = reshape((/3,3,2,2,2,2,0,0,0/),(/3,3/))

          test2_loc_borders(:,:,1) = reshape((/1,1,2,3/),(/2,2/))
          test2_loc_borders(:,:,2) = reshape((/1,1,2,3/),(/2,2/))
          test2_loc_borders(:,:,3) = reshape((/1,1,3,2/),(/2,2/))
          test2_loc_borders(:,:,4) = reshape((/1,1,3,2/),(/2,2/))

          test2_array(:,:,1,1) = reshape((/2,2,0,2,3,0,3,3,0/),(/3,3/))
          test2_array(:,:,2,1) = reshape((/2,2,0,3,2,0,3,3,0/),(/3,3/))
          test2_array(:,:,1,2) = reshape((/3,3,0,2,3,0,2,2,0/),(/3,3/))
          test2_array(:,:,2,2) = reshape((/3,3,0,3,2,0,2,2,0/),(/3,3/))
          test2_array(:,:,1,3) = reshape((/2,2,3,2,3,3,0,0,0/),(/3,3/))
          test2_array(:,:,2,3) = reshape((/2,3,3,2,2,3,0,0,0/),(/3,3/))
          test2_array(:,:,1,4) = reshape((/3,2,2,3,3,2,0,0,0/),(/3,3/))
          test2_array(:,:,2,4) = reshape((/3,3,2,3,2,2,0,0,0/),(/3,3/))


          test_validated = .true.


          do k=1,4
             do l=1,2

                array1 = reshape((/
     $               0,0,0,0,0,0,0,0,0/),(/3,3/))

                array2 =  reshape((/
     $               0,0,0,0,0,0,0,0,0/),(/3,3/))


                ! output
                call get_edge_crenel_id_param(
     $               bc_sections(:,k),
     $               side(l),
     $               grdpts_ex_borders,
     $               loc_borders1,
     $               array1,
     $               loc_borders2,
     $               array2)
                
                ! validation
                test_loc = is_int_matrix_validated(
     $               grdpts_ex_borders,
     $               test_grdpts_ex_borders(:,:,l,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''grdpts_ex_borders('',2I2,'') failed'')', k,l
                end if

                 test_loc = is_int_matrix_validated(
     $               loc_borders1,
     $               test1_loc_borders(:,:,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''loc_borders1('',2I2,'') failed'')', k,l
                end if

                test_loc = is_int_matrix_validated(
     $               array1,
     $               test1_array(:,:,l,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''array1('',2I2,'') failed'')', k,l
                end if

                test_loc = is_int_matrix_validated(
     $               loc_borders2,
     $               test2_loc_borders(:,:,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''loc_borders2('',2I2,'') failed'')', k,l
                end if

                test_loc = is_int_matrix_validated(
     $               array2,
     $               test2_array(:,:,l,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''array2('',2I2,'') failed'')', k,l
                end if

             end do
          end do

        end function test_get_edge_crenel_id_param


      end program test_bf_layer_bc_sections_icr
