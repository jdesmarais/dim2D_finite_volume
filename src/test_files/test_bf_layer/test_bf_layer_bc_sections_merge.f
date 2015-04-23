      program test_bf_layer_bc_sections_merge

        use bf_layer_bc_sections_merge_module, only :
     $       reallocate_bc_sections_for_merge,
     $       update_corner_for_merge,
     $       update_anticorner_for_merge,
     $       test_grdpts_id_config,
     $       get_edge_test_param,
     $       get_anticorner_test_param,
     $       get_extent_bc_section_edge

        use check_data_module, only :
     $       is_int_validated,
     $       is_int_vector_validated,
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
     $       W_overlap,
     $       interior_pt,
     $       bc_interior_pt,
     $       bc_pt

        use parameters_constant, only :
     $       left, right

        use parameters_kind, only :
     $       ikind


        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        test_loc = test_get_extent_bc_section_edge(detailled)
        test_validated=test_validated.and.test_loc
        print '(''test_get_extent_bc_section_edge: '',L1)', test_loc
        print '()'


        test_loc = test_get_anticorner_test_param(detailled)
        test_validated=test_validated.and.test_loc
        print '(''test_get_anticorner_test_param: '',L1)', test_loc
        print '()'


        test_loc = test_get_edge_test_param(detailled)
        test_validated=test_validated.and.test_loc
        print '(''test_get_edge_test_param: '',L1)', test_loc
        print '()'


        test_loc = test_test_grdpts_id_config(detailled)
        test_validated=test_validated.and.test_loc
        print '(''test_test_grdpts_id_config: '',L1)', test_loc
        print '()'


        test_loc = test_update_anticorner_for_merge(detailled)
        test_validated=test_validated.and.test_loc
        print '(''test_update_anticorner_for_merge: '',L1)', test_loc
        print '()'


        test_loc = test_update_corner_for_merge(detailled)
        test_validated=test_validated.and.test_loc
        print '(''test_update_corner_for_merge: '',L1)', test_loc
        print '()'


        test_loc = test_reallocate_bc_sections_for_merge(detailled)
        test_validated=test_validated.and.test_loc
        print '(''test_reallocate_bc_sections_for_merge: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_reallocate_bc_sections_for_merge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(:,:), allocatable :: bc_sections
          integer                              :: nb_ele_removed
          integer, dimension(5,1)              :: test_bc_sections


          test_validated = .true.

          nb_ele_removed = 0


          ! input
          allocate(bc_sections(5,3))

          bc_sections = reshape((/
     $         no_bc_procedure_type, 1,2, no_overlap, no_overlap,
     $         NW_edge_type, 1,3, no_overlap, no_overlap,
     $         no_bc_procedure_type, 1,3, no_overlap, no_overlap/),
     $         (/5,3/))

          test_bc_sections(:,1) = bc_sections(:,2)

          ! output
          call reallocate_bc_sections_for_merge(
     $         bc_sections,2)

          ! validation
          test_validated = is_int_matrix_validated(
     $         bc_sections,
     $         test_bc_sections,
     $         detailled)

        end function test_reallocate_bc_sections_for_merge

        
        function test_update_corner_for_merge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(5,3) :: bc_sections
          integer, dimension(5,3) :: test_bc_sections

          integer                      :: anticorner_type
          integer(ikind), dimension(2) :: anticorner_position

          integer                      :: corner_type
          integer(ikind), dimension(2) :: corner_position
          integer                      :: corner_overlap

          integer :: nb_ele_removed
          logical :: test_loc


          test_validated = .true.

          nb_ele_removed = 0


          ! input
          bc_sections = reshape((/
     $         NE_edge_type, 1,2, no_overlap, no_overlap,
     $         NW_edge_type, 1,3, no_overlap, no_overlap,
     $         SE_corner_type, 1,3, no_overlap, no_overlap/),
     $         (/5,3/))

          test_bc_sections = bc_sections

          anticorner_type     = NE_edge_type
          anticorner_position = [1,3]

          corner_type = SW_corner_type
          corner_position = [1,3]
          corner_overlap = S_overlap

          ! output
          call update_corner_for_merge(
     $         anticorner_type,
     $         anticorner_position,
     $         corner_type,
     $         corner_position,
     $         corner_overlap,
     $         bc_sections,
     $         nb_ele_removed)

          ! validation
          test_loc = is_int_matrix_validated(
     $         bc_sections,
     $         test_bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''update_corner_for_merge(1) failed'')'
          end if

          test_loc = nb_ele_removed.eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_ele_removed(1) failed'')'
          end if
          

          ! input
          anticorner_type       = NW_edge_type
          corner_type           = SE_corner_type
          test_bc_sections(1,2) = no_bc_procedure_type
          test_bc_sections(5,3) = S_overlap

          ! output
          call update_corner_for_merge(
     $         anticorner_type,
     $         anticorner_position,
     $         corner_type,
     $         corner_position,
     $         corner_overlap,
     $         bc_sections,
     $         nb_ele_removed)

          ! validation
          test_loc = is_int_matrix_validated(
     $         bc_sections,
     $         test_bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''update_anticorner_for_merge(2) failed'')'
          end if

          test_loc = nb_ele_removed.eq.1
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_ele_removed(2) failed'')'
          end if

        end function test_update_corner_for_merge
        

        function test_update_anticorner_for_merge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(5,3) :: bc_sections
          integer, dimension(5,3) :: test_bc_sections

          integer                      :: anticorner_type
          integer(ikind), dimension(2) :: anticorner_position

          integer :: nb_ele_removed
          logical :: test_loc


          test_validated = .true.

          nb_ele_removed = 0


          ! input
          bc_sections = reshape((/
     $         NE_edge_type, 1,2, no_overlap, no_overlap,
     $         NW_edge_type, 1,3, no_overlap, no_overlap,
     $         SE_corner_type, 1,3, no_overlap, no_overlap/),
     $         (/5,3/))

          test_bc_sections = bc_sections

          anticorner_type     = NE_edge_type
          anticorner_position = [1,3]

          ! output
          call update_anticorner_for_merge(
     $         anticorner_type,
     $         anticorner_position,
     $         bc_sections,
     $         nb_ele_removed)

          ! validation
          test_loc = is_int_matrix_validated(
     $         bc_sections,
     $         test_bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''update_anticorner_for_merge(1) failed'')'
          end if

          test_loc = nb_ele_removed.eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_ele_removed(1) failed'')'
          end if
          

          ! input
          anticorner_type     = NW_edge_type
          test_bc_sections(1,2) = no_bc_procedure_type

          ! output
          call update_anticorner_for_merge(
     $         anticorner_type,
     $         anticorner_position,
     $         bc_sections,
     $         nb_ele_removed)

          ! validation
          test_loc = is_int_matrix_validated(
     $         bc_sections,
     $         test_bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''update_anticorner_for_merge(2) failed'')'
          end if

          test_loc = nb_ele_removed.eq.1
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_ele_removed(2) failed'')'
          end if

        end function test_update_anticorner_for_merge


        function test_test_grdpts_id_config(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer, dimension(12)  :: config_original
          integer, dimension(12)  :: config_modified
          integer, dimension(3,4) :: config
          integer, dimension(3,4) :: test_config

          integer :: k
          logical :: test_loc


          test_validated = .true.


          config_original = (/(k,k=1,12)/)
          config = reshape(config_original,(/3,4/))
          test_config = config

          test_loc = test_grdpts_id_config(config,test_config)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test('',I2,'') failed'')', 0
          end if

          do k=1,12
             
             config_modified    = config_original
             config_modified(k) = 0

             test_config = reshape(config_modified,(/3,4/))

             test_loc = .not.test_grdpts_id_config(config,test_config)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_test_grdpts_id_config


        function test_get_edge_test_param(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          logical       , dimension(2)   :: side
          integer(ikind), dimension(5,4) :: bc_sections

          integer(ikind), dimension(2,2,2,4) :: test_grdpts_ex_borders
          integer(ikind), dimension(2,2,4)   :: test_merge_loc_borders
          integer(ikind), dimension(4,4,2,4) :: test_merge_array
          integer(ikind), dimension(2,2,4)   :: test_over_loc_borders
          integer(ikind), dimension(4,4,2,4) :: test_over_array
          integer       , dimension(2,4)     :: test_merge_anticorner_type
          integer       , dimension(2,2,4)   :: test_merge_anticorner_position
          integer       , dimension(2,4)     :: test_over_corner_type
          integer       , dimension(2,2,4)   :: test_over_corner_position
          integer       , dimension(4)       :: test_over_corner_overlap
          integer(ikind), dimension(3,2,4)   :: test_edge_new_position

          integer(ikind), dimension(2,2) :: grdpts_ex_borders
          integer(ikind), dimension(2,2) :: merge_loc_borders
          integer(ikind), dimension(4,4) :: merge_array
          integer(ikind), dimension(2,2) :: over_loc_borders
          integer(ikind), dimension(4,4) :: over_array
          integer(ikind)                 :: merge_anticorner_type
          integer(ikind), dimension(2)   :: merge_anticorner_position
          integer(ikind)                 :: over_corner_type
          integer(ikind), dimension(2)   :: over_corner_position
          integer                        :: over_corner_overlap
          integer(ikind), dimension(3)   :: edge_new_position

          integer :: k,l
          logical :: test_loc


          side = [left,right]

          bc_sections(:,1) = [N_edge_type,1,2,5,no_overlap]
          bc_sections(:,2) = [S_edge_type,1,2,5,no_overlap]
          bc_sections(:,3) = [E_edge_type,2,1,5,no_overlap]
          bc_sections(:,4) = [W_edge_type,2,1,5,no_overlap]


          test_grdpts_ex_borders(:,:,1,1) = reshape((/-1, 2,0,5/),(/2,2/))
          test_grdpts_ex_borders(:,:,2,1) = reshape((/ 6, 2,7,5/),(/2,2/))
          test_grdpts_ex_borders(:,:,1,2) = reshape((/-1, 0,0,3/),(/2,2/))
          test_grdpts_ex_borders(:,:,2,2) = reshape((/ 6, 0,7,3/),(/2,2/))
          test_grdpts_ex_borders(:,:,1,3) = reshape((/ 2,-1,5,0/),(/2,2/))
          test_grdpts_ex_borders(:,:,2,3) = reshape((/ 2, 6,5,7/),(/2,2/))
          test_grdpts_ex_borders(:,:,1,4) = reshape((/ 0,-1,3,0/),(/2,2/))
          test_grdpts_ex_borders(:,:,2,4) = reshape((/ 0, 6,3,7/),(/2,2/))

          test_merge_loc_borders(:,:,1) = reshape((/1,1,2,4/),(/2,2/))
          test_merge_loc_borders(:,:,2) = reshape((/1,1,2,4/),(/2,2/))
          test_merge_loc_borders(:,:,3) = reshape((/1,1,4,2/),(/2,2/))
          test_merge_loc_borders(:,:,4) = reshape((/1,1,4,2/),(/2,2/))

          test_merge_array(:,:,1,1) = reshape((/2,2,0,0,2,3,0,0,2,3,0,0,3,3,0,0/),(/4,4/))
          test_merge_array(:,:,2,1) = reshape((/2,2,0,0,3,2,0,0,3,2,0,0,3,3,0,0/),(/4,4/))
          test_merge_array(:,:,1,2) = reshape((/3,3,0,0,2,3,0,0,2,3,0,0,2,2,0,0/),(/4,4/))
          test_merge_array(:,:,2,2) = reshape((/3,3,0,0,3,2,0,0,3,2,0,0,2,2,0,0/),(/4,4/))
          test_merge_array(:,:,1,3) = reshape((/2,2,2,3,2,3,3,3,0,0,0,0,0,0,0,0/),(/4,4/))
          test_merge_array(:,:,2,3) = reshape((/2,3,3,3,2,2,2,3,0,0,0,0,0,0,0,0/),(/4,4/))
          test_merge_array(:,:,1,4) = reshape((/3,2,2,2,3,3,3,2,0,0,0,0,0,0,0,0/),(/4,4/))
          test_merge_array(:,:,2,4) = reshape((/3,3,3,2,3,2,2,2,0,0,0,0,0,0,0,0/),(/4,4/))

          test_over_loc_borders(:,:,1) = reshape((/1,1,2,3/),(/2,2/))
          test_over_loc_borders(:,:,2) = reshape((/1,2,2,4/),(/2,2/))
          test_over_loc_borders(:,:,3) = reshape((/1,1,3,2/),(/2,2/))
          test_over_loc_borders(:,:,4) = reshape((/2,1,4,2/),(/2,2/))

          test_over_array(:,:,1,1) = reshape((/2,2,0,0,2,3,0,0,3,3,0,0,0,0,0,0/),(/4,4/))
          test_over_array(:,:,2,1) = reshape((/2,2,0,0,3,2,0,0,3,3,0,0,0,0,0,0/),(/4,4/))
          test_over_array(:,:,1,2) = reshape((/0,0,0,0,3,3,0,0,2,3,0,0,2,2,0,0/),(/4,4/))
          test_over_array(:,:,2,2) = reshape((/0,0,0,0,3,3,0,0,3,2,0,0,2,2,0,0/),(/4,4/))
          test_over_array(:,:,1,3) = reshape((/2,2,3,0,2,3,3,0,0,0,0,0,0,0,0,0/),(/4,4/))
          test_over_array(:,:,2,3) = reshape((/2,3,3,0,2,2,3,0,0,0,0,0,0,0,0,0/),(/4,4/))
          test_over_array(:,:,1,4) = reshape((/0,3,2,2,0,3,3,2,0,0,0,0,0,0,0,0/),(/4,4/))
          test_over_array(:,:,2,4) = reshape((/0,3,3,2,0,3,2,2,0,0,0,0,0,0,0,0/),(/4,4/))

          test_merge_anticorner_type(:,1) = [NE_edge_type,NW_edge_type]
          test_merge_anticorner_type(:,2) = [SE_edge_type,SW_edge_type]
          test_merge_anticorner_type(:,3) = [NE_edge_type,SE_edge_type]
          test_merge_anticorner_type(:,4) = [NW_edge_type,SW_edge_type]

          test_merge_anticorner_position(:,1,1) = [-1, 2]
          test_merge_anticorner_position(:,2,1) = [ 6, 2]
          test_merge_anticorner_position(:,1,2) = [-1, 2]
          test_merge_anticorner_position(:,2,2) = [ 6, 2]
          test_merge_anticorner_position(:,1,3) = [ 2,-1]
          test_merge_anticorner_position(:,2,3) = [ 2, 6]
          test_merge_anticorner_position(:,1,4) = [ 2,-1]
          test_merge_anticorner_position(:,2,4) = [ 2, 6]

          test_over_corner_type(1,1) = NE_corner_type
          test_over_corner_type(2,1) = NW_corner_type
          test_over_corner_type(1,2) = SE_corner_type
          test_over_corner_type(2,2) = SW_corner_type
          test_over_corner_type(1,3) = NE_corner_type
          test_over_corner_type(2,3) = SE_corner_type
          test_over_corner_type(1,4) = NW_corner_type
          test_over_corner_type(2,4) = SW_corner_type

          test_over_corner_position(:,1,1) = [-1, 3]
          test_over_corner_position(:,2,1) = [ 6, 3]
          test_over_corner_position(:,1,2) = [-1, 1]
          test_over_corner_position(:,2,2) = [ 6, 1]
          test_over_corner_position(:,1,3) = [ 3,-1]
          test_over_corner_position(:,2,3) = [ 3, 6]
          test_over_corner_position(:,1,4) = [ 1,-1]
          test_over_corner_position(:,2,4) = [ 1, 6]

          test_over_corner_overlap = [S_overlap,N_overlap,W_overlap,E_overlap]

          test_edge_new_position(:,1,1) = [-1, 2, 5]
          test_edge_new_position(:,2,1) = [ 1, 2, 7]
          test_edge_new_position(:,1,2) = [-1, 2, 5]
          test_edge_new_position(:,2,2) = [ 1, 2, 7]
          test_edge_new_position(:,1,3) = [ 2,-1, 5]
          test_edge_new_position(:,2,3) = [ 2, 1, 7]
          test_edge_new_position(:,1,4) = [ 2,-1, 5]
          test_edge_new_position(:,2,4) = [ 2, 1, 7]

          
          test_validated = .true.


          do k=1,4
             do l=1,2

                merge_array = reshape((/
     $               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/),(/4,4/))

                over_array =  reshape((/
     $               0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/),(/4,4/))

                ! output
                call get_edge_test_param(
     $               bc_sections(:,k),
     $               side(l),
     $               
     $               grdpts_ex_borders,
     $               
     $               merge_loc_borders,
     $               merge_array,
     $               
     $               over_loc_borders,
     $               over_array,
     $               
     $               merge_anticorner_type,
     $               merge_anticorner_position,
     $               
     $               over_corner_type,
     $               over_corner_position,
     $               over_corner_overlap,
     $               
     $               edge_new_position)

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
     $               merge_loc_borders,
     $               test_merge_loc_borders(:,:,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''merge_loc_borders('',2I2,'') failed'')', k,l
                end if

                test_loc = is_int_matrix_validated(
     $               merge_array,
     $               test_merge_array(:,:,l,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''merge_array('',2I2,'') failed'')', k,l
                end if

                test_loc = is_int_matrix_validated(
     $               over_loc_borders,
     $               test_over_loc_borders(:,:,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''over_loc_borders('',2I2,'') failed'')', k,l
                end if

                test_loc = is_int_matrix_validated(
     $               over_array,
     $               test_over_array(:,:,l,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''over_array('',2I2,'') failed'')', k,l
                end if

                test_loc = is_int_validated(
     $               merge_anticorner_type,
     $               test_merge_anticorner_type(l,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''merge_anticorner_type('',2I2,'') failed'')', k,l
                end if

                test_loc = is_int_vector_validated(
     $               merge_anticorner_position,
     $               test_merge_anticorner_position(:,l,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''merge_anticorner_position('',2I2,'') failed'')', k,l
                end if

                test_loc = is_int_validated(
     $               over_corner_type,
     $               test_over_corner_type(l,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''over_corner_type('',2I2,'') failed'')', k,l
                end if

                test_loc = is_int_vector_validated(
     $               over_corner_position,
     $               test_over_corner_position(:,l,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''over_corner_position('',2I2,'') failed'')', k,l
                end if
                
                test_loc = is_int_validated(
     $               over_corner_overlap,
     $               test_over_corner_overlap(k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''over_corner_overlap('',2I2,'') failed'')', k,l
                end if

                test_loc = is_int_vector_validated(
     $               edge_new_position,
     $               test_edge_new_position(:,l,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''edge_new_position('',2I2,'') failed'')', k,l
                end if

             end do
          end do


        end function test_get_edge_test_param


        function test_get_anticorner_test_param(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          logical       , dimension(2)   :: side
          integer(ikind), dimension(5,4) :: bc_sections

          integer(ikind), dimension(2,2,2,4) :: test_grdpts_ex_borders
          integer(ikind), dimension(2,2,4)   :: test1_loc_borders
          integer(ikind), dimension(2,2,2,4) :: test1_array
          integer(ikind), dimension(2,2,2,4) :: test2_loc_borders
          integer(ikind), dimension(2,2,2,4) :: test2_array
          integer(ikind), dimension(4)       :: test_new_anticorner_type

          integer(ikind), dimension(2,2) :: grdpts_ex_borders
          integer(ikind), dimension(2,2) :: loc_borders1
          integer(ikind), dimension(2,2) :: array1
          integer(ikind), dimension(2,2) :: loc_borders2
          integer(ikind), dimension(2,2) :: array2
          integer(ikind)                 :: new_anticorner_type

          integer :: k,l
          logical :: test_loc


          side = [left,right]

          bc_sections(:,1) = [NW_edge_type,1,2,no_overlap,no_overlap]
          bc_sections(:,2) = [NE_edge_type,1,2,no_overlap,no_overlap]
          bc_sections(:,3) = [SW_edge_type,1,2,no_overlap,no_overlap]
          bc_sections(:,4) = [SE_edge_type,1,2,no_overlap,no_overlap]


          test_grdpts_ex_borders(:,:,1,1) = reshape((/-1, 2,0,3/),(/2,2/))
          test_grdpts_ex_borders(:,:,2,1) = reshape((/ 1, 4,2,5/),(/2,2/))
          test_grdpts_ex_borders(:,:,1,2) = reshape((/ 3, 2,4,3/),(/2,2/))
          test_grdpts_ex_borders(:,:,2,2) = reshape((/ 1, 4,2,5/),(/2,2/))
          test_grdpts_ex_borders(:,:,1,3) = reshape((/-1, 2,0,3/),(/2,2/))
          test_grdpts_ex_borders(:,:,2,3) = reshape((/ 1, 0,2,1/),(/2,2/))
          test_grdpts_ex_borders(:,:,1,4) = reshape((/ 3, 2,4,3/),(/2,2/))
          test_grdpts_ex_borders(:,:,2,4) = reshape((/ 1, 0,2,1/),(/2,2/))

          test1_loc_borders(:,:,1) = reshape((/1,1,2,2/),(/2,2/))
          test1_loc_borders(:,:,2) = reshape((/1,1,2,2/),(/2,2/))
          test1_loc_borders(:,:,3) = reshape((/1,1,2,2/),(/2,2/))
          test1_loc_borders(:,:,4) = reshape((/1,1,2,2/),(/2,2/))

          test1_array(:,:,1,1) = reshape((/3,2,3,3/),(/2,2/))
          test1_array(:,:,2,1) = reshape((/3,2,3,3/),(/2,2/))
          test1_array(:,:,1,2) = reshape((/2,3,3,3/),(/2,2/))
          test1_array(:,:,2,2) = reshape((/2,3,3,3/),(/2,2/))
          test1_array(:,:,1,3) = reshape((/3,3,3,2/),(/2,2/))
          test1_array(:,:,2,3) = reshape((/3,3,3,2/),(/2,2/))
          test1_array(:,:,1,4) = reshape((/3,3,2,3/),(/2,2/))
          test1_array(:,:,2,4) = reshape((/3,3,2,3/),(/2,2/))

          test2_loc_borders(:,:,1,1) = reshape((/2,1,2,2/),(/2,2/))
          test2_loc_borders(:,:,2,1) = reshape((/1,1,2,1/),(/2,2/))
          test2_loc_borders(:,:,1,2) = reshape((/1,1,1,2/),(/2,2/))
          test2_loc_borders(:,:,2,2) = reshape((/1,1,2,1/),(/2,2/))
          test2_loc_borders(:,:,1,3) = reshape((/2,1,2,2/),(/2,2/))
          test2_loc_borders(:,:,2,3) = reshape((/1,2,2,2/),(/2,2/))
          test2_loc_borders(:,:,1,4) = reshape((/1,1,1,2/),(/2,2/))
          test2_loc_borders(:,:,2,4) = reshape((/1,2,2,2/),(/2,2/))

          test2_array(:,:,1,1) = reshape((/0,3,0,3/),(/2,2/))
          test2_array(:,:,2,1) = reshape((/3,3,0,0/),(/2,2/))
          test2_array(:,:,1,2) = reshape((/3,0,3,0/),(/2,2/))
          test2_array(:,:,2,2) = reshape((/3,3,0,0/),(/2,2/))
          test2_array(:,:,1,3) = reshape((/0,3,0,3/),(/2,2/))
          test2_array(:,:,2,3) = reshape((/0,0,3,3/),(/2,2/))
          test2_array(:,:,1,4) = reshape((/3,0,3,0/),(/2,2/))
          test2_array(:,:,2,4) = reshape((/0,0,3,3/),(/2,2/))

          test_new_anticorner_type(1) = NW_corner_type
          test_new_anticorner_type(2) = NE_corner_type
          test_new_anticorner_type(3) = SW_corner_type
          test_new_anticorner_type(4) = SE_corner_type

          
          test_validated = .true.


          do k=1,4
             do l=1,2

                array1 = reshape((/0,0,0,0/),(/2,2/))

                array2 = reshape((/0,0,0,0/),(/2,2/))

                ! output
                call get_anticorner_test_param(
     $               bc_sections(:,k),
     $               side(l),
     $               
     $               grdpts_ex_borders,
     $               
     $               loc_borders1,
     $               array1,
     $               
     $               loc_borders2,
     $               array2,
     $               
     $               new_anticorner_type)

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
     $               test2_loc_borders(:,:,l,k),
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

                test_loc = is_int_validated(
     $               new_anticorner_type,
     $               test_new_anticorner_type(k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''new_anticorner_type('',2I2,'') failed'')', k,l
                end if

             end do
          end do

        end function test_get_anticorner_test_param


        function test_get_extent_bc_section_edge(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer       , parameter             :: nb_tests = 4
          integer(ikind), dimension(5,nb_tests) :: test_bc_sections
          integer(ikind), dimension(nb_tests)   :: test_extent

          integer :: k
          logical :: test_loc


          test_bc_sections = reshape((/
     $         N_edge_type,1,2,4, no_overlap,
     $         S_edge_type,1,2,3, no_overlap,
     $         E_edge_type,7,1,9, no_overlap,
     $         W_edge_type,1,5,9, no_overlap/),
     $         (/5,4/))

          test_extent = [4,3,9,5]

          
          test_validated = .true.


          do k=1, nb_tests

             test_loc = is_int_validated(
     $            get_extent_bc_section_edge(test_bc_sections(:,k)),
     $            test_extent(k),
     $            detailled)

             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test( '',I2,'' ) failed'')', k
             end if

          end do

        end function test_get_extent_bc_section_edge

      end program test_bf_layer_bc_sections_merge
