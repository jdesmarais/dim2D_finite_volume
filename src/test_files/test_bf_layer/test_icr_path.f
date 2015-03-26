      program test_icr_path
      
        use bf_interface_coords_class, only :
     $       bf_interface_coords

        use bf_sublayer_class, only :
     $       bf_sublayer

        use dim2d_parameters, only :
     $       cv_r

        use icr_path_class, only :
     $       icr_path

        use check_data_module, only :
     $       is_real_validated,
     $       is_real_vector_validated,
     $       is_int_vector_validated,
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W

        use parameters_constant, only :
     $       N,S,E,
     $       obc_eigenqties_bc,
     $       obc_eigenqties_lin

        use parameters_input, only :
     $       nx,ny,ne,
     $       obc_eigenqties_strategy

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

        
        test_loc = test_reinitialize(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_reinitialize: '',L1)', test_loc
        print '()'


        test_loc = test_add_pt_to_path(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_add_pt_to_path: '',L1)', test_loc
        print '()'

        
        test_loc = test_share_update_operations_with(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_share_update_operations_with: '',L1)', test_loc
        print '()'        


        test_loc = test_update_alignment(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_update_alignment: '',L1)', test_loc
        print '()'


        test_loc = test_are_pts_in_same_mainlayer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_are_pts_in_same_mainlayer: '',L1)', test_loc
        print '()'


        test_loc = test_are_pts_in_same_path(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_are_pts_in_same_path: '',L1)', test_loc
        print '()'


        test_loc = test_stage_first_grdpt_for_update(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_stage_first_grdpt_for_update: '',L1)', test_loc
        print '()'


        test_loc = test_stage_next_grdpt_for_update(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_stage_next_grdpt_for_update: '',L1)', test_loc
        print '()'


        test_loc = test_stage(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_stage: '',L1)', test_loc
        print '()'


        test_loc = test_should_bf_layers_be_merged(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_should_bf_layers_be_merged: '',L1)', test_loc
        print '()'


        test_loc = test_update_allocation_bf_layer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_update_allocation_bf_layer: '',L1)', test_loc
        print '()'


        test_loc = test_commit(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_commit: '',L1)', test_loc
        print '()'

        
        test_loc = test_merge(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_merge: '',L1)', test_loc
        print '()'


        test_loc = test_remove(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_remove: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_ini(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(icr_path) :: icr_path_used


          test_validated = .true.
          

          call icr_path_used%ini()

          test_loc = .not.associated(icr_path_used%matching_sublayer)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''matching sublayer failed'')'
          end if

          test_loc = icr_path_used%nb_pts.eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_pts failed'')'
          end if

          test_loc = icr_path_used%ends.eqv.(.false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''ends failed'')'
          end if

        end function test_ini


        function test_reinitialize(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(icr_path) :: icr_path_used


          test_validated = .true.
          

          call icr_path_used%reinitialize()

          test_loc = .not.associated(icr_path_used%matching_sublayer)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''matching sublayer failed'')'
          end if

          test_loc = icr_path_used%nb_pts.eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_pts failed'')'
          end if

          test_loc = icr_path_used%ends.eqv.(.false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''ends failed'')'
          end if

        end function test_reinitialize


        function test_add_pt_to_path(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(icr_path)                  :: icr_path_used
          integer(ikind), dimension(2,11) :: test_gen_coords

          integer :: i,j
          logical :: test_loc


          test_validated = .true.


          !input
          test_gen_coords = reshape((/
     $         ((2*(j-1) + (i-1), i=1,2),j=1,11)/),
     $         (/2,11/))

          call icr_path_used%ini()


          !output
          do j=1,size(test_gen_coords,2)
             call icr_path_used%add_pt_to_path(test_gen_coords(:,j))
          end do


          !validation
          test_loc = icr_path_used%nb_pts.eq.11
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_pts failed'')'
          end if

          if(test_loc) then
             test_loc = is_int_matrix_validated(
     $            icr_path_used%pts(:,1:icr_path_used%nb_pts),
     $            test_gen_coords,
     $            detailled)
             test_validated=test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''gen_coords failed'')'
             end if

          end if

        end function test_add_pt_to_path


        function test_share_update_operations_with(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(icr_path) :: path1
          type(icr_path) :: path2
          logical        :: test_share
          logical        :: share

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,11

             !input
             call ini_paths_should_merge(k,path1,path2,test_share)

             !output
             share = path1%share_update_operations_with(path2)

             !validation
             test_loc = share.eqv.test_share
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')',k
             end if
        
          end do

        end function test_share_update_operations_with      

        
        subroutine ini_paths_should_merge(test_id,path1,path2,test_share)

          implicit none

          integer       , intent(in)    :: test_id
          type(icr_path), intent(inout) :: path1
          type(icr_path), intent(inout) :: path2
          logical       , intent(out)   :: test_share


          select case(test_id)

            case(1)
               call path1%ini()
               call path2%ini()
               path1%nb_pts= 1
               test_share = .false.

            case(2)
               call path1%ini()
               call path2%ini()
               path2%nb_pts= 1
               test_share = .false.

            case(3)
               call path1%ini()
               call path2%ini()
               path1%nb_pts= 2
               path2%nb_pts= 1
               
               path1%mainlayer_id = N
               path2%mainlayer_id = S

               test_share = .false.

            case(4)
               call path1%ini()
               call path2%ini()
               path1%nb_pts= 2
               path2%nb_pts= 1
               
               path1%mainlayer_id = N
               path2%mainlayer_id = N

               path1%alignment = reshape((/
     $              align_W+2, align_N, align_W+3, align_N/),
     $              (/2,2/))

               path2%alignment = reshape((/
     $              align_W+10, align_N, align_W+11, align_N/),
     $              (/2,2/))

               test_share = .false.

            case(5)
               path1%alignment = reshape((/
     $              align_W+3, align_N, align_W+4, align_N/),
     $              (/2,2/))
               test_share = .true.

            case(6)
               path1%alignment = reshape((/
     $              align_W+17, align_N, align_W+18, align_N/),
     $              (/2,2/))
               test_share = .true.

            case(7)
               path1%alignment = reshape((/
     $              align_W+18, align_N, align_W+19, align_N/),
     $              (/2,2/))
               test_share = .false.

            case(8)
               call path1%ini()
               call path2%ini()
               path1%nb_pts= 2
               path2%nb_pts= 1
               
               path1%mainlayer_id = E
               path2%mainlayer_id = E

               path1%alignment = reshape((/
     $              align_E, align_S+1, align_E, align_S+2/),
     $              (/2,2/))

               path2%alignment = reshape((/
     $              align_E, align_S+9, align_E, align_S+10/),
     $              (/2,2/))

               test_share = .false.

            case(9)
               path1%alignment = reshape((/
     $              align_E, align_S+2, align_E, align_S+3/),
     $              (/2,2/))
               test_share = .true.

            case(10)
               path1%alignment = reshape((/
     $              align_E, align_S+16, align_E, align_S+17/),
     $              (/2,2/))
               test_share = .true.

            case(11)
               path1%alignment = reshape((/
     $              align_E, align_S+17, align_E, align_S+18/),
     $              (/2,2/))
               test_share = .false.

          end select

        end subroutine ini_paths_should_merge


        function is_path_validated(path1,path2,detailled)
     $     result(test_validated)

          implicit none

          type(icr_path), intent(in) :: path1
          type(icr_path), intent(in) :: path2
          logical       , intent(in) :: detailled
          logical                    :: test_validated

          logical :: test_loc


          test_validated = .true.


          !nb_pts
          test_loc = path1%nb_pts.eq.path2%nb_pts
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_pts: '',I2,''->'',I2)', 
     $            path1%nb_pts,
     $            path2%nb_pts
          end if

          if(path1%nb_pts.gt.0) then

             !mainlayer_id
             test_loc = path1%mainlayer_id.eq.path2%mainlayer_id
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''mainlayer_id: '',I2,''->'',I2)', 
     $               path1%mainlayer_id,
     $               path2%mainlayer_id
             end if
   
             !alignment
             test_loc = is_int_matrix_validated(
     $            path1%alignment,
     $            path2%alignment,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''alignment failed'')'
             end if
   
             !matching sublayer
             if(
     $            associated(path1%matching_sublayer).and.
     $            associated(path2%matching_sublayer)) then
             
                test_loc = associated(
     $               path1%matching_sublayer,
     $               path2%matching_sublayer)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''matching_sublayer failed'')'
                end if
   
             end if

             !pts
             if(path1%nb_pts.eq.path2%nb_pts) then
                test_loc = is_int_matrix_validated(
     $               path1%pts(:,1:path1%nb_pts),
     $               path2%pts(:,1:path2%nb_pts),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''pts failed'')'
                end if
             end if

          end if

        end function is_path_validated


        function test_update_alignment(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(icr_path)                          :: icr_path_used
          integer       , parameter               :: nb_tests = 4
          integer(ikind), dimension(2,nb_tests)   :: gen_coords
          integer(ikind), dimension(2,2,nb_tests) :: test_alignment

          integer :: k


          test_validated = .true.


          !input
          !------------------------------------------------------------
          !no update from the exchange of neighbors
          test_alignment(:,:,1) = reshape((/
     $         align_W+5,align_N,align_W+5,align_N/),
     $         (/2,2/))

          !no update from the gen_coords
          gen_coords(:,2) = [align_W+5,align_N]

          test_alignment(:,:,2) = reshape((/
     $         align_W+5,align_N,align_W+5,align_N/),
     $         (/2,2/))

          !update from the gen_coords
          gen_coords(:,3) = [align_W+6,align_N+1]
          
          test_alignment(:,:,3) = reshape((/
     $         align_W+5,align_N,align_W+6,align_N+1/),
     $         (/2,2/))

          !update from the exchange of neighbors
          gen_coords(:,4) = [align_W+4,align_N+2]
          
          test_alignment(:,:,4) = reshape((/
     $         align_W+1,align_N,align_W+6,align_N+2/),
     $         (/2,2/))


          !tests
          !------------------------------------------------------------
          !test w/o optional argument passed
          call icr_path_used%ini()
          icr_path_used%mainlayer_id = N
          icr_path_used%alignment = reshape((/
     $         align_W+5,align_N,align_W+5,align_N/),
     $         (/2,2/))

          call icr_path_used%update_alignment()

          test_loc = is_int_matrix_validated(
     $         icr_path_used%alignment,
     $         test_alignment(:,:,1),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test('',I2,'') failed'')', 1
          end if


          do k=2, nb_tests

             !output
             call icr_path_used%update_alignment(gen_coords=gen_coords(:,k))

             !validation
             test_loc = is_int_matrix_validated(
     $            icr_path_used%alignment,
     $            test_alignment(:,:,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_update_alignment


        function test_are_pts_in_same_mainlayer(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(icr_path)                 :: icr_path_used
          integer , parameter            :: nb_tests=2
          integer, dimension(2,nb_tests) :: gen_coords
          logical, dimension(nb_tests)   :: test_same_mainlayer
          logical                        :: same_mainlayer

          integer :: k
          logical :: test_loc


          test_validated = .true.


          call icr_path_used%ini()
          icr_path_used%mainlayer_id = N
          gen_coords(:,1) = [align_E,align_N-1]
          gen_coords(:,2) = [align_E,align_N]
          test_same_mainlayer(1) = .false.
          test_same_mainlayer(2) = .true.

          do k=1,nb_tests

             !output
             same_mainlayer = icr_path_used%are_pts_in_same_mainlayer(
     $            gen_coords(:,k))

             !validation
             test_loc = same_mainlayer.eqv.test_same_mainlayer(k)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_are_pts_in_same_mainlayer


        function test_are_pts_in_same_path(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(icr_path)                 :: icr_path_used
          integer , parameter            :: nb_tests=2
          integer, dimension(2,nb_tests) :: gen_coords
          logical, dimension(nb_tests)   :: test_same_path
          integer                        :: tolerance
          logical                        :: same_path

          integer :: k
          logical :: test_loc


          test_validated = .true.


          call icr_path_used%ini()
          icr_path_used%mainlayer_id = N
          icr_path_used%alignment = reshape((/
     $         align_E-3,align_N,align_E+1,align_N/),
     $         (/2,2/))

          gen_coords(:,1) = [align_E+3,align_N]
          gen_coords(:,2) = [align_E+2,align_N]
          test_same_path(1) = .false.
          test_same_path(2) = .true.
          tolerance = 1


          do k=1,nb_tests

             !output
             same_path= icr_path_used%are_pts_in_same_path(
     $            gen_coords(:,k),
     $            tolerance)

             !validation
             test_loc = same_path.eqv.test_same_path(k)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'',1) failed'')', k
             end if

          end do


          call icr_path_used%ini()
          icr_path_used%mainlayer_id = E
          icr_path_used%alignment = reshape((/
     $         align_E,align_S,align_E,align_S+3/),
     $         (/2,2/))

          gen_coords(:,1) = [align_E,align_S+5]
          gen_coords(:,2) = [align_E,align_S+4]
          test_same_path(1) = .false.
          test_same_path(2) = .true.
          tolerance = 1


          do k=1,nb_tests

             !output
             same_path= icr_path_used%are_pts_in_same_path(
     $            gen_coords(:,k),
     $            tolerance)

             !validation
             test_loc = same_path.eqv.test_same_path(k)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'',2) failed'')', k
             end if

          end do

        end function test_are_pts_in_same_path


        function test_stage_first_grdpt_for_update(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(icr_path)                          :: icr_path_used
          type(bf_interface_coords)               :: bf_interface_used
          integer       , parameter               :: nb_tests = 4
          integer(ikind), dimension(2,nb_tests)   :: gen_coords
          integer       , dimension(nb_tests)     :: tolerance
          integer(ikind), dimension(2,2,nb_tests) :: test_alignment
          logical       , dimension(nb_tests)     :: test_matching_sublayer
          

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          integer(ikind), dimension(2,2)   :: alignment_tmp
          type(bf_sublayer), pointer       :: bf_sublayer_N_ptr

          integer :: k
          logical :: test_loc


          test_validated = .true.


          !input
          !------------------------------------------------------------
          call bf_interface_used%ini(interior_x_map,interior_y_map)

          alignment_tmp = reshape((/
     $         align_W+11,align_N,align_W+11,align_N/),
     $         (/2,2/))

          bf_sublayer_N_ptr => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         alignment_tmp)


          !test inputs
          !------------------------------------------------------------
          !no matching sublayer
          gen_coords(:,1) = [align_W+17,align_N]

          tolerance(1) = 1
          
          test_alignment(:,:,1) = reshape((/
     $         align_W+17,align_N,align_W+17,align_N/),
     $         (/2,2/))

          test_matching_sublayer(1) = .false.


          !matching sublayer, no influence from the neighbor update
          gen_coords(:,2) = [align_W+14,align_N]
          
          tolerance(2) = 1

          test_alignment(:,:,2) = reshape((/
     $         align_W+11,align_N,align_W+14,align_N/),
     $         (/2,2/))

          test_matching_sublayer(2) = .true.


          !no matching sublayer, influence from the neighbor update
          gen_coords(:,3) = [align_W+4,align_N]
          
          tolerance(3) = 1

          test_alignment(:,:,3) = reshape((/
     $         align_W+1,align_N,align_W+4,align_N/),
     $         (/2,2/))

          test_matching_sublayer(3) = .false.


          !matching sublayer, influence from the neighbor update
          gen_coords(:,4) = [align_W+4,align_N]

          tolerance(4) = 5

          test_alignment(:,:,4) = reshape((/
     $         align_W+1,align_N,align_W+11,align_N/),
     $         (/2,2/))

          test_matching_sublayer(4) = .true.


          !tests
          !------------------------------------------------------------
          do k=1,2

             !input
             call icr_path_used%ini()

             !output
             call icr_path_used%stage_first_grdpt_for_update(
     $            gen_coords(:,k),
     $            bf_interface_used,
     $            tolerance(k))

             !validation
             !------------------------------------------------------------
             !mainlayer_id
             test_loc = icr_path_used%mainlayer_id.eq.N
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''mainlayer id('',I2,'') failed'')',k
             end if

             !matching sublayer
             if(test_matching_sublayer(k)) then
                test_loc = associated(
     $               icr_path_used%matching_sublayer,
     $               bf_sublayer_N_ptr)
             else
                test_loc = .not.associated(icr_path_used%matching_sublayer)
             end if
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''matching sulayer('',I2,'') failed'')',k
             end if

             !alignment
             test_loc = is_int_matrix_validated(
     $            icr_path_used%alignment,
     $            test_alignment(:,:,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''alignment('',I2,'') failed'')',k
             end if

             !pt in path
             test_loc = is_int_vector_validated(
     $            icr_path_used%pts(:,1),
     $            gen_coords(:,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''pts('',I2,'') failed'')',k
             end if

          end do

        end function test_stage_first_grdpt_for_update


        function test_stage_next_grdpt_for_update(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(icr_path)                          :: icr_path_used
          type(bf_interface_coords)               :: bf_interface_used
          integer       , parameter               :: nb_tests = 5

          integer(ikind), dimension(2)            :: gen_coords_ini
          integer                                 :: tolerance_pts_same_sublayer
          integer                                 :: tolerance_pts_same_path

          integer(ikind), dimension(2,nb_tests)   :: gen_coords
          integer(ikind), dimension(2,2,nb_tests) :: test_alignment
          integer       , dimension(nb_tests)     :: test_matching_sublayer
          logical       , dimension(nb_tests)     :: test_ends

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          integer(ikind), dimension(2,2)   :: alignment_tmp
          type(bf_sublayer), pointer       :: bf_sublayer_N1_ptr
          type(bf_sublayer), pointer       :: bf_sublayer_N2_ptr

          integer :: k
          logical :: test_loc


          test_validated = .true.


          !input
          !------------------------------------------------------------
          call bf_interface_used%ini(interior_x_map,interior_y_map)

          tolerance_pts_same_path     = 2
          tolerance_pts_same_sublayer = 1

          gen_coords_ini = [align_W+5,align_N]


          !not same mainlayer          
          gen_coords(:,1)           = [align_W+5,align_N-1]
          test_ends(1)              = .true.
          test_matching_sublayer(1) = 0
          test_alignment(:,:,1)     = reshape((/
     $                                    align_W+5,align_N,
     $                                    align_W+5,align_N/),
     $                                (/2,2/))

          gen_coords(:,2)           = [align_W+13,align_N]
          test_ends(2)              = .true.
          test_matching_sublayer(2) = 0
          test_alignment(:,:,2)     = reshape((/
     $                                    align_W+5,align_N,
     $                                    align_W+5,align_N/),
     $                                (/2,2/))

          gen_coords(:,3)           = [align_W+7,align_N]
          test_ends(3)              = .false.
          test_matching_sublayer(3) = 0
          test_alignment(:,:,3)     = reshape((/
     $                                    align_W+5 ,align_N,
     $                                    align_W+7,align_N/),
     $                                (/2,2/))

          gen_coords(:,4)           = [align_W+7,align_N]
          test_ends(4)              = .false.
          test_matching_sublayer(4) = 2
          test_alignment(:,:,4)     = reshape((/
     $                                    align_W+5 ,align_N,
     $                                    align_W+10,align_N/),
     $                                (/2,2/))

          gen_coords(:,5)           = [align_W+7,align_N]
          test_ends(5)              = .false.
          test_matching_sublayer(5) = 1
          test_alignment(:,:,5)     = reshape((/
     $                                    align_W+1 ,align_N,
     $                                    align_W+7,align_N/),
     $                                (/2,2/))

          

          do k=1,nb_tests

             !input
             if(k.eq.4) then
                alignment_tmp = reshape((/
     $               align_W+10,align_N,align_W+10,align_N/),
     $               (/2,2/))
                
                bf_sublayer_N2_ptr => bf_interface_used%allocate_sublayer(
     $               N,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               alignment_tmp)
             end if

             if(k.eq.5) then
                alignment_tmp = reshape((/
     $               align_W+4,align_N,align_W+5,align_N/),
     $               (/2,2/))
                
                bf_sublayer_N1_ptr => bf_interface_used%allocate_sublayer(
     $               N,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               alignment_tmp)
             end if

             call icr_path_used%ini()
             call icr_path_used%stage_first_grdpt_for_update(
     $            gen_coords_ini,
     $            bf_interface_used,
     $            tolerance_pts_same_sublayer)

             
             call icr_path_used%stage_next_grdpt_for_update(
     $            gen_coords(:,k),
     $            bf_interface_used,
     $            tolerance_pts_same_path,
     $            tolerance_pts_same_sublayer)


             !validation
             !------------------------------------------------------------
             !ends
             test_loc = icr_path_used%ends.eqv.test_ends(k)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test ends('',I2,'') failed'')', k
             end if

             if(.not.icr_path_used%ends) then
             !matching sublayer
                select case(test_matching_sublayer(k))
                  case(0)
                     test_loc = .not.associated(icr_path_used%matching_sublayer)
                  case(1)
                     test_loc = associated(icr_path_used%matching_sublayer,bf_sublayer_N1_ptr)
                  case(2)
                     test_loc = associated(icr_path_used%matching_sublayer,bf_sublayer_N2_ptr)
                end select
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test matching_sublayer('',I2,'') failed'')',k
                end if

             !alignment
                test_loc = is_int_matrix_validated(
     $               icr_path_used%alignment,
     $               test_alignment(:,:,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test alignment('',I2,'') failed'')',k
                end if
                
             !pts
                test_loc = is_int_matrix_validated(
     $               icr_path_used%pts(:,1:icr_path_used%nb_pts),
     $               reshape((/
     $                 gen_coords_ini(1), gen_coords_ini(2),
     $                 gen_coords(1,k)  , gen_coords(2,k)/),
     $                 (/2,2/)),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test pts('',I2,'') failed'')',k
                end if

             end if

          end do

        end function test_stage_next_grdpt_for_update


        function test_stage(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(icr_path)                          :: icr_path_used
          type(bf_interface_coords)               :: bf_interface_used
          integer       , parameter               :: nb_tests = 5

          integer(ikind), dimension(2)            :: gen_coords_ini
          integer                                 :: tolerance_pts_same_sublayer
          integer                                 :: tolerance_pts_same_path

          integer(ikind), dimension(2,nb_tests)   :: gen_coords
          integer(ikind), dimension(2,2,nb_tests) :: test_alignment
          integer       , dimension(nb_tests)     :: test_matching_sublayer
          logical       , dimension(nb_tests)     :: test_ends

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          integer(ikind), dimension(2,2)   :: alignment_tmp
          type(bf_sublayer), pointer       :: bf_sublayer_N1_ptr
          type(bf_sublayer), pointer       :: bf_sublayer_N2_ptr

          integer :: k
          logical :: test_loc


          test_validated = .true.


          !input
          !------------------------------------------------------------
          call bf_interface_used%ini(interior_x_map,interior_y_map)

          tolerance_pts_same_path     = 2
          tolerance_pts_same_sublayer = 1

          gen_coords_ini = [align_W+5,align_N]


          !not same mainlayer          
          gen_coords(:,1)           = [align_W+5,align_N-1]
          test_ends(1)              = .true.
          test_matching_sublayer(1) = 0
          test_alignment(:,:,1)     = reshape((/
     $                                    align_W+5,align_N,
     $                                    align_W+5,align_N/),
     $                                (/2,2/))

          gen_coords(:,2)           = [align_W+12,align_N]
          test_ends(2)              = .true.
          test_matching_sublayer(2) = 0
          test_alignment(:,:,2)     = reshape((/
     $                                    align_W+5,align_N,
     $                                    align_W+5,align_N/),
     $                                (/2,2/))

          gen_coords(:,3)           = [align_W+11,align_N]
          test_ends(3)              = .false.
          test_matching_sublayer(3) = 0
          test_alignment(:,:,3)     = reshape((/
     $                                    align_W+5 ,align_N,
     $                                    align_W+11,align_N/),
     $                                (/2,2/))

          gen_coords(:,4)           = [align_W+11,align_N]
          test_ends(4)              = .false.
          test_matching_sublayer(4) = 2
          test_alignment(:,:,4)     = reshape((/
     $                                    align_W+5 ,align_N,
     $                                    align_W+17,align_N/),
     $                                (/2,2/))

          gen_coords(:,5)           = [align_W+11,align_N]
          test_ends(5)              = .false.
          test_matching_sublayer(5) = 1
          test_alignment(:,:,5)     = reshape((/
     $                                    align_W+1 ,align_N,
     $                                    align_W+11,align_N/),
     $                                (/2,2/))

          

          do k=1,nb_tests

             !input
             if(k.eq.4) then
                alignment_tmp = reshape((/
     $               align_W+17,align_N,align_W+17,align_N/),
     $               (/2,2/))
                
                bf_sublayer_N2_ptr => bf_interface_used%allocate_sublayer(
     $               N,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               alignment_tmp)
             end if

             if(k.eq.5) then
                alignment_tmp = reshape((/
     $               align_W+4,align_N,align_W+4,align_N/),
     $               (/2,2/))
                
                bf_sublayer_N1_ptr => bf_interface_used%allocate_sublayer(
     $               N,
     $               interior_x_map,
     $               interior_y_map,
     $               interior_nodes,
     $               alignment_tmp)
             end if

             call icr_path_used%ini()
             call icr_path_used%stage(gen_coords_ini,bf_interface_used)
             call icr_path_used%stage(gen_coords(:,k),bf_interface_used)


             !validation
             !------------------------------------------------------------
             !ends
             test_loc = icr_path_used%ends.eqv.test_ends(k)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test ends('',I2,'') failed'')', k
             end if

             if(.not.icr_path_used%ends) then
             !matching sublayer
                select case(test_matching_sublayer(k))
                  case(0)
                     test_loc = .not.associated(icr_path_used%matching_sublayer)
                  case(1)
                     test_loc = associated(icr_path_used%matching_sublayer,bf_sublayer_N1_ptr)
                  case(2)
                     test_loc = associated(icr_path_used%matching_sublayer,bf_sublayer_N2_ptr)
                end select
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test matching_sublayer('',I2,'') failed'')',k
                end if

             !alignment
                test_loc = is_int_matrix_validated(
     $               icr_path_used%alignment,
     $               test_alignment(:,:,k),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test alignment('',I2,'') failed'')',k
                end if
                
             !pts
                test_loc = is_int_matrix_validated(
     $               icr_path_used%pts(:,1:icr_path_used%nb_pts),
     $               reshape((/
     $                 gen_coords_ini(1), gen_coords_ini(2),
     $                 gen_coords(1,k)  , gen_coords(2,k)/),
     $                 (/2,2/)),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test pts('',I2,'') failed'')',k
                end if

             end if

          end do

        end function test_stage


       function test_should_bf_layers_be_merged(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(icr_path)                   :: icr_path_used
          integer , parameter              :: nb_tests=2
          integer, dimension(2,2,nb_tests) :: test_alignment
          logical, dimension(nb_tests)     :: test_merge_layers
          logical                          :: merge_layers

          type(bf_sublayer) :: bf_sublayer_ptr

          integer :: k
          logical :: test_loc


          test_validated = .true.


          bf_sublayer_ptr%localization = N
          bf_sublayer_ptr%alignment = reshape((/
     $         align_W+5,align_N,align_W+5,align_N/),
     $         (/2,2/))


          test_alignment(:,:,1) = reshape((/
     $         align_W+11,align_N,align_W+11,align_N/),
     $         (/2,2/))
          test_merge_layers(1) = .true.

          test_alignment(:,:,2) = reshape((/
     $         align_W+12,align_N,align_W+12,align_N/),
     $         (/2,2/))
          test_merge_layers(2) = .false.


          do k=1,nb_tests

             !input
             call icr_path_used%ini()
             icr_path_used%mainlayer_id = N
             icr_path_used%alignment = test_alignment(:,:,k)

             !output
             merge_layers = icr_path_used%should_bf_layers_be_merged(
     $            bf_sublayer_ptr)

             !validation
             test_loc = merge_layers.eqv.test_merge_layers(k)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_should_bf_layers_be_merged


        function test_update_allocation_bf_layer(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(icr_path)            :: icr_path_used
          type(bf_interface_coords) :: bf_interface_used

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          type(bf_sublayer), pointer :: modified_sublayer

          logical :: test_loc


          test_validated = .true.


          call bf_interface_used%ini(interior_x_map,interior_y_map)


          !test allocation
          !------------------------------------------------------------

          !input
          call icr_path_used%ini()
          call icr_path_used%stage([align_W+5,align_N],bf_interface_used)

          !output
          modified_sublayer => icr_path_used%update_allocation_bf_layer(
     $         bf_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes)
          
          !validation
          test_loc = is_int_matrix_validated(
     $         modified_sublayer%get_alignment_tab(),
     $         reshape((/
     $            align_W+5,align_N,align_W+5,align_N/),
     $            (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test allocation failed'')'
          end if


          !test reallocation w/o neighbor
          !------------------------------------------------------------
          !input
          call icr_path_used%reinitialize()
          call icr_path_used%stage([align_W+6,align_N],bf_interface_used)

          !output
          modified_sublayer => icr_path_used%update_allocation_bf_layer(
     $         bf_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes)
          
          !validation
          test_loc = is_int_matrix_validated(
     $         modified_sublayer%get_alignment_tab(),
     $         reshape((/
     $            align_W+5,align_N,align_W+6,align_N/),
     $            (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test reallocation w/o neighbor failed'')'
          end if


          !test reallocation w/ neighbor too far
          !------------------------------------------------------------
          !input
          call icr_path_used%reinitialize()
          call icr_path_used%stage([align_W+20,align_N],bf_interface_used)

          modified_sublayer => icr_path_used%update_allocation_bf_layer(
     $         bf_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes)

          call icr_path_used%reinitialize()
          call icr_path_used%stage([align_W+7,align_N],bf_interface_used)

          modified_sublayer => icr_path_used%update_allocation_bf_layer(
     $         bf_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes)
          
          !validation
          test_loc = is_int_matrix_validated(
     $         modified_sublayer%get_alignment_tab(),
     $         reshape((/
     $            align_W+5,align_N,align_W+7,align_N/),
     $            (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test reallocation w/ neighbor too far failed'')'
          end if


          !test merge
          !------------------------------------------------------------
          !input
          call icr_path_used%reinitialize()
          call icr_path_used%stage([align_W+13,align_N],bf_interface_used)
          call icr_path_used%stage([align_W+14,align_N],bf_interface_used)

          modified_sublayer => icr_path_used%update_allocation_bf_layer(
     $         bf_interface_used,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes)

          !validation
          test_loc = is_int_matrix_validated(
     $         modified_sublayer%get_alignment_tab(),
     $         reshape((/
     $            align_W+5,align_N,align_W+20,align_N/),
     $            (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test merge failed'')'
          end if

        end function test_update_allocation_bf_layer


        function test_commit(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(icr_path)            :: icr_path_used
          type(bf_interface_coords) :: bf_interface_used
          type(pmodel_eq)           :: p_model

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes0
          real(rkind), dimension(nx,ny,ne) :: interior_nodes1

          real(rkind) :: t
          real(rkind) :: dt

          type(bf_sublayer), pointer :: new_sublayer

          integer(ikind)             :: i,j
          real(rkind), dimension(ne) :: test_newgrdpt

          logical        :: test_loc


          test_validated = .true.


          !test description:
          !------------------------------------------------------------
          !    create a new East buffer layer for one new
          !interior_pt. We check the configuration of the
          !grdpts_id as well as the value of one new bc_pt
          !------------------------------------------------------------

          !input
          !------------------------------------------------------------
          call bf_interface_used%ini(interior_x_map,interior_y_map)

          call p_model%initial_conditions%ini_far_field()

          t  = 0.0d0
          dt = 0.25d0

          interior_x_map = (/ ((i-1)*1.00d0,i=1,nx) /)
          interior_y_map = (/ ((j-1)*0.25d0,j=1,ny) /)

          interior_nodes0(align_E-1:align_E+1,
     $                    align_S+7:align_S+9,:)
     $         = reshape((/
     $         1.48d0, 1.30d0, 1.35d0,
     $         1.26d0, 1.45d0, 1.40d0,
     $         1.46d0, 1.27d0, 1.47d0,
     $         
     $         0.128d0, 0.127d0, 0.142d0,
     $         1.138d0, 0.148d0, 0.132d0,
     $         0.146d0, 0.143d0, 0.145d0,
     $         
     $         0.0050d0, 0.020d0, 0.060d0,
     $         0.0025d0, 0.001d0, 0.015d0,
     $         0.0100d0, 0.002d0, 0.050d0,
     $         
     $         4.88d0, 4.870d0, 4.855d0,
     $         4.85d0, 4.865d0, 4.845d0,
     $         4.89d0, 4.870d0, 4.860d0/),
     $         (/3,3,ne/))

          interior_nodes1(align_E-1:align_E+1,
     $                    align_S+7:align_S+9,:)
     $         = reshape((/
     $         1.50d0, 1.455d0, 1.48d0,
     $         1.20d0, 1.350d0, 1.25d0,
     $         1.49d0, 1.250d0, 1.40d0,
     $         
     $         0.128d0, 0.450d0, 0.135d0,
     $         0.148d0, 0.150d0, 0.122d0,
     $         0.142d0, 1.152d0, 0.236d0,
     $         
     $         0.006d0, 0.0600d0, 0.020d0,
     $         0.000d0, 0.0028d0, 0.035d0,
     $         0.020d0, 0.0030d0, 0.040d0,
     $         
     $         4.876d0, 4.825d0, 4.862d0,
     $         4.890d0, 4.871d0, 4.892d0,
     $         4.865d0, 4.757d0, 4.895d0/),
     $         (/3,3,ne/))

          select case(obc_eigenqties_strategy)

            case(obc_eigenqties_bc)
               test_newgrdpt = [
     $              1.22383078395524d0,
     $              0.39531842478603d0,
     $              -0.21050217290879d0,
     $              4.19684181018688d0]
               
            case(obc_eigenqties_lin)
               test_newgrdpt = [
     $              1.21167555521982d0,
     $              0.35901827468671d0,
     $              -0.20475732388332d0,
     $              4.20002914561351d0]
          
            case default
               print '(''test_bf_newgrdpt_prim'')'
               print '(''test_sym_compute_newgrdpt_x'')'
               print '(''obc_eigenqties_strategy not recognized'')'
               stop ''
               
          end select


          !output
          !------------------------------------------------------------
          call icr_path_used%ini()
          call icr_path_used%stage([align_E,align_S+8],bf_interface_used)
          call icr_path_used%commit(
     $         bf_interface_used,
     $         p_model,
     $         t,
     $         dt,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes0,
     $         interior_nodes1)
          call icr_path_used%reinitialize()


          !validation
          !------------------------------------------------------------
          new_sublayer => bf_interface_used%mainlayer_pointers(E)%get_head_sublayer()

          !verification of the alignment
          test_loc = is_int_matrix_validated(
     $         new_sublayer%get_alignment_tab(),
     $         reshape((/
     $            align_E, align_S+8, align_E, align_S+8/),
     $            (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment failed'')'
          end if

          !verification of the configuration of the grdpts_id
          test_loc = is_int_matrix_validated(
     $         new_sublayer%grdpts_id,
     $         reshape((/
     $            1,1,2,3,3,
     $            1,1,2,2,3,
     $            1,1,1,2,3,
     $            1,1,2,2,3,
     $            1,1,2,3,3/),
     $            (/5,5/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test configuration grdpts_id failed'')'
          end if

          !verification of the value of one of the new grid-points
          test_loc = is_real_vector_validated(
     $         new_sublayer%nodes(5,3,:),
     $         test_newgrdpt,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test newgrdpt failed'')'
          end if

        end function test_commit


        function test_merge(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_sublayer), pointer :: bf_sublayer_ptr
          type(icr_path)             :: path1
          type(icr_path)             :: path2
          type(icr_path)             :: path_test
          logical                    :: test_check_for_merge

          integer :: k
          logical :: test_loc

          test_validated = .true.


          allocate(bf_sublayer_ptr)

          do k=1,4

             !input
             call get_param_test_merge(
     $            k,
     $            bf_sublayer_ptr,
     $            path1,
     $            path2,
     $            path_test,
     $            test_check_for_merge)

             !output
             call path1%merge(path2,check_for_merge=test_check_for_merge)

             !validation
             test_loc = is_path_validated(
     $            path1,
     $            path_test,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if
             
          end do

        end function test_merge

        
        subroutine get_param_test_merge(
     $     test_id,
     $     bf_sublayer_ptr,
     $     path1,
     $     path2,
     $     path_test,
     $     test_check_for_merge)

          implicit none

          integer                   , intent(in)    :: test_id
          type(bf_sublayer), pointer, intent(in)    :: bf_sublayer_ptr
          type(icr_path)            , intent(inout) :: path1
          type(icr_path)            , intent(inout) :: path2
          type(icr_path)            , intent(inout) :: path_test
          logical                   , intent(out)   :: test_check_for_merge


          call path1%remove()
          call path2%remove()
          call path_test%remove()
          

          !path1
          !------------------------------------------------------------
          call path1%ini()
          path1%mainlayer_id = N
          path1%nb_pts = 2
          allocate(path1%pts(2,4))
          path1%pts(:,1) = [1,2]
          path1%pts(:,2) = [3,4]

          path1%alignment = reshape((/
     $         align_W+1,align_N,align_W+2,align_N/),
     $         (/2,2/))

          
          !path2
          !------------------------------------------------------------
          call path2%ini()
          path2%mainlayer_id = N
          path2%nb_pts = 3
          allocate(path2%pts(2,5))
          path2%pts(:,1) = [5,6]
          path2%pts(:,2) = [7,8]
          path2%pts(:,3) = [9,10]

          path2%alignment = reshape((/
     $         align_W+8,align_N,align_W+9,align_N/),
     $         (/2,2/))


          !path_test
          !------------------------------------------------------------
          call path_test%ini()
          path_test%mainlayer_id = N
          path_test%nb_pts = 5
          allocate(path_test%pts(2,10))
          path_test%pts(:,1) = [1,2]
          path_test%pts(:,2) = [3,4]
          path_test%pts(:,3) = [5,6]
          path_test%pts(:,4) = [7,8]
          path_test%pts(:,5) = [9,10]

          path_test%alignment = reshape((/
     $         align_W+1,align_N,align_W+9,align_N/),
     $         (/2,2/))


          test_check_for_merge = .true.


          select case(test_id)


            !both paths should not be merged
            !------------------------------------------------------------
            case(1)               

               path2%mainlayer_id = S
               path2%alignment = reshape((/
     $              align_E,align_S,align_E+1,align_S/),
     $              (/2,2/))

               test_check_for_merge = .false.

               path_test%mainlayer_id = N
               path_test%alignment = reshape((/
     $              align_W+1,align_S,align_E+1,align_N/),
     $              (/2,2/))


            !both paths can be merged: no matching sublayers
            !------------------------------------------------------------
            case(2)
               

            !both paths can be merged: matching sublayer on 1
            !------------------------------------------------------------
            case(3)
               path1%matching_sublayer     => bf_sublayer_ptr
               path_test%matching_sublayer => bf_sublayer_ptr


            !both paths can be merged: matching sublayer on 2
            !------------------------------------------------------------
            case(4)
               path2%matching_sublayer     => bf_sublayer_ptr
               path_test%matching_sublayer => bf_sublayer_ptr
               
          end select

        end subroutine get_param_test_merge


        function test_remove(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          type(icr_path)      :: path_used
          

          call path_used%ini()
          
          path_used%nb_pts = 3
          allocate(path_used%pts(2,10))
          
          call path_used%remove()

          test_validated = .not.allocated(path_used%pts)
          if(detailled.and.(.not.test_validated)) then
             print '(''deallocation failed'')'
          end if


        end function test_remove


        subroutine check_inputs()

          implicit none

          type(pmodel_eq)            :: p_model
          real(rkind), dimension(ne) :: far_field

          if(ne.ne.4) then
             stop 'the test requires ne=4: DIM2D model'
          end if

          if(.not.is_real_validated(cv_r,2.5d0,.false.)) then
             stop 'the test requires c_v_r=2.5'
          end if

          call p_model%initial_conditions%ini_far_field()

          far_field = p_model%get_far_field(0.0d0,1.0d0,1.0d0)

          if(.not.is_real_vector_validated(
     $         far_field,
     $         [1.46510213931996d0,0.146510214d0,0.0d0,2.84673289046992d0],
     $         .true.)) then
             print '(''the test requires p_model%get_far_field(t,x,y)='')'
             print '(''[1.465102139d0,0.14651021d0,0.0d0,2.84673289d0]'')'
             print '()'
             print '(''T0 should be 0.95'')'
             print '(''flow_direction should be x-direction'')'
             print '(''ic_choice should be newgrdpt_test'')'
             print '()'
             stop ''
             
          end if

        end subroutine check_inputs

      end program test_icr_path
