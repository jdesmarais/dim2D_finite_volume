      program test_bf_interface_dyn

        use bf_interface_dyn_class, only :
     $     bf_interface_dyn

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_real_matrix3D_validated,
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       align_N,align_S,
     $       align_E,align_W

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled      = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_allocate_sublayer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_allocate_sublayer: '',L1)', test_loc
        print '()'


        test_loc = test_reallocate_sublayer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_reallocate_sublayer: '',L1)', test_loc
        print '()'


        test_loc = test_merge_sublayer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_merge_sublayer: '',L1)', test_loc
        print '()'


        test_loc = test_remove_sublayer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_remove_sublayer: '',L1)', test_loc
        print '()'


        test_loc = test_uniformize_mainlayer_interfaces(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_uniformize_mainlayer_interfaces: '',L1)', test_loc
        print '()'


        contains


        function test_allocate_sublayer(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,5
             
             test_loc = perform_test_allocate_sublayer(k,detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
                print '()'
             end if

          end do

        end function test_allocate_sublayer


        function test_reallocate_sublayer(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,3
             
             test_loc = perform_test_reallocate_sublayer(2*(k-1)+1,detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
                print '()'
             end if

          end do

        end function test_reallocate_sublayer


        function test_merge_sublayer(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,3
             
             test_loc = perform_test_merge_sublayer(2*(k-1)+1,detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
                print '()'
             end if

          end do

        end function test_merge_sublayer


        function test_remove_sublayer(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,5
             
             test_loc = perform_test_remove_sublayer(k,detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
                print '()'
             end if

          end do

        end function test_remove_sublayer


        function test_uniformize_mainlayer_interfaces(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_dyn)         :: bf_interface_used
          integer(ikind), dimension(2,2) :: bf_alignment_tmp

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          integer    , dimension(2,2,4,2)  :: test_alignment

          type(bf_sublayer), pointer :: added_sublayer
          type(bf_sublayer), pointer :: sublayer_tested

          integer(ikind) :: i,j
          integer        :: k


          test_validated = .true.


          !input
          !------------------------------------------------------------
          interior_nodes = reshape((/
     $         (((200*(k-1)+20*(j-1)+(i-1),i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          call bf_interface_used%ini(interior_x_map,interior_y_map)

          !NW_interface(N)
          bf_alignment_tmp = reshape((/
     $         align_W+4, align_N, align_W+10, align_N+1/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          !NE_interface(N)
          bf_alignment_tmp = reshape((/
     $         align_E-10, align_N, align_E-4, align_N+1/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          !SW_interface(W)
          bf_alignment_tmp = reshape((/
     $         align_W-5, align_S+1, align_W, align_S+2/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         W,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          !NW_interface(W)
          bf_alignment_tmp = reshape((/
     $         align_W-7, align_N-2, align_W, align_N-1/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         W,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          !SE_interface(E)/NE_interface(E)
          bf_alignment_tmp = reshape((/
     $         align_E, align_S+4, align_E+3, align_N-4/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         E,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          !SW_interface(S)
          bf_alignment_tmp = reshape((/
     $         align_W-10, align_S-1, align_W+5, align_S/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         S,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          !SE_interface(S)
          bf_alignment_tmp = reshape((/
     $         align_E-9, align_S-1, align_E+5, align_S/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         S,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)


          !output
          !------------------------------------------------------------
          call bf_interface_used%uniformize_mainlayer_interfaces(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes)


          !validation
          !------------------------------------------------------------
          test_alignment(:,:,N,1) = reshape((/
     $         align_W-7, align_N, align_W+10, align_N+1/),
     $         (/2,2/))
          test_alignment(:,:,N,2) = reshape((/
     $         align_E-10, align_N, align_E+5, align_N+1/),
     $         (/2,2/))

          test_alignment(:,:,S,1) = reshape((/
     $         align_W-10, align_S-1, align_W+5, align_S/),
     $         (/2,2/))
          test_alignment(:,:,S,2) = reshape((/
     $         align_E-9, align_S-1, align_E+5, align_S/),
     $         (/2,2/))

          test_alignment(:,:,E,1) = reshape((/
     $         align_E, align_S+1, align_E+5, align_N-1/),
     $         (/2,2/))

          test_alignment(:,:,W,1) = reshape((/
     $         align_W-10, align_S+1, align_W, align_S+2/),
     $         (/2,2/))
          test_alignment(:,:,W,2) = reshape((/
     $         align_W-7, align_N-2, align_W, align_N-1/),
     $         (/2,2/))
          
          !NW_interface(N)
          sublayer_tested => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         sublayer_tested%alignment,
     $         test_alignment(:,:,N,1),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test NW_interface_N failed'')'
          end if

          !NE_interface(N)
          sublayer_tested => bf_interface_used%mainlayer_pointers(N)%ptr%get_tail_sublayer()
          test_loc = is_int_matrix_validated(
     $         sublayer_tested%alignment,
     $         test_alignment(:,:,N,2),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test NE_interface_N failed'')'
          end if

          !SW_interface(S)
          sublayer_tested => bf_interface_used%mainlayer_pointers(S)%ptr%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         sublayer_tested%alignment,
     $         test_alignment(:,:,S,1),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test SW_interface_S failed'')'
          end if

          !SE_interface(S)
          sublayer_tested => bf_interface_used%mainlayer_pointers(S)%ptr%get_tail_sublayer()
          test_loc = is_int_matrix_validated(
     $         sublayer_tested%alignment,
     $         test_alignment(:,:,S,2),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test NE_interface_N failed'')'
          end if


          !SW_interface(W)
          sublayer_tested => bf_interface_used%mainlayer_pointers(W)%ptr%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         sublayer_tested%alignment,
     $         test_alignment(:,:,W,1),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test SW_interface_W failed'')'
          end if

          !NW_interface(W)
          sublayer_tested => bf_interface_used%mainlayer_pointers(W)%ptr%get_tail_sublayer()
          test_loc = is_int_matrix_validated(
     $         sublayer_tested%alignment,
     $         test_alignment(:,:,W,2),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test NW_interface_W failed'')'
          end if


          !SE_interface(E)/NE_interface(E)
          sublayer_tested => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         sublayer_tested%alignment,
     $         test_alignment(:,:,E,1),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test SE_interface_E/NE_interface_E failed'')'
          end if

        end function test_uniformize_mainlayer_interfaces


        function perform_test_allocate_sublayer(test_id,detailled)
     $     result(test_validated)

          implicit none

          integer, intent(in) :: test_id
          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_dyn)                     :: bf_interface_used
          type(bf_sublayer), pointer                 :: bf_sublayer_tested
          integer(ikind), dimension(2,2)             :: bf_alignment_test
          real(rkind), dimension(:,:,:), allocatable :: nodes_test
          

          !input+output
          call ini_bf_interface_for_test(
     $         test_id,
     $         bf_interface_used)


          !param for validation
          call get_param_test_allocate(
     $         test_id,
     $         bf_interface_used,
     $         bf_sublayer_tested,
     $         bf_alignment_test,
     $         nodes_test)


          !validation
          test_validated = is_test_allocation_validated(
     $         bf_sublayer_tested,
     $         bf_alignment_test,
     $         nodes_test,
     $         detailled)

        end function perform_test_allocate_sublayer


        function perform_test_reallocate_sublayer(test_id,detailled)
     $     result(test_validated)

          implicit none

          integer, intent(in) :: test_id
          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_dyn)                     :: bf_interface_used
          integer(ikind), dimension(2,2)             :: bf_alignment
          type(bf_sublayer), pointer                 :: bf_sublayer_tested
          integer(ikind), dimension(2,2)             :: bf_alignment_test
          real(rkind), dimension(:,:,:), allocatable :: nodes_test

          real(rkind), dimension(nx)                 :: interior_x_map
          real(rkind), dimension(ny)                 :: interior_y_map
          real(rkind), dimension(nx,ny,ne)           :: interior_nodes
          
          integer(ikind) :: i,j
          integer        :: k


          !input
          call ini_bf_interface_for_test(
     $         test_id,
     $         bf_interface_used)


          !param for validation
          call get_param_test_reallocate(
     $         test_id,
     $         bf_interface_used,
     $         bf_alignment,
     $         bf_sublayer_tested,
     $         bf_alignment_test,
     $         nodes_test)


          !output
          interior_nodes = reshape((/
     $         (((200*(k-1)+20*(j-1)+(i-1),i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          call bf_interface_used%reallocate_sublayer(
     $         bf_sublayer_tested,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)


          !validation
          test_validated = is_test_allocation_validated(
     $         bf_sublayer_tested,
     $         bf_alignment_test,
     $         nodes_test,
     $         detailled)

        end function perform_test_reallocate_sublayer


        function perform_test_merge_sublayer(test_id,detailled)
     $     result(test_validated)

          implicit none

          integer, intent(in) :: test_id
          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_dyn)                     :: bf_interface_used
          integer(ikind), dimension(2,2)             :: bf_alignment
          type(bf_sublayer), pointer                 :: bf_sublayer1
          type(bf_sublayer), pointer                 :: bf_sublayer2
          type(bf_sublayer), pointer                 :: merged_sublayer
          integer(ikind), dimension(2,2)             :: bf_alignment_test
          real(rkind), dimension(:,:,:), allocatable :: nodes_test

          real(rkind), dimension(nx)                 :: interior_x_map
          real(rkind), dimension(ny)                 :: interior_y_map
          real(rkind), dimension(nx,ny,ne)           :: interior_nodes
          
          integer(ikind) :: i,j
          integer        :: k


          !input
          call ini_bf_interface_for_test(
     $         test_id,
     $         bf_interface_used)


          !param for validation
          call get_param_test_merge(
     $         test_id,
     $         bf_interface_used,
     $         bf_alignment,
     $         bf_sublayer1,
     $         bf_sublayer2,
     $         bf_alignment_test,
     $         nodes_test)


          !output
          interior_nodes = reshape((/
     $         (((200*(k-1)+20*(j-1)+(i-1),i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          merged_sublayer => bf_interface_used%merge_sublayers(
     $         bf_sublayer1,
     $         bf_sublayer2,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)


          !validation
          test_validated = is_test_allocation_validated(
     $         merged_sublayer,
     $         bf_alignment_test,
     $         nodes_test,
     $         detailled)

        end function perform_test_merge_sublayer


        function perform_test_remove_sublayer(test_id,detailled)
     $     result(test_validated)

          implicit none

          integer, intent(in) :: test_id
          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_dyn)     :: bf_interface_used
          type(bf_sublayer), pointer :: bf_removed_ptr
          type(bf_sublayer), pointer :: bf_sublayer_tested
          logical :: test_loc

          
          test_validated = .true.


          !input
          call ini_bf_interface_for_test(
     $         test_id,
     $         bf_interface_used)


          !output
          bf_removed_ptr => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
          call bf_interface_used%remove_sublayer(bf_removed_ptr)
          

          !validation
          test_loc = .not.associated(bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer())
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test sublayer not associated failed'')'
          end if

          test_loc = .not.associated(bf_interface_used%mainlayer_interfaces%NW_interface_N_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test NW_interface_N_ptr not associated failed'')'
          end if

          test_loc = .not.associated(bf_interface_used%mainlayer_interfaces%NE_interface_N_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test NE_interface_N_ptr not associated failed'')'
          end if

          if(associated(bf_interface_used%mainlayer_interfaces%NW_interface_W_ptr)) then
             bf_sublayer_tested => bf_interface_used%mainlayer_interfaces%NW_interface_W_ptr
             test_loc = .not.(bf_sublayer_tested%can_exchange_with_neighbor2())
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test NW_interface_W_ptr not exchange with neighbor2 failed'')'
             end if
          end if

          if(associated(bf_interface_used%mainlayer_interfaces%NE_interface_E_ptr)) then
             bf_sublayer_tested => bf_interface_used%mainlayer_interfaces%NE_interface_E_ptr
             test_loc = .not.(bf_sublayer_tested%can_exchange_with_neighbor2())
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test NE_interface_E_ptr not exchange with neighbor2 failed'')'
             end if
          end if

        end function perform_test_remove_sublayer


        subroutine ini_bf_interface_for_test(
     $     test_id,
     $     bf_interface_used)

          implicit none

          integer               , intent(in)    :: test_id
          type(bf_interface_dyn), intent(inout) :: bf_interface_used


          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          type(bf_sublayer), pointer :: added_sublayer

          integer(ikind) :: i,j
          integer        :: k

          integer(ikind), dimension(2,2) :: bf_alignment


          interior_nodes = reshape((/
     $         (((200*(k-1)+20*(j-1)+(i-1),i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          call bf_interface_used%ini(
     $         interior_x_map,
     $         interior_y_map)


          select case(test_id)
            case(1)

               bf_alignment = reshape((/
     $              align_W+5,align_N,align_W+10,align_N+1/),
     $              (/2,2/))

               added_sublayer => bf_interface_used%allocate_sublayer(
     $              N,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment)

            case(2)

               bf_alignment = reshape((/
     $              align_W+4,align_N,align_W+9,align_N+1/),
     $              (/2,2/))
               
               added_sublayer => bf_interface_used%allocate_sublayer(
     $              N,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment)

            case(3)

               bf_alignment = reshape((/
     $              align_W-7,align_N-2,align_W,align_N-1/),
     $              (/2,2/))

               added_sublayer => bf_interface_used%allocate_sublayer(
     $              W,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment)

               added_sublayer%nodes = reshape((/
     $              (((200*(k-1)+20*(align_N-5+j-1)+(align_W-10+i-1),i=1,12),j=1,6),k=1,ne)/),
     $              (/12,6,ne/))

               bf_alignment = reshape((/
     $              align_W+4,align_N,align_W+9,align_N+1/),
     $              (/2,2/))

               added_sublayer => bf_interface_used%allocate_sublayer(
     $              N,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment)


            case(4)

               bf_alignment = reshape((/
     $              align_E-9,align_N,align_E-4,align_N+1/),
     $              (/2,2/))

               added_sublayer => bf_interface_used%allocate_sublayer(
     $              N,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment)

            case(5)

               bf_alignment = reshape((/
     $              align_E,align_N-2,align_E+7,align_N-1/),
     $              (/2,2/))

               added_sublayer => bf_interface_used%allocate_sublayer(
     $              E,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment)

               added_sublayer%nodes = reshape((/
     $              (((200*(k-1)+20*(align_N-5+j-1)+(align_E-3+i-1),i=1,12),j=1,6),k=1,ne)/),
     $              (/12,6,ne/))


               bf_alignment = reshape((/
     $              align_E-9,align_N,align_E-4,align_N+1/),
     $              (/2,2/))
               
               added_sublayer => bf_interface_used%allocate_sublayer(
     $              N,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment)
               
            end select

        end subroutine ini_bf_interface_for_test


        subroutine get_param_test_allocate(
     $     test_id,
     $     bf_interface_used,
     $     bf_sublayer_tested,
     $     bf_alignment_test,
     $     nodes_test)

          implicit none

          integer                                   , intent(in)  :: test_id
          type(bf_interface_dyn)                    , intent(in)  :: bf_interface_used
          type(bf_sublayer), pointer                , intent(out) :: bf_sublayer_tested
          integer(ikind), dimension(2,2)            , intent(out) :: bf_alignment_test
          real(rkind), dimension(:,:,:), allocatable, intent(out) :: nodes_test
          

          integer(ikind) :: i,j
          integer        :: k


          bf_sublayer_tested => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()

          select case(test_id)
            case(1)

               bf_alignment_test = reshape((/
     $              align_W+5,align_N,align_W+10,align_N+1/),
     $              (/2,2/))
               
               allocate(nodes_test(10,6,ne))
               nodes_test = reshape((/
     $              (((200*(k-1)+20*(align_N-3+j-1)+(align_W+2+i-1),i=1,10),j=1,6),k=1,ne)/),
     $              (/10,6,ne/))


            case(2)

               bf_alignment_test = reshape((/
     $              align_W+1,align_N,align_W+9,align_N+1/),
     $              (/2,2/))
               
               allocate(nodes_test(13,6,ne))
               nodes_test = reshape((/
     $              (((200*(k-1)+20*(align_N-3+j-1)+(align_W-2+i-1),i=1,13),j=1,6),k=1,ne)/),
     $              (/13,6,ne/))

            case(3)

               bf_alignment_test = reshape((/
     $              align_W-7,align_N,align_W+9,align_N+1/),
     $              (/2,2/))
               
               allocate(nodes_test(21,6,ne))
               nodes_test = reshape((/
     $              (((200*(k-1)+20*(align_N-3+j-1)+(align_W-10+i-1),i=1,21),j=1,6),k=1,ne)/),
     $              (/21,6,ne/))

            case(4)

               bf_alignment_test = reshape((/
     $              align_E-9,align_N,align_E-1,align_N+1/),
     $              (/2,2/))
               
               allocate(nodes_test(13,6,ne))
               nodes_test = reshape((/
     $              (((200*(k-1)+20*(align_N-3+j-1)+(align_E-12+i-1),i=1,13),j=1,6),k=1,ne)/),
     $              (/13,6,ne/))

            case(5)
               
               bf_alignment_test = reshape((/
     $              align_E-9,align_N,align_E+7,align_N+1/),
     $              (/2,2/))
               
               allocate(nodes_test(21,6,ne))
               nodes_test = reshape((/
     $              (((200*(k-1)+20*(align_N-3+j-1)+(align_E-12+i-1),i=1,21),j=1,6),k=1,ne)/),
     $              (/21,6,ne/))

          end select

        end subroutine get_param_test_allocate


        subroutine get_param_test_reallocate(
     $     test_id,
     $     bf_interface_used,
     $     bf_alignment,
     $     bf_sublayer_tested,
     $     bf_alignment_test,
     $     nodes_test)

          implicit none

          integer                                   , intent(in)     :: test_id
          type(bf_interface_dyn)                    , intent(inout)  :: bf_interface_used
          integer(ikind), dimension(2,2)            , intent(out)    :: bf_alignment
          type(bf_sublayer), pointer                , intent(out)    :: bf_sublayer_tested
          integer(ikind), dimension(2,2)            , intent(out)    :: bf_alignment_test
          real(rkind), dimension(:,:,:), allocatable, intent(out)    :: nodes_test
          

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          
          type(bf_sublayer), pointer :: added_sublayer

          integer(ikind), dimension(2,2) :: bf_alignment_tmp

          integer(ikind) :: i,j
          integer        :: k


          interior_nodes = reshape((/
     $         (((200*(k-1)+20*(j-1)+(i-1),i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))


          select case(test_id)
            case(1)

               bf_alignment = reshape((/
     $              align_W+4,align_N,align_E-4,align_N+1/),
     $              (/2,2/))

               bf_alignment_test = reshape((/
     $              align_W+1,align_N,align_E-1,align_N+1/),
     $              (/2,2/))
               
               allocate(nodes_test(nx,6,ne))
               nodes_test = reshape((/
     $              (((200*(k-1)+20*(align_N-3+j-1)+(i-1),i=1,nx),j=1,6),k=1,ne)/),
     $              (/nx,6,ne/))


            case(3)

               added_sublayer => bf_interface_used%mainlayer_pointers(W)%ptr%get_head_sublayer()

               bf_alignment_tmp = reshape((/
     $                  align_W-9,align_N-2,align_W,align_N-1/),
     $                  (/2,2/))

               call bf_interface_used%reallocate_sublayer(
     $              added_sublayer,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               added_sublayer%nodes = reshape((/
     $              (((200*(k-1)+20*(align_N-5+j-1)+(align_W-12+i-1),i=1,14),j=1,6),k=1,ne)/),
     $              (/14,6,ne/))

               bf_alignment = reshape((/
     $              align_W-8,align_N,align_W+9,align_N+1/),
     $              (/2,2/))

               bf_alignment_test = reshape((/
     $              align_W-9,align_N,align_W+9,align_N+1/),
     $              (/2,2/))
               
               allocate(nodes_test(23,6,ne))
               nodes_test = reshape((/
     $              (((200*(k-1)+20*(align_N-3+j-1)+(align_W-12+i-1),i=1,23),j=1,6),k=1,ne)/),
     $              (/23,6,ne/))

            case(5)

               added_sublayer => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()

               bf_alignment_tmp = reshape((/
     $                  align_E,align_N-2,align_E+9,align_N-1/),
     $                  (/2,2/))

               call bf_interface_used%reallocate_sublayer(
     $              added_sublayer,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)
               
               added_sublayer%nodes = reshape((/
     $              (((200*(k-1)+20*(align_N-5+j-1)+(align_E-3+i-1),i=1,14),j=1,6),k=1,ne)/),
     $              (/14,6,ne/))

               bf_alignment = reshape((/
     $              align_E-9,align_N,align_E+8,align_N+1/),
     $              (/2,2/))

               bf_alignment_test = reshape((/
     $              align_E-9,align_N,align_E+9,align_N+1/),
     $              (/2,2/))
               
               allocate(nodes_test(23,6,ne))
               nodes_test = reshape((/
     $              (((200*(k-1)+20*(align_N-3+j-1)+(align_E-12+i-1),i=1,23),j=1,6),k=1,ne)/),
     $              (/23,6,ne/))

          end select

          bf_sublayer_tested => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()

        end subroutine get_param_test_reallocate


        subroutine get_param_test_merge(
     $     test_id,
     $     bf_interface_used,
     $     bf_alignment,
     $     bf_sublayer1,
     $     bf_sublayer2,
     $     bf_alignment_test,
     $     nodes_test)

          implicit none

          integer                                   , intent(in)     :: test_id
          type(bf_interface_dyn)                    , intent(inout)  :: bf_interface_used
          integer(ikind), dimension(2,2)            , intent(out)    :: bf_alignment
          type(bf_sublayer), pointer                , intent(out)    :: bf_sublayer1
          type(bf_sublayer), pointer                , intent(out)    :: bf_sublayer2
          integer(ikind), dimension(2,2)            , intent(out)    :: bf_alignment_test
          real(rkind), dimension(:,:,:), allocatable, intent(out)    :: nodes_test
          

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          
          type(bf_sublayer), pointer :: added_sublayer

          integer(ikind), dimension(2,2) :: bf_alignment_tmp

          integer(ikind) :: i,j
          integer        :: k


          interior_nodes = reshape((/
     $         (((200*(k-1)+20*(j-1)+(i-1),i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))


          
          select case(test_id)
            case(1)

               bf_alignment_tmp = reshape((/
     $              align_E-10, align_N, align_E-5, align_N+1/),
     $              (/2,2/))

               added_sublayer => bf_interface_used%allocate_sublayer(
     $              N,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               bf_alignment = reshape((/
     $              align_W+4,align_N,align_E-4,align_N+1/),
     $              (/2,2/))

               bf_alignment_test = reshape((/
     $              align_W+1,align_N,align_E-1,align_N+1/),
     $              (/2,2/))
               
               allocate(nodes_test(nx,6,ne))
               nodes_test = reshape((/
     $              (((200*(k-1)+20*(align_N-3+j-1)+(i-1),i=1,nx),j=1,6),k=1,ne)/),
     $              (/nx,6,ne/))


            case(3)

               bf_alignment_tmp = reshape((/
     $              align_E-10, align_N, align_E-5, align_N+1/),
     $              (/2,2/))

               added_sublayer => bf_interface_used%allocate_sublayer(
     $              N,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               bf_alignment = reshape((/
     $              align_W-7,align_N,align_E-4,align_N+1/),
     $              (/2,2/))

               bf_alignment_test = reshape((/
     $              align_W-7,align_N,align_E-1,align_N+1/),
     $              (/2,2/))
               
               allocate(nodes_test(nx+8,6,ne))
               nodes_test = reshape((/
     $              (((200*(k-1)+20*(align_N-3+j-1)+(align_W-10+i-1),i=1,nx+8),j=1,6),k=1,ne)/),
     $              (/nx+8,6,ne/))

            case(5)

               bf_alignment_tmp = reshape((/
     $              align_W+5, align_N, align_W+10, align_N+1/),
     $              (/2,2/))

               added_sublayer => bf_interface_used%allocate_sublayer(
     $              N,
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              bf_alignment_tmp)

               bf_alignment = reshape((/
     $              align_W+4,align_N,align_E+7,align_N+1/),
     $              (/2,2/))

               bf_alignment_test = reshape((/
     $              align_W+1,align_N,align_E+7,align_N+1/),
     $              (/2,2/))
               
               allocate(nodes_test(nx+8,6,ne))
               nodes_test = reshape((/
     $              (((200*(k-1)+20*(align_N-3+j-1)+(i-1),i=1,nx+8),j=1,6),k=1,ne)/),
     $              (/nx+8,6,ne/))

          end select

          bf_sublayer1 => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
          bf_sublayer2 => bf_interface_used%mainlayer_pointers(N)%ptr%get_tail_sublayer()


        end subroutine get_param_test_merge


        function is_test_allocation_validated(
     $     bf_sublayer_tested,
     $     bf_alignment_test,
     $     bf_nodes_test,
     $     detailled)
     $     result(test_validated)

          implicit none

          type(bf_sublayer)               , intent(in) :: bf_sublayer_tested
          integer(ikind), dimension(2,2)  , intent(in) :: bf_alignment_test
          real(rkind)   , dimension(:,:,:), intent(in) :: bf_nodes_test
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated


          logical :: test_loc


          test_validated = .true.


          test_loc = is_int_matrix_validated(
     $         bf_sublayer_tested%alignment,
     $         bf_alignment_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test alignment failed'')'
          end if

          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_tested%nodes(:,1:4,:),
     $         bf_nodes_test(:,1:4,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nodes failed'')'
          end if        


        end function is_test_allocation_validated


        subroutine check_inputs()

          implicit none

          if(.not.(
     $         (nx.eq.30).and.
     $         (ny.eq.35))) then
             
             print '(''the test requires:'')'
             print '(''  - nx=30'')'
             print '(''  - ny=35'')'
             stop ''

          end if


        end subroutine check_inputs

      end program test_bf_interface_dyn
