      !------------------------------------------------------------       
      !     __N1__    N2    __N3_
      !    |     x|  |  |  |     |
      !    |______|  |__|  |_____|
      !     _   _____________   _
      ! W2 | | |             | |x| E2
      !    |_| |             | |_|
      !     _  |  interior   |  _
      ! W1 | | |             | | | E1
      !    |x| |_____________| |_|
      !     ______    __    _____
      !    |      |  |  |  |     |
      !    |______|  |_x|  |_____|
      !       S1      S2      S3
      !
      ! x: the cross indicates where the obc are undermined
      !    preventing the removal of the corresponding buffer
      !    layer
      !------------------------------------------------------------
      program test_bf_interface

        use bf_interface_class, only :
     $       bf_interface

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W,
     $       interior_pt,
     $       search_dcr

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne

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


        test_loc = test_remove_inactivated_bf_layers(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_remove_inactivated_bf_layers: '',L1)', test_loc
        print '()'


        contains


        function test_remove_inactivated_bf_layers(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          type(bf_interface)               :: bf_interface_used
          type(bf_sublayer), pointer       :: bf_sublayer_ptr
          type(pmodel_eq)                  :: p_model
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          type(bf_sublayer), pointer :: N1_sublayer_ptr
          type(bf_sublayer), pointer :: N3_sublayer_ptr
          type(bf_sublayer), pointer :: S1_sublayer_ptr
          type(bf_sublayer), pointer :: E2_sublayer_ptr
          type(bf_sublayer), pointer :: W1_sublayer_ptr
          type(bf_sublayer), pointer :: W2_sublayer_ptr

          integer(ikind) :: i,j
          logical :: remain_status


          logical :: test_loc
          


          test_validated = .true.


          !       
          !     __N1__    N2    __N3_
          !    |     x|  |  |  |     |
          !    |______|  |__|  |_____|
          !     _   _____________   _
          ! W2 | | |             | |x| E2
          !    |_| |             | |_|
          !     _  |  interior   |  _
          ! W1 | | |             | | | E1
          !    |x| |_____________| |_|
          !     ______    __    _____
          !    |      |  |  |  |     |
          !    |______|  |_x|  |_____|
          !       S1      S2      S3
          
          ! we compute the remain status of the
          ! buffer layers to checked whether the
          ! initialization was correct
          !
          ! we ask the bf_interface_used to remove 
          ! inactivated buffer layers:
          !  - S3: since S3 and E1 are inactivated
          !  - N2: since N2 is inactivated
          !============================================================
          !input
          !============================================================
          call ini_bf_interface_for_tests(bf_interface_used)

          interior_nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,nx),j=1,ny)/), (/nx,ny/))
          
          !check remain status
          !------------------------------------------------------------
          !N1
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%get_head_sublayer()
          remain_status = bf_sublayer_ptr%should_remain(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)

          if(.not.remain_status) then
             print '(''remain_status(N1) failed'')'
          end if

          !N2
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%get_head_sublayer()
          bf_sublayer_ptr => bf_sublayer_ptr%get_next()
          remain_status = bf_sublayer_ptr%should_remain(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)

          if(remain_status) then
             print '(''remain_status(N2) failed'')'
          end if

          !N3
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%get_tail_sublayer()
          remain_status = bf_sublayer_ptr%should_remain(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)

          if(remain_status) then
             print '(''remain_status(N3) failed'')'
          end if


          !S1
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%get_head_sublayer()
          remain_status = bf_sublayer_ptr%should_remain(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)

          if(remain_status) then
             print '(''remain_status(S1) failed'')'
          end if

          !S2
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%get_head_sublayer()
          bf_sublayer_ptr => bf_sublayer_ptr%get_next()
          remain_status = bf_sublayer_ptr%should_remain(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)

          if(.not.remain_status) then
             print '(''remain_status(S2) failed'')'
          end if

          !S3
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%get_tail_sublayer()
          remain_status = bf_sublayer_ptr%should_remain(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)

          if(remain_status) then
             print '(''remain_status(S3) failed'')'
          end if          

          !E1
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%get_head_sublayer()
          remain_status = bf_sublayer_ptr%should_remain(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)

          if(remain_status) then
             print '(''remain_status(E1) failed'')'
          end if

          !E2
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%get_tail_sublayer()
          remain_status = bf_sublayer_ptr%should_remain(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)

          if(.not.remain_status) then
             print '(''remain_status(E2) failed'')'
          end if

          !W1
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%get_head_sublayer()
          remain_status = bf_sublayer_ptr%should_remain(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)

          if(.not.remain_status) then
             print '(''remain_status(W1) failed'')'
          end if

          !W2
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%get_tail_sublayer()
          remain_status = bf_sublayer_ptr%should_remain(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)

          if(remain_status) then
             print '(''remain_status(W2) failed'')'
          end if


          !============================================================
          !output
          !============================================================
          call bf_interface_used%remove_inactivated_bf_layers(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)


          !============================================================
          !validation
          !============================================================

          !check the number of sublayers in the main layers
          !------------------------------------------------------------
          test_loc = bf_interface_used%mainlayer_pointers(N)%get_nb_sublayers().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_sublayers N failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(S)%get_nb_sublayers().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_sublayers S failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(E)%get_nb_sublayers().eq.1
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_sublayers E failed'')'
          end if

          test_loc = bf_interface_used%mainlayer_pointers(W)%get_nb_sublayers().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_sublayers W failed'')'
          end if


          !check the alignment of the remaining buffer layers
          !------------------------------------------------------------
          !N1
          !............................................................
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%get_alignment_tab(),
     $         reshape((/
     $         align_W+1,align_N,align_W+5,align_N+1/),
     $         (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test N1 failed'')'
          end if

          !N3
          !............................................................
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%get_tail_sublayer()
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%get_alignment_tab(),
     $         reshape((/
     $         align_E-5,align_N,align_E-1,align_N+1/),
     $         (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test N3 failed'')'
          end if

          !S1
          !............................................................
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%get_alignment_tab(),
     $         reshape((/
     $         align_W+1,align_S-1,align_W+5,align_S/),
     $         (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test S1 failed'')'
          end if

          !S2
          !............................................................
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%get_tail_sublayer()
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%get_alignment_tab(),
     $         reshape((/
     $         align_W+12,align_S-1,align_W+12,align_S/),
     $         (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test S2 failed'')'
          end if

          !W1
          !............................................................
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%get_alignment_tab(),
     $         reshape((/
     $         align_W-4,align_S+1,align_W,align_S+6/),
     $         (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test W1 failed'')'
          end if

          !W2
          !............................................................
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%get_tail_sublayer()
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%get_alignment_tab(),
     $         reshape((/
     $         align_W-4,align_N-6,align_W,align_N-1/),
     $         (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test W2 failed'')'
          end if

          !E2
          !............................................................
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%get_alignment_tab(),
     $         reshape((/
     $         align_E,align_N-6,align_E+4,align_N-1/),
     $         (/2,2/)),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test E2 failed'')'
          end if

          ! check the references to the buffer layers at the interface
          ! between main layers
          !------------------------------------------------------------
          N1_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%get_head_sublayer()
          N3_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%get_tail_sublayer()

          S1_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%get_head_sublayer()

          W1_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%get_head_sublayer()
          W2_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%get_tail_sublayer()

          E2_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%get_tail_sublayer()
          
          !NE_interface
          !............................................................
          test_loc = associated(bf_interface_used%mainlayer_interfaces%NE_interface_N_ptr,
     $                          N3_sublayer_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NE_interface_N_ptr failed'')'
          end if

          test_loc = associated(bf_interface_used%mainlayer_interfaces%NE_interface_E_ptr,
     $                          E2_sublayer_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NE_interface_E_ptr failed'')'
          end if

          !NW_interface
          !............................................................
          test_loc = associated(bf_interface_used%mainlayer_interfaces%NW_interface_N_ptr,
     $                          N1_sublayer_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NW_interface_N_ptr failed'')'
          end if

          test_loc = associated(bf_interface_used%mainlayer_interfaces%NW_interface_W_ptr,
     $                          W2_sublayer_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NW_interface_W_ptr failed'')'
          end if

          !SE_interface
          !............................................................
          test_loc = .not.associated(bf_interface_used%mainlayer_interfaces%SE_interface_S_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SE_interface_S_ptr failed'')'
          end if

          test_loc = .not.associated(bf_interface_used%mainlayer_interfaces%SE_interface_E_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SE_interface_E_ptr failed'')'
          end if

          !SW_interface
          !............................................................
          test_loc = associated(bf_interface_used%mainlayer_interfaces%SW_interface_S_ptr,
     $                          S1_sublayer_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SW_interface_S_ptr failed'')'
          end if

          test_loc = associated(bf_interface_used%mainlayer_interfaces%SW_interface_W_ptr,
     $                          W1_sublayer_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SW_interface_W_ptr failed'')'
          end if


        end function test_remove_inactivated_bf_layers


        subroutine ini_bf_interface_for_tests(bf_interface_used)

          implicit none

          type(bf_interface), intent(inout) :: bf_interface_used

          real(rkind)   , dimension(nx)       :: interior_x_map
          real(rkind)   , dimension(ny)       :: interior_y_map
          real(rkind)   , dimension(nx,ny,ne) :: interior_nodes
          integer(ikind), dimension(2,2)      :: bf_alignment


          type(bf_sublayer), pointer :: added_sublayer
          integer(ikind) :: i,j


          interior_nodes(:,:,1) = reshape(
     $         (/((1.0d0, i=1,nx),j=1,ny)/),
     $         (/nx,ny/))


          !first North buffer layer: cannot be removed
          !------------------------------------------------------------
          bf_alignment = reshape((/
     $         align_W+1,align_N,align_W+5,align_N+1/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          added_sublayer%grdpts_id = reshape(
     $         (/((interior_pt, i=1,9), j=1,6)/),
     $         (/9,6/))

          added_sublayer%nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,9),j=1,6)/),(/9,6/))
          added_sublayer%nodes(8,5,1) = -1.0d0



          !second North buffer layer: can be removed
          !------------------------------------------------------------
          bf_alignment = reshape((/
     $         align_W+12,align_N,align_W+12,align_N+1/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          added_sublayer%grdpts_id = reshape(
     $         (/((interior_pt, i=1,6), j=1,6)/),
     $         (/6,6/))

          added_sublayer%nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,6),j=1,6)/),(/6,6/))


          !third North buffer layer: cannot be removed b/c of E neighbor
          !------------------------------------------------------------
          bf_alignment = reshape((/
     $         align_E-5,align_N,align_E-1,align_N+1/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          added_sublayer%grdpts_id = reshape(
     $         (/((interior_pt, i=1,6), j=1,6)/),
     $         (/6,6/))

          added_sublayer%nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,6),j=1,6)/),(/6,6/))


          !first South buffer layer: cannot be removed b/c of W neighbor
          !------------------------------------------------------------
          bf_alignment = reshape((/
     $         align_W+1,align_S-1,align_W+5,align_S/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         S,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          added_sublayer%grdpts_id = reshape(
     $         (/((interior_pt, i=1,9), j=1,6)/),
     $         (/9,6/))

          added_sublayer%nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,9),j=1,6)/),(/9,6/))


          !second South buffer layer: cannot be removed
          !------------------------------------------------------------
          bf_alignment = reshape((/
     $         align_W+12,align_S-1,align_W+12,align_S/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         S,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          added_sublayer%grdpts_id = reshape(
     $         (/((interior_pt, i=1,6), j=1,6)/),
     $         (/6,6/))

          added_sublayer%nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,6),j=1,6)/),(/6,6/))
          added_sublayer%nodes(3,2,1) = -1.0d0


          !third South buffer layer: can be removed
          !------------------------------------------------------------
          bf_alignment = reshape((/
     $         align_E-5,align_S-1,align_E-1,align_S/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         S,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          added_sublayer%grdpts_id = reshape(
     $         (/((interior_pt, i=1,6), j=1,6)/),
     $         (/6,6/))

          added_sublayer%nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,6),j=1,6)/),(/6,6/))


          !first West buffer layer: cannot be removed
          !------------------------------------------------------------
          bf_alignment = reshape((/
     $         align_W-4,align_S+1,align_W,align_S+6/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         W,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          added_sublayer%grdpts_id = reshape(
     $         (/((interior_pt, i=1,9), j=1,10)/),
     $         (/9,10/))

          added_sublayer%nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,9),j=1,10)/),(/9,10/))
          added_sublayer%nodes(4,3,1) = -1.0d0


          !second West buffer layer: cannot be removed b/c of N neighbor
          !------------------------------------------------------------
          bf_alignment = reshape((/
     $         align_W-4,align_N-6,align_W,align_N-1/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         W,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          added_sublayer%grdpts_id = reshape(
     $         (/((interior_pt, i=1,9), j=1,10)/),
     $         (/9,10/))

          added_sublayer%nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,9),j=1,10)/),(/9,10/))
          added_sublayer%nodes(2,3,1) = -1.0d0


          !first East buffer layer: can be removed
          !------------------------------------------------------------
          bf_alignment = reshape((/
     $         align_E,align_S+1,align_E+4,align_S+6/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         E,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          added_sublayer%grdpts_id = reshape(
     $         (/((interior_pt, i=1,9), j=1,10)/),
     $         (/9,10/))

          added_sublayer%nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,9),j=1,10)/),(/9,10/))


          !second East buffer layer: cannot be removed
          !------------------------------------------------------------
          bf_alignment = reshape((/
     $         align_E,align_N-6,align_E+4,align_N-1/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         E,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment)

          added_sublayer%grdpts_id = reshape(
     $         (/((interior_pt, i=1,9), j=1,10)/),
     $         (/9,10/))

          added_sublayer%nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,9),j=1,10)/),(/9,10/))
          added_sublayer%nodes(5,3,1) = -1.0d0

        end subroutine ini_bf_interface_for_tests


        subroutine check_inputs()

          implicit none

          type(pmodel_eq)                :: p_model
          real(rkind), dimension(3)      :: x_map
          real(rkind), dimension(3)      :: y_map
          real(rkind), dimension(3,3,ne) :: nodes
          
          logical :: test_loc

          if(search_dcr.ne.4) then

             print '(''the test requires: '')'
             print '(''search_dcr=4'')'
             stop ''
             
          end if

          nodes(2,2,1) = -1.0d0
          test_loc = p_model%are_openbc_undermined(x_map,y_map,nodes)
          if(.not.test_loc) then
             print '(''the test requires: '')'
             print '(''openbc_undermined if nodes(2,2,1)<0'')'
             stop ''
          end if

          
          nodes(2,2,1) = 1.0d0
          test_loc = .not.p_model%are_openbc_undermined(x_map,y_map,nodes)
          if(.not.test_loc) then
             print '(''the test requires: '')'
             print '(''.not.openbc_undermined if nodes(2,2,1)>0'')'
             stop ''
          end if

          if(((align_W+14).gt.(align_E-9)).or.((align_S+8).gt.(align_N-10))) then
              print '(''the test requires: '')'
              print '(''   - alignW+14<align_E-8'')'
              print '(''   - align_S+8<align_N-9'')'
              print '()'
              stop ''
           end if

        end subroutine check_inputs

      end program test_bf_interface
