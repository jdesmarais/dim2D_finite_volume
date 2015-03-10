      program test_bf_interface_sync

        use bf_interface_sync_class, only :
     $       bf_interface_sync

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_real_matrix3D_validated

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W,
     $       interior_pt,
     $       NW_interface_type,
     $       NE_interface_type,
     $       SW_interface_type,
     $       SE_interface_type

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

        
        detailled = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_sync_nodes(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_sync_nodes: '',L1)', test_loc
        print '()'


        contains

        function test_sync_nodes(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_interface_sync)          :: bf_interface_used
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          real(rkind), dimension(nx,ny,ne) :: interior_nodes_test

          real(rkind), dimension(12,6,ne)  :: N_1_nodes_test
          real(rkind), dimension(12,6,ne)  :: N_2_nodes_test
          real(rkind), dimension(12,6,ne)  :: S_1_nodes_test
          real(rkind), dimension(12,6,ne)  :: S_2_nodes_test
          real(rkind), dimension( 9,6,ne)  :: E_1_nodes_test
          real(rkind), dimension( 9,6,ne)  :: E_2_nodes_test
          real(rkind), dimension( 9,6,ne)  :: W_1_nodes_test
          real(rkind), dimension( 9,6,ne)  :: W_2_nodes_test

          real(rkind), dimension(12,6,ne)  :: N_1_nodes
          real(rkind), dimension(12,6,ne)  :: N_2_nodes
          real(rkind), dimension(12,6,ne)  :: S_1_nodes
          real(rkind), dimension(12,6,ne)  :: S_2_nodes
          real(rkind), dimension( 9,6,ne)  :: E_1_nodes
          real(rkind), dimension( 9,6,ne)  :: E_2_nodes
          real(rkind), dimension( 9,6,ne)  :: W_1_nodes
          real(rkind), dimension( 9,6,ne)  :: W_2_nodes

          type(bf_sublayer), pointer :: current_sublayer

          integer(ikind) :: i,j
          integer        :: k


          test_validated = .true.


          !input
          !interior nodes initialization
          !------------------------------------------------------------
          interior_nodes_test = reshape((/
     $         (((200*(k-1) + 20*(j-1) + (i-1),i=1,nx),j=1,ny),k=1,ne)/),
     $         (/nx,ny,ne/))

          interior_nodes = interior_nodes_test
          
          interior_nodes(1:7,1:2,:) = reshape((/
     $         (((-99.0d0, i=1,7),j=1,2),k=1,ne) /),
     $         (/7,2,ne/))

          interior_nodes(14:20,1:2,:) = reshape((/
     $         (((-99.0d0, i=1,7),j=1,2),k=1,ne) /),
     $         (/7,2,ne/))

          interior_nodes(1:7,24:25,:) = reshape((/
     $         (((-99.0d0, i=1,7),j=1,2),k=1,ne) /),
     $         (/7,2,ne/))

          interior_nodes(14:20,24:25,:) = reshape((/
     $         (((-99.0d0, i=1,7),j=1,2),k=1,ne) /),
     $         (/7,2,ne/))

          interior_nodes(1:2,1:6,:) = reshape((/
     $         (((-99.0d0, i=1,2),j=1,6),k=1,ne) /),
     $         (/2,6,ne/))

          interior_nodes(1:2,20:25,:) = reshape((/
     $         (((-99.0d0, i=1,2),j=1,6),k=1,ne) /),
     $         (/2,6,ne/))

          interior_nodes(19:20,1:6,:) = reshape((/
     $         (((-99.0d0, i=1,2),j=1,6),k=1,ne) /),
     $         (/2,6,ne/))

          interior_nodes(19:20,20:25,:) = reshape((/
     $         (((-99.0d0, i=1,2),j=1,6),k=1,ne) /),
     $         (/2,6,ne/))
          

          !initialization of the buffer layers nodes
          !------------------------------------------------------------
          N_1_nodes_test = reshape((/
     $         (((200*(k-1) + 20*(align_N-3+j-1) + (align_W-7+i-1),i=1,12),j=1,6),k=1,ne)/),
     $         (/12,6,ne/))

          N_2_nodes_test = reshape((/
     $         (((200*(k-1) + 20*(align_N-3+j-1) + (align_E-6+i-1),i=1,12),j=1,6),k=1,ne)/),
     $         (/12,6,ne/))

          S_1_nodes_test = reshape((/
     $         (((200*(k-1) + 20*(align_S-4+j-1) + (align_W-7+i-1),i=1,12),j=1,6),k=1,ne)/),
     $         (/12,6,ne/))

          S_2_nodes_test = reshape((/
     $         (((200*(k-1) + 20*(align_S-4+j-1) + (align_E-6+i-1),i=1,12),j=1,6),k=1,ne)/),
     $         (/12,6,ne/))

          E_1_nodes_test = reshape((/
     $         (((200*(k-1) + 20*(align_S-2+j-1) + (align_E-3+i-1),i=1,9),j=1,6),k=1,ne)/),
     $         (/9,6,ne/))

          E_2_nodes_test = reshape((/
     $         (((200*(k-1) + 20*(align_N-5+j-1) + (align_E-3+i-1),i=1,9),j=1,6),k=1,ne)/),
     $         (/9,6,ne/))

          W_1_nodes_test = reshape((/
     $         (((200*(k-1) + 20*(align_S-2+j-1) + (align_W-7+i-1),i=1,9),j=1,6),k=1,ne)/),
     $         (/9,6,ne/))

          W_2_nodes_test = reshape((/
     $         (((200*(k-1) + 20*(align_N-5+j-1) + (align_W-7+i-1),i=1,9),j=1,6),k=1,ne)/),
     $         (/9,6,ne/))

          N_1_nodes = N_1_nodes_test
          N_2_nodes = N_2_nodes_test
          S_1_nodes = S_1_nodes_test
          S_2_nodes = S_2_nodes_test

          E_1_nodes = E_1_nodes_test
          E_2_nodes = E_2_nodes_test
          W_1_nodes = W_1_nodes_test
          W_2_nodes = W_2_nodes_test


          N_1_nodes(:,1:2,:) = reshape((/
     $         (((-99.0d0,i=1,12),j=1,2),k=1,ne)/),(/12,2,ne/))

          N_2_nodes(:,1:2,:) = reshape((/
     $         (((-99.0d0,i=1,12),j=1,2),k=1,ne)/),(/12,2,ne/))

          S_1_nodes(:,5:6,:) = reshape((/
     $         (((-99.0d0,i=1,12),j=1,2),k=1,ne)/),(/12,2,ne/))

          S_2_nodes(:,5:6,:) = reshape((/
     $         (((-99.0d0,i=1,12),j=1,2),k=1,ne)/),(/12,2,ne/))

          E_1_nodes(:,1:2,:) = reshape((/
     $         (((-99.0d0,i=1,9),j=1,2),k=1,ne)/),(/9,2,ne/))
          
          E_1_nodes(1:2,:,:) = reshape((/
     $         (((-99.0d0,i=1,2),j=1,6),k=1,ne)/),(/2,6,ne/))

          E_2_nodes(:,5:6,:) = reshape((/
     $         (((-99.0d0,i=1,9),j=1,2),k=1,ne)/),(/9,2,ne/))

          E_2_nodes(1:2,:,:) = reshape((/
     $         (((-99.0d0,i=1,2),j=1,6),k=1,ne)/),(/2,6,ne/))

          W_1_nodes(:,1:2,:) = reshape((/
     $         (((-99.0d0,i=1,9),j=1,2),k=1,ne)/),(/9,2,ne/))
          
          W_1_nodes(8:9,:,:) = reshape((/
     $         (((-99.0d0,i=1,2),j=1,6),k=1,ne)/),(/2,6,ne/))

          W_2_nodes(:,5:6,:) = reshape((/
     $         (((-99.0d0,i=1,9),j=1,2),k=1,ne)/),(/9,2,ne/))

          W_2_nodes(8:9,:,:) = reshape((/
     $         (((-99.0d0,i=1,2),j=1,6),k=1,ne)/),(/2,6,ne/))


          !set the buffer layers inside the bf_interface
          !------------------------------------------------------------
          call bf_interface_used%ini(interior_x_map,interior_y_map)


          !two North layers
          call bf_interface_used%mainlayer_pointers(N)%ini_mainlayer(N)
          current_sublayer => bf_interface_used%mainlayer_pointers(N)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_W-4, align_N, align_W+3, align_N+1/),
     $            (/2,2/)))
          current_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,12),j=1,6)/),(/12,6/))
          current_sublayer%nodes = N_1_nodes

          call bf_interface_used%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $         NW_interface_type, current_sublayer)

          current_sublayer => bf_interface_used%mainlayer_pointers(N)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_E-3, align_N, align_E+4, align_N+1/),
     $            (/2,2/)))
          current_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,12),j=1,6)/),(/12,6/))
          current_sublayer%nodes = N_2_nodes

          call bf_interface_used%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $         NE_interface_type, current_sublayer)


          !two south layers
          call bf_interface_used%mainlayer_pointers(S)%ini_mainlayer(S)
          current_sublayer => bf_interface_used%mainlayer_pointers(S)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_W-4, align_S-1, align_W+3, align_S/),
     $            (/2,2/)))
          current_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,12),j=1,6)/),(/12,6/))
          current_sublayer%nodes = S_1_nodes

          call bf_interface_used%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $         SW_interface_type, current_sublayer)


          current_sublayer => bf_interface_used%mainlayer_pointers(S)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_E-3, align_S-1, align_E+4, align_S/),
     $            (/2,2/)))
          current_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,12),j=1,6)/),(/12,6/))
          current_sublayer%nodes = S_2_nodes

          call bf_interface_used%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $         SE_interface_type, current_sublayer)


          !two east layers
          call bf_interface_used%mainlayer_pointers(E)%ini_mainlayer(E)
          current_sublayer => bf_interface_used%mainlayer_pointers(E)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_E, align_S+1, align_E+4, align_S+2/),
     $            (/2,2/)))
          current_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,9),j=1,6)/),(/9,6/))
          current_sublayer%nodes = E_1_nodes

          call bf_interface_used%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $         SE_interface_type, current_sublayer)


          current_sublayer => bf_interface_used%mainlayer_pointers(E)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_E, align_N-2, align_E+4, align_N-1/),
     $            (/2,2/)))
          current_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,9),j=1,6)/),(/9,6/))
          current_sublayer%nodes = E_2_nodes

          call bf_interface_used%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $         NE_interface_type, current_sublayer)


          !two west layers
          call bf_interface_used%mainlayer_pointers(W)%ini_mainlayer(W)
          current_sublayer => bf_interface_used%mainlayer_pointers(W)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_W-4, align_S+1, align_W, align_S+2/),
     $            (/2,2/)))
          current_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,9),j=1,6)/),(/9,6/))
          current_sublayer%nodes = W_1_nodes
          
          call bf_interface_used%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $         SW_interface_type, current_sublayer)


          current_sublayer => bf_interface_used%mainlayer_pointers(W)%ptr%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $            align_W-4, align_N-2, align_W, align_N-1/),
     $            (/2,2/)))
          current_sublayer%grdpts_id = reshape((/
     $         ((interior_pt,i=1,9),j=1,6)/),(/9,6/))
          current_sublayer%nodes = W_2_nodes

          call bf_interface_used%mainlayer_interfaces%set_mainlayer_interface_bf_layer(
     $         NW_interface_type, current_sublayer)


          !output
          call bf_interface_used%sync_nodes(interior_nodes)


          !validation
          !interior_nodes validation
          test_loc = is_real_matrix3D_validated(
     $         interior_nodes,
     $         interior_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''interior_nodes synchrnization failed'')'
          end if

          !north layers validation
          current_sublayer => bf_interface_used%mainlayer_interfaces%NW_interface_N_ptr
          test_loc = is_real_matrix3D_validated(
     $         current_sublayer%nodes,
     $         N_1_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''N_1_nodes synchronization failed'')'
          end if

          current_sublayer => bf_interface_used%mainlayer_interfaces%NE_interface_N_ptr
          test_loc = is_real_matrix3D_validated(
     $         current_sublayer%nodes,
     $         N_2_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''N_2_nodes synchronization failed'')'
          end if

          !south layers validation
          current_sublayer => bf_interface_used%mainlayer_interfaces%SW_interface_S_ptr
          test_loc = is_real_matrix3D_validated(
     $         current_sublayer%nodes,
     $         S_1_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''S_1_nodes synchronization failed'')'
          end if

          current_sublayer => bf_interface_used%mainlayer_interfaces%SE_interface_S_ptr
          test_loc = is_real_matrix3D_validated(
     $         current_sublayer%nodes,
     $         S_2_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''S_2_nodes synchronization failed'')'
          end if

          !east layers validation
          current_sublayer => bf_interface_used%mainlayer_interfaces%SE_interface_E_ptr
          test_loc = is_real_matrix3D_validated(
     $         current_sublayer%nodes,
     $         E_1_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''E_1_nodes synchronization failed'')'
          end if

          current_sublayer => bf_interface_used%mainlayer_interfaces%NE_interface_E_ptr
          test_loc = is_real_matrix3D_validated(
     $         current_sublayer%nodes,
     $         E_2_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''E_2_nodes synchrnization failed'')'
          end if

          !west layers validation
          current_sublayer => bf_interface_used%mainlayer_interfaces%SW_interface_W_ptr
          test_loc = is_real_matrix3D_validated(
     $         current_sublayer%nodes,
     $         W_1_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''W_1_nodes synchronization failed'')'
          end if

          current_sublayer => bf_interface_used%mainlayer_interfaces%NW_interface_W_ptr
          test_loc = is_real_matrix3D_validated(
     $         current_sublayer%nodes,
     $         W_2_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''W_2_nodes synchrnization failed'')'
          end if

        end function test_sync_nodes


        subroutine check_inputs()

          if(.not.(
     $       (nx.eq.20).and.
     $       (ny.eq.25))) then
             
             print '(''the test requires:'')'
             print '(''nx=20'')'
             print '(''ny=25'')'
             stop ''

          end if

        end subroutine check_inputs

      end program test_bf_interface_sync
