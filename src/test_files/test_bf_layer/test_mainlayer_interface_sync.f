      program test_mainlayer_interface_sync

        use bf_layer_errors_module, only :
     $     error_mainlayer_interface_type,
     $     error_mainlayer_interface_incompatible

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_real_matrix3D_validated

        use mainlayer_interface_sync_class, only :
     $     mainlayer_interface_sync

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W,
     $       
     $       NW_interface_type, NE_interface_type,
     $       SW_interface_type, SE_interface_type

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


        test_loc = test_sync_nodes_at_mainlayer_interfaces(detailled)
        test_validated = test_validated.and.test_loc
        print '(''sync_nodes_at_mainlayer_interfaces: '',L1)', test_loc
        print '()'


        contains


        function test_sync_nodes_at_mainlayer_interfaces(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mainlayer_interface_sync)  :: mainlayer_interface_used
          type(bf_sublayer), pointer      :: bf_sublayer_added

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          real(rkind), dimension(10,6,ne) :: NE_interface_N_nodes
          real(rkind), dimension(10,6,ne) :: NW_interface_N_nodes
          real(rkind), dimension(10,6,ne) :: SE_interface_S_nodes
          real(rkind), dimension(10,6,ne) :: SW_interface_S_nodes

          real(rkind), dimension(8,6,ne)  :: NE_interface_E_nodes
          real(rkind), dimension(8,6,ne)  :: NW_interface_W_nodes
          real(rkind), dimension(8,6,ne)  :: SE_interface_E_nodes
          real(rkind), dimension(8,6,ne)  :: SW_interface_W_nodes

          real(rkind), dimension(10,6,ne) :: NE_interface_N_nodes_test
          real(rkind), dimension(10,6,ne) :: NW_interface_N_nodes_test
          real(rkind), dimension(10,6,ne) :: SE_interface_S_nodes_test
          real(rkind), dimension(10,6,ne) :: SW_interface_S_nodes_test

          real(rkind), dimension(8,6,ne)  :: NE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne)  :: NW_interface_W_nodes_test
          real(rkind), dimension(8,6,ne)  :: SE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne)  :: SW_interface_W_nodes_test

          integer(ikind) :: i,j
          integer        :: k


          test_validated = .true.


          !input
          !------------------------------------------------------------
          call mainlayer_interface_used%ini()

          !North and South layers
          NW_interface_N_nodes_test = reshape((/
     $         (((100*(k-1)+10*(align_N-3+j-1)+(align_W-6+i-1),i=1,10),j=1,6),k=1,ne)/),
     $         (/10,6,ne/))

          NE_interface_N_nodes_test = reshape((/
     $         (((100*(k-1)+10*(align_N-3+j-1)+(align_E-5+i-1),i=1,10),j=1,6),k=1,ne)/),
     $         (/10,6,ne/))

          SW_interface_S_nodes_test = reshape((/
     $         (((100*(k-1)+10*(align_S-4+j-1)+(align_W-6+i-1),i=1,10),j=1,6),k=1,ne)/),
     $         (/10,6,ne/))

          SE_interface_S_nodes_test = reshape((/
     $         (((100*(k-1)+10*(align_S-4+j-1)+(align_E-5+i-1),i=1,10),j=1,6),k=1,ne)/),
     $         (/10,6,ne/))

          !East and West layers
          NW_interface_W_nodes_test = reshape((/
     $         (((100*(k-1)+10*(align_N-5+j-1)+(align_W-6+i-1),i=1,8),j=1,6),k=1,ne)/),
     $         (/8,6,ne/))

          NE_interface_E_nodes_test = reshape((/
     $         (((100*(k-1)+10*(align_N-5+j-1)+(align_E-3+i-1),i=1,8),j=1,6),k=1,ne)/),
     $         (/8,6,ne/))

          SW_interface_W_nodes_test = reshape((/
     $         (((100*(k-1)+10*(align_S-2+j-1)+(align_W-6+i-1),i=1,8),j=1,6),k=1,ne)/),
     $         (/8,6,ne/))

          SE_interface_E_nodes_test = reshape((/
     $         (((100*(k-1)+10*(align_S-2+j-1)+(align_E-3+i-1),i=1,8),j=1,6),k=1,ne)/),
     $         (/8,6,ne/))


          !removing the layers exchanged for the test
          NW_interface_N_nodes = NW_interface_N_nodes_test
          NW_interface_W_nodes = NW_interface_W_nodes_test
          NE_interface_N_nodes = NE_interface_N_nodes_test
          NE_interface_E_nodes = NE_interface_E_nodes_test
          SW_interface_S_nodes = SW_interface_S_nodes_test
          SW_interface_W_nodes = SW_interface_W_nodes_test
          SE_interface_S_nodes = SE_interface_S_nodes_test
          SE_interface_E_nodes = SE_interface_E_nodes_test

          
          !North layers
          NW_interface_N_nodes(:,1:2,:) = reshape((/
     $         (((-99.0d0,i=1,8),j=1,2),k=1,ne)/),
     $         (/10,2,ne/))

          NE_interface_N_nodes(:,1:2,:) = reshape((/
     $         (((-99.0d0,i=1,8),j=1,2),k=1,ne)/),
     $         (/10,2,ne/))
          
          !East and West layers with North
          NW_interface_W_nodes(:,5:6,:)
     $         = reshape((/
     $         (((-99.0d0,i=1,8),j=5,6),k=1,ne)/),
     $         (/8,2,ne/))

          NE_interface_E_nodes(:,5:6,:)
     $         = reshape((/
     $         (((-99.0d0,i=1,8),j=5,6),k=1,ne)/),
     $         (/8,2,ne/))

          !South layers
          SW_interface_S_nodes(:,5:6,:)
     $         = reshape((/
     $         (((-99.0d0,i=1,10),j=5,6),k=1,ne)/),
     $         (/10,2,ne/))

          SE_interface_S_nodes(:,5:6,:)
     $         = reshape((/
     $         (((-99.0d0,i=1,10),j=5,6),k=1,ne)/),
     $         (/10,2,ne/))

          !East and West layers with S
          SW_interface_W_nodes(:,1:2,:)
     $         = reshape((/
     $         (((-99.0d0,i=1,8),j=1,2),k=1,ne)/),
     $         (/8,2,ne/))

          SE_interface_E_nodes(:,1:2,:)
     $         = reshape((/
     $         (((-99.0d0,i=1,8),j=1,2),k=1,ne)/),
     $         (/8,2,ne/))


          !NW_interface(N)
          allocate(bf_sublayer_added)
          call bf_sublayer_added%ini(N)
          call bf_sublayer_added%allocate_bf_layer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape(
     $           (/align_W-3, align_N, align_W+2, align_N+1/),
     $           (/2,2/))
     $         )
          bf_sublayer_added%nodes = NW_interface_N_nodes
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         NW_interface_type,
     $         bf_sublayer_added)


          !NW_interface(W)
          allocate(bf_sublayer_added)
          call bf_sublayer_added%ini(W)
          call bf_sublayer_added%allocate_bf_layer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape(
     $           (/align_W-3, align_N-2, align_W, align_N-1/),
     $           (/2,2/))
     $         )
          bf_sublayer_added%nodes = NW_interface_W_nodes
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         NW_interface_type,
     $         bf_sublayer_added)


          !NE_interface(N)
          allocate(bf_sublayer_added)
          call bf_sublayer_added%ini(N)
          call bf_sublayer_added%allocate_bf_layer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape(
     $           (/align_E-2, align_N, align_E+3, align_N+1/),
     $           (/2,2/))
     $         )
          bf_sublayer_added%nodes = NE_interface_N_nodes
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         NE_interface_type,
     $         bf_sublayer_added)

          
          !NE_interface(E)
          allocate(bf_sublayer_added)
          call bf_sublayer_added%ini(E)
          call bf_sublayer_added%allocate_bf_layer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape(
     $           (/align_E, align_N-2, align_E+3, align_N-1/),
     $           (/2,2/))
     $         )
          bf_sublayer_added%nodes = NE_interface_E_nodes
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         NE_interface_type,
     $         bf_sublayer_added)


          !SW_interface(S)
          allocate(bf_sublayer_added)
          call bf_sublayer_added%ini(S)
          call bf_sublayer_added%allocate_bf_layer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape(
     $           (/align_W-3, align_S-1, align_W+2, align_S/),
     $           (/2,2/))
     $         )
          bf_sublayer_added%nodes = SW_interface_S_nodes
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SW_interface_type,
     $         bf_sublayer_added)


          !SW_interface(W)
          allocate(bf_sublayer_added)
          call bf_sublayer_added%ini(W)
          call bf_sublayer_added%allocate_bf_layer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape(
     $           (/align_W-3, align_S+1, align_W, align_S+2/),
     $           (/2,2/))
     $         )
          bf_sublayer_added%nodes = SW_interface_W_nodes
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SW_interface_type,
     $         bf_sublayer_added)


          !SE_interface(S)
          allocate(bf_sublayer_added)
          call bf_sublayer_added%ini(S)
          call bf_sublayer_added%allocate_bf_layer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape(
     $           (/align_E-2, align_S-1, align_E+3, align_S/),
     $           (/2,2/))
     $         )
          bf_sublayer_added%nodes = SE_interface_S_nodes
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         bf_sublayer_added)


          !SE_interface(E)
          allocate(bf_sublayer_added)
          call bf_sublayer_added%ini(E)
          call bf_sublayer_added%allocate_bf_layer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape(
     $           (/align_E, align_S+1, align_E+3, align_S+2/),
     $           (/2,2/))
     $         )
          bf_sublayer_added%nodes = SE_interface_E_nodes
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         bf_sublayer_added)


          !output
          !------------------------------------------------------------
          call mainlayer_interface_used%sync_nodes_at_mainlayer_interfaces()


          !validation
          !------------------------------------------------------------
          !NW_interface(N)
          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%NW_interface_N_ptr%nodes(1:8,:,:),
     $         NW_interface_N_nodes_test(1:8,:,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NW_interface(N) failed'')'
          end if


          !NW_interface(W)
          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%NW_interface_W_ptr%nodes(1:8,:,:),
     $         NW_interface_W_nodes_test(1:8,:,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NW_interface(W) failed'')'
          end if


          !NE_interface(N)
          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%NE_interface_N_ptr%nodes(3:10,:,:),
     $         NE_interface_N_nodes_test(3:10,:,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NE_interface(N) failed'')'
          end if


          !NE_interface(E)
          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%NE_interface_E_ptr%nodes(1:8,:,:),
     $         NE_interface_E_nodes_test(1:8,:,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NE_interface(E) failed'')'
          end if


          !SW_interface(S)
          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%SW_interface_S_ptr%nodes(1:8,:,:),
     $         SW_interface_S_nodes_test(1:8,:,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SW_interface(S) failed'')'
          end if


          !SW_interface(W)
          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%SW_interface_W_ptr%nodes(1:8,:,:),
     $         SW_interface_W_nodes_test(1:8,:,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SW_interface(W) failed'')'
          end if


          !SE_interface(S)
          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%SE_interface_S_ptr%nodes(3:10,:,:),
     $         SE_interface_S_nodes_test(3:10,:,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SE_interface(S) failed'')'
          end if


          !SE_interface(E)
          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%SE_interface_E_ptr%nodes,
     $         SE_interface_E_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SE_interface(E) failed'')'
          end if


        end function test_sync_nodes_at_mainlayer_interfaces

      end program test_mainlayer_interface_sync
