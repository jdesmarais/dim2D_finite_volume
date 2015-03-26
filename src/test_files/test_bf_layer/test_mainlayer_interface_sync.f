      program test_mainlayer_interface_sync

        use bf_layer_errors_module, only :
     $     error_mainlayer_interface_type,
     $     error_mainlayer_interface_incompatible

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_int_matrix_validated,
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


        call check_inputs()


        test_loc = test_sync_nodes_at_mainlayer_interfaces(detailled)
        test_validated = test_validated.and.test_loc
        print '(''sync_nodes_at_mainlayer_interfaces: '',L1)', test_loc
        print '()'


        test_loc = test_copy_from_neighbors(detailled)
        test_validated = test_validated.and.test_loc
        print '(''copy_from_mainlayer_neighbors: '',L1)', test_loc
        print '()'


        test_loc = test_copy_to_neighbors(detailled)
        test_validated = test_validated.and.test_loc
        print '(''copy_to_mainlayer_neighbors: '',L1)', test_loc
        print '()'
        

        print '(''test_validated: '',L1)', test_validated

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


        function test_copy_from_neighbors(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          type(mainlayer_interface_sync)  :: mainlayer_interface_used

          !nodes for the tests
          real(rkind), dimension(10,6,ne) :: NE_interface_N_nodes_test
          real(rkind), dimension(10,6,ne) :: NW_interface_N_nodes_test
          real(rkind), dimension(10,6,ne) :: SE_interface_S_nodes_test
          real(rkind), dimension(10,6,ne) :: SW_interface_S_nodes_test

          real(rkind), dimension(8,6,ne)  :: NE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne)  :: NW_interface_W_nodes_test
          real(rkind), dimension(8,6,ne)  :: SE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne)  :: SW_interface_W_nodes_test

          integer, dimension(10,6)        :: NE_interface_N_grdpts_id_test
          integer, dimension(10,6)        :: NW_interface_N_grdpts_id_test
          integer, dimension(10,6)        :: SE_interface_S_grdpts_id_test
          integer, dimension(10,6)        :: SW_interface_S_grdpts_id_test
                                          
          integer, dimension(8,6)         :: NE_interface_E_grdpts_id_test
          integer, dimension(8,6)         :: NW_interface_W_grdpts_id_test
          integer, dimension(8,6)         :: SE_interface_E_grdpts_id_test
          integer, dimension(8,6)         :: SW_interface_W_grdpts_id_test

          logical                         :: test_loc
          integer                         :: k
          
          test_validated = .true.


          !input
          !------------------------------------------------------------
          call ini_test_copy_from_to_mainlayer_neighbors(
     $         mainlayer_interface_used,
     $         
     $         NE_interface_N_nodes_test,
     $         NW_interface_N_nodes_test,
     $         SE_interface_S_nodes_test,
     $         SW_interface_S_nodes_test,
     $         
     $         NE_interface_E_nodes_test,
     $         NW_interface_W_nodes_test,
     $         SE_interface_E_nodes_test,
     $         SW_interface_W_nodes_test,
     $         
     $         NE_interface_N_grdpts_id_test,
     $         NW_interface_N_grdpts_id_test,
     $         SE_interface_S_grdpts_id_test,
     $         SW_interface_S_grdpts_id_test,
     $         
     $         NE_interface_E_grdpts_id_test,
     $         NW_interface_W_grdpts_id_test,
     $         SE_interface_E_grdpts_id_test,
     $         SW_interface_W_grdpts_id_test)


          !perform the tests
          !------------------------------------------------------------
          do k=1,8

             test_loc = perform_test_copy_from_neighbors(
     $            k,detailled,
     $            
     $            mainlayer_interface_used,
     $            
     $            NE_interface_N_nodes_test,
     $            NW_interface_N_nodes_test,
     $            SE_interface_S_nodes_test,
     $            SW_interface_S_nodes_test,
     $            
     $            NE_interface_E_nodes_test,
     $            NW_interface_W_nodes_test,
     $            SE_interface_E_nodes_test,
     $            SW_interface_W_nodes_test,
     $            
     $            NE_interface_N_grdpts_id_test,
     $            NW_interface_N_grdpts_id_test,
     $            SE_interface_S_grdpts_id_test,
     $            SW_interface_S_grdpts_id_test,
     $            
     $            NE_interface_E_grdpts_id_test,
     $            NW_interface_W_grdpts_id_test,
     $            SE_interface_E_grdpts_id_test,
     $            SW_interface_W_grdpts_id_test)

             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_copy_from_neighbors


        function test_copy_to_neighbors(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          type(mainlayer_interface_sync)  :: mainlayer_interface_used

          !nodes for the tests
          real(rkind), dimension(10,6,ne) :: NE_interface_N_nodes_test
          real(rkind), dimension(10,6,ne) :: NW_interface_N_nodes_test
          real(rkind), dimension(10,6,ne) :: SE_interface_S_nodes_test
          real(rkind), dimension(10,6,ne) :: SW_interface_S_nodes_test

          real(rkind), dimension(8,6,ne)  :: NE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne)  :: NW_interface_W_nodes_test
          real(rkind), dimension(8,6,ne)  :: SE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne)  :: SW_interface_W_nodes_test

          integer, dimension(10,6)        :: NE_interface_N_grdpts_id_test
          integer, dimension(10,6)        :: NW_interface_N_grdpts_id_test
          integer, dimension(10,6)        :: SE_interface_S_grdpts_id_test
          integer, dimension(10,6)        :: SW_interface_S_grdpts_id_test
                                          
          integer, dimension(8,6)         :: NE_interface_E_grdpts_id_test
          integer, dimension(8,6)         :: NW_interface_W_grdpts_id_test
          integer, dimension(8,6)         :: SE_interface_E_grdpts_id_test
          integer, dimension(8,6)         :: SW_interface_W_grdpts_id_test

          logical                         :: test_loc
          integer                         :: k
          
          test_validated = .true.


          !input
          !------------------------------------------------------------
          call ini_test_copy_from_to_mainlayer_neighbors(
     $         mainlayer_interface_used,
     $         
     $         NE_interface_N_nodes_test,
     $         NW_interface_N_nodes_test,
     $         SE_interface_S_nodes_test,
     $         SW_interface_S_nodes_test,
     $         
     $         NE_interface_E_nodes_test,
     $         NW_interface_W_nodes_test,
     $         SE_interface_E_nodes_test,
     $         SW_interface_W_nodes_test,
     $         
     $         NE_interface_N_grdpts_id_test,
     $         NW_interface_N_grdpts_id_test,
     $         SE_interface_S_grdpts_id_test,
     $         SW_interface_S_grdpts_id_test,
     $         
     $         NE_interface_E_grdpts_id_test,
     $         NW_interface_W_grdpts_id_test,
     $         SE_interface_E_grdpts_id_test,
     $         SW_interface_W_grdpts_id_test)


          !perform the tests
          !------------------------------------------------------------
          do k=1,8

             test_loc = perform_test_copy_to_neighbors(
     $            k,detailled,
     $            
     $            mainlayer_interface_used,
     $            
     $            NE_interface_N_nodes_test,
     $            NW_interface_N_nodes_test,
     $            SE_interface_S_nodes_test,
     $            SW_interface_S_nodes_test,
     $            
     $            NE_interface_E_nodes_test,
     $            NW_interface_W_nodes_test,
     $            SE_interface_E_nodes_test,
     $            SW_interface_W_nodes_test,
     $            
     $            NE_interface_N_grdpts_id_test,
     $            NW_interface_N_grdpts_id_test,
     $            SE_interface_S_grdpts_id_test,
     $            SW_interface_S_grdpts_id_test,
     $            
     $            NE_interface_E_grdpts_id_test,
     $            NW_interface_W_grdpts_id_test,
     $            SE_interface_E_grdpts_id_test,
     $            SW_interface_W_grdpts_id_test)

             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_copy_to_neighbors

        
        subroutine ini_test_copy_from_to_mainlayer_neighbors(
     $     mainlayer_interface_used,
     $         
     $     NE_interface_N_nodes_test,
     $     NW_interface_N_nodes_test,
     $     SE_interface_S_nodes_test,
     $     SW_interface_S_nodes_test,
     $     
     $     NE_interface_E_nodes_test,
     $     NW_interface_W_nodes_test,
     $     SE_interface_E_nodes_test,
     $     SW_interface_W_nodes_test,
     $     
     $     NE_interface_N_grdpts_id_test,
     $     NW_interface_N_grdpts_id_test,
     $     SE_interface_S_grdpts_id_test,
     $     SW_interface_S_grdpts_id_test,
     $     
     $     NE_interface_E_grdpts_id_test,
     $     NW_interface_W_grdpts_id_test,
     $     SE_interface_E_grdpts_id_test,
     $     SW_interface_W_grdpts_id_test)

          implicit none

          type(mainlayer_interface_sync) , intent(inout) :: mainlayer_interface_used

          real(rkind), dimension(10,6,ne), intent(out)   :: NE_interface_N_nodes_test
          real(rkind), dimension(10,6,ne), intent(out)   :: NW_interface_N_nodes_test
          real(rkind), dimension(10,6,ne), intent(out)   :: SE_interface_S_nodes_test
          real(rkind), dimension(10,6,ne), intent(out)   :: SW_interface_S_nodes_test

          real(rkind), dimension(8,6,ne) , intent(out)   :: NE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne) , intent(out)   :: NW_interface_W_nodes_test
          real(rkind), dimension(8,6,ne) , intent(out)   :: SE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne) , intent(out)   :: SW_interface_W_nodes_test

          integer, dimension(10,6)       , intent(out)   :: NE_interface_N_grdpts_id_test
          integer, dimension(10,6)       , intent(out)   :: NW_interface_N_grdpts_id_test
          integer, dimension(10,6)       , intent(out)   :: SE_interface_S_grdpts_id_test
          integer, dimension(10,6)       , intent(out)   :: SW_interface_S_grdpts_id_test
                                                         
          integer, dimension(8,6)        , intent(out)   :: NE_interface_E_grdpts_id_test
          integer, dimension(8,6)        , intent(out)   :: NW_interface_W_grdpts_id_test
          integer, dimension(8,6)        , intent(out)   :: SE_interface_E_grdpts_id_test
          integer, dimension(8,6)        , intent(out)   :: SW_interface_W_grdpts_id_test


          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          type(bf_sublayer), pointer :: bf_sublayer_added

          integer(ikind) :: i,j
          integer        :: k

          
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

          NW_interface_N_grdpts_id_test = reshape((/
     $         ((10*(align_N-3+j-1)+(align_W-6+i-1),i=1,10),j=1,6)/),
     $         (/10,6/))

          NE_interface_N_grdpts_id_test = reshape((/
     $         ((10*(align_N-3+j-1)+(align_E-5+i-1),i=1,10),j=1,6)/),
     $         (/10,6/))

          SW_interface_S_grdpts_id_test = reshape((/
     $         ((10*(align_S-4+j-1)+(align_W-6+i-1),i=1,10),j=1,6)/),
     $         (/10,6/))

          SE_interface_S_grdpts_id_test = reshape((/
     $         ((10*(align_S-4+j-1)+(align_E-5+i-1),i=1,10),j=1,6)/),
     $         (/10,6/))


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

          NW_interface_W_grdpts_id_test = reshape((/
     $         ((10*(align_N-5+j-1)+(align_W-6+i-1),i=1,8),j=1,6)/),
     $         (/8,6/))

          NE_interface_E_grdpts_id_test = reshape((/
     $         ((10*(align_N-5+j-1)+(align_E-3+i-1),i=1,8),j=1,6)/),
     $         (/8,6/))

          SW_interface_W_grdpts_id_test = reshape((/
     $         ((10*(align_S-2+j-1)+(align_W-6+i-1),i=1,8),j=1,6)/),
     $         (/8,6/))

          SE_interface_E_grdpts_id_test = reshape((/
     $         ((10*(align_S-2+j-1)+(align_E-3+i-1),i=1,8),j=1,6)/),
     $         (/8,6/))


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
          bf_sublayer_added%grdpts_id = NW_interface_N_grdpts_id_test
          bf_sublayer_added%nodes     = NW_interface_N_nodes_test
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
          bf_sublayer_added%grdpts_id = NW_interface_W_grdpts_id_test
          bf_sublayer_added%nodes     = NW_interface_W_nodes_test
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
          bf_sublayer_added%grdpts_id = NE_interface_N_grdpts_id_test
          bf_sublayer_added%nodes     = NE_interface_N_nodes_test
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
          bf_sublayer_added%grdpts_id = NE_interface_E_grdpts_id_test
          bf_sublayer_added%nodes     = NE_interface_E_nodes_test
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
          bf_sublayer_added%grdpts_id = SW_interface_S_grdpts_id_test
          bf_sublayer_added%nodes     = SW_interface_S_nodes_test
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
          bf_sublayer_added%grdpts_id = SW_interface_W_grdpts_id_test
          bf_sublayer_added%nodes     = SW_interface_W_nodes_test
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
          bf_sublayer_added%grdpts_id = SE_interface_S_grdpts_id_test
          bf_sublayer_added%nodes     = SE_interface_S_nodes_test
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
          bf_sublayer_added%grdpts_id = SE_interface_E_grdpts_id_test
          bf_sublayer_added%nodes     = SE_interface_E_nodes_test
          call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $         SE_interface_type,
     $         bf_sublayer_added)

        end subroutine ini_test_copy_from_to_mainlayer_neighbors


        function perform_test_copy_from_neighbors(
     $     test_id, detailled,
     $            
     $     mainlayer_interface_used,
     $     
     $     NE_interface_N_nodes_test,
     $     NW_interface_N_nodes_test,
     $     SE_interface_S_nodes_test,
     $     SW_interface_S_nodes_test,
     $     
     $     NE_interface_E_nodes_test,
     $     NW_interface_W_nodes_test,
     $     SE_interface_E_nodes_test,
     $     SW_interface_W_nodes_test,
     $     
     $     NE_interface_N_grdpts_id_test,
     $     NW_interface_N_grdpts_id_test,
     $     SE_interface_S_grdpts_id_test,
     $     SW_interface_S_grdpts_id_test,
     $     
     $     NE_interface_E_grdpts_id_test,
     $     NW_interface_W_grdpts_id_test,
     $     SE_interface_E_grdpts_id_test,
     $     SW_interface_W_grdpts_id_test)
     $     result(test_validated)


          implicit none

          integer                        , intent(in)    :: test_id
          logical                        , intent(in)    :: detailled
          logical                                        :: test_validated

          type(mainlayer_interface_sync) , intent(inout) :: mainlayer_interface_used

          real(rkind), dimension(10,6,ne), intent(in)    :: NE_interface_N_nodes_test
          real(rkind), dimension(10,6,ne), intent(in)    :: NW_interface_N_nodes_test
          real(rkind), dimension(10,6,ne), intent(in)    :: SE_interface_S_nodes_test
          real(rkind), dimension(10,6,ne), intent(in)    :: SW_interface_S_nodes_test
                                                         
          real(rkind), dimension(8,6,ne) , intent(in)    :: NE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne) , intent(in)    :: NW_interface_W_nodes_test
          real(rkind), dimension(8,6,ne) , intent(in)    :: SE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne) , intent(in)    :: SW_interface_W_nodes_test
                                                         
          integer, dimension(10,6)       , intent(in)    :: NE_interface_N_grdpts_id_test
          integer, dimension(10,6)       , intent(in)    :: NW_interface_N_grdpts_id_test
          integer, dimension(10,6)       , intent(in)    :: SE_interface_S_grdpts_id_test
          integer, dimension(10,6)       , intent(in)    :: SW_interface_S_grdpts_id_test
                                                         
          integer, dimension(8,6)        , intent(in)    :: NE_interface_E_grdpts_id_test
          integer, dimension(8,6)        , intent(in)    :: NW_interface_W_grdpts_id_test
          integer, dimension(8,6)        , intent(in)    :: SE_interface_E_grdpts_id_test
          integer, dimension(8,6)        , intent(in)    :: SW_interface_W_grdpts_id_test

          type(bf_sublayer), pointer :: current_sublayer


          !input
          call get_param_for_test_copy_from_neighbors(
     $         test_id, mainlayer_interface_used,
     $         current_sublayer)

          !output
          call mainlayer_interface_used%copy_from_mainlayer_neighbors(current_sublayer)

          !validation
          test_validated = is_test_copy_from_to_neighbor_validated(
     $         mainlayer_interface_used,
     $         NE_interface_N_nodes_test,
     $         NW_interface_N_nodes_test,
     $         SE_interface_S_nodes_test,
     $         SW_interface_S_nodes_test,
     $         
     $         NE_interface_E_nodes_test,
     $         NW_interface_W_nodes_test,
     $         SE_interface_E_nodes_test,
     $         SW_interface_W_nodes_test,
     $         
     $         NE_interface_N_grdpts_id_test,
     $         NW_interface_N_grdpts_id_test,
     $         SE_interface_S_grdpts_id_test,
     $         SW_interface_S_grdpts_id_test,
     $         
     $         NE_interface_E_grdpts_id_test,
     $         NW_interface_W_grdpts_id_test,
     $         SE_interface_E_grdpts_id_test,
     $         SW_interface_W_grdpts_id_test,
     $         
     $         detailled)

        end function perform_test_copy_from_neighbors


        function perform_test_copy_to_neighbors(
     $     test_id, detailled,
     $            
     $     mainlayer_interface_used,
     $     
     $     NE_interface_N_nodes_test,
     $     NW_interface_N_nodes_test,
     $     SE_interface_S_nodes_test,
     $     SW_interface_S_nodes_test,
     $     
     $     NE_interface_E_nodes_test,
     $     NW_interface_W_nodes_test,
     $     SE_interface_E_nodes_test,
     $     SW_interface_W_nodes_test,
     $     
     $     NE_interface_N_grdpts_id_test,
     $     NW_interface_N_grdpts_id_test,
     $     SE_interface_S_grdpts_id_test,
     $     SW_interface_S_grdpts_id_test,
     $     
     $     NE_interface_E_grdpts_id_test,
     $     NW_interface_W_grdpts_id_test,
     $     SE_interface_E_grdpts_id_test,
     $     SW_interface_W_grdpts_id_test)
     $     result(test_validated)


          implicit none

          integer                        , intent(in)    :: test_id
          logical                        , intent(in)    :: detailled
          logical                                        :: test_validated

          type(mainlayer_interface_sync) , intent(inout) :: mainlayer_interface_used

          real(rkind), dimension(10,6,ne), intent(in)    :: NE_interface_N_nodes_test
          real(rkind), dimension(10,6,ne), intent(in)    :: NW_interface_N_nodes_test
          real(rkind), dimension(10,6,ne), intent(in)    :: SE_interface_S_nodes_test
          real(rkind), dimension(10,6,ne), intent(in)    :: SW_interface_S_nodes_test
                                                         
          real(rkind), dimension(8,6,ne) , intent(in)    :: NE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne) , intent(in)    :: NW_interface_W_nodes_test
          real(rkind), dimension(8,6,ne) , intent(in)    :: SE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne) , intent(in)    :: SW_interface_W_nodes_test
                                                         
          integer, dimension(10,6)       , intent(in)    :: NE_interface_N_grdpts_id_test
          integer, dimension(10,6)       , intent(in)    :: NW_interface_N_grdpts_id_test
          integer, dimension(10,6)       , intent(in)    :: SE_interface_S_grdpts_id_test
          integer, dimension(10,6)       , intent(in)    :: SW_interface_S_grdpts_id_test
                                                         
          integer, dimension(8,6)        , intent(in)    :: NE_interface_E_grdpts_id_test
          integer, dimension(8,6)        , intent(in)    :: NW_interface_W_grdpts_id_test
          integer, dimension(8,6)        , intent(in)    :: SE_interface_E_grdpts_id_test
          integer, dimension(8,6)        , intent(in)    :: SW_interface_W_grdpts_id_test

          type(bf_sublayer), pointer :: current_sublayer


          !input
          call get_param_for_test_copy_to_neighbors(
     $         test_id, mainlayer_interface_used,
     $         current_sublayer)

          !output
          call mainlayer_interface_used%copy_to_mainlayer_neighbors(current_sublayer)

          !validation
          test_validated = is_test_copy_from_to_neighbor_validated(
     $         mainlayer_interface_used,
     $         NE_interface_N_nodes_test,
     $         NW_interface_N_nodes_test,
     $         SE_interface_S_nodes_test,
     $         SW_interface_S_nodes_test,
     $         
     $         NE_interface_E_nodes_test,
     $         NW_interface_W_nodes_test,
     $         SE_interface_E_nodes_test,
     $         SW_interface_W_nodes_test,
     $         
     $         NE_interface_N_grdpts_id_test,
     $         NW_interface_N_grdpts_id_test,
     $         SE_interface_S_grdpts_id_test,
     $         SW_interface_S_grdpts_id_test,
     $         
     $         NE_interface_E_grdpts_id_test,
     $         NW_interface_W_grdpts_id_test,
     $         SE_interface_E_grdpts_id_test,
     $         SW_interface_W_grdpts_id_test,
     $         
     $         detailled)

        end function perform_test_copy_to_neighbors


        subroutine get_param_for_test_copy_from_neighbors(
     $     test_id,
     $     mainlayer_interface_used,
     $     current_sublayer)

          implicit none

          integer                       , intent(in)    :: test_id
          type(mainlayer_interface_sync), intent(inout) :: mainlayer_interface_used
          type(bf_sublayer), pointer    , intent(out)   :: current_sublayer


          integer(ikind) :: i,j
          integer        :: k

          select case(test_id)

            case(1)
               mainlayer_interface_used%NW_interface_N_ptr%grdpts_id(1:8,1:4) =
     $              reshape((/ ((-99,i=1,8),j=1,4)/), (/8,4/))

               mainlayer_interface_used%NW_interface_N_ptr%nodes(1:8,1:4,:) =
     $              reshape((/(((-99.0d0,i=1,8),j=1,4),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%NW_interface_N_ptr

            case(2)
               mainlayer_interface_used%NE_interface_N_ptr%grdpts_id(3:10,1:4) =
     $              reshape((/ ((-99,i=3,10),j=1,4)/), (/8,4/))

               mainlayer_interface_used%NE_interface_N_ptr%nodes(3:10,1:4,:) =
     $              reshape((/(((-99.0d0,i=3,10),j=1,4),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%NE_interface_N_ptr

            case(3)
               mainlayer_interface_used%SW_interface_S_ptr%grdpts_id(1:8,3:6) =
     $              reshape((/ ((-99,i=1,8),j=3,6)/), (/8,4/))

               mainlayer_interface_used%SW_interface_S_ptr%nodes(1:8,3:6,:) =
     $              reshape((/(((-99.0d0,i=1,8),j=3,6),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%SW_interface_S_ptr

            case(4)
               mainlayer_interface_used%SE_interface_S_ptr%grdpts_id(3:10,3:6) =
     $              reshape((/ ((-99,i=3,10),j=3,6)/), (/8,4/))

               mainlayer_interface_used%SE_interface_S_ptr%nodes(3:10,3:6,:) =
     $              reshape((/(((-99.0d0,i=3,10),j=3,6),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%SE_interface_S_ptr

            case(5)
               mainlayer_interface_used%SE_interface_E_ptr%grdpts_id(1:8,1:4) =
     $              reshape((/ ((-99,i=1,8),j=1,4)/), (/8,4/))

               mainlayer_interface_used%SE_interface_E_ptr%nodes(1:8,1:4,:) =
     $              reshape((/(((-99.0d0,i=1,8),j=1,4),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%SE_interface_E_ptr

            case(6)
               mainlayer_interface_used%NE_interface_E_ptr%grdpts_id(1:8,3:6) =
     $              reshape((/ ((-99,i=1,8),j=1,4)/), (/8,4/))

               mainlayer_interface_used%NE_interface_E_ptr%nodes(1:8,3:6,:) =
     $              reshape((/(((-99.0d0,i=1,8),j=1,4),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%NE_interface_E_ptr

            case(7)
               mainlayer_interface_used%SW_interface_W_ptr%grdpts_id(1:8,1:4) =
     $              reshape((/ ((-99,i=1,8),j=1,4)/), (/8,4/))

               mainlayer_interface_used%SW_interface_W_ptr%nodes(1:8,1:4,:) =
     $              reshape((/(((-99.0d0,i=1,8),j=1,4),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%SW_interface_W_ptr

            case(8)
               mainlayer_interface_used%NW_interface_W_ptr%grdpts_id(1:8,3:6) =
     $              reshape((/ ((-99,i=1,8),j=1,4)/), (/8,4/))

               mainlayer_interface_used%NW_interface_W_ptr%nodes(1:8,3:6,:) =
     $              reshape((/(((-99.0d0,i=1,8),j=1,4),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%NW_interface_W_ptr

          end select

        end subroutine get_param_for_test_copy_from_neighbors


        subroutine get_param_for_test_copy_to_neighbors(
     $     test_id,
     $     mainlayer_interface_used,
     $     current_sublayer)

          implicit none

          integer                       , intent(in)    :: test_id
          type(mainlayer_interface_sync), intent(inout) :: mainlayer_interface_used
          type(bf_sublayer), pointer    , intent(out)   :: current_sublayer


          integer(ikind) :: i,j
          integer        :: k

          select case(test_id)

            case(1)
               mainlayer_interface_used%NW_interface_W_ptr%grdpts_id(1:8,3:6) =
     $              reshape((/ ((-99,i=1,8),j=1,4)/), (/8,4/))

               mainlayer_interface_used%NW_interface_W_ptr%nodes(1:8,3:6,:) =
     $              reshape((/(((-99.0d0,i=1,8),j=1,4),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%NW_interface_N_ptr

            case(2)
               mainlayer_interface_used%NE_interface_E_ptr%grdpts_id(1:8,3:6) =
     $              reshape((/ ((-99,i=1,8),j=1,4)/), (/8,4/))

               mainlayer_interface_used%NE_interface_E_ptr%nodes(1:8,3:6,:) =
     $              reshape((/(((-99.0d0,i=1,8),j=1,4),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%NE_interface_N_ptr

            case(3)
               mainlayer_interface_used%SW_interface_W_ptr%grdpts_id(1:8,1:4) =
     $              reshape((/ ((-99,i=1,8),j=1,4)/), (/8,4/))

               mainlayer_interface_used%SW_interface_W_ptr%nodes(1:8,1:4,:) =
     $              reshape((/(((-99.0d0,i=1,8),j=1,4),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%SW_interface_S_ptr

            case(4)
               mainlayer_interface_used%SE_interface_E_ptr%grdpts_id(1:8,1:4) =
     $              reshape((/ ((-99,i=1,8),j=1,4)/), (/8,4/))

               mainlayer_interface_used%SE_interface_E_ptr%nodes(1:8,1:4,:) =
     $              reshape((/(((-99.0d0,i=1,8),j=1,4),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%SE_interface_S_ptr

            case(5)
               mainlayer_interface_used%SE_interface_S_ptr%grdpts_id(3:10,3:6) =
     $              reshape((/ ((-99,i=3,10),j=3,6)/), (/8,4/))

               mainlayer_interface_used%SE_interface_S_ptr%nodes(3:10,3:6,:) =
     $              reshape((/(((-99.0d0,i=3,10),j=3,6),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%SE_interface_E_ptr

            case(6)
               mainlayer_interface_used%NE_interface_N_ptr%grdpts_id(3:10,1:4) =
     $              reshape((/ ((-99,i=3,10),j=1,4)/), (/8,4/))

               mainlayer_interface_used%NE_interface_N_ptr%nodes(3:10,1:4,:) =
     $              reshape((/(((-99.0d0,i=3,10),j=1,4),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%NE_interface_E_ptr

            case(7)
               mainlayer_interface_used%SW_interface_S_ptr%grdpts_id(1:8,3:6) =
     $              reshape((/ ((-99,i=1,8),j=3,6)/), (/8,4/))

               mainlayer_interface_used%SW_interface_S_ptr%nodes(1:8,3:6,:) =
     $              reshape((/(((-99.0d0,i=1,8),j=3,6),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%SW_interface_W_ptr

            case(8)
               mainlayer_interface_used%NW_interface_N_ptr%grdpts_id(1:8,1:4) =
     $              reshape((/ ((-99,i=1,8),j=1,4)/), (/8,4/))

               mainlayer_interface_used%NW_interface_N_ptr%nodes(1:8,1:4,:) =
     $              reshape((/(((-99.0d0,i=1,8),j=1,4),k=1,ne)/), (/8,4,ne/))

               current_sublayer => mainlayer_interface_used%NW_interface_W_ptr

          end select

        end subroutine get_param_for_test_copy_to_neighbors


        function is_test_copy_from_to_neighbor_validated(
     $         mainlayer_interface_used,
     $     
     $         NE_interface_N_nodes_test,
     $         NW_interface_N_nodes_test,
     $         SE_interface_S_nodes_test,
     $         SW_interface_S_nodes_test,
     $         
     $         NE_interface_E_nodes_test,
     $         NW_interface_W_nodes_test,
     $         SE_interface_E_nodes_test,
     $         SW_interface_W_nodes_test,
     $         
     $         NE_interface_N_grdpts_id_test,
     $         NW_interface_N_grdpts_id_test,
     $         SE_interface_S_grdpts_id_test,
     $         SW_interface_S_grdpts_id_test,
     $         
     $         NE_interface_E_grdpts_id_test,
     $         NW_interface_W_grdpts_id_test,
     $         SE_interface_E_grdpts_id_test,
     $         SW_interface_W_grdpts_id_test,
     $     
     $         detailled)
     $     result(test_validated)

          implicit none

          type(mainlayer_interface_sync) , intent(in) :: mainlayer_interface_used
          real(rkind), dimension(10,6,ne), intent(in) :: NE_interface_N_nodes_test
          real(rkind), dimension(10,6,ne), intent(in) :: NW_interface_N_nodes_test
          real(rkind), dimension(10,6,ne), intent(in) :: SE_interface_S_nodes_test
          real(rkind), dimension(10,6,ne), intent(in) :: SW_interface_S_nodes_test

          real(rkind), dimension(8,6,ne) , intent(in) :: NE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne) , intent(in) :: NW_interface_W_nodes_test
          real(rkind), dimension(8,6,ne) , intent(in) :: SE_interface_E_nodes_test
          real(rkind), dimension(8,6,ne) , intent(in) :: SW_interface_W_nodes_test

          integer, dimension(10,6)       , intent(in) :: NE_interface_N_grdpts_id_test
          integer, dimension(10,6)       , intent(in) :: NW_interface_N_grdpts_id_test
          integer, dimension(10,6)       , intent(in) :: SE_interface_S_grdpts_id_test
          integer, dimension(10,6)       , intent(in) :: SW_interface_S_grdpts_id_test

          integer, dimension(8,6)        , intent(in) :: NE_interface_E_grdpts_id_test
          integer, dimension(8,6)        , intent(in) :: NW_interface_W_grdpts_id_test
          integer, dimension(8,6)        , intent(in) :: SE_interface_E_grdpts_id_test
          integer, dimension(8,6)        , intent(in) :: SW_interface_W_grdpts_id_test

          logical                        , intent(in) :: detailled

          logical                                     :: test_validated


          logical :: test_loc


          test_validated = .true.


          !grdpts_id
          !------------------------------------------------------------
          !north
          test_loc = is_int_matrix_validated(
     $         mainlayer_interface_used%NW_interface_N_ptr%grdpts_id,
     $         NW_interface_N_grdpts_id_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NW_interface_N_grdpts failed'')'
          end if

          test_loc = is_int_matrix_validated(
     $         mainlayer_interface_used%NE_interface_N_ptr%grdpts_id,
     $         NE_interface_N_grdpts_id_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NE_interface_N_grdpts failed'')'
          end if

          !south
          test_loc = is_int_matrix_validated(
     $         mainlayer_interface_used%SW_interface_S_ptr%grdpts_id,
     $         SW_interface_S_grdpts_id_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SW_interface_S_grdpts failed'')'
          end if

          test_loc = is_int_matrix_validated(
     $         mainlayer_interface_used%SE_interface_S_ptr%grdpts_id,
     $         SE_interface_S_grdpts_id_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SE_interface_S_grdpts failed'')'
          end if

          !west
          test_loc = is_int_matrix_validated(
     $         mainlayer_interface_used%NW_interface_W_ptr%grdpts_id,
     $         NW_interface_W_grdpts_id_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NW_interface_W_grdpts failed'')'
          end if

          test_loc = is_int_matrix_validated(
     $         mainlayer_interface_used%SW_interface_W_ptr%grdpts_id,
     $         SW_interface_W_grdpts_id_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SW_interface_W_grdpts failed'')'
          end if

          !east
          test_loc = is_int_matrix_validated(
     $         mainlayer_interface_used%NE_interface_E_ptr%grdpts_id,
     $         NE_interface_E_grdpts_id_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NE_interface_E_grdpts failed'')'
          end if

          test_loc = is_int_matrix_validated(
     $         mainlayer_interface_used%SE_interface_E_ptr%grdpts_id,
     $         SE_interface_E_grdpts_id_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SE_interface_E_grdpts failed'')'
          end if


          !nodes
          !------------------------------------------------------------
          !north
          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%NW_interface_N_ptr%nodes,
     $         NW_interface_N_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NW_interface_N_nodes failed'')'
          end if

          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%NE_interface_N_ptr%nodes,
     $         NE_interface_N_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NE_interface_N_nodes failed'')'
          end if

          !south
          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%SW_interface_S_ptr%nodes,
     $         SW_interface_S_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SW_interface_S_nodes failed'')'
          end if

          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%SE_interface_S_ptr%nodes,
     $         SE_interface_S_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SE_interface_S_nodes failed'')'
          end if

          !west
          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%NW_interface_W_ptr%nodes,
     $         NW_interface_W_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NW_interface_W_nodes failed'')'
          end if

          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%SW_interface_W_ptr%nodes,
     $         SW_interface_W_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SW_interface_W_nodes failed'')'
          end if

          !east
          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%NE_interface_E_ptr%nodes,
     $         NE_interface_E_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''NE_interface_E_nodes failed'')'
          end if

          test_loc = is_real_matrix3D_validated(
     $         mainlayer_interface_used%SE_interface_E_ptr%nodes,
     $         SE_interface_E_nodes_test,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''SE_interface_E_nodes failed'')'
          end if


        end function is_test_copy_from_to_neighbor_validated


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


      end program test_mainlayer_interface_sync
