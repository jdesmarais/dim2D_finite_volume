      program test_mainlayer_interface_dyn

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_int_vector_validated,
     $       is_int_matrix_validated

        use mainlayer_interface_dyn_class, only :
     $       mainlayer_interface_dyn

        use parameters_bf_layer, only :
     $       align_N, align_S,
     $       align_E, align_W,
     $       
     $       NE_interface_type, NW_interface_type,
     $       SE_interface_type, SW_interface_type

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


        test_loc = test_update_alignment_and_sync_properties(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_update_alignment_and_sync_properties: '',L1)', test_loc
        print '()'


        test_loc = test_uniformize_alignment_for_mainlayer_interface(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_uniformize_alignment_for_mainlayer_interface: '',L1)', test_loc
        print '()'


        contains


        function test_update_alignment_and_sync_properties(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(mainlayer_interface_dyn) :: mainlayer_interface_used

          integer, dimension(2,2,4,2) :: test_alignment
          integer, dimension(2,2,4,2) :: test_alignment_after
          integer, dimension(4,4,2)   :: test_sync_after

          logical               :: can_exchange_with_neighbor1
          logical               :: can_exchange_with_neighbor2
          integer               :: mainlayer_interface_type1
          integer               :: mainlayer_interface_type2
          integer, dimension(4) :: test_sync

          integer :: k,l
          logical :: test_loc

          
          test_validated = .true.


          test_alignment(:,:,N,1) = reshape((/ align_W+5, align_N  , align_E-5, align_N+1 /),(/2,2/))
          test_alignment(:,:,N,2) = reshape((/ align_W+4, align_N  , align_E-4, align_N+1 /),(/2,2/))
          test_alignment(:,:,S,1) = reshape((/ align_W+5, align_S-1, align_E-5, align_S   /),(/2,2/))
          test_alignment(:,:,S,2) = reshape((/ align_W+4, align_S-1, align_E-4, align_S   /),(/2,2/))
          test_alignment(:,:,E,1) = reshape((/ align_W-5, align_S+5, align_W  , align_N-5 /),(/2,2/))
          test_alignment(:,:,E,2) = reshape((/ align_W-5, align_S+4, align_W  , align_N-4 /),(/2,2/))
          test_alignment(:,:,W,1) = reshape((/ align_E  , align_S+5, align_E+5, align_N-5 /),(/2,2/))
          test_alignment(:,:,W,2) = reshape((/ align_E  , align_S+4, align_E+5, align_N-4 /),(/2,2/))

          test_alignment_after(:,:,N,1) = reshape((/ align_W+5, align_N  , align_E-5, align_N+1 /),(/2,2/))
          test_alignment_after(:,:,N,2) = reshape((/ align_W+1, align_N  , align_E-1, align_N+1 /),(/2,2/))
          test_alignment_after(:,:,S,1) = reshape((/ align_W+5, align_S-1, align_E-5, align_S   /),(/2,2/))
          test_alignment_after(:,:,S,2) = reshape((/ align_W+1, align_S-1, align_E-1, align_S   /),(/2,2/))
          test_alignment_after(:,:,E,1) = reshape((/ align_W-5, align_S+5, align_W  , align_N-5 /),(/2,2/))
          test_alignment_after(:,:,E,2) = reshape((/ align_W-5, align_S+1, align_W  , align_N-1 /),(/2,2/))
          test_alignment_after(:,:,W,1) = reshape((/ align_E  , align_S+5, align_E+5, align_N-5 /),(/2,2/))
          test_alignment_after(:,:,W,2) = reshape((/ align_E  , align_S+1, align_E+5, align_N-1 /),(/2,2/))

          test_sync_after(:,N,1) = [0,0,0,0]
          test_sync_after(:,N,2) = [1,NW_interface_type,1,NE_interface_type]
          test_sync_after(:,S,1) = [0,0,0,0]
          test_sync_after(:,S,2) = [1,SW_interface_type,1,SE_interface_type]
          test_sync_after(:,E,1) = [0,0,0,0]
          test_sync_after(:,E,2) = [1,SE_interface_type,1,NE_interface_type]
          test_sync_after(:,W,1) = [0,0,0,0]
          test_sync_after(:,W,2) = [1,SW_interface_type,1,NW_interface_type]


          !loop over the main layers
          do k=1,4
             !loop over the neighbors
             do l=1,2

                !output
                call mainlayer_interface_used%update_alignment_and_sync_properties(
     $               k,
     $               test_alignment(:,:,k,l),
     $               can_exchange_with_neighbor1,
     $               mainlayer_interface_type1,
     $               can_exchange_with_neighbor2,
     $               mainlayer_interface_type2)

                test_loc = is_int_matrix_validated(
     $               test_alignment(:,:,k,l),
     $               test_alignment_after(:,:,k,l),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''alignment('',2I2,'') failed'')',k,l
                end if

                if(can_exchange_with_neighbor1) then 
                   test_sync(1) = 1
                else
                   test_sync(1) = 0
                end if
                test_sync(2) = mainlayer_interface_type1                

                if(can_exchange_with_neighbor2) then 
                   test_sync(3) = 1
                else
                   test_sync(3) = 0
                end if
                test_sync(4) = mainlayer_interface_type2

                test_loc = is_int_vector_validated(
     $               test_sync,
     $               test_sync_after(:,k,l),
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''sync('',2I2,'') failed'')',k,l
                end if                

             end do
          end do


          !test with existing neighbors
          test_loc = perform_test_with_N_neighbors(
     $         detailled)
          test_validated = test_validated.and.test_loc

          test_loc = perform_test_with_S_neighbors(
     $         detailled)
          test_validated = test_validated.and.test_loc

          test_loc = perform_test_with_E_neighbors(
     $         detailled)
          test_validated = test_validated.and.test_loc

          test_loc = perform_test_with_W_neighbors(
     $         detailled)
          test_validated = test_validated.and.test_loc

       end function test_update_alignment_and_sync_properties


       function perform_test_with_N_neighbors(detailled)
     $     result(test_validated)

         implicit none
         
         logical, intent(in) :: detailled
         logical             :: test_validated

         type(mainlayer_interface_dyn) :: mainlayer_interface_used

         logical               :: can_exchange_with_neighbor1
         logical               :: can_exchange_with_neighbor2
         integer               :: mainlayer_interface_type1
         integer               :: mainlayer_interface_type2

         integer, dimension(2,2) :: alignment_in
         integer, dimension(2,2) :: alignment_out
         integer, dimension(4)   :: sync_in
         integer, dimension(4)   :: sync_out

         logical :: test_loc

         test_validated = .true.
          
         call mainlayer_interface_used%ini()
         
         allocate(mainlayer_interface_used%NW_interface_W_ptr)
         mainlayer_interface_used%NW_interface_W_ptr%alignment = reshape(
     $        (/align_W-5,align_N-3,align_W,align_N-1/),
     $        (/2,2/))
         
         allocate(mainlayer_interface_used%NE_interface_E_ptr)
         mainlayer_interface_used%NE_interface_E_ptr%alignment = reshape(
     $        (/align_E,align_N-3,align_E+5,align_N-1/),
     $        (/2,2/))
         
         alignment_in  = reshape(
     $        (/align_W+4, align_N, align_E-4, align_N+1/),
     $        (/2,2/))
         alignment_out = reshape(
     $        (/align_W-5, align_N, align_E+5, align_N+1/),
     $        (/2,2/))
         call mainlayer_interface_used%update_alignment_and_sync_properties(
     $        N,
     $        alignment_in,
     $        can_exchange_with_neighbor1,
     $        mainlayer_interface_type1,
     $        can_exchange_with_neighbor2,
     $        mainlayer_interface_type2)

         sync_out = [1,NW_interface_type,1,NE_interface_type]

         test_loc = is_int_matrix_validated(
     $        alignment_in,
     $        alignment_out,
     $        detailled)
         test_validated = test_validated.and.test_loc
         if(detailled.and.(.not.test_loc)) then
            print '(''alignment('',2I2,'') failed'')',N,3
         end if

         if(can_exchange_with_neighbor1) then 
            sync_in(1) = 1
         else
            sync_in(1) = 0
         end if
         sync_in(2) = mainlayer_interface_type1                
         if(can_exchange_with_neighbor2) then 
            sync_in(3) = 1
         else
            sync_in(3) = 0
         end if
         sync_in(4) = mainlayer_interface_type2

         test_loc = is_int_vector_validated(
     $        sync_in,
     $        sync_out,
     $        detailled)
         test_validated = test_validated.and.test_loc
         if(detailled.and.(.not.test_loc)) then
            print '(''sync('',2I2,'') failed'')',N,3
         end if

        end function perform_test_with_N_neighbors


        function perform_test_with_S_neighbors(detailled)
     $     result(test_validated)

         implicit none
         
         logical, intent(in) :: detailled
         logical             :: test_validated

         type(mainlayer_interface_dyn) :: mainlayer_interface_used

         logical               :: can_exchange_with_neighbor1
         logical               :: can_exchange_with_neighbor2
         integer               :: mainlayer_interface_type1
         integer               :: mainlayer_interface_type2

         integer, dimension(2,2) :: alignment_in
         integer, dimension(2,2) :: alignment_out
         integer, dimension(4)   :: sync_in
         integer, dimension(4)   :: sync_out

         logical :: test_loc

         test_validated = .true.
          
         call mainlayer_interface_used%ini()
         
         allocate(mainlayer_interface_used%SW_interface_W_ptr)
         mainlayer_interface_used%SW_interface_W_ptr%alignment = reshape(
     $        (/align_W-5,align_S+1,align_W,align_S+3/),
     $        (/2,2/))
         
         allocate(mainlayer_interface_used%SE_interface_E_ptr)
         mainlayer_interface_used%SE_interface_E_ptr%alignment = reshape(
     $        (/align_E,align_S+1,align_E+5,align_S+3/),
     $        (/2,2/))
         
         alignment_in  = reshape(
     $        (/align_W+4, align_S-1, align_E-4, align_S/),
     $        (/2,2/))
         alignment_out = reshape(
     $        (/align_W-5, align_S-1, align_E+5, align_S/),
     $        (/2,2/))
         call mainlayer_interface_used%update_alignment_and_sync_properties(
     $        S,
     $        alignment_in,
     $        can_exchange_with_neighbor1,
     $        mainlayer_interface_type1,
     $        can_exchange_with_neighbor2,
     $        mainlayer_interface_type2)

         sync_out = [1,SW_interface_type,1,SE_interface_type]

         test_loc = is_int_matrix_validated(
     $        alignment_in,
     $        alignment_out,
     $        detailled)
         test_validated = test_validated.and.test_loc
         if(detailled.and.(.not.test_loc)) then
            print '(''alignment('',2I2,'') failed'')',S,3
         end if

         if(can_exchange_with_neighbor1) then 
            sync_in(1) = 1
         else
            sync_in(1) = 0
         end if
         sync_in(2) = mainlayer_interface_type1                
         if(can_exchange_with_neighbor2) then 
            sync_in(3) = 1
         else
            sync_in(3) = 0
         end if
         sync_in(4) = mainlayer_interface_type2

         test_loc = is_int_vector_validated(
     $        sync_in,
     $        sync_out,
     $        detailled)
         test_validated = test_validated.and.test_loc
         if(detailled.and.(.not.test_loc)) then
            print '(''sync('',2I2,'') failed'')',S,3
         end if

        end function perform_test_with_S_neighbors


        function perform_test_with_E_neighbors(detailled)
     $     result(test_validated)

         implicit none
         
         logical, intent(in) :: detailled
         logical             :: test_validated

         type(mainlayer_interface_dyn) :: mainlayer_interface_used

         logical               :: can_exchange_with_neighbor1
         logical               :: can_exchange_with_neighbor2
         integer               :: mainlayer_interface_type1
         integer               :: mainlayer_interface_type2

         integer, dimension(2,2) :: alignment_in
         integer, dimension(2,2) :: alignment_out
         integer, dimension(4)   :: sync_in
         integer, dimension(4)   :: sync_out

         logical :: test_loc

         test_validated = .true.
          
         call mainlayer_interface_used%ini()
         
         allocate(mainlayer_interface_used%NE_interface_N_ptr)
         mainlayer_interface_used%NE_interface_N_ptr%alignment = reshape(
     $        (/align_E-6,align_N,align_E+10,align_N+1/),
     $        (/2,2/))
         
         allocate(mainlayer_interface_used%SE_interface_S_ptr)
         mainlayer_interface_used%SE_interface_S_ptr%alignment = reshape(
     $        (/align_E-6,align_S-1,align_E+15,align_S/),
     $        (/2,2/))
         
         alignment_in  = reshape(
     $        (/align_E, align_S+4, align_E+5 , align_N-4/),
     $        (/2,2/))
         alignment_out = reshape(
     $        (/align_E, align_S+1, align_E+15, align_N-1/),
     $        (/2,2/))

         call mainlayer_interface_used%update_alignment_and_sync_properties(
     $        E,
     $        alignment_in,
     $        can_exchange_with_neighbor1,
     $        mainlayer_interface_type1,
     $        can_exchange_with_neighbor2,
     $        mainlayer_interface_type2)

         sync_out = [1,SE_interface_type,1,NE_interface_type]

         test_loc = is_int_matrix_validated(
     $        alignment_in,
     $        alignment_out,
     $        detailled)
         test_validated = test_validated.and.test_loc
         if(detailled.and.(.not.test_loc)) then
            print '(''alignment('',2I2,'') failed'')',E,3
         end if

         if(can_exchange_with_neighbor1) then 
            sync_in(1) = 1
         else
            sync_in(1) = 0
         end if
         sync_in(2) = mainlayer_interface_type1                
         if(can_exchange_with_neighbor2) then 
            sync_in(3) = 1
         else
            sync_in(3) = 0
         end if
         sync_in(4) = mainlayer_interface_type2

         test_loc = is_int_vector_validated(
     $        sync_in,
     $        sync_out,
     $        detailled)
         test_validated = test_validated.and.test_loc
         if(detailled.and.(.not.test_loc)) then
            print '(''sync('',2I2,'') failed'')',E,3
         end if

        end function perform_test_with_E_neighbors


        function perform_test_with_W_neighbors(detailled)
     $     result(test_validated)

         implicit none
         
         logical, intent(in) :: detailled
         logical             :: test_validated

         type(mainlayer_interface_dyn) :: mainlayer_interface_used

         logical               :: can_exchange_with_neighbor1
         logical               :: can_exchange_with_neighbor2
         integer               :: mainlayer_interface_type1
         integer               :: mainlayer_interface_type2

         integer, dimension(2,2) :: alignment_in
         integer, dimension(2,2) :: alignment_out
         integer, dimension(4)   :: sync_in
         integer, dimension(4)   :: sync_out

         logical :: test_loc

         test_validated = .true.
          
         call mainlayer_interface_used%ini()
         
         allocate(mainlayer_interface_used%NW_interface_N_ptr)
         mainlayer_interface_used%NW_interface_N_ptr%alignment = reshape(
     $        (/align_W-10,align_N,align_W+6,align_N+1/),
     $        (/2,2/))
         
         allocate(mainlayer_interface_used%SW_interface_S_ptr)
         mainlayer_interface_used%SW_interface_S_ptr%alignment = reshape(
     $        (/align_W-15,align_S-1,align_W+6,align_S/),
     $        (/2,2/))
         
         alignment_in  = reshape(
     $        (/align_W-5 , align_S+4, align_W, align_N-4/),
     $        (/2,2/))
         alignment_out = reshape(
     $        (/align_W-15, align_S+1, align_W, align_N-1/),
     $        (/2,2/))

         call mainlayer_interface_used%update_alignment_and_sync_properties(
     $        W,
     $        alignment_in,
     $        can_exchange_with_neighbor1,
     $        mainlayer_interface_type1,
     $        can_exchange_with_neighbor2,
     $        mainlayer_interface_type2)

         sync_out = [1,SW_interface_type,1,NW_interface_type]

         test_loc = is_int_matrix_validated(
     $        alignment_in,
     $        alignment_out,
     $        detailled)
         test_validated = test_validated.and.test_loc
         if(detailled.and.(.not.test_loc)) then
            print '(''alignment('',2I2,'') failed'')',E,3
         end if

         if(can_exchange_with_neighbor1) then 
            sync_in(1) = 1
         else
            sync_in(1) = 0
         end if
         sync_in(2) = mainlayer_interface_type1                
         if(can_exchange_with_neighbor2) then 
            sync_in(3) = 1
         else
            sync_in(3) = 0
         end if
         sync_in(4) = mainlayer_interface_type2

         test_loc = is_int_vector_validated(
     $        sync_in,
     $        sync_out,
     $        detailled)
         test_validated = test_validated.and.test_loc
         if(detailled.and.(.not.test_loc)) then
            print '(''sync('',2I2,'') failed'')',W,3
         end if

        end function perform_test_with_W_neighbors


        function test_uniformize_alignment_for_mainlayer_interface(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,20

             test_loc = perform_test_uniformize(k,detailled)
             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_uniformize_alignment_for_mainlayer_interface


        function perform_test_uniformize(test_id,detailled)
     $     result(test_validated)

          implicit none

          integer, intent(in) :: test_id
          logical, intent(in) :: detailled
          logical             :: test_validated

          type(mainlayer_interface_dyn)  :: mainlayer_interface_used
          type(bf_sublayer), pointer     :: bf_sublayer_tested
          logical                        :: test_should_be_reallocated
          integer(ikind), dimension(2,2) :: test_bf_alignment
          logical                        :: should_be_reallocated
          integer(ikind), dimension(2,2) :: bf_alignment

          logical :: test_loc


          test_validated = .true.


          !input
          call get_param_test_uniformize(
     $         test_id,mainlayer_interface_used,
     $         bf_sublayer_tested,
     $         test_should_be_reallocated,
     $         test_bf_alignment)

          !output
          bf_alignment = mainlayer_interface_used%uniformize_alignment_for_mainlayer_interface(
     $         bf_sublayer_tested,
     $         should_be_reallocated)

          !validation
          test_loc = should_be_reallocated.eqv.test_should_be_reallocated
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''should_be_reallocated failed'')'
          end if

          test_loc = is_int_matrix_validated(
     $         bf_alignment,
     $         test_bf_alignment,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''bf_alignment failed'')'
          end if

        end function perform_test_uniformize


        subroutine get_param_test_uniformize(
     $     test_id,
     $     mainlayer_interface_used,
     $     bf_sublayer_tested,
     $     test_should_be_reallocated,
     $     test_bf_alignment)

          implicit none

          integer                       , intent(in)    :: test_id
          type(mainlayer_interface_dyn) , intent(inout) :: mainlayer_interface_used
          type(bf_sublayer), pointer    , intent(out)   :: bf_sublayer_tested
          logical                       , intent(out)   :: test_should_be_reallocated
          integer(ikind), dimension(2,2), intent(out)   :: test_bf_alignment


          type(bf_sublayer), pointer       :: bf_sublayer_ptr
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          


          call mainlayer_interface_used%ini()


          select case(test_id)
            case(1)

               !north buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(N)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-5,align_N,align_W+6,align_N+1/),
     $                (/2,2/)))               

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               test_should_be_reallocated = .false.
               test_bf_alignment = reshape(
     $              (/align_W-5,align_N,align_W+6,align_N+1/),
     $              (/2,2/))

            case(2)

               !north buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(N)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-5,align_N,align_W+6,align_N+1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !west buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(W)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-10,align_N-2,align_W,align_N-1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_W-10,align_N,align_W+6,align_N+1/),
     $              (/2,2/))

            case(3)

               !north buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(N)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-5,align_N,align_W+6,align_N+1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !west buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(W)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-10,align_S+1,align_W,align_N-1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr)

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr)

               !south buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(S)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-15,align_S-1,align_W+6,align_S/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_W-15,align_N,align_W+6,align_N+1/),
     $              (/2,2/))

            case(4)

               !north buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(N)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E-6,align_N,align_E+5,align_N+1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               test_should_be_reallocated = .false.
               test_bf_alignment = reshape(
     $              (/align_E-6,align_N,align_E+5,align_N+1/),
     $              (/2,2/))

            case(5)

               !north buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(N)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E-6,align_N,align_E+5,align_N+1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               test_should_be_reallocated = .false.
               test_bf_alignment = reshape(
     $              (/align_E-6,align_N,align_E+5,align_N+1/),
     $              (/2,2/))

               !east buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(E)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E,align_N-2,align_E+10,align_N-1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_E-6,align_N,align_E+10,align_N+1/),
     $              (/2,2/))

            case(6)

               !north buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(N)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E-6,align_N,align_E+5,align_N+1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !west buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(E)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E,align_S+1,align_E+10,align_N-1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr)

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr)

               !south buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(S)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E-6,align_S-1,align_E+15,align_S/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_E-6,align_N,align_E+15,align_N+1/),
     $              (/2,2/))

            case(7)

               !south buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(S)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-5,align_S-1,align_W+6,align_S/),
     $                (/2,2/)))               

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               test_should_be_reallocated = .false.
               test_bf_alignment = reshape(
     $              (/align_W-5,align_S-1,align_W+6,align_S/),
     $              (/2,2/))

            case(8)

               !south buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(S)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-5,align_S-1,align_W+6,align_S/),
     $                (/2,2/)))               

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !west buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(W)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-10,align_S+1,align_W,align_S+2/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_W-10,align_S-1,align_W+6,align_S/),
     $              (/2,2/))

            case(9)

               !south buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(S)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-5,align_S-1,align_W+6,align_S/),
     $                (/2,2/)))               

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !west buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(W)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-10,align_S+1,align_W,align_N-1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr)

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr)

               !north buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(N)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-15,align_N,align_W+6,align_N+1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_W-15,align_S-1,align_W+6,align_S/),
     $              (/2,2/))

            case(10)

               !south buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(S)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E-6,align_S-1,align_E+5,align_S/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               test_should_be_reallocated = .false.
               test_bf_alignment = reshape(
     $              (/align_E-6,align_S-1,align_E+5,align_S/),
     $              (/2,2/))

            case(11)

               !south buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(S)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E-6,align_S-1,align_E+5,align_S/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !east buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(E)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E,align_S+1,align_E+10,align_S+2/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_E-6,align_S-1,align_E+10,align_S/),
     $              (/2,2/))

            case(12)

               !south buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(S)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E-6,align_S-1,align_E+5,align_S/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !east buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(E)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E,align_S+1,align_E+10,align_N-1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr)

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr)

               !north buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(N)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E-6,align_N,align_E+15,align_N+1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_E-6,align_S-1,align_E+15,align_S/),
     $              (/2,2/))


            case(13)

               !east buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(E)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E,align_S+1,align_E+5,align_S+2/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               test_should_be_reallocated = .false.
               test_bf_alignment = reshape(
     $              (/align_E,align_S+1,align_E+5,align_S+2/),
     $              (/2,2/))

            case(14)

               !east buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(E)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E,align_S+1,align_E+5,align_S+2/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !south buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(S)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E-5,align_S-1,align_E+10,align_S/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_E,align_S+1,align_E+10,align_S+2/),
     $              (/2,2/))

            case(15)

               !east buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(E)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E,align_N-2,align_E+5,align_N-1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !north buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(N)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E-5,align_N,align_E+10,align_N+1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_E,align_N-2,align_E+10,align_N-1/),
     $              (/2,2/))

            case(16)

               !east buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(E)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E,align_S+1,align_E+5,align_N-1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr)

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !south buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(S)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E-5,align_S-1,align_E+15,align_S/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SE_interface_type,
     $              bf_sublayer_ptr)

               !north buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(N)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_E-5,align_N,align_E+10,align_N+1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NE_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_E,align_S+1,align_E+15,align_N-1/),
     $              (/2,2/))

            case(17)

               !west buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(W)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-5,align_S+1,align_W,align_S+2/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               test_should_be_reallocated = .false.
               test_bf_alignment = reshape(
     $              (/align_W-5,align_S+1,align_W,align_S+2/),
     $              (/2,2/))

            case(18)

               !west buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(W)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-5,align_S+1,align_W,align_S+2/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !south buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(S)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-10,align_S-1,align_W+5,align_S/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_W-10,align_S+1,align_W,align_S+2/),
     $              (/2,2/))

            case(19)

               !west buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(W)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-5,align_N-2,align_W,align_N-1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !north buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(N)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-10,align_N,align_W+5,align_N+1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_W-10,align_N-2,align_W,align_N-1/),
     $              (/2,2/))

            case(20)

               !west buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(W)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-5,align_S+1,align_W,align_N-1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr)

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr)

               bf_sublayer_tested => bf_sublayer_ptr

               !south buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(S)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-15,align_S-1,align_W+5,align_S/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              SW_interface_type,
     $              bf_sublayer_ptr)

               !north buffer layer
               allocate(bf_sublayer_ptr)
               call bf_sublayer_ptr%ini(N)
               call bf_sublayer_ptr%allocate_bf_layer(
     $              interior_x_map,
     $              interior_y_map,
     $              interior_nodes,
     $              reshape(
     $                (/align_W-10,align_N,align_W+5,align_N+1/),
     $                (/2,2/)))

               call mainlayer_interface_used%set_mainlayer_interface_bf_layer(
     $              NW_interface_type,
     $              bf_sublayer_ptr)

               test_should_be_reallocated = .true.
               test_bf_alignment = reshape(
     $              (/align_W-15,align_S+1,align_W,align_N-1/),
     $              (/2,2/))

            end select

        end subroutine get_param_test_uniformize


      end program test_mainlayer_interface_dyn
