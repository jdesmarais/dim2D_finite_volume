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
      program test_dcr_interface

        use dcr_interface_class, only :
     $       dcr_interface

        use bf_interface_icr_class, only :
     $       bf_interface_icr

        use bf_sublayer_class, only :
     $       bf_sublayer

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


        test_loc = test_ini(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_ini: '',L1)', test_loc
        print '()'


        test_loc = test_not_in_no_check_list(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_not_in_no_check_list: '',L1)', test_loc
        print '()'

        
        test_loc = test_prevent_neighbor_removal(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_prevent_neighbor_removal: '',L1)', test_loc
        print '()'


        test_loc = test_stage(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_stage: '',L1)', test_loc
        print '()'


        test_loc = test_check_if_neighbors_remain(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_check_if_neighbors_remain: '',L1)', test_loc
        print '()'


        test_loc = test_finalize_domain_decrease(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_finalize_domain_decrease: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_ini(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(dcr_interface)    :: dcr_interface_used
          type(bf_interface_icr) :: bf_interface_used
          integer, dimension(4)  :: test_size

          integer :: k
          logical :: test_loc


          test_validated =.true.


          !intialize the bf_interface
          !------------------------------------------------------------
          ! - 3 North buffer layers
          ! - 3 South buffer layers
          ! - 2 East buffer layers
          ! - 2 West buffer layers
          !------------------------------------------------------------
          call ini_bf_interface_for_tests(bf_interface_used)
          

          !input
          !------------------------------------------------------------
          test_size(N) = 3
          test_size(S) = 3
          test_size(E) = 2
          test_size(W) = 2

          
          !output
          !------------------------------------------------------------
          call dcr_interface_used%ini(bf_interface_used)


          !validation
          !------------------------------------------------------------
          !check that the number of buffer layers in each main layer
          !matches the size for the references in the lists
          !no_check_list and double_check_list
          do k=1,4

             test_loc = allocated(dcr_interface_used%no_check_list(k)%list)
             test_validated = test_validated.and.test_loc
             if(test_loc) then
                test_loc = size(dcr_interface_used%no_check_list(k)%list,1).eq.test_size(k)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test size(1,'',I2,'') failed'')',k
                end if
             else
                if(detailled) then
                   print '(''test allocation(1,'',I2,'') failed'')',k
                end if
             end if

             test_loc = allocated(dcr_interface_used%no_check_list(k)%list)
             test_validated = test_validated.and.test_loc
             if(test_loc) then
                test_loc = size(dcr_interface_used%no_check_list(k)%list,1).eq.test_size(k)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test size(2,'',I2,'') failed'')',k
                end if
             else
                if(detailled) then
                   print '(''test allocation(2,'',I2,'') failed'')',k
                end if
             end if

          end do
          
        end function test_ini


        
        subroutine ini_bf_interface_for_tests(bf_interface_used)

          implicit none

          type(bf_interface_icr), intent(inout) :: bf_interface_used

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


        function test_not_in_no_check_list(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(dcr_interface)        :: dcr_interface_used
          type(bf_interface_icr)     :: bf_interface_used
          type(bf_sublayer), pointer :: bf_sublayer1_ptr
          type(bf_sublayer), pointer :: bf_sublayer2_ptr
          type(bf_sublayer), pointer :: bf_sublayer3_ptr

          logical :: test_loc


          test_validated = .true.


          !input
          !------------------------------------------------------------
          call ini_bf_interface_for_tests(bf_interface_used)
          call dcr_interface_used%ini(bf_interface_used)

          allocate(bf_sublayer1_ptr)
          allocate(bf_sublayer2_ptr)
          allocate(bf_sublayer3_ptr)
          
          call dcr_interface_used%no_check_list(N)%add_ele(bf_sublayer1_ptr)
          call dcr_interface_used%no_check_list(N)%add_ele(bf_sublayer2_ptr)
          

          !output+validation
          !------------------------------------------------------------
          test_loc = .not.(dcr_interface_used%not_in_no_check_list(N,bf_sublayer1_ptr))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''not_in_no_check_list(1) failed'')'
          end if

          test_loc = .not.(dcr_interface_used%not_in_no_check_list(N,bf_sublayer2_ptr))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''not_in_no_check_list(2) failed'')'
          end if

          test_loc = dcr_interface_used%not_in_no_check_list(S,bf_sublayer2_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''not_in_no_check_list(3) failed'')'
          end if

          test_loc = dcr_interface_used%not_in_no_check_list(N,bf_sublayer3_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''not_in_no_check_list(4) failed'')'
          end if


          !cleaning
          !------------------------------------------------------------
          deallocate(bf_sublayer1_ptr)
          deallocate(bf_sublayer2_ptr)
          deallocate(bf_sublayer3_ptr)


        end function test_not_in_no_check_list


        function test_prevent_neighbor_removal(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          type(dcr_interface)        :: dcr_interface_used
          type(bf_interface_icr)     :: bf_interface_used
          type(bf_sublayer), pointer :: bf_sublayer_ptr

          logical :: test_loc

          test_validated = .true.


          !input
          !------------------------------------------------------------
          call ini_bf_interface_for_tests(bf_interface_used)
          call dcr_interface_used%ini(bf_interface_used)

          
          !output
          !------------------------------------------------------------
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%get_tail_sublayer()
          call bf_sublayer_ptr%set_remain_status(.false.)

          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%get_tail_sublayer()
          call bf_sublayer_ptr%set_remain_status(.false.)

          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%get_head_sublayer()
          call dcr_interface_used%prevent_neighbor_removal(
     $         bf_interface_used,
     $         N,bf_sublayer_ptr)

          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%get_tail_sublayer()
          call dcr_interface_used%prevent_neighbor_removal(
     $         bf_interface_used,
     $         N,bf_sublayer_ptr)


          !validation
          !------------------------------------------------------------
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%get_tail_sublayer()
          test_loc = bf_sublayer_ptr%get_remain_status()
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''removal status(W) failed'')'
          end if

          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%get_tail_sublayer()
          test_loc = bf_sublayer_ptr%get_remain_status()
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''removal status(E) failed'')'
          end if

        end function test_prevent_neighbor_removal


        function test_stage(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          type(dcr_interface)        :: dcr_interface_used
          type(bf_interface_icr)     :: bf_interface_used
          type(bf_sublayer), pointer :: bf_sublayer_ptr

          logical :: test_loc

          test_validated = .true.


          !input
          !------------------------------------------------------------
          call ini_bf_interface_for_tests(bf_interface_used)
          call dcr_interface_used%ini(bf_interface_used)


          !test removal of a buffer layer that does not depends on
          !neighbors
          !============================================================

          !output
          !------------------------------------------------------------
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%get_head_sublayer()
          bf_sublayer_ptr => bf_sublayer_ptr%get_next()

          call dcr_interface_used%stage(
     $         bf_interface_used,
     $         N,
     $         bf_sublayer_ptr)


          !validation
          !------------------------------------------------------------
          test_loc = bf_interface_used%mainlayer_pointers(N)%get_nb_sublayers().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''remove w/o neighbors failed'')'
          end if


          !test removal of a buffer layer that depends on neighbors
          !============================================================
          !output
          !------------------------------------------------------------
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%get_head_sublayer()
          bf_sublayer_ptr => bf_sublayer_ptr%get_next()

          call dcr_interface_used%stage(
     $         bf_interface_used,
     $         W,
     $         bf_sublayer_ptr)


          !validation
          !------------------------------------------------------------
          test_loc = bf_interface_used%mainlayer_pointers(W)%get_nb_sublayers().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''remove w/ neighbors, nb_sublayers failed'')'
          end if

          !veirfy that double_check contains bf_sublayer_ptr
          test_loc = .not.(dcr_interface_used%double_check_list(W)%does_not_contain(bf_sublayer_ptr))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''remove w/ neighbors, double_check failed'')'
          end if


        end function test_stage


        function test_check_if_neighbors_remain(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          type(dcr_interface)              :: dcr_interface_used
          type(bf_interface_icr)           :: bf_interface_used
          type(bf_sublayer), pointer       :: bf_sublayer_ptr
          type(pmodel_eq)                  :: p_model
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

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
          
          ! we will compute the remain status of the sublayers
          ! north 1 and 3 and south 1 and 3
          !
          ! ask the buffer layers W2 and E2 whether they have
          ! neighbors that remain: for W2, N1 is remaining; for E2,
          ! N3 does not remain
          !
          ! ask the buffer layers W1 and E1 whether they have
          ! neighbors that remain: for W1, S1 is not remaining; for E1,
          ! S1 does not remain

          !input
          !------------------------------------------------------------
          call ini_bf_interface_for_tests(bf_interface_used)
          call dcr_interface_used%ini(bf_interface_used)

          interior_nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,nx),j=1,ny)/), (/nx,ny/))
          
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


          !output+validation
          !------------------------------------------------------------
          !W1
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%get_head_sublayer()
          test_loc = .not.(dcr_interface_used%check_if_neighbors_remain(
     $         bf_interface_used,
     $         W,
     $         bf_sublayer_ptr))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test W1 failed'')'
          end if

          !W2
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%get_tail_sublayer()
          test_loc = dcr_interface_used%check_if_neighbors_remain(
     $         bf_interface_used,
     $         W,
     $         bf_sublayer_ptr)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test W2 failed'')'
          end if

          !E1
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%get_head_sublayer()
          test_loc = .not.(dcr_interface_used%check_if_neighbors_remain(
     $         bf_interface_used,
     $         E,
     $         bf_sublayer_ptr))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test E1 failed'')'
          end if

          !E2
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%get_tail_sublayer()
          test_loc = .not.(dcr_interface_used%check_if_neighbors_remain(
     $         bf_interface_used,
     $         E,
     $         bf_sublayer_ptr))
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test E2 failed'')'
          end if

        end function test_check_if_neighbors_remain


        function test_finalize_domain_decrease(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          type(dcr_interface)              :: dcr_interface_used
          type(bf_interface_icr)           :: bf_interface_used
          type(bf_sublayer), pointer       :: bf_sublayer_ptr
          type(pmodel_eq)                  :: p_model
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

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
          ! buffer layers S3 and N1
          !
          ! we stage the buffer layers W2 and E1 for
          ! removal, b/c of their neighbor dependencies
          ! they should be put in the double_check_list
          !
          ! we then finalize_domain_decrease, S3 should be
          ! removed since E1 an be removed and W2 should
          ! remain as N1 cannot be removed
          !============================================================
          !input
          !============================================================
          call ini_bf_interface_for_tests(bf_interface_used)
          call dcr_interface_used%ini(bf_interface_used)

          interior_nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,nx),j=1,ny)/), (/nx,ny/))
          
          !compute remain status
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

          !stage W2 and E1 for removal
          !------------------------------------------------------------
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%get_tail_sublayer()
          call dcr_interface_used%stage(bf_interface_used,W,bf_sublayer_ptr)

          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%get_head_sublayer()
          call dcr_interface_used%stage(bf_interface_used,E,bf_sublayer_ptr)
        

          !output+validation
          !============================================================
          call dcr_interface_used%finalize_domain_decrease(bf_interface_used)

          !W2
          test_loc = bf_interface_used%mainlayer_pointers(W)%get_nb_sublayers().eq.2
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test W2 failed'')'
          end if

          !E1
          test_loc = bf_interface_used%mainlayer_pointers(E)%get_nb_sublayers().eq.1
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test E1 failed'')'
          end if

          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%get_head_sublayer()
          test_loc = bf_sublayer_ptr%can_exchange_with_neighbor2()
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test removal E1 failed'')'
          end if
          

        end function test_finalize_domain_decrease


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

      end program test_dcr_interface
