      program test_bf_mainlayer_basic

        use bf_mainlayer_basic_class, only :
     $       bf_mainlayer_basic

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       align_E,align_W,
     $       align_N,align_S

        use parameters_constant, only :
     $       N,S,E,W,
     $       x_direction,
     $       y_direction

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_input, only :
     $       rkind

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        test_loc = test_ini(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_ini: '',L1)', test_loc
        print '()'


        test_loc = test_add_sublayer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_add_sublayer: '',L1)', test_loc
        print '()'


        test_loc = test_remove_sublayer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_remove_sublayer: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_ini(detailled)
     $       result(test_validated)

          implicit none
          
          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_mainlayer_basic) :: bf_mainlayer_used
          logical :: test_loc


          test_validated = .true.

          call bf_mainlayer_used%ini(N)

          test_loc = bf_mainlayer_used%mainlayer_id.eq.N
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''mainlayer_id failed'')'
          end if

          test_loc = bf_mainlayer_used%nb_sublayers.eq.0
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_sublayers failed'')'
          end if
          
          test_loc = .not.associated(bf_mainlayer_used%head_sublayer)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''head_sublayer failed'')'
          end if

          test_loc = .not.associated(bf_mainlayer_used%tail_sublayer)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''tail_sublayer failed'')'
          end if

        end function test_ini


        function test_add_sublayer(detailled)
     $       result(test_validated)

          implicit none
          
          logical, intent(in) :: detailled
          logical             :: test_validated
          
          integer :: k,l
          logical :: test_loc

          test_validated = .true.

          do l=1,4
             do k=1,6
                test_loc = perform_test_add_sublayer(l,k,detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test('',I2,'') failed'')'
                end if
             end do
          end do

        end function test_add_sublayer


        function perform_test_add_sublayer(
     $     mainlayer_id,
     $     test_id,
     $     detailled)
     $     result(test_validated)

          integer, intent(in) :: mainlayer_id
          integer, intent(in) :: test_id
          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_mainlayer_basic)  :: bf_mainlayer_used

          integer                   :: nb_sublayers_test
          integer, dimension(2,2,3) :: bf_alignment_added
          integer, dimension(2,2,3) :: bf_alignment_tested

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          type(bf_sublayer), pointer :: current_sublayer

          integer :: k

          !input
          call bf_mainlayer_used%ini(mainlayer_id)

          call get_param_test_add_sublayer(
     $         mainlayer_id,
     $         test_id,
     $         nb_sublayers_test,
     $         bf_alignment_added,
     $         bf_alignment_tested)

          !output
          do k=1, nb_sublayers_test

             current_sublayer => bf_mainlayer_used%add_sublayer(
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            bf_alignment_added(:,:,k))

          end do

          !validation
          test_validated = test_content_bf_mainlayer(
     $         bf_mainlayer_used,
     $         nb_sublayers_test,
     $         bf_alignment_tested,
     $         detailled)

        end function perform_test_add_sublayer

      
        function test_content_bf_mainlayer(
     $     bf_mainlayer_used,
     $     nb_sublayers_test,
     $     bf_alignment_tested,
     $     detailled)
     $     result(test_validated)

          implicit none

          type(bf_mainlayer_basic) , intent(in) :: bf_mainlayer_used
          integer                  , intent(in) :: nb_sublayers_test
          integer, dimension(:,:,:), intent(in) :: bf_alignment_tested
          logical                  , intent(in) :: detailled
          logical                               :: test_validated

          type(bf_sublayer), pointer :: current_sublayer
          integer                    :: nb_sublayers
          integer, dimension(2,2)    :: bf_alignment
          integer                    :: k
          
          logical :: test_loc

          
          test_validated = .true.


          !nb_sublayers
          nb_sublayers = bf_mainlayer_used%get_nb_sublayers()
          test_loc = nb_sublayers.eq.nb_sublayers_test
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nb_sublayers failed'')'
          end if

          !alignment
          if(nb_sublayers>0) then

             current_sublayer => bf_mainlayer_used%get_head_sublayer()
             do k=1, nb_sublayers

                bf_alignment = current_sublayer%get_alignment_tab()

                test_loc = is_int_matrix_validated(
     $               bf_alignment,
     $               bf_alignment_tested(:,:,k),
     $               detailled)
                test_validated = test_validated.and.test_loc

                if(detailled.and.(.not.test_loc)) then
                   print '(''test failed for sublayer '',I2)',k
                end if

                current_sublayer => current_sublayer%get_next()
             end do          

          end if

        end function test_content_bf_mainlayer


        subroutine get_param_test_add_sublayer(
     $     mainlayer_id,
     $     test_id,
     $     nb_sublayers,
     $     bf_alignment_added,
     $     bf_alignment_tested)

          implicit none

          integer                  , intent(in)  :: mainlayer_id
          integer                  , intent(in)  :: test_id
          integer                  , intent(out) :: nb_sublayers
          integer, dimension(2,2,3), intent(out) :: bf_alignment_added
          integer, dimension(2,2,3), intent(out) :: bf_alignment_tested

          integer :: k

          select case(test_id)
            case(1)
               nb_sublayers = 1
               bf_alignment_added(:,:,1) = reshape((/
     $              1,1,3,3/),
     $              (/2,2/))
               bf_alignment_tested = bf_alignment_added

            case(2)
               nb_sublayers = 2
               bf_alignment_added = reshape((/
     $              1,1,3,3,
     $              4,4,5,5,
     $              0,0,0,0/),
     $              (/2,2,3/))
               bf_alignment_tested = bf_alignment_added

            case(3)
               nb_sublayers = 2
               bf_alignment_added = reshape((/
     $              4,4,5,5,
     $              1,1,3,3,
     $              0,0,0,0/),
     $              (/2,2,3/))
               bf_alignment_tested = reshape((/
     $              1,1,3,3,
     $              4,4,5,5,
     $              0,0,0,0/),
     $              (/2,2,3/))

            case(4)
               nb_sublayers = 3
               bf_alignment_added = reshape((/
     $              1,1,3,3,
     $              4,4,5,5,
     $              6,6,7,7/),
     $              (/2,2,3/))
               bf_alignment_tested = bf_alignment_added

            case(5)
               nb_sublayers = 3
               bf_alignment_added = reshape((/
     $              1,1,3,3,
     $              6,6,7,7,
     $              4,4,5,5/),
     $              (/2,2,3/))
               bf_alignment_tested = reshape((/
     $              1,1,3,3,
     $              4,4,5,5,
     $              6,6,7,7/),
     $              (/2,2,3/))

            case(6)
               nb_sublayers = 3
               bf_alignment_added = reshape((/
     $              6,6,7,7,
     $              4,4,5,5,
     $              1,1,3,3/),
     $              (/2,2,3/))
               bf_alignment_tested = reshape((/
     $              1,1,3,3,
     $              4,4,5,5,
     $              6,6,7,7/),
     $              (/2,2,3/))

          end select

          select case(mainlayer_id)

            case(N)
               do k=1, nb_sublayers
                  bf_alignment_added(2,1,k) = align_N
                  bf_alignment_added(2,2,k) = align_N+10
                  bf_alignment_tested(2,1,k) = align_N
                  bf_alignment_tested(2,2,k) = align_N+10
               end do

            case(S)
               do k=1, nb_sublayers
                  bf_alignment_added(2,1,k) = align_S-10
                  bf_alignment_added(2,2,k) = align_S
                  bf_alignment_tested(2,1,k) = align_S-10
                  bf_alignment_tested(2,2,k) = align_S
               end do

            case(W)
               do k=1, nb_sublayers
                  bf_alignment_added(1,1,k) = align_W-10
                  bf_alignment_added(1,2,k) = align_W
                  bf_alignment_tested(1,1,k) = align_W-10
                  bf_alignment_tested(1,2,k) = align_W
               end do

            case(E)
               do k=1, nb_sublayers
                  bf_alignment_added(1,1,k) = align_E
                  bf_alignment_added(1,2,k) = align_E+10
                  bf_alignment_tested(1,1,k) = align_E
                  bf_alignment_tested(1,2,k) = align_E+10
               end do

          end select

        end subroutine get_param_test_add_sublayer


        function test_remove_sublayer(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated
          
          integer :: k
          logical :: test_loc   

          test_validated = .true.

          do k=1,6

             test_loc = perform_test_remove_sublayer(k,detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_remove_sublayer


        function perform_test_remove_sublayer(test_id,detailled)
     $     result(test_validated)

          implicit none

          integer, intent(in) :: test_id
          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_mainlayer_basic)  :: bf_mainlayer_used
          integer, dimension(2,2,3) :: bf_alignment1
          integer, dimension(2,2,2) :: bf_alignment2
          integer, dimension(3)     :: remove_order

          integer, dimension(2,2,3) :: bf_alignment_added
          
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          type(bf_sublayer), pointer :: bf_sublayer_ptr1
          type(bf_sublayer), pointer :: bf_sublayer_ptr2
          type(bf_sublayer), pointer :: bf_sublayer_ptr3


          test_validated = .true.

          
          bf_alignment_added = reshape((/
     $              1,align_N,3,align_N,
     $              4,align_N,5,align_N,
     $              6,align_N,7,align_N/),
     $              (/2,2,3/))

          call bf_mainlayer_used%ini(N)

          bf_sublayer_ptr1 => bf_mainlayer_used%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_added(:,:,1))

          bf_sublayer_ptr2 => bf_mainlayer_used%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_added(:,:,2))

          bf_sublayer_ptr3 => bf_mainlayer_used%add_sublayer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_added(:,:,3))

          call get_param_test_remove_sublayer(
     $         test_id,
     $         bf_alignment1,
     $         bf_alignment2,
     $         remove_order)

          !remove first
          select case(remove_order(1))
            case(1)
               call bf_mainlayer_used%remove_sublayer(bf_sublayer_ptr1)
            case(2)
               call bf_mainlayer_used%remove_sublayer(bf_sublayer_ptr2)
            case(3)
               call bf_mainlayer_used%remove_sublayer(bf_sublayer_ptr3)
          end select
               
          !test what is left
          test_loc = test_content_bf_mainlayer(
     $         bf_mainlayer_used,
     $         2,
     $         bf_alignment1,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''remove first failed'')'
          end if


          !remove second
          select case(remove_order(2))
            case(1)
               call bf_mainlayer_used%remove_sublayer(bf_sublayer_ptr1)
            case(2)
               call bf_mainlayer_used%remove_sublayer(bf_sublayer_ptr2)
            case(3)
               call bf_mainlayer_used%remove_sublayer(bf_sublayer_ptr3)
          end select
          
          !test what is left
          test_loc = test_content_bf_mainlayer(
     $         bf_mainlayer_used,
     $         1,
     $         bf_alignment2,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''remove second failed'')'
          end if

          !remove third
          select case(remove_order(3))
            case(1)
               call bf_mainlayer_used%remove_sublayer(bf_sublayer_ptr1)
            case(2)
               call bf_mainlayer_used%remove_sublayer(bf_sublayer_ptr2)
            case(3)
               call bf_mainlayer_used%remove_sublayer(bf_sublayer_ptr3)
          end select
          
          !test what is left
          test_loc = test_content_bf_mainlayer(
     $         bf_mainlayer_used,
     $         0,
     $         bf_alignment2,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''remove third failed'')'
          end if

        end function perform_test_remove_sublayer


        subroutine get_param_test_remove_sublayer(
     $     test_id,
     $     bf_alignment1,
     $     bf_alignment2,
     $     remove_order)

          implicit none

          integer                  , intent(in)  :: test_id
          integer, dimension(2,2,2), intent(out) :: bf_alignment1
          integer, dimension(2,2,1), intent(out) :: bf_alignment2
          integer, dimension(3)    , intent(out) :: remove_order

          select case(test_id)

            case(1)
               bf_alignment1 = reshape((/
     $              4,align_N,5,align_N,
     $              6,align_N,7,align_N/),
     $              (/2,2,2/))

               bf_alignment2 = reshape((/
     $              6,align_N,7,align_N/),
     $              (/2,2,1/))

               remove_order = [1,2,3]

            case(2)
               bf_alignment1 = reshape((/
     $              4,align_N,5,align_N,
     $              6,align_N,7,align_N/),
     $              (/2,2,2/))

               bf_alignment2 = reshape((/
     $              4,align_N,5,align_N/),
     $              (/2,2,1/))

               remove_order = [1,3,2]

            case(3)
               bf_alignment1 = reshape((/
     $              1,align_N,3,align_N,
     $              6,align_N,7,align_N/),
     $              (/2,2,2/))

               bf_alignment2 = reshape((/
     $              6,align_N,7,align_N/),
     $              (/2,2,1/))

               remove_order = [2,1,3]

            case(4)
               bf_alignment1 = reshape((/
     $              1,align_N,3,align_N,
     $              6,align_N,7,align_N/),
     $              (/2,2,2/))

               bf_alignment2 = reshape((/
     $              1,align_N,3,align_N/),
     $              (/2,2,1/))

               remove_order = [2,3,1]

            case(5)
               bf_alignment1 = reshape((/
     $              1,align_N,3,align_N,
     $              4,align_N,5,align_N/),
     $              (/2,2,2/))

               bf_alignment2 = reshape((/
     $              4,align_N,5,align_N/),
     $              (/2,2,1/))

               remove_order = [3,1,2]

            case(6)
               bf_alignment1 = reshape((/
     $              1,align_N,3,align_N,
     $              4,align_N,5,align_N/),
     $              (/2,2,2/))

               bf_alignment2 = reshape((/
     $              1,align_N,3,align_N/),
     $              (/2,2,1/))

               remove_order = [3,2,1]

          end select

        end subroutine get_param_test_remove_sublayer

      end program test_bf_mainlayer_basic
