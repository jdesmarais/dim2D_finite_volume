      program test_bf_layer_grdpts_id_update

        use bf_layer_grdpts_id_update_class, only :
     $     bf_layer_grdpts_id_update

        use check_data_module, only :
     $       is_int_matrix_validated,
     $       is_boolean_matrix_validated

        use parameters_bf_layer, only :
     $       interior_pt,
     $       bc_interior_pt,
     $       bc_pt,
     $       no_pt,
     $       BF_SUCCESS

        use parameters_constant, only :
     $       E

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled      = .true.
        test_validated = .true.


        test_loc = test_check_grdpts_id_pt(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_check_grdpts_id_pt: '',L1)', test_loc
        print '()'


        test_loc = test_update_grdpts_id(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_update_grdpts_id: '',L1)', test_loc
        print '()'        


        test_loc = test_finalize_update_grdpts_id(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_finalize_update_grdpts_id_around_new_interior_pt: '',L1)', test_loc
        print '()'


        test_loc = test_detect_bc_interior_pt_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_detect_bc_interior_pt_crenel: '',L1)', test_loc
        print '()'


        test_loc = test_can_interior_crenel_be_curbed(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_can_interior_crenel_be_curbed: '',L1)', test_loc
        print '()'


        test_loc = test_detect_and_curb_bc_pt_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_detect_and_curb_bc_pt_crenel: '',L1)', test_loc
        print '()'

        
        print '(''test_validated: '',L1)', test_validated

        contains

        function test_detect_and_curb_bc_pt_crenel(detailled)
     $       result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_layer_grdpts_id_update) :: bf_layer_used
          integer                         :: test_i
          integer                         :: test_j
          integer, dimension(5,5)         :: test_grdpts_id
          logical                         :: ierror
          logical                         :: test_ierror

          logical :: test_loc
          integer :: k


          test_validated = .true.
          

          do k=1,3

             !input
             call get_test_param_detect_crenel(
     $            k,
     $            bf_layer_used,
     $            test_i,
     $            test_j,
     $            test_grdpts_id,
     $            test_ierror)


             !output
             call bf_layer_used%detect_and_curb_bc_pt_crenel(
     $            test_i,test_j,ierror)
             

             !validation
             test_loc = ierror.eqv.test_ierror
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') ierror failed'')',k
             end if

             if(ierror.eqv.BF_SUCCESS) then
                test_loc = is_int_matrix_validated(
     $               bf_layer_used%grdpts_id,
     $               test_grdpts_id,
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test('',I2,'') grdpts_id failed'')'
                end if
             end if

          end do


        end function test_detect_and_curb_bc_pt_crenel


        subroutine get_test_param_detect_crenel(
     $     test_id,
     $     bf_layer_used,
     $     test_i,
     $     test_j,
     $     test_grdpts_id,
     $     test_ierror)

          implicit none

          integer                        , intent(in)    :: test_id
          type(bf_layer_grdpts_id_update), intent(inout) :: bf_layer_used
          integer                        , intent(out)   :: test_i
          integer                        , intent(out)   :: test_j
          integer, dimension(5,5)        , intent(out)   :: test_grdpts_id
          logical                        , intent(out)   :: test_ierror


          test_grdpts_id = reshape((/
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3,
     $              1,1,1,2,3/),
     $              (/5,5/))

          select case(test_id)
            case(1)
               bf_layer_used%localization=E
               allocate(bf_layer_used%grdpts_id(5,5))
               bf_layer_used%grdpts_id = reshape((/
     $              1,1,1,2,3,
     $              1,1,2,2,3,
     $              1,1,2,3,3,
     $              1,1,2,2,3,
     $              1,1,1,2,3/),
     $              (/5,5/))

               call bf_layer_used%set_neighbor1_share(.false.)
               call bf_layer_used%set_neighbor2_share(.false.)
               test_i = 4
               test_j = 3               
               test_ierror = BF_SUCCESS

            case(2)
               test_i = 5
               test_j = 2
               test_ierror = BF_SUCCESS

            case(3)
               call bf_layer_used%set_neighbor1_share(.true.)
               test_i = 5
               test_j = 1
               test_ierror = .not.BF_SUCCESS

          end select

        end subroutine get_test_param_detect_crenel


        function test_finalize_update_grdpts_id(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_layer_grdpts_id_update) :: bf_layer_used
          integer                         :: i_prev
          integer                         :: j_prev
          integer                         :: i
          integer                         :: j,k
          integer, dimension(11,11)       :: test_grdpts_id


          test_validated = .true.


          allocate(bf_layer_used%grdpts_id(11,11))

          i_prev = 6
          j_prev = 6

          do j=j_prev-4, j_prev+4
             do i=i_prev-4, i_prev+4

                !input
                call ini_grdpts_id(
     $               bf_layer_used%grdpts_id,
     $               i_prev,j_prev,
     $               i,j,
     $               test_grdpts_id)

                !output
                call bf_layer_used%finalize_update_grdpts_id_around_new_interior_pt(
     $               i_prev, j_prev,
     $               i,j)                               


                !validation
                test_loc = is_int_matrix_validated(
     $               bf_layer_used%grdpts_id,
     $               test_grdpts_id,
     $               detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test('',2I3,'') failed'')', i-i_prev,j-j_prev

                   print '(''grdpts_id'')'
                   do k=1,11
                      print '(11I2)', bf_layer_used%grdpts_id(:,11-(k-1))
                   end do
                   print '()'
                   
                   print '(''test_grdpts_id'')'
                   do k=1,11
                      print '(11I2)', test_grdpts_id(:,11-(k-1))
                   end do
                   print '()' 
                end if                

             end do
          end do

        end function test_finalize_update_grdpts_id


        subroutine ini_grdpts_id(
     $     grdpts_id,
     $     i_prev,j_prev,
     $     i_center,j_center,
     $     test_grdpts_id)

          implicit none

          integer, dimension(11,11), intent(out) :: grdpts_id
          integer                  , intent(in)  :: i_prev
          integer                  , intent(in)  :: j_prev
          integer                  , intent(in)  :: i_center
          integer                  , intent(in)  :: j_center
          integer, dimension(11,11), intent(out) :: test_grdpts_id


          integer :: i,j


          test_grdpts_id = reshape((/
     $         ((bc_pt, i=1,11),j=1,11)/),
     $         (/11,11/))

          do j=j_prev-2,j_prev+2
             do i=i_prev-2,i_prev+2
                test_grdpts_id(i,j) = no_pt
             end do
          end do

          do j=j_prev-1,j_prev+1
             do i=i_prev-1,i_prev+1
                test_grdpts_id(i,j) = interior_pt
             end do
          end do

          do j=j_center-1,j_center+1
             do i=i_center-1,i_center+1
                if(test_grdpts_id(i,j).eq.no_pt) then
                   test_grdpts_id(i,j) = bc_interior_pt
                end if
             end do
          end do

          do j=j_prev-2,j_prev+2
             do i=i_prev-2,i_prev+2
                if(test_grdpts_id(i,j).eq.no_pt) then
                   test_grdpts_id(i,j) = bc_pt
                end if
             end do
          end do


          grdpts_id = reshape((/
     $         ((bc_pt, i=1,11),j=1,11)/),
     $         (/11,11/))

          do j=j_prev-1,j_prev+1
             do i=i_prev-1,i_prev+1
                grdpts_id(i,j) = interior_pt
             end do
          end do

        end subroutine ini_grdpts_id


        function test_update_grdpts_id(detailled)
     $     result(test_validated)

          implicit none
          
          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_grdpts_id_update) :: bf_layer_used
          logical , dimension(5,5)        :: compute_newgrdpt
          integer , dimension(5,5)        :: test_grdpts_id
          logical , dimension(5,5)        :: test_compute_newgrdpt

          integer :: i,j,k
          logical :: test_loc


          test_validated = .true.


          allocate(bf_layer_used%grdpts_id(5,5))

          do k=1,3

             !input
             call ini_test_update_grdpts_id(
     $            k,
     $            bf_layer_used%grdpts_id,
     $            test_grdpts_id,
     $            test_compute_newgrdpt)

             !output
             do j=1,5
                do i=1,5
                   compute_newgrdpt(i,j) =
     $                  bf_layer_used%update_grdpt_id_next_to_new_interior_pt(
     $                  i,j,3,3)
                end do
             end do

             !validation
             test_loc = is_int_matrix_validated(
     $            bf_layer_used%grdpts_id,
     $            test_grdpts_id,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') grdpts_id failed'')'
             end if

             test_loc = is_boolean_matrix_validated(
     $            compute_newgrdpt,
     $            test_compute_newgrdpt,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') compute_newgrdpt failed'')'
             end if             

          end do

        end function test_update_grdpts_id


        subroutine ini_test_update_grdpts_id(
     $     test_id,
     $     grdpts_id,
     $     test_grdpts_id,
     $     test_compute_newgrdpt)

          implicit none

          integer                , intent(in)  :: test_id
          integer, dimension(5,5), intent(out) :: grdpts_id
          integer, dimension(5,5), intent(out) :: test_grdpts_id
          logical, dimension(5,5), intent(out) :: test_compute_newgrdpt

          integer :: i,j


          select case(test_id)
            case(1)
               grdpts_id = reshape((/
     $              ((no_pt,i=1,5),j=1,5)/),
     $              (/5,5/))

               test_grdpts_id = reshape((/
     $              3,3,3,3,3,
     $              3,2,2,2,3,
     $              3,2,2,2,3,
     $              3,2,2,2,3,
     $              3,3,3,3,3/),
     $              (/5,5/))

               test_compute_newgrdpt = reshape((/
     $              ((.true.,i=1,5),j=1,5)/),
     $              (/5,5/))

            case(2)
               grdpts_id = reshape((/
     $              ((bc_pt,i=1,5),j=1,5)/),
     $              (/5,5/))

               test_grdpts_id = reshape((/
     $              3,3,3,3,3,
     $              3,2,2,2,3,
     $              3,2,2,2,3,
     $              3,2,2,2,3,
     $              3,3,3,3,3/),
     $              (/5,5/))

               test_compute_newgrdpt = reshape((/
     $              ((.false.,i=1,5),j=1,5)/),
     $              (/5,5/))

            case(3)
               grdpts_id = reshape((/
     $              ((bc_interior_pt,i=1,5),j=1,5)/),
     $              (/5,5/))

               test_grdpts_id = grdpts_id

               test_compute_newgrdpt = reshape((/
     $              ((.false.,i=1,5),j=1,5)/),
     $              (/5,5/))

          end select

        end subroutine ini_test_update_grdpts_id


        function test_check_grdpts_id_pt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_layer_grdpts_id_update) :: bf_layer_used

          allocate(bf_layer_used%grdpts_id(1,1))

          bf_layer_used%grdpts_id(1,1) = bc_pt

          test_validated = bf_layer_used%check_grdpts_id_pt(1,1,bc_pt)
          if(detailled.and.(.not.test_validated)) then
             print '(I2,'' -> '',I2)', bf_layer_used%grdpts_id(1,1), bc_pt
          end if

        end function test_check_grdpts_id_pt


        function test_detect_bc_interior_pt_crenel(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_grdpts_id_update) :: bf_layer_used
          
          integer :: test_i
          integer :: test_j
          
          logical :: test_detect
          logical :: test_ierror

          logical :: detect
          logical :: ierror

          integer :: k


          test_validated = .true.


          do k=1,3
          
             !input
             call get_param_detect(
     $            k,
     $            bf_layer_used,
     $            test_i, test_j,
     $            test_detect,
     $            test_ierror)

             !output
             detect = bf_layer_used%detect_bc_interior_pt_crenel(
     $            test_i,test_j,ierror)

             !validation
             test_loc = ierror.eqv.test_ierror
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test ierror('',I2,'') failed'')', k
             end if

             if(ierror.eqv.BF_SUCCESS) then
                test_loc = detect.eqv.test_detect
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test detect('',I2,'') failed'')', k
                end if
             end if

          end do

        end function test_detect_bc_interior_pt_crenel


        subroutine get_param_detect(
     $     test_id,
     $     bf_layer_used,
     $     test_i, test_j,
     $     test_detect,
     $     test_ierror)

          implicit none

          integer                        , intent(in)    :: test_id
          type(bf_layer_grdpts_id_update), intent(inout) :: bf_layer_used
          integer                        , intent(out)   :: test_i
          integer                        , intent(out)   :: test_j
          logical                        , intent(out)   :: test_detect
          logical                        , intent(out)   :: test_ierror


          select case(test_id)
            case(1)
               bf_layer_used%localization = E
               allocate(bf_layer_used%grdpts_id(5,5))
               bf_layer_used%grdpts_id = reshape((/
     $              1,1,2,3,0,
     $              1,1,2,3,0,
     $              1,1,2,3,0,
     $              1,1,2,3,0,
     $              1,1,2,3,0/),
     $              (/5,5/))

               test_i = 3
               test_j = 3
               
               test_detect = .false.
               test_ierror = BF_SUCCESS

            case(2)
               bf_layer_used%localization = E
               call bf_layer_used%set_neighbor1_share(.false.)
               call bf_layer_used%set_neighbor2_share(.false.)

               test_i = 5
               test_j = 3

               test_detect = .false.
               test_ierror = BF_SUCCESS

            case(3)
               call bf_layer_used%set_neighbor1_share(.true.)

               test_i = 3
               test_j = 1

               test_ierror = .not.BF_SUCCESS

          end select

        end subroutine get_param_detect


        function test_can_interior_crenel_be_curbed(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer_grdpts_id_update) :: bf_layer_used
          
          integer :: test_i
          integer :: test_j
          
          logical :: test_can_be_curbed
          logical :: test_ierror

          logical :: can_be_curbed
          logical :: ierror

          integer :: k


          test_validated = .true.


          do k=1,2
          
             !input
             call get_param_curb(
     $            k,
     $            bf_layer_used,
     $            test_i, test_j,
     $            test_can_be_curbed,
     $            test_ierror)

             !output
             can_be_curbed = bf_layer_used%can_bc_interior_pt_crenel_be_curbed(
     $            test_i,test_j,ierror)

             !validation
             test_loc = ierror.eqv.test_ierror
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test ierror('',I2,'') failed'')', k
             end if

             if(ierror.eqv.BF_SUCCESS) then
                test_loc = can_be_curbed.eqv.test_can_be_curbed
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test can_be_curbed('',I2,'') failed'')', k
                end if
             end if

          end do

        end function test_can_interior_crenel_be_curbed


        subroutine get_param_curb(
     $     test_id,
     $     bf_layer_used,
     $     test_i, test_j,
     $     test_can_be_curbed,
     $     test_ierror)

          implicit none

          integer                        , intent(in)    :: test_id
          type(bf_layer_grdpts_id_update), intent(inout) :: bf_layer_used
          integer                        , intent(out)   :: test_i
          integer                        , intent(out)   :: test_j
          logical                        , intent(out)   :: test_can_be_curbed
          logical                        , intent(out)   :: test_ierror


          select case(test_id)
            case(1)
               bf_layer_used%localization = E
               allocate(bf_layer_used%grdpts_id(5,5))
               bf_layer_used%grdpts_id = reshape((/
     $              1,1,2,3,0,
     $              1,1,2,3,0,
     $              1,1,2,3,0,
     $              1,1,2,3,0,
     $              1,1,2,3,0/),
     $              (/5,5/))

               test_i = 3
               test_j = 3
               
               test_can_be_curbed = .false.
               test_ierror = BF_SUCCESS

            case(2)
               test_i = 4
               test_j = 3

               test_can_be_curbed = .false.
               test_ierror = .not.BF_SUCCESS

          end select

        end subroutine get_param_curb

      end program test_bf_layer_grdpts_id_update
