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
     $       no_pt

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

        
        contains

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

      end program test_bf_layer_grdpts_id_update
