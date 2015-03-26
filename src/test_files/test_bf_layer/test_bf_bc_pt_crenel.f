      program test_bf_bc_pt_crenel

        use bf_bc_pt_crenel_module, only :
     $     check_if_bc_pt_crenel,
     $     remove_bc_pt_crenel

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       no_pt,
     $       bc_pt

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled      = .true.
        test_validated = .true.


        test_loc = test_check_if_bc_pt_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_check_if_bc_pt_crenel: '',L1)', test_loc
        print '()'


        test_loc = test_remove_bc_pt_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_remove_bc_pt_crenel: '',L1)', test_loc
        print '()'

        print '(''test_validated: '',L1)', test_validated


        contains

        function test_remove_bc_pt_crenel(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(3,3) :: grdpts_id
          integer, dimension(3,3) :: test_grdpts_id

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,8

             !input
             call get_test_param_remove(
     $            k,
     $            grdpts_id,
     $            test_grdpts_id)

             !output
             call remove_bc_pt_crenel(
     $            grdpts_id,2,2)

             !validation
             test_loc = is_int_matrix_validated(
     $            grdpts_id,
     $            test_grdpts_id,
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')',k
             end if


          end do

        end function test_remove_bc_pt_crenel


        subroutine get_test_param_remove(
     $     test_id,
     $     grdpts_id,
     $     test_grdpts_id)

          implicit none

          integer                , intent(in)  :: test_id
          integer, dimension(3,3), intent(out) :: grdpts_id
          integer, dimension(3,3), intent(out) :: test_grdpts_id


          select case(test_id)
            case(1)
               grdpts_id = reshape((/
     $              3,2,2,
     $              3,3,2,
     $              3,3,3/),
     $              (/3,3/))

               test_grdpts_id = reshape((/
     $              3,2,1,
     $              3,2,2,
     $              3,3,3/),
     $              (/3,3/))

            case(2)
               grdpts_id = reshape((/
     $              2,2,3,
     $              2,3,3,
     $              3,3,3/),
     $              (/3,3/))

               test_grdpts_id = reshape((/
     $              1,2,3,
     $              2,2,3,
     $              3,3,3/),
     $              (/3,3/))

            case(3)
               grdpts_id = reshape((/
     $              3,3,3,
     $              3,3,2,
     $              3,2,2/),
     $              (/3,3/))

               test_grdpts_id = reshape((/
     $              3,3,3,
     $              3,2,2,
     $              3,2,1/),
     $              (/3,3/))

            case(4)
               grdpts_id = reshape((/
     $              3,3,3,
     $              2,3,3,
     $              2,2,3/),
     $              (/3,3/))

               test_grdpts_id = reshape((/
     $              3,3,3,
     $              2,2,3,
     $              1,2,3/),
     $              (/3,3/))

            case(5)
               grdpts_id = reshape((/
     $              2,2,3,
     $              2,3,3,
     $              2,2,3/),
     $              (/3,3/))

               test_grdpts_id = reshape((/
     $              1,2,3,
     $              1,2,3,
     $              1,2,3/),
     $              (/3,3/))

            case(6)
               grdpts_id = reshape((/
     $              3,2,2,
     $              3,3,2,
     $              3,2,2/),
     $              (/3,3/))

               test_grdpts_id = reshape((/
     $              3,2,1,
     $              3,2,1,
     $              3,2,1/),
     $              (/3,3/))

            case(7)
               grdpts_id = reshape((/
     $              2,2,2,
     $              2,3,2,
     $              3,3,3/),
     $              (/3,3/))

               test_grdpts_id = reshape((/
     $              1,1,1,
     $              2,2,2,
     $              3,3,3/),
     $              (/3,3/))

            case(8)
               grdpts_id = reshape((/
     $              3,3,3,
     $              2,3,2,
     $              2,2,2/),
     $              (/3,3/))

               test_grdpts_id = reshape((/
     $              3,3,3,
     $              2,2,2,
     $              1,1,1/),
     $              (/3,3/))

          end select

        end subroutine get_test_param_remove


        function test_check_if_bc_pt_crenel(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(3,3) :: grdpts_id
          integer, dimension(8)   :: grdpts_id_inserted

          integer :: i,j
          integer :: k


          test_validated = .true.

          
          grdpts_id = reshape((/
     $            ((bc_pt,i=1,3),j=1,3)/),
     $            (/3,3/))

          test_loc = check_if_bc_pt_crenel(grdpts_id,2,2)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test('',I2,'') failed'')',0
          end if


          do k=1,8
             
             !input
             grdpts_id_inserted = (/ (bc_pt, i=1,8) /)
             grdpts_id_inserted(k) = no_pt

             grdpts_id = reshape((/
     $            ((bc_pt,i=1,3),j=1,3)/),
     $            (/3,3/))

             grdpts_id(1:3,1) = grdpts_id_inserted(1:3)
             grdpts_id(1,2)   = grdpts_id_inserted(4)
             grdpts_id(3,2)   = grdpts_id_inserted(5)
             grdpts_id(1:3,3) = grdpts_id_inserted(6:8)

             !output+validation
             test_loc = .not.check_if_bc_pt_crenel(grdpts_id,2,2)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')',k
             end if

          end do

        end function test_check_if_bc_pt_crenel

      end program test_bf_bc_pt_crenel
