      program test_bf_bc_interior_pt_crenel

        use bf_bc_interior_pt_crenel_module, only :
     $     check_if_bc_interior_pt_crenel

        use parameters_bf_layer, only :
     $       bc_interior_pt,
     $       bc_pt

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled      = .true.
        test_validated = .true.


        test_loc = test_check_if_bc_interior_pt_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_check_if_bc_interior_pt_crenel: '',L1)', test_loc
        print '()'


        contains


        function test_check_if_bc_interior_pt_crenel(detailled)
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
     $            ((bc_interior_pt,i=1,3),j=1,3)/),
     $            (/3,3/))

          test_loc = check_if_bc_interior_pt_crenel(grdpts_id,2,2)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test('',I2,'') failed'')',0
          end if


          do k=1,8
             
             !input
             grdpts_id_inserted = (/ (bc_interior_pt, i=1,8) /)
             grdpts_id_inserted(k) = bc_pt

             grdpts_id = reshape((/
     $            ((bc_interior_pt,i=1,3),j=1,3)/),
     $            (/3,3/))

             grdpts_id(1:3,1) = grdpts_id_inserted(1:3)
             grdpts_id(1,2)   = grdpts_id_inserted(4)
             grdpts_id(3,2)   = grdpts_id_inserted(5)
             grdpts_id(1:3,3) = grdpts_id_inserted(6:8)

             !output+validation
             test_loc = .not.check_if_bc_interior_pt_crenel(grdpts_id,2,2)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')',k
             end if

          end do

        end function test_check_if_bc_interior_pt_crenel

      end program test_bf_bc_interior_pt_crenel
