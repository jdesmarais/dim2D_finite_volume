      program test_bf_bc_crenel

        use bf_bc_crenel_module, only :
     $     is_temp_array_needed_for_bc_crenel,
     $     detect_bc_double_crenel,
     $     curb_bc_double_crenel,
     $     detect_and_curb_bc_single_crenel

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       N,S,E,W,
     $       bc_pt,bc_interior_pt,interior_pt,no_pt

        use parameters_kind, only :
     $       ikind, rkind

        implicit none


        logical :: test_loc
        logical :: test_validated
        logical :: detailled


        detailled = .true.
        test_validated = .true.

        
        test_loc = test_is_temp_array_needed_for_bc_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_is_temp_array_needed_for_bc_crenel: '',L1)', test_loc
        print '()'


        test_loc = test_detect_bc_double_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_detect_bc_double_crenel: '',L1)', test_loc
        print '()'


        test_loc = test_curb_bc_double_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_curb_bc_double_crenel: '',L1)', test_loc
        print '()'


        test_loc = test_curb_bc_single_crenel(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_curb_bc_single_crenel: '',L1)', test_loc
        print '()'

        contains


        function test_curb_bc_single_crenel(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          integer                 :: k
          integer, dimension(3,3) :: bf_grdpts_id
          integer, dimension(3,3) :: test_bf_grdpts_id
          logical                 :: bc_pt_crenel_exists
          logical                 :: test_bc_pt_crenel_exists
          logical                 :: test_loc
          

          test_validated = .true.


          do k=1,8

             call get_test_curb_bc_single_crenel_param(
     $            k,
     $            bf_grdpts_id,
     $            test_bf_grdpts_id,
     $            test_bc_pt_crenel_exists)

             bc_pt_crenel_exists = detect_and_curb_bc_single_crenel(
     $            [2,2],
     $            [3,3],
     $            bf_grdpts_id)

             test_loc = is_int_matrix_validated(
     $            bf_grdpts_id,
     $            test_bf_grdpts_id,
     $            detailled=detailled).and.
     $            (bc_pt_crenel_exists.eqv.test_bc_pt_crenel_exists)

             test_validated = test_validated.and.test_loc

             if(detailled) then
                if(.not.test_loc) then
                   print '(''test '',I2,'' failed'')', k
                   print '(''crenel_exists: ''L1, '' -> '',L1)',
     $                  bc_pt_crenel_exists, test_bc_pt_crenel_exists
                   print *, bf_grdpts_id, ' -> ', test_bf_grdpts_id
                   print '()'
                else
                   print '(''test '',I2,'' validated'')', k
                end if                
             end if

          end do

        end function test_curb_bc_single_crenel
        

        subroutine get_test_curb_bc_single_crenel_param(
     $     test_id,
     $     bf_grdpts_id,
     $     test_bf_grdpts_id,
     $     bc_pt_crenel_exists)

          implicit none

          integer                , intent(in)  :: test_id
          integer, dimension(3,3), intent(out) :: bf_grdpts_id
          integer, dimension(3,3), intent(out) :: test_bf_grdpts_id
          logical                , intent(out) :: bc_pt_crenel_exists


          select case(test_id)
            case(1)

               bf_grdpts_id = reshape((/
     $              bc_interior_pt, bc_interior_pt, bc_pt,
     $              bc_pt         , bc_pt         , bc_pt,
     $              no_pt         , no_pt         , no_pt/),
     $              (/3,3/))

               test_bf_grdpts_id = bf_grdpts_id
               bc_pt_crenel_exists = .false.

            case(2)

               bf_grdpts_id = reshape((/
     $              bc_pt         , bc_pt         , no_pt,
     $              bc_interior_pt, bc_pt         , no_pt,
     $              bc_interior_pt, bc_pt         , no_pt/),
     $              (/3,3/))

               test_bf_grdpts_id = bf_grdpts_id
               bc_pt_crenel_exists = .false.

            case(3)

               bf_grdpts_id = reshape((/
     $              no_pt         , no_pt         , no_pt,
     $              bc_pt         , bc_pt         , bc_pt,
     $              bc_pt         , bc_interior_pt, bc_interior_pt/),
     $              (/3,3/))

               test_bf_grdpts_id = bf_grdpts_id
               bc_pt_crenel_exists = .false.

            case(4)

               bf_grdpts_id = reshape((/
     $              no_pt, bc_pt, bc_interior_pt,
     $              no_pt, bc_pt, bc_interior_pt,
     $              no_pt, bc_pt, bc_pt/),
     $              (/3,3/))

               test_bf_grdpts_id = bf_grdpts_id
               bc_pt_crenel_exists = .false.

            case(5)
               
               bf_grdpts_id = reshape((/
     $              bc_pt         , bc_pt         , bc_pt,
     $              bc_interior_pt, bc_pt         , bc_interior_pt,
     $              bc_interior_pt, bc_interior_pt, bc_interior_pt/),
     $              (/3,3/))

               test_bf_grdpts_id = bf_grdpts_id
               test_bf_grdpts_id(2,2) = bc_interior_pt
               bc_pt_crenel_exists = .true.

            case(6)
               
               bf_grdpts_id = reshape((/
     $              bc_pt, bc_interior_pt, bc_interior_pt,
     $              bc_pt, bc_pt         , bc_interior_pt,
     $              bc_pt, bc_interior_pt, bc_interior_pt/),
     $              (/3,3/))

               test_bf_grdpts_id = bf_grdpts_id
               test_bf_grdpts_id(2,2) = bc_interior_pt
               bc_pt_crenel_exists = .true.

            case(7)
               
               bf_grdpts_id = reshape((/
     $              bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $              bc_interior_pt, bc_pt         , bc_interior_pt,
     $              bc_pt         , bc_pt         , bc_pt/),
     $              (/3,3/))

               test_bf_grdpts_id = bf_grdpts_id
               test_bf_grdpts_id(2,2) = bc_interior_pt
               bc_pt_crenel_exists = .true.

            case(8)
               
               bf_grdpts_id = reshape((/
     $              bc_interior_pt, bc_interior_pt, bc_pt,
     $              bc_interior_pt, bc_pt         , bc_pt,
     $              bc_interior_pt, bc_interior_pt, bc_pt/),
     $              (/3,3/))

               test_bf_grdpts_id = bf_grdpts_id
               test_bf_grdpts_id(2,2) = bc_interior_pt
               bc_pt_crenel_exists = .true.

            case default
               print '(''test_bf_bc_crenel'')'
               print '(''get_test_curb_bc_single_crenel_param'')'
               print '(''test no implemented: '',I2)', test_id

          end select             

        end subroutine get_test_curb_bc_single_crenel_param

        
        function test_curb_bc_double_crenel(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          
          integer                 :: k
          integer, dimension(2)   :: bc_pt_double_crenel_coords
          integer, dimension(4,4) :: bf_grdpts_id
          integer, dimension(4,4) :: test_bf_grdpts_id
          logical                 :: test_loc
          

          test_validated = .true.


          do k=1,4

             call get_test_curb_bc_double_crenel_param(
     $            k,
     $            bc_pt_double_crenel_coords,
     $            bf_grdpts_id,
     $            test_bf_grdpts_id)

             call curb_bc_double_crenel(
     $            bc_pt_double_crenel_coords,
     $            [4,4],
     $            bf_grdpts_id)

             test_loc = is_int_matrix_validated(
     $            bf_grdpts_id,
     $            test_bf_grdpts_id,
     $            detailled=detailled)

             test_validated = test_validated.and.test_loc

             if(detailled) then
                if(.not.test_loc) then
                   print '(''test '',I2,'' failed'')', k
                   print *, bf_grdpts_id, ' -> ', test_bf_grdpts_id
                   print '()'
                else
                   print '(''test '',I2,'' validated'')', k
                end if                
             end if

          end do

        end function test_curb_bc_double_crenel


        subroutine get_test_curb_bc_double_crenel_param(
     $     test_id,
     $     bc_pt_double_crenel_coords,
     $     bf_grdpts_id,
     $     test_bf_grdpts_id)

          implicit none

          integer                , intent(in)  :: test_id
          integer, dimension(2)  , intent(out) :: bc_pt_double_crenel_coords
          integer, dimension(4,4), intent(out) :: bf_grdpts_id
          integer, dimension(4,4), intent(out) :: test_bf_grdpts_id


          select case(test_id)
            case(1)

               bc_pt_double_crenel_coords = [3,2]

               bf_grdpts_id = reshape((/
     $              interior_pt, bc_interior_pt, bc_interior_pt, bc_pt,
     $              interior_pt, bc_interior_pt, bc_pt, bc_pt,
     $              interior_pt, bc_interior_pt, bc_pt, bc_pt,
     $              interior_pt, bc_interior_pt, bc_interior_pt, bc_pt/),
     $              (/4,4/))

               test_bf_grdpts_id = reshape((/
     $              interior_pt, interior_pt, bc_interior_pt, bc_pt,
     $              interior_pt, interior_pt, bc_interior_pt, bc_pt,
     $              interior_pt, interior_pt, bc_interior_pt, bc_pt,
     $              interior_pt, interior_pt, bc_interior_pt, bc_pt/),
     $              (/4,4/))

            case(2)

               bc_pt_double_crenel_coords = [1,2]

               bf_grdpts_id = reshape((/
     $              bc_pt, bc_interior_pt, bc_interior_pt, interior_pt,
     $              bc_pt, bc_pt         , bc_interior_pt, interior_pt,
     $              bc_pt, bc_pt         , bc_interior_pt, interior_pt,
     $              bc_pt, bc_interior_pt, bc_interior_pt, interior_pt
     $              /),
     $              (/4,4/))

               test_bf_grdpts_id = reshape((/
     $              bc_pt, bc_interior_pt, interior_pt, interior_pt,
     $              bc_pt, bc_interior_pt, interior_pt, interior_pt,
     $              bc_pt, bc_interior_pt, interior_pt, interior_pt,
     $              bc_pt, bc_interior_pt, interior_pt, interior_pt/),
     $              (/4,4/))

            case(3)

               bc_pt_double_crenel_coords = [2,1]

               bf_grdpts_id = reshape((/
     $              bc_pt         , bc_pt         , bc_pt         , bc_pt         ,
     $              bc_interior_pt, bc_pt         , bc_pt         , bc_interior_pt,
     $              bc_interior_pt, bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $              interior_pt   , interior_pt   , interior_pt   , interior_pt   
     $              /),
     $              (/4,4/))

               test_bf_grdpts_id = reshape((/
     $              bc_pt         , bc_pt         , bc_pt         , bc_pt         ,
     $              bc_interior_pt, bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $              interior_pt   , interior_pt,    interior_pt   , interior_pt,
     $              interior_pt   , interior_pt   , interior_pt   , interior_pt   
     $              /),
     $              (/4,4/))

            case(4)

               bc_pt_double_crenel_coords = [2,3]

               bf_grdpts_id = reshape((/
     $              interior_pt   , interior_pt   , interior_pt   , interior_pt   ,
     $              bc_interior_pt, bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $              bc_interior_pt, bc_pt         , bc_pt         , bc_interior_pt,
     $              bc_pt         , bc_pt         , bc_pt         , bc_pt
     $              /),
     $              (/4,4/))

               test_bf_grdpts_id = reshape((/
     $              interior_pt   , interior_pt   , interior_pt   , interior_pt,   
     $              interior_pt   , interior_pt,    interior_pt   , interior_pt,
     $              bc_interior_pt, bc_interior_pt, bc_interior_pt, bc_interior_pt,
     $              bc_pt         , bc_pt         , bc_pt         , bc_pt
     $              /),
     $              (/4,4/))

            case default
               print '(''test_bf_grdpts_id_bc_double_crenel'')'
               print '(''get_test_curb_bc_double_crenel_param'')'
               print '(''test not implemneted: '',I2)', test_id
               stop ''
               
          end select

        end subroutine get_test_curb_bc_double_crenel_param
        

        function test_detect_bc_double_crenel(detailled)
     $       result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer                              :: k
          integer, dimension(2)                :: bf_sizes
          integer, dimension(2)                :: cpt_local_coords
          integer, dimension(:,:), allocatable :: bf_grdpts_id
          integer, dimension(2)                :: test_bc_pt_double_crenel_coords
          logical                              :: test_bc_pt_double_crenel_exists

          integer, dimension(2)                :: bc_pt_double_crenel_coords
          logical                              :: bc_pt_double_crenel_exists


          test_validated = .true.

          do k=1, 22

             call get_test_detect_bc_double_crenel_param(
     $            k,
     $            cpt_local_coords,
     $            bf_sizes,
     $            bf_grdpts_id,
     $            test_bc_pt_double_crenel_coords,
     $            test_bc_pt_double_crenel_exists)

             bc_pt_double_crenel_exists = detect_bc_double_crenel(
     $            cpt_local_coords,
     $            bf_sizes,
     $            bf_grdpts_id,
     $            bc_pt_double_crenel_coords)

             if(test_bc_pt_double_crenel_exists) then
                test_loc = 
     $               (bc_pt_double_crenel_exists.eqv.test_bc_pt_double_crenel_exists).and.
     $               (bc_pt_double_crenel_coords(1).eq.test_bc_pt_double_crenel_coords(1)).and.
     $               (bc_pt_double_crenel_coords(2).eq.test_bc_pt_double_crenel_coords(2))

             else
                test_loc = 
     $               (bc_pt_double_crenel_exists.eqv.test_bc_pt_double_crenel_exists)

             end if

             if(detailled) then
                if(.not.test_loc) then
                   print '(''test '',I2,'' failed'')', k

                   print '(''bc_pt_exists: '',L2,'' -> '',L2)',
     $                  bc_pt_double_crenel_exists,
     $                  test_bc_pt_double_crenel_exists

                   if(test_bc_pt_double_crenel_exists) then
                      print '(''crenel_coords(1): '',I2,'' -> '',I2)',
     $                     bc_pt_double_crenel_coords(1), test_bc_pt_double_crenel_coords(1)
                      print '(''crenel_coords(2): '',I2,'' -> '',I2)',
     $                     bc_pt_double_crenel_coords(2), test_bc_pt_double_crenel_coords(2)
                      print '()'

                   end if

                else
                   print '(''test '',I2,'' validated'')', k
                end if
             end if

             test_validated = test_validated.and.test_loc

             deallocate(bf_grdpts_id)

          end do

        end function test_detect_bc_double_crenel


        subroutine get_test_detect_bc_double_crenel_param(
     $     test_id,
     $     cpt_local_coords,
     $     bf_sizes,
     $     bf_grdpts_id,
     $     test_bc_pt_double_crenel_coords,
     $     test_bc_pt_double_crenel_exists)
        
          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(2)               , intent(out) :: cpt_local_coords
          integer, dimension(2)               , intent(out) :: bf_sizes
          integer, dimension(:,:), allocatable, intent(out) :: bf_grdpts_id
          integer, dimension(2)               , intent(out) :: test_bc_pt_double_crenel_coords
          logical                             , intent(out) :: test_bc_pt_double_crenel_exists


          test_bc_pt_double_crenel_coords = [0,0]

          select case(test_id)
            case(1)
               allocate(bf_grdpts_id(3,3))
               bf_grdpts_id = reshape((/
     $              bc_interior_pt,bc_pt,no_pt,
     $              bc_interior_pt,bc_pt,no_pt,
     $              bc_interior_pt,bc_pt,no_pt/),
     $              (/3,3/))
               test_bc_pt_double_crenel_exists = .false.
               cpt_local_coords = [2,2]
               
            case(2)
               allocate(bf_grdpts_id(3,3))
               bf_grdpts_id = reshape((/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt/),
     $              (/3,3/))
               test_bc_pt_double_crenel_exists = .false.
               cpt_local_coords = [2,2]

            case(3)
               allocate(bf_grdpts_id(3,3))
               bf_grdpts_id = reshape((/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_pt,no_pt/),
     $              (/3,3/))
               test_bc_pt_double_crenel_exists = .false.
               cpt_local_coords = [2,2]

            case(4)
               allocate(bf_grdpts_id(3,3))
               bf_grdpts_id = reshape((/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              no_pt,bc_pt,bc_interior_pt/),
     $              (/3,3/))
               test_bc_pt_double_crenel_exists = .false.
               cpt_local_coords = [2,2]

            case(5)
               allocate(bf_grdpts_id(3,3))
               bf_grdpts_id = reshape((/
     $              no_pt,bc_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt/),
     $              (/3,3/))
               test_bc_pt_double_crenel_exists = .false.
               cpt_local_coords = [2,2]

            case(6)
               allocate(bf_grdpts_id(3,3))
               bf_grdpts_id = reshape((/
     $              bc_interior_pt,bc_pt,no_pt,
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt/),
     $              (/3,3/))
               test_bc_pt_double_crenel_exists = .false.
               cpt_local_coords = [2,2]

            case(7)
               allocate(bf_grdpts_id(3,3))
               bf_grdpts_id = reshape((/
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt/),
     $              (/3,3/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [1,1]
               cpt_local_coords = [2,2]

            case(8)
               allocate(bf_grdpts_id(2,2))
               bf_grdpts_id = reshape((/
     $              bc_pt,bc_pt,
     $              bc_pt,bc_pt/),
     $              (/2,2/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [1,1]
               cpt_local_coords = [2,2]

            case(9)
               allocate(bf_grdpts_id(2,3))
               bf_grdpts_id = reshape((/
     $              bc_pt,bc_pt,
     $              bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt/),
     $              (/2,3/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [1,1]
               cpt_local_coords = [2,2]

            case(10)
               allocate(bf_grdpts_id(3,2))
               bf_grdpts_id = reshape((/
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt/),
     $              (/3,2/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [1,1]
               cpt_local_coords = [2,2]

            case(11)
               allocate(bf_grdpts_id(3,3))
               bf_grdpts_id = reshape((/
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt/),
     $              (/3,3/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [2,1]
               cpt_local_coords = [2,2]

            case(12)
               allocate(bf_grdpts_id(2,2))
               bf_grdpts_id = reshape((/
     $              bc_pt,bc_pt,
     $              bc_pt,bc_pt/),
     $              (/2,2/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [1,1]
               cpt_local_coords = [1,2]

            case(13)
               allocate(bf_grdpts_id(2,3))
               bf_grdpts_id = reshape((/
     $              bc_pt,bc_pt,
     $              bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt/),
     $              (/2,3/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [1,1]
               cpt_local_coords = [1,2]

            case(14)
               allocate(bf_grdpts_id(3,2))
               bf_grdpts_id = reshape((/
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_pt,bc_pt/),
     $              (/3,2/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [2,1]
               cpt_local_coords = [2,2]

            case(15)
               allocate(bf_grdpts_id(3,3))
               bf_grdpts_id = reshape((/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt/),
     $              (/3,3/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [1,2]
               cpt_local_coords = [2,2]

            case(16)
               allocate(bf_grdpts_id(2,2))
               bf_grdpts_id = reshape((/
     $              bc_pt,bc_pt,
     $              bc_pt,bc_pt/),
     $              (/2,2/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [1,1]
               cpt_local_coords = [2,1]

            case(17)
               allocate(bf_grdpts_id(2,3))
               bf_grdpts_id = reshape((/
     $              bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,
     $              bc_pt,bc_pt/),
     $              (/2,3/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [1,2]
               cpt_local_coords = [2,2]

            case(18)
               allocate(bf_grdpts_id(3,2))
               bf_grdpts_id = reshape((/
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt/),
     $              (/3,2/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [1,1]
               cpt_local_coords = [2,1]

            case(19)
               allocate(bf_grdpts_id(3,3))
               bf_grdpts_id = reshape((/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_pt,bc_pt/),
     $              (/3,3/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [2,2]
               cpt_local_coords = [2,2]

            case(20)
               allocate(bf_grdpts_id(2,2))
               bf_grdpts_id = reshape((/
     $              bc_pt,bc_pt,
     $              bc_pt,bc_pt/),
     $              (/2,2/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [1,1]
               cpt_local_coords = [1,1]

            case(21)
               allocate(bf_grdpts_id(2,3))
               bf_grdpts_id = reshape((/
     $              bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,
     $              bc_pt,bc_pt/),
     $              (/2,3/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [1,2]
               cpt_local_coords = [1,2]

            case(22)
               allocate(bf_grdpts_id(3,2))
               bf_grdpts_id = reshape((/
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_pt,bc_pt/),
     $              (/3,2/))
               test_bc_pt_double_crenel_exists = .true.
               test_bc_pt_double_crenel_coords = [2,1]
               cpt_local_coords = [2,1]

            case default
               print '(''test_bf_grdpts_id_bc_double_crenel'')'
               print '(''get_test_detect_bc_double_crenel_param'')'
               print '(''test not implemented: '',I2)', test_id

          end select

          bf_sizes = [size(bf_grdpts_id,1),size(bf_grdpts_id,2)]

        end subroutine get_test_detect_bc_double_crenel_param


        function test_is_temp_array_needed_for_bc_crenel(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer               :: k
          integer               :: bf_localization
          integer, dimension(2) :: bf_sizes
          integer, dimension(2) :: cpt_local_coords
          logical               :: test_temp_array_needed
          logical               :: temp_array_needed
          logical               :: test_loc


          test_validated = .true.


          do k=0,47

             call get_test_array_needed_param(
     $            k,
     $            bf_localization,
     $            bf_sizes,
     $            cpt_local_coords,
     $            test_temp_array_needed)

             temp_array_needed = is_temp_array_needed_for_bc_crenel(
     $            bf_localization,
     $            bf_sizes,
     $            cpt_local_coords)

             test_loc = temp_array_needed.eqv.test_temp_array_needed

             if(detailled) then
                if(.not.test_loc) then
                   print '(''test '',I2,''failed'')', k
                   print '(L1,'' -> '',L1)', temp_array_needed, test_temp_array_needed
                   print '(''bf_localization:  '',I2)', bf_localization
                   print '(''bf_sizes       :  '',2I3)', bf_sizes
                   print '(''cpt_local_coords: '',2I3)', cpt_local_coords
                   print ''
                else
                   print '(''test '',I2,'' validated'')', k
                end if
             end if

             test_validated = test_validated.and.test_loc

          end do

        end function test_is_temp_array_needed_for_bc_crenel


        subroutine get_test_array_needed_param(
     $     test_id,
     $     bf_localization,
     $     bf_sizes,
     $     cpt_local_coords,
     $     test_temp_array_needed)

          implicit none

          integer              , intent(in)  :: test_id
          integer              , intent(out) :: bf_localization
          integer, dimension(2), intent(out) :: bf_sizes
          integer, dimension(2), intent(out) :: cpt_local_coords
          logical              , intent(out) :: test_temp_array_needed

          integer :: test_id_cpt_coords
          integer :: test_id_localization

          bf_sizes = [10,8]

          test_id_cpt_coords   = mod(test_id,12)
          test_id_localization = (test_id-mod(test_id,12))/12


          select case(test_id_cpt_coords)
            case(0)
               cpt_local_coords(1)=3
               cpt_local_coords(2)=3
            case(1)
               cpt_local_coords(1)=3
               cpt_local_coords(2)=2
            case(2)
               cpt_local_coords(1)=2
               cpt_local_coords(2)=3

            case(3)
               cpt_local_coords(1)=8
               cpt_local_coords(2)=3
            case(4)
               cpt_local_coords(1)=8
               cpt_local_coords(2)=2
            case(5)
               cpt_local_coords(1)=9
               cpt_local_coords(2)=3

            case(6)
               cpt_local_coords(1)=3
               cpt_local_coords(2)=6
            case(7)
               cpt_local_coords(1)=3
               cpt_local_coords(2)=7
            case(8)
               cpt_local_coords(1)=2
               cpt_local_coords(2)=6

            case(9)
               cpt_local_coords(1)=8
               cpt_local_coords(2)=6
            case(10)
               cpt_local_coords(1)=8
               cpt_local_coords(2)=7
            case(11)
               cpt_local_coords(1)=9
               cpt_local_coords(2)=6

            case default
               print '(''test_bf_grdpts_id_bc_double_crenel'')'
               print '(''get_test_array_needed_param'')'
               print '(''test_cpt_coords not implemented'',I2)', test_id_cpt_coords
               stop '' 

          end select


          select case(test_id_localization)
            case(0)
               bf_localization = N

               select case(test_id_cpt_coords)
                 case(0,3,6,7,9,10)
                    test_temp_array_needed=.false.
                 case default
                    test_temp_array_needed=.true.
               end select

            case(1)
               bf_localization = S

               select case(test_id_cpt_coords)
                 case(0,1,3,4,6,9)
                    test_temp_array_needed=.false.
                 case default
                    test_temp_array_needed=.true.
               end select

            case(2)
               bf_localization = E

               select case(test_id_cpt_coords)
                 case(0,3,5,6,9,11)
                    test_temp_array_needed=.false.
                 case default
                    test_temp_array_needed=.true.
               end select

            case(3)
               bf_localization = W
               
               select case(test_id_cpt_coords)
                 case(0,2,3,6,8,9)
                    test_temp_array_needed=.false.
                 case default
                    test_temp_array_needed=.true.
               end select
               
            case default
               print '(''test_bf_grdpts_id_bc_double_crenel'')'
               print '(''get_test_array_needed_param'')'
               print '(''test_localization not implemented'',I2)', test_id_localization
               stop ''

          end select

        end subroutine get_test_array_needed_param

      end program test_bf_bc_crenel
