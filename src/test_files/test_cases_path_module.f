      module test_cases_path_module

        use bf_layer_path_class, only : bf_layer_path
        use parameters_constant, only : N,S,E,W,N_E,N_W,S_E,S_W
        use parameters_input   , only : nx,ny,bc_size

        implicit none

        private
        public :: ini_path


        contains


        subroutine ini_path(
     $       current_path,
     $       test_case_id,
     $       corner_id,
     $       corner_order,
     $       bf_corner_distance)

          implicit none

          type(bf_layer_path)          , intent(inout) :: current_path
          integer                      , intent(in)    :: test_case_id
          integer            , optional, intent(in)    :: corner_id
          integer            , optional, intent(in)    :: corner_order
          integer            , optional, intent(in)    :: bf_corner_distance


          integer :: corner_id_i, corner_order_i, bf_corner_distance_i


          !initialize the optional argments
          if(present(corner_id)) then
             corner_id_i = corner_id
          else
             corner_id_i = N
          end if

          if(present(corner_order)) then
             corner_order_i = corner_order
          else
             corner_order_i = 1
          end if

          if(present(bf_corner_distance)) then
             bf_corner_distance_i = bf_corner_distance
          else
             bf_corner_distance_i = 3
          end if


          !select the test case
          select case(test_case_id)
            case(0)
               call ini_path_testcase0(
     $              current_path,
     $              corner_id,
     $              corner_order,
     $              bf_corner_distance)
            case default
               print '(''test_cases_path_module'')'
               print '(''ini_path'')'
               print '(''test case ID not recognized'')'
               print '(''test_case_id: '',I2)', test_case_id
               stop 'change test_case_id'
          end select

        end subroutine ini_path



        subroutine ini_path_testcase0(
     $       current_path,
     $       corner_id,
     $       corner_order,
     $       bf_corner_distance)

          implicit none

          type(bf_layer_path), intent(inout) :: current_path
          integer            , intent(in)    :: corner_id
          integer            , intent(in)    :: corner_order
          integer            , intent(in)    :: bf_corner_distance


          if((corner_order.ne.1).and.(corner_order.ne.2)) then
             stop 'corner_order not recognized'
          end if

          current_path%ends             = .true.
          current_path%ends_with_corner = .true.
          current_path%corner_id        = corner_id

          select case(corner_id)

            case(N_E)
               if(corner_order.eq.1) then
                  current_path%mainlayer = N
                  current_path%alignment(1,1) = nx-bc_size-bf_corner_distance-4
                  current_path%alignment(1,2) = nx-bc_size-bf_corner_distance
                  current_path%alignment(2,1) = ny-1
                  current_path%alignment(2,2) = ny-1            
               else
                  current_path%mainlayer = E
                  current_path%alignment(1,1) = nx-1
                  current_path%alignment(1,2) = nx-1
                  current_path%alignment(2,1) = ny-bc_size-bf_corner_distance-4
                  current_path%alignment(2,2) = ny-bc_size-bf_corner_distance
               end if

            case(N_W)
               if(corner_order.eq.2) then
                  current_path%mainlayer = N
                  current_path%alignment(1,1) = 1+bc_size+bf_corner_distance
                  current_path%alignment(1,2) = 1+bc_size+bf_corner_distance+4
                  current_path%alignment(2,1) = ny-1
                  current_path%alignment(2,2) = ny-1
               else
                  current_path%mainlayer = W
                  current_path%alignment(1,1) = 1
                  current_path%alignment(1,2) = 1
                  current_path%alignment(2,1) = ny-bc_size-bf_corner_distance-4
                  current_path%alignment(2,2) = ny-bc_size-bf_corner_distance
               end if

            case(S_E)
               if(corner_order.eq.2) then
                  current_path%mainlayer = S
                  current_path%alignment(1,1) = nx-bc_size-bf_corner_distance-4
                  current_path%alignment(1,2) = nx-bc_size-bf_corner_distance
                  current_path%alignment(2,1) = 1
                  current_path%alignment(2,2) = 1
               else
                  current_path%mainlayer = E
                  current_path%alignment(1,1) = nx-1
                  current_path%alignment(1,2) = nx-1
                  current_path%alignment(2,1) = 1+bc_size+bf_corner_distance
                  current_path%alignment(2,2) = 1+bc_size+bf_corner_distance+4
               end if

            case(S_W)
               if(corner_order.eq.1) then
                  current_path%mainlayer = S
                  current_path%alignment(1,1) = 1+bc_size+bf_corner_distance
                  current_path%alignment(1,2) = 1+bc_size+bf_corner_distance+4
                  current_path%alignment(2,1) = 1
                  current_path%alignment(2,2) = 1
               else
                  current_path%mainlayer = W
                  current_path%alignment(1,1) = 1
                  current_path%alignment(1,2) = 1
                  current_path%alignment(2,1) = 1+bc_size+bf_corner_distance
                  current_path%alignment(2,2) = 1+bc_size+bf_corner_distance+4
               end if

            case default
               print '(''test_cases_path_module'')'
               print '(''ini_path'')'
               stop 'corner_id not recognized'
          end select 
          

        end subroutine ini_path_testcase0

      end module test_cases_path_module
