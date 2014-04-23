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
          

        end subroutine ini_path

      end module test_cases_path_module
