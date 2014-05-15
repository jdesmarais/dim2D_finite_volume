      !< subroutines shared by bf_layer_class and bf_layer_path_class
      !> needed to implement corretcly the update of sublayers at the
      !> corner if a corner point is activated
      module bf_corner_module

        use parameters_input   , only : nx,ny,bc_size
        use parameters_constant, only : N,S,E,W,N_E,N_W,S_E,S_W,
     $                                  x_direction, y_direction,
     $                                  min_border, max_border
        use parameters_kind    , only : ikind

        implicit none

        private
        public :: is_alignment_compatible_with_corner,
     $            get_default_alignment_and_neighbors


        contains


        !< investigate whether the grid point needed by the new
        !> corner buffer layer can be included in an existing
        !> buffer layer or be included in the next path leading
        !> to a buffer layer modification
        function is_alignment_compatible_with_corner(
     $       corner_id, mainlayer_id,
     $       alignment, neighbors,
     $       new_alignment, new_neighbors)
     $       result(compatible)

          implicit none

          integer                       , intent(in)  :: corner_id
          integer                       , intent(in)  :: mainlayer_id
          integer(ikind), dimension(2,2), intent(in)  :: alignment
          logical       , dimension(4)  , intent(in)  :: neighbors
          integer(ikind), dimension(2,2), intent(out) :: new_alignment
          logical       , dimension(4)  , intent(out) :: new_neighbors
          logical                                     :: compatible

          integer :: direction
          integer :: n_direction
          integer :: border


          !we investigate whether the grid point needed by the new
          !corner buffer layer can be included in an existing
          !buffer layer or not
          !if the buffer layer investigated is N,S, the condition
          !relies on i_min or i_max
          !if the buffer layer investigated is E,W, the condition
          !relies on j_min or j_max
          select case(mainlayer_id)

            case(N,S)

               direction = x_direction

               select case(corner_id)
                 case(N_E,S_E)

                    n_direction = nx-bc_size
                    border      = max_border

                    compatible = is_compatible(
     $                   alignment, direction, border, n_direction)

                    if(compatible) then

                       new_neighbors(N)=neighbors(N)
                       new_neighbors(S)=neighbors(S)
                       new_neighbors(E)=.true.
                       new_neighbors(W)=neighbors(W)

                       new_alignment(1,1) = alignment(1,1)
                       new_alignment(2,1) = alignment(2,1)
                       new_alignment(1,2) = nx - bc_size
                       new_alignment(2,2) = alignment(2,2)

                    end if

                 case(N_W,S_W)
                    n_direction = bc_size+1
                    border      = min_border

                    compatible = is_compatible(
     $                   alignment, direction, border, n_direction)

                    if(compatible) then

                       new_neighbors(N)=neighbors(N)
                       new_neighbors(S)=neighbors(S)
                       new_neighbors(E)=neighbors(E)
                       new_neighbors(W)=.true.

                       new_alignment(1,1) = bc_size+1
                       new_alignment(2,1) = alignment(2,1)
                       new_alignment(1,2) = alignment(1,2)
                       new_alignment(2,2) = alignment(2,2)

                    end if
                 case default
                    call compatibility_issue(corner_id, mainlayer_id)
               end select

            case(E,W)

               direction = y_direction

               select case(corner_id)

                 case(N_W,N_E)
                    n_direction = ny-bc_size
                    border      = max_border

                    compatible = is_compatible(
     $                   alignment, direction, border, n_direction)

                    if(compatible) then

                       new_neighbors(N)=.true.
                       new_neighbors(S)=neighbors(S)
                       new_neighbors(E)=neighbors(E)
                       new_neighbors(W)=neighbors(W)

                       new_alignment(1,1) = alignment(1,1)
                       new_alignment(2,1) = alignment(2,1)
                       new_alignment(1,2) = alignment(1,2)
                       new_alignment(2,2) = ny-bc_size

                    end if

                 case(S_E,S_W)
                    n_direction = bc_size+1
                    border      = min_border

                    compatible = is_compatible(
     $                   alignment, direction, border, n_direction)

                    if(compatible) then

                       new_neighbors(N)=neighbors(N)
                       new_neighbors(S)=.true.
                       new_neighbors(E)=neighbors(E)
                       new_neighbors(W)=neighbors(W)

                       new_alignment(1,1) = alignment(1,1)
                       new_alignment(2,1) = bc_size+1
                       new_alignment(1,2) = alignment(1,2)
                       new_alignment(2,2) = alignment(2,2)

                    end if

                 case default
                    call compatibility_issue(corner_id, mainlayer_id)

               end select

            case default
               call compatibility_issue(corner_id, mainlayer_id)

          end select

        end function is_alignment_compatible_with_corner


        !< check is the alignment is compatible with the corner identified
        !> by if direction and its position
        function is_compatible(alignment, direction, border, n_direction)

          implicit none

          integer(ikind), dimension(2,2), intent(in) :: alignment
          integer                       , intent(in) :: direction
          integer                       , intent(in) :: border
          integer(ikind)                , intent(in) :: n_direction
          logical                                    :: is_compatible

          is_compatible = abs(alignment(direction,border)-n_direction)
     $         .le.(2*bc_size+1)

        end function is_compatible


        subroutine get_default_alignment_and_neighbors(
     $     corner_id, mainlayer_id,
     $     new_alignment, new_neighbors)

          implicit none

          integer                       , intent(in)  :: corner_id
          integer                       , intent(in)  :: mainlayer_id
          integer(ikind), dimension(2,2), intent(out) :: new_alignment
          logical       , dimension(4)  , intent(out) :: new_neighbors


          select case(corner_id)

            case(N_E)
               select case(mainlayer_id)
                 case(N)
                    new_alignment(1,1) = nx-bc_size
                    new_alignment(1,2) = nx-bc_size
                    new_alignment(2,1) = ny+1
                    new_alignment(2,2) = ny+1
                    
                    new_neighbors(N) = .false.
                    new_neighbors(S) = .false.
                    new_neighbors(E) = .true.
                    new_neighbors(W) = .false.
                 case(E)
                    new_alignment(1,1) = nx+1
                    new_alignment(1,2) = nx+1
                    new_alignment(2,1) = ny-bc_size
                    new_alignment(2,2) = ny-bc_size
                    
                    new_neighbors(N) = .true.
                    new_neighbors(S) = .false.
                    new_neighbors(E) = .false.
                    new_neighbors(W) = .false.
                 case default
                    call compatibility_issue(corner_id, mainlayer_id)
                 end select

            case(N_W)
               select case(mainlayer_id)
                 case(N)
                    new_alignment(1,1) = bc_size+1
                    new_alignment(1,2) = bc_size+1
                    new_alignment(2,1) = ny+1
                    new_alignment(2,2) = ny+1
                    
                    new_neighbors(N) = .false.
                    new_neighbors(S) = .false.
                    new_neighbors(E) = .false.
                    new_neighbors(W) = .true.
                 case(W)
                    new_alignment(1,1) = 0
                    new_alignment(1,2) = 0
                    new_alignment(2,1) = ny-bc_size
                    new_alignment(2,2) = ny-bc_size
                    
                    new_neighbors(N) = .true.
                    new_neighbors(S) = .false.
                    new_neighbors(E) = .false.
                    new_neighbors(W) = .false.
                 case default
                    call compatibility_issue(corner_id, mainlayer_id)
                 end select

            case(S_E)
               select case(mainlayer_id)
                 case(S)
                    new_alignment(1,1) = nx-bc_size
                    new_alignment(1,2) = nx-bc_size
                    new_alignment(2,1) = 0
                    new_alignment(2,2) = 0
                    
                    new_neighbors(N) = .false.
                    new_neighbors(S) = .false.
                    new_neighbors(E) = .true.
                    new_neighbors(W) = .false.
                 case(E)
                    new_alignment(1,1) = nx+1
                    new_alignment(1,2) = nx+1
                    new_alignment(2,1) = bc_size+1
                    new_alignment(2,2) = bc_size+1
                    
                    new_neighbors(N) = .false.
                    new_neighbors(S) = .true.
                    new_neighbors(E) = .false.
                    new_neighbors(W) = .false.
                 case default
                    call compatibility_issue(corner_id, mainlayer_id)
                 end select

            case(S_W)
               select case(mainlayer_id)
                 case(S)
                    new_alignment(1,1) = bc_size+1
                    new_alignment(1,2) = bc_size+1
                    new_alignment(2,1) = 0
                    new_alignment(2,2) = 0
                    
                    new_neighbors(N) = .false.
                    new_neighbors(S) = .false.
                    new_neighbors(E) = .false.
                    new_neighbors(W) = .true.
                 case(W)
                    new_alignment(1,1) = 0
                    new_alignment(1,2) = 0
                    new_alignment(2,1) = bc_size+1
                    new_alignment(2,2) = bc_size+1
                    
                    new_neighbors(N) = .false.
                    new_neighbors(S) = .true.
                    new_neighbors(E) = .false.
                    new_neighbors(W) = .false.
                 case default
                    call compatibility_issue(corner_id, mainlayer_id)
                 end select

            case default
               print '(''bf_corner_module'')'
               print '(''get_default_alignment'')'
               print '(''corner not recognized'')'
               print '(''corner_id: '',I2)', corner_id
               stop 'change corner_id'
          end select

        end subroutine get_default_alignment_and_neighbors



        subroutine compatibility_issue(corner_id, mainlayer_id)
        
          implicit none

          integer, intent(in) :: corner_id
          integer, intent(in) :: mainlayer_id

          print '(''bf_corner_module'')'
          print '(''is_alignment_compatible_with_corner'')'
          print '(''corner_id incompatible with mainlayer_id'')'
          print '(''corner_id: '', I2)', corner_id
          print '(''mainlayer_id: '', I2)', mainlayer_id
          stop 'change corner or mainlayer_id'

        end subroutine compatibility_issue

      end module bf_corner_module
