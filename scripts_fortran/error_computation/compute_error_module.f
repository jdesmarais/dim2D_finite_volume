      module compute_error_module

        use compare_type_module, only :
     $     compare_rkind

        use netcdf

        use parameters_cst, only :
     $       NOT_SUCCESS,
     $       SUCCESS

        use parameters_kind, only :
     $       rkind,
     $       ikind

        implicit none

        private
        public ::
     $       compute_relative_error,
     $       get_index_coord

        contains


        !var1: small_domain
        !var2: large_domain (reference)
        !compute relative error
        subroutine compute_relative_error(
     $       var1,
     $       var2,
     $       error,
     $       error_max,
     $       x_error_max,
     $       y_error_max,
     $       x_map,
     $       y_map,
     $       ierror)

          implicit none

          real(rkind), dimension(:,:,:)                       , intent(in)  :: var1
          real(rkind), dimension(:,:,:)                       , intent(in)  :: var2
          real(rkind), dimension(:,:,:), allocatable          , intent(out) :: error
          real(rkind), dimension(:)    , allocatable, optional, intent(out) :: error_max
          real(rkind), dimension(:)    , allocatable, optional, intent(out) :: x_error_max
          real(rkind), dimension(:)    , allocatable, optional, intent(out) :: y_error_max
          real(rkind), dimension(:)                 , optional, intent(in)  :: x_map
          real(rkind), dimension(:)                 , optional, intent(in)  :: y_map
          integer                                   , optional, intent(out) :: ierror


          integer                                :: ierror_op
          integer(ikind)                         :: size_x
          integer(ikind)                         :: size_y
          integer(ikind)                         :: size_e
          integer(ikind)                         :: i,j
          integer                                :: k
          logical                                :: test_inputs
          real(rkind), dimension(:), allocatable :: error_max_op
          logical                                :: compute_xy_error_max
          

          ierror_op = NOT_SUCCESS


          !1) verify that var1 and var2 have the same shape
          size_x = size(var1,1)
          size_y = size(var1,2)
          size_e = size(var1,3)
          test_inputs =
     $         (size_x.eq.(size(var2,1))).and.
     $         (size_y.eq.(size(var2,2))).and.
     $         (size_e.eq.(size(var2,3)))

          if(.not.test_inputs) then

             if(size_x.ne.size(var2,1)) then
                print '(''size_x do not match'')'
                ierror = NOT_SUCCESS
             end if

             if(size_y.ne.size(var2,2)) then
                print '(''size_y do not match'')'
                ierror = NOT_SUCCESS
             end if

             if(size_e.ne.size(var2,3)) then
                print '(''size_e do not match'')'
                ierror = NOT_SUCCESS
             end if

          end if


          !2) verify that the x_map, y_map have the same shape
          !   as size_x and size_y
          if(present(x_map).and.present(y_map)) then
             test_inputs =
     $            (size(x_map,1).eq.size_x).and.
     $            (size(y_map,1).eq.size_y)

             if(size_x.ne.size(x_map,1)) then
                print '(''size_x do not match for x_map'')'
                ierror = NOT_SUCCESS
             end if

             if(size_y.ne.size(y_map,1)) then
                print '(''size_y do not match for y_map'')'
                ierror = NOT_SUCCESS
             end if

          end if

          compute_xy_error_max =
     $         present(x_map).and.
     $         present(y_map).and.
     $         present(x_error_max).and.
     $         present(y_error_max)


          !2) allocate the error arrays
          allocate(error(size_x,size_y,size_e))
          allocate(error_max_op(size_e))

          if(compute_xy_error_max) then
             allocate(x_error_max(size_e))
             allocate(y_error_max(size_e))
          end if


          !3) fill the error array
          do k=1, size_e

             error_max_op(k) = abs(var1(1,1,k)-var2(1,1,k))/abs(var2(1,1,k))

             if(compute_xy_error_max) then
                x_error_max(k)  = x_map(1)
                y_error_max(k)  = y_map(1)
             end if

             do j=1, size_y
                do i=1, size_x

                   error(i,j,k) = abs(var1(i,j,k)-var2(i,j,k))/abs(var2(i,j,k))

                   if(error(i,j,k).gt.error_max_op(k)) then
                      error_max_op(k) = error(i,j,k)

                      if(compute_xy_error_max) then
                         x_error_max(k) = x_map(i)
                         y_error_max(k) = y_map(j)
                      end if

                   end if

                end do
             end do
          end do

          ierror_op = SUCCESS


          !4) create the outputs
          if(present(error_max)) then
             call MOVE_ALLOC(error_max_op,error_max)
          else
             deallocate(error_max_op)
          end if

          if(present(ierror)) then
             ierror = ierror_op
          end if

        end subroutine compute_relative_error


        !get the index 'i' corresponding to 'coord' in 'map'
        !such that map[i]=coord
        function get_index_coord(
     $       coord,
     $       coord_map,
     $       side,
     $       detailled,
     $       ierror)
     $       result(index)
        
          implicit none

          real(rkind)              , intent(in)  :: coord
          real(rkind), dimension(:), intent(in)  :: coord_map
          logical    , optional    , intent(in)  :: side
          logical    , optional    , intent(in)  :: detailled
          integer    , optional    , intent(out) :: ierror
          integer                                :: index

          logical :: side_op
          logical :: detailled_op
          integer :: ierror_op

          integer :: i

          ierror_op = NOT_SUCCESS


          if(present(side)) then
             side_op = side
          else
             side_op = .true.
          end if

          
          if(present(detailled)) then
             detailled_op = detailled
          else
             detailled_op = .false.
          end if


          !start from the right side
          if(side_op) then

             do i=size(coord_map,1), 1, -1

                if(compare_rkind(
     $               coord,
     $               coord_map(i),
     $               detailled=detailled_op)
     $            ) then

                   index = i
                   ierror_op = SUCCESS
                   exit
                end if

             end do
             

          !start from the left side
          else

             do i=1, size(coord_map,1), +1

                if(compare_rkind(coord,coord_map(i),detailled=detailled_op)) then
                   index = i
                   ierror_op = SUCCESS
                   exit
                end if

             end do

          end if


          !set ierror
          if(present(ierror)) then
             ierror = ierror_op
          end if

        end function get_index_coord

      end module compute_error_module
