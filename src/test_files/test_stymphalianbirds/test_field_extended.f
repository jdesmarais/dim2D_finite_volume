      program test_field_extended

        use field_class, only :
     $       field

        use field_extended_class, only :
     $       field_extended

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        type(field)          :: field_used
        type(field_extended) :: field_extended_used

        logical              :: detailled
        logical              :: test_validated


        !initialized the fields
        call field_used%ini()
        call field_used%apply_bc_on_nodes()

        call field_extended_used%ini()
        call field_extended_used%apply_bc_on_nodes()


        !test whether the computation of the time
        !derivatives in the same in the field_used
        !and the field_extended when no buffer layer
        !are present
        detailled = .true.
        test_validated = test_compute_time_dev(
     $       field_used,
     $       field_extended_used,
     $       detailled)
        print '(''test_compute_time_dev: '',L1)', test_validated

        contains


        function test_compute_time_dev(
     $       field_used,
     $       field_extended_used,
     $       detailled)
     $       result(test_validated)

          implicit none

          class(field)         , intent(in) :: field_used
          class(field_extended), intent(in) :: field_extended_used
          logical              , intent(in) :: detailled
          logical                           :: test_validated


          real(rkind), dimension(nx,ny,ne) :: timedev_f
          real(rkind), dimension(nx,ny,ne) :: timedev_f_ext
          
          timedev_f     = field_used%compute_time_dev()
          timedev_f_ext = field_extended_used%compute_time_dev()

          test_validated = compare_arrays(
     $         timedev_f,
     $         timedev_f_ext,
     $         detailled)

        end function test_compute_time_dev


        function compare_arrays(
     $     timedev_f,
     $     timedev_f_ext,
     $     detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(in) :: timedev_f
          real(rkind), dimension(nx,ny,ne), intent(in) :: timedev_f_ext
          logical                         , intent(in) :: detailled
          logical                                      :: test_validated

          integer(ikind) :: i,j
          integer        :: k

          test_validated = .true.

          do k=1, ne
             do j=1, ny
                do i=1, nx
                   
                   test_validated = test_validated.and.is_test_validated(
     $                  timedev_f(i,j,k),
     $                  timedev_f_ext(i,j,k),
     $                  detailled)

                end do
             end do
          end do

        end function compare_arrays


        function is_test_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical    , intent(in) :: detailled
          logical                 :: test_validated

                       
          test_validated=abs(
     $         int(var*10000.)-
     $         sign(int(abs(cst*10000.)),int(cst*10000.))).le.1

          if((.not.test_validated).and.detailled) then
             print *, int(var*1e5)
             print *, int(cst*1e5)
          end if
          
        end function is_test_validated

      end program test_field_extended
