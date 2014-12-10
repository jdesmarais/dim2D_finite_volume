      module check_data_module

        use parameters_kind, only :
     $     rkind

        implicit none


        private
        public ::
     $       is_test_validated,
     $       is_vector_validated,
     $       is_matrix_validated

        
        contains


        function is_test_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical                 :: detailled
          logical                 :: test_validated

          if(detailled) then
             print *, nint(var*1e5)
             print *, nint(cst*1e5)
          end if
          
          test_validated=abs(
     $         nint(var*1e5)-
     $         nint(cst*1e5)).le.1
          
        end function is_test_validated


        function is_vector_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), dimension(:), intent(in) :: var
          real(rkind), dimension(:), intent(in) :: cst
          logical                  , intent(in) :: detailled
          logical                               :: test_validated

          logical :: test_loc
          integer :: i

          test_validated = .true.

          do i=1,size(var,1)
             test_loc = is_test_validated(var(i),cst(i),.false.)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''['',I2'']:'',F8.3,'' -> '',F8.3)', 
     $               i, var(i), cst(i)
             end if
          end do

        end function is_vector_validated


        function is_matrix_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), dimension(:,:), intent(in) :: var
          real(rkind), dimension(:,:), intent(in) :: cst
          logical                    , intent(in) :: detailled
          logical                                 :: test_validated

          logical :: test_loc
          integer :: i,j

          test_validated = .true.

          do j=1,size(var,2)
             do i=1,size(var,1)
                test_loc = is_test_validated(var(i,j),cst(i,j),.false.)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''['',2I2'']:'',F8.3,'' -> '',F8.3)', 
     $                  i,j,
     $                  var(i,j), cst(i,j)
                end if
             end do
          end do

        end function is_matrix_validated

      end module check_data_module
