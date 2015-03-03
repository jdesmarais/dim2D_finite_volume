      module check_data_module

        use parameters_kind, only :
     $     rkind

        implicit none


        private
        public ::
     $       is_real_validated,
     $       is_real_vector_validated,
     $       is_real_matrix_validated,
     $       is_real_matrix3D_validated,
     $       is_int_vector_validated,
     $       is_int_matrix_validated,
     $       is_int_matrix3D_validated

        
        contains


        function is_real_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical                 :: detailled
          logical                 :: test_validated

          test_validated = abs(var-cst)<1e-10

          if(detailled.and.(.not.test_validated)) then
             print *, var
             print *, cst
          end if
          
        end function is_real_validated


        function is_real_vector_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), dimension(:), intent(in) :: var
          real(rkind), dimension(:), intent(in) :: cst
          logical                  , intent(in) :: detailled
          logical                               :: test_validated

          logical :: test_loc
          integer :: i

          test_validated = .true.

          if(size(var,1).eq.size(cst,1)) then

             do i=1,size(var,1)
                test_loc = is_real_validated(var(i),cst(i),detailled)
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''['',I2'']:'',F8.3,'' -> '',F8.3)', 
     $                  i, var(i), cst(i)
                end if
             end do

          else

             test_validated = .false.
             print '(''sizes do not match'')'
             print '(''  - size_x : '',I2,'' -> '',I2)', size(var,1), size(cst,1)
             print '()'

          end if

        end function is_real_vector_validated


        function is_real_matrix_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), dimension(:,:), intent(in) :: var
          real(rkind), dimension(:,:), intent(in) :: cst
          logical    , optional      , intent(in) :: detailled
          logical                                 :: test_validated

          logical :: test_loc
          integer :: i,j
          logical :: detailled_op

          if(present(detailled)) then
             detailled_op = detailled
          else
             detailled_op = .false.
          end if

          test_validated = .true.

          if((size(var,1).eq.size(cst,1)).and.
     $         (size(var,2).eq.size(cst,2))) then

             do j=1,size(var,2)
                do i=1,size(var,1)
                   test_loc = is_real_validated(var(i,j),cst(i,j),detailled)
                   test_validated = test_validated.and.test_loc
                   if(detailled_op.and.(.not.test_loc)) then
                      print '(''['',2I2'']:'',F8.3,'' -> '',F8.3)', 
     $                     i,j,
     $                     var(i,j), cst(i,j)
                   end if
                end do
             end do

          else

             test_validated = .false.
             print '(''sizes do not match'')'
             print '(''  - size_x : '',I2,'' -> '',I2)', size(var,1), size(cst,1)
             print '(''  - size_y : '',I2,'' -> '',I2)', size(var,2), size(cst,2)
             print '()'

          end if

        end function is_real_matrix_validated


        function is_real_matrix3D_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: var
          real(rkind), dimension(:,:,:), intent(in) :: cst
          logical    , optional        , intent(in) :: detailled
          logical                                   :: test_validated

          logical :: test_loc
          integer :: i,j,k
          logical :: detailled_op

          if(present(detailled)) then
             detailled_op = detailled
          else
             detailled_op = .false.
          end if

          test_validated = .true.

          if((size(var,1).eq.size(cst,1)).and.
     $       (size(var,2).eq.size(cst,2)).and.
     $       (size(var,3).eq.size(cst,3))) then
             
             do k=1,size(var,3)
                do j=1,size(var,2)
                   do i=1,size(var,1)
                      test_loc = is_real_validated(var(i,j,k),cst(i,j,k),detailled)
                      test_validated = test_validated.and.test_loc
                      if(detailled_op.and.(.not.test_loc)) then
                         print '(''['',3I3'']:'',F8.3,'' -> '',F8.3)', 
     $                        i,j,k,
     $                        var(i,j,k), cst(i,j,k)
                      end if
                   end do
                end do
             end do

          else

             test_validated = .false.
             print '(''sizes do not match'')'
             print '(''  - size_x : '',I2,'' -> '',I2)', size(var,1), size(cst,1)
             print '(''  - size_y : '',I2,'' -> '',I2)', size(var,2), size(cst,2)
             print '(''  - size_z : '',I2,'' -> '',I2)', size(var,3), size(cst,3)
             print '()'

          end if

        end function is_real_matrix3D_validated


        function is_int_vector_validated(
     $     int_vector,
     $     int_vector_cst,
     $     detailled)
     $     result(test_validated)

          implicit none

          integer, dimension(:), intent(in) :: int_vector
          integer, dimension(:), intent(in) :: int_vector_cst
          logical, optional    , intent(in) :: detailled
          logical                           :: test_validated


          integer :: i
          logical :: test_loc
          logical :: detailled_op


          if(present(detailled)) then
             detailled_op = detailled
          else
             detailled_op = .false.
          end if


          test_validated = .true.


          if(size(int_vector,1).eq.size(int_vector_cst,1)) then

             do i=1, size(int_vector)

                test_loc = int_vector(i).eq.int_vector_cst(i)
                test_validated = test_validated.and.test_loc

                if(detailled_op.and.(.not.test_loc)) then

                   print '(''['',I2,'']: '',I5, '' -> '',I5)',
     $                  i, int_vector(i), int_vector_cst(i)

                end if

             end do

          else

             test_validated = .false.
             print '(''sizes do not match'')'
             print '(''  - size_x : '',I2,'' -> '',I2)', size(int_vector,1), size(int_vector_cst,1)
             print '()'

          end if

        end function is_int_vector_validated


        function is_int_matrix_validated(
     $     int_matrix,
     $     int_matrix_cst,
     $     detailled)
     $     result(test_validated)

          implicit none

          integer, dimension(:,:), intent(in) :: int_matrix
          integer, dimension(:,:), intent(in) :: int_matrix_cst
          logical, optional      , intent(in) :: detailled
          logical                             :: test_validated


          integer :: i,j
          logical :: test_loc
          logical :: detailled_op


          if(present(detailled)) then
             detailled_op = detailled
          else
             detailled_op = .false.
          end if

          test_validated = .true.

          if((size(int_matrix,1).eq.size(int_matrix_cst,1)).and.
     $       (size(int_matrix,2).eq.size(int_matrix_cst,2))) then

             do j=1, size(int_matrix,2)
                do i=1, size(int_matrix,1)
                   
                   test_loc = int_matrix(i,j).eq.int_matrix_cst(i,j)
                   
                   if(detailled_op.and.(.not.test_loc)) then
                      print '(''['',2I3'']:'',I5,'' -> '',I5)', 
     $                     i,j,
     $                     int_matrix(i,j),
     $                     int_matrix_cst(i,j)
                   end if
                
                   test_validated = test_validated.and.test_loc
                   
                end do
             end do

          else

             test_validated = .false.
             print '(''sizes do not match'')'
             print '(''  - size_x : '',I2,'' -> '',I2)', size(int_matrix,1), size(int_matrix_cst,1)
             print '(''  - size_y : '',I2,'' -> '',I2)', size(int_matrix,2), size(int_matrix_cst,2)
             print '()'
          end if

        end function is_int_matrix_validated


        function is_int_matrix3D_validated(
     $     int_matrix,
     $     int_matrix_cst,
     $     detailled)
     $     result(test_validated)

          implicit none

          integer, dimension(:,:,:), intent(in) :: int_matrix
          integer, dimension(:,:,:), intent(in) :: int_matrix_cst
          logical, optional        , intent(in) :: detailled
          logical                               :: test_validated


          integer :: i,j,k
          logical :: test_loc
          logical :: detailled_op


          if(present(detailled)) then
             detailled_op = detailled
          else
             detailled_op = .false.
          end if

          test_validated = .true.

          if((size(int_matrix,1).eq.size(int_matrix_cst,1)).and.
     $         (size(int_matrix,2).eq.size(int_matrix_cst,2)).and.
     $         (size(int_matrix,3).eq.size(int_matrix_cst,3))) then

             do k=1, size(int_matrix,3)
                do j=1, size(int_matrix,2)
                   do i=1, size(int_matrix,1)
                      
                      test_loc = int_matrix(i,j,k).eq.int_matrix_cst(i,j,k)
                      
                      if(detailled_op.and.(.not.test_loc)) then
                         print '(''['',3I3'']:'',I5,'' -> '',I5)',
     $                        i,j,k,
     $                        int_matrix(i,j,k),
     $                        int_matrix_cst(i,j,k)
                      end if
                      
                      test_validated = test_validated.and.test_loc
                      
                   end do
                end do
             end do

          else

             test_validated = .false.
             print '(''sizes do not match'')'
             print '(''  - size_x : '',I2,'' -> '',I2)', size(int_matrix,1), size(int_matrix_cst,1)
             print '(''  - size_y : '',I2,'' -> '',I2)', size(int_matrix,2), size(int_matrix_cst,2)
             print '(''  - size_z : '',I2,'' -> '',I2)', size(int_matrix,3), size(int_matrix_cst,3)
             print '()'

          end if

        end function is_int_matrix3D_validated

      end module check_data_module
