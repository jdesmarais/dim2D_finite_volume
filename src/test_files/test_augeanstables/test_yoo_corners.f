      program test_yoo_corners

        use lodi_corner_inflow_inflow_class, only :
     $     get_lodi_A_inflow_inflow

        use parameters_kind, only :
     $       rkind

        implicit none

        logical :: test_validated
        logical :: detailled


        print '(''test get_lodi_A'')'
        print '(''---------------------------------------'')'

        detailled = .false.

        test_validated = test_get_lodi_A(detailled)

        print '()'


        contains


        function is_test_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical    , intent(in) :: detailled
          logical                 :: test_validated

          if(detailled) then
             print *, int(var*1e5)
             print *, int(cst*1e5)
          end if
          
          test_validated=abs(
     $         int(var*10000.)-
     $         sign(int(abs(cst*10000.)),int(cst*10000.))).le.1
          
        end function is_test_validated


        function is_matrix_validated(matrix,test_data,detailled)
     $     result(test_validated)
        
          implicit none

          real(rkind), dimension(:,:), intent(in) :: matrix
          real(rkind), dimension(:,:), intent(in) :: test_data
          logical                    , intent(in) :: detailled
          logical                                 :: test_validated

          integer :: i,j
          logical :: test_loc

          test_validated = .true.

          do j=1, size(matrix,2)
             do i=1, size(matrix,1)
                
                test_loc = is_test_validated(
     $               matrix(i,j),
     $               test_data(i,j),
     $               detailled)

                test_validated = test_validated.and.test_loc
                if(detailled) then
                   print '(''('',2I2,''):'',L2)',i,j,test_loc
                end if

             end do
          end do

        end function is_matrix_validated


        function test_get_lodi_A(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(6,6) :: lodi_A

          real(rkind), dimension(6,6) :: test_data
          real(rkind)                 :: md
          real(rkind)                 :: c
          real(rkind)                 :: tau
          integer                     :: sign_x
          integer                     :: sign_y

          !inputs
          md     =  0.3d0
          c      =  0.5d0
          tau    =  0.2d0
          sign_x =  1
          sign_y = -1

          !fct tested
          lodi_A = get_lodi_A_inflow_inflow(
     $         md,c,sign_x,sign_y,tau)
          
          !data output expected
          test_data(1,1) =  1.71957672d0
          test_data(2,1) =  0.0d0
          test_data(3,1) = -3.527336861d0
          test_data(4,1) =  0.423280423d0
          test_data(5,1) =  0.0d0
          test_data(6,1) =  5.996472663d0

          test_data(1,2) =  0.0d0
          test_data(2,2) =  2.777777778d0
          test_data(3,2) =  0.0d0
          test_data(4,2) =  0.0d0
          test_data(5,2) =  -2.22222222d0
          test_data(6,2) =  0.0d0

          test_data(1,3) = -0.158730159d0
          test_data(2,3) =  0.0d0
          test_data(3,3) =  2.248677249d0
          test_data(4,3) = -0.26984127d0
          test_data(5,3) =  0.0d0
          test_data(6,3) = -1.322751323d0

          test_data(1,4) =  0.423280423d0
          test_data(2,4) =  0.0d0
          test_data(3,4) = -5.996472663d0
          test_data(4,4) =  1.71957672d0
          test_data(5,4) =  0.0d0
          test_data(6,4) =  3.527336861d0

          test_data(1,5) =  0.0d0
          test_data(2,5) = -2.222222222d0
          test_data(3,5) =  0.0d0
          test_data(4,5) =  0.0d0
          test_data(5,5) =  2.777777778d0
          test_data(6,5) =  0.0d0

          test_data(1,6) =  0.26984127d0
          test_data(2,6) =  0.0d0
          test_data(3,6) = -1.322751323d0
          test_data(4,6) =  0.158730159d0
          test_data(5,6) =  0.0d0
          test_data(6,6) =  2.248677249d0

          !test the output
          test_validated = is_matrix_validated(
     $         lodi_A,
     $         test_data,
     $         detailled)

          if(.not.detailled) then
             print '(''test_get_lodi_A: '',L2)', test_validated
          end if

        end function test_get_lodi_A

      end program test_yoo_corners
