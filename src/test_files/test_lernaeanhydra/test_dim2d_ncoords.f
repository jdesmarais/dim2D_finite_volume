      program test_dim2d_ncoords

        use dim2d_ncoords_module, only :
     $       compute_n1_eigenvalues_dim2d,
     $       compute_n2_eigenvalues_dim2d,
     $       compute_n1_lefteigenvector_dim2d,
     $       compute_n1_righteigenvector_dim2d,
     $       compute_n2_lefteigenvector_dim2d,
     $       compute_n2_righteigenvector_dim2d,
     $       compute_n1_transM_dim2d,
     $       compute_n2_transM_dim2d

        use parameters_input, only :
     $       ne
        
        use parameters_kind, only :
     $       ikind,
     $       rkind


        implicit none


        logical :: test_loc
        logical :: test_validated
        logical :: detailled

        detailled = .true.


        test_loc = test_compute_n1_eigenvalues(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_n1_eigenvalues: '',L1)', test_loc
        print '()'


        test_loc = test_compute_n2_eigenvalues(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_n2_eigenvalues: '',L1)', test_loc
        print '()'


        test_loc = test_compute_n1_lefteigenvector(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_n1_lefteigenvector: '',L1)', test_loc
        print '()'


        test_loc = test_compute_n1_righteigenvector(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_n1_righteigenvector: '',L1)', test_loc
        print '()'


        test_loc = test_compute_n2_lefteigenvector(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_n2_lefteigenvector: '',L1)', test_loc
        print '()'


        test_loc = test_compute_n2_righteigenvector(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_n2_righteigenvector: '',L1)', test_loc
        print '()'


        test_loc = test_compute_n1_transM(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_n1_transM: '',L1)', test_loc
        print '()'


        test_loc = test_compute_n2_transM(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_n2_transM: '',L1)', test_loc
        print '()'

        
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


        function is_matrix_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), dimension(ne,ne), intent(in) :: var
          real(rkind), dimension(ne,ne), intent(in) :: cst
          logical                      , intent(in) :: detailled
          logical                                   :: test_validated

          logical :: test_loc
          integer :: i,j

          test_validated = .true.

          do j=1,ne
             do i=1,ne
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


        function test_compute_n1_eigenvalues(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          logical :: test_loc
          integer :: i

          real(rkind), dimension(ne) :: nodes
          real(rkind), dimension(ne) :: eigenvalues
          real(rkind), dimension(ne) :: test_eigenvalues

          test_validated = .true.


          nodes            = [4.2d0,5.8d0,2.6d0,2.2d0]
          eigenvalues      = compute_n1_eigenvalues_dim2d(nodes)
          test_eigenvalues = [0.538748024d0,0.538748024,-3.550921501d0,4.628417549d0]


          do i=1,ne

             test_loc = is_test_validated(
     $            eigenvalues(i),
     $            test_eigenvalues(i),
     $            .false.)
             if(detailled.and.(.not.test_loc)) then
                print '(''['',I2,'']: '',F10.3,'' -> '', F10.3)',
     $               i,
     $               eigenvalues(i),
     $               test_eigenvalues(i)
             end if

             test_validated = test_validated.and.test_loc

          end do

        end function test_compute_n1_eigenvalues


        function test_compute_n2_eigenvalues(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          logical :: test_loc
          integer :: i

          real(rkind), dimension(ne) :: nodes
          real(rkind), dimension(ne) :: eigenvalues
          real(rkind), dimension(ne) :: test_eigenvalues

          test_validated = .true.


          nodes            = [4.2d0,5.8d0,2.6d0,2.2d0]
          eigenvalues      = compute_n2_eigenvalues_dim2d(nodes)
          test_eigenvalues = [1.414213562d0,1.414213562d0,-2.675455963d0,5.503883088d0]


          do i=1,ne

             test_loc = is_test_validated(
     $            eigenvalues(i),
     $            test_eigenvalues(i),
     $            .false.)
             if(detailled.and.(.not.test_loc)) then
                print '(''['',I2,'']: '',F10.3,'' -> '', F10.3)',
     $               i,
     $               eigenvalues(i),
     $               test_eigenvalues(i)
             end if

             test_validated = test_validated.and.test_loc

          end do

        end function test_compute_n2_eigenvalues


        function test_compute_n1_lefteigenvector(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: eigenvect
          real(rkind), dimension(ne,ne) :: test_eigenvect

          test_validated = .true.


          nodes            = [4.2d0,5.8d0,2.6d0,2.2d0]

          eigenvect        = compute_n1_lefteigenvector_dim2d(nodes)

          test_eigenvect = reshape((/
     $         -0.238095238d0,  0.119047619d0,  0.119047619d0,  0.0d0,
     $          1.575156930d0, -0.082566195d0, -0.037012432d0,  0.059789314d0,
     $         -3.708213258d0, -0.755440337d0,	1.755440337d0, -0.5d0,
     $         -5.911514633d0,  2.136392718d0, -1.136392718d0, -0.5d0/),
     $         (/4,4/))

          test_validated = is_matrix_validated(
     $         eigenvect,
     $         test_eigenvect,
     $         detailled)

        end function test_compute_n1_lefteigenvector


        function test_compute_n1_righteigenvector(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: eigenvect
          real(rkind), dimension(ne,ne) :: test_eigenvect

          test_validated = .true.


          nodes            = [4.2d0,5.8d0,2.6d0,2.2d0]

          eigenvect        = compute_n1_righteigenvector_dim2d(nodes)

          test_eigenvect = reshape((/
     $         0.0d0,  1.0d0        ,  0.059789314d0,  0.059789314d0,
     $         4.2d0,  1.380952381d0, -0.090334519d0,  0.255466909d0,
     $         4.2d0,  0.619047619d0,  0.209913146d0, -0.135888282d0,
     $         8.4d0, -7.329478458d0, -1.569958365d0, -1.306490610d0/),
     $         (/4,4/))

          test_validated = is_matrix_validated(
     $         eigenvect,
     $         test_eigenvect,
     $         detailled)

        end function test_compute_n1_righteigenvector


        function test_compute_n2_lefteigenvector(detailled)
     $     result(test_validated)
        
          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: eigenvect
          real(rkind), dimension(ne,ne) :: test_eigenvect

          test_validated = .true.


          nodes            = [4.2d0,5.8d0,2.6d0,2.2d0]

          eigenvect        = compute_n2_lefteigenvector_dim2d(nodes)

          test_eigenvect = reshape((/
     $         0.090702948d0, -0.119047619d0,  0.119047619d0,  0.0d0,
     $         1.575156930d0, -0.082566195d0, -0.037012432d0,  0.059789314d0,
     $        -1.918030891d0, -0.755440337d0, -1.136392718d0, -0.5d0,
     $        -7.701697000d0,  2.136392718d0,  1.755440337d0, -0.5d0/),
     $         (/4,4/))

          test_validated = is_matrix_validated(
     $         eigenvect,
     $         test_eigenvect,
     $         detailled)

        end function test_compute_n2_lefteigenvector


        function test_compute_n2_righteigenvector(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: eigenvect
          real(rkind), dimension(ne,ne) :: test_eigenvect

          test_validated = .true.


          nodes            = [4.2d0,5.8d0,2.6d0,2.2d0]

          eigenvect        = compute_n2_righteigenvector_dim2d(nodes)

          test_eigenvect = reshape((/
     $         0.0d0,  1.0d0        ,  0.059789314d0,  0.059789314d0,
     $        -4.2d0,  1.380952381d0, -0.090334519d0,  0.255466909d0,
     $         4.2d0,  0.619047619d0, -0.135888282d0,  0.209913146d0,
     $        -3.2d0, -7.329478458d0, -1.784025916d0, -1.092423060d0/),
     $         (/4,4/))

          test_validated = is_matrix_validated(
     $         eigenvect,
     $         test_eigenvect,
     $         detailled)

        end function test_compute_n2_righteigenvector


        function test_compute_n1_transM(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: eigenvect
          real(rkind), dimension(ne,ne) :: test_eigenvect

          test_validated = .true.


          nodes            = [4.2d0,5.8d0,2.6d0,2.2d0]

          eigenvect        = compute_n1_transM_dim2d(nodes)

          test_eigenvect = reshape((/
     $         5.55112e-17  ,  0.707106781d0,  0.707106781d0,  0.0d0        ,
     $        -8.755136411d0,  3.367175149d0,  1.414213562d0, -0.707106781d0,
     $        -7.677640364d0,  1.414213562d0,  2.289679101d0, -0.707106781d0,
     $        20.414381220d0,-15.056403850d0,-16.133899890d0,  0.0d0/),
     $         (/4,4/))

          test_validated = is_matrix_validated(
     $         eigenvect,
     $         test_eigenvect,
     $         detailled)

        end function test_compute_n1_transM


        function test_compute_n2_transM(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(ne)    :: nodes
          real(rkind), dimension(ne,ne) :: eigenvect
          real(rkind), dimension(ne,ne) :: test_eigenvect

          test_validated = .true.


          nodes            = [4.2d0,5.8d0,2.6d0,2.2d0]

          eigenvect        = compute_n2_transM_dim2d(nodes)

          test_eigenvect = reshape((/
     $         5.55112E-17  ,  0.707106781d0, -0.707106781d0,  0.0d0        ,
     $        -7.546160191d0,  2.491709610d0, -0.538748024d0, -0.707106781d0,
     $         6.468664144d0, -0.538748024d0, -0.336717515d0,  0.707106781d0,
     $         7.776907130d0,-16.265380070d0, 17.342876110d0,  0.0d0/),
     $         (/4,4/))

          test_validated = is_matrix_validated(
     $         eigenvect,
     $         test_eigenvect,
     $         detailled)

        end function test_compute_n2_transM

      end program test_dim2d_ncoords
