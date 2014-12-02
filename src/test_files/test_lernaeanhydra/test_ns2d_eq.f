      !test the procedure implemented in ns2d/pmodel_eq_class
      program test_ns2d_eq
      
        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq


        implicit none


        logical :: detailled
        logical :: test_validated


        detailled = .false.


        test_validated = test_compute_x_eigenvalues(detailled)
        print '(''test_compute_x_eigenvalues: '',L1)', test_validated
        print '()'

        test_validated = test_compute_y_eigenvalues(detailled)
        print '(''test_compute_y_eigenvalues: '',L1)', test_validated
        print '()'

        test_validated = test_compute_x_lefteigenvector(detailled)
        print '(''test_compute_x_lefteigenvector: '',L1)', test_validated
        print '()'

        test_validated = test_compute_x_righteigenvector(detailled)
        print '(''test_compute_x_righteigenvector: '',L1)', test_validated
        print '()'

        test_validated = test_compute_y_lefteigenvector(detailled)
        print '(''test_compute_y_lefteigenvector: '',L1)', test_validated
        print '()'

        test_validated = test_compute_y_righteigenvector(detailled)
        print '(''test_compute_y_righteigenvector: '',L1)', test_validated
        print '()'

        test_validated = test_compute_x_transM(detailled)
        print '(''test_compute_x_transM: '',L1)', test_validated
        print '()'

        test_validated = test_compute_y_transM(detailled)
        print '(''test_compute_y_transM: '',L1)', test_validated
        print '()'

        detailled = .false.

        test_validated = test_compute_n1_lefteigenvector(detailled)
        print '(''test_compute_n1_lefteigenvector: '',L1)', test_validated
        print '()'

        test_validated = test_compute_n1_righteigenvector(detailled)
        print '(''test_compute_n1_righteigenvector: '',L1)', test_validated
        print '()'

        test_validated = test_compute_n2_lefteigenvector(detailled)
        print '(''test_compute_n2_lefteigenvector: '',L1)', test_validated
        print '()'

        test_validated = test_compute_n2_righteigenvector(detailled)
        print '(''test_compute_n2_righteigenvector: '',L1)', test_validated
        print '()'

        test_validated = test_compute_n1_transM(detailled)
        print '(''test_compute_n1_transM: '',L1)', test_validated
        print '()'

        test_validated = test_compute_n2_transM(detailled)
        print '(''test_compute_n2_transM: '',L1)', test_validated
        print '()'


        detailled = .true.

        contains

        function test_compute_x_eigenvalues(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne) :: nodes_test
          real(rkind), dimension(ne) :: eigenvalues_test
          real(rkind), dimension(ne) :: eigenvalues
          type(pmodel_eq)            :: p_model
          integer                    :: k
          logical                    :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvalues_test = [0.316326531d0, 0.316326531d0, -0.316123056d0, 0.948776118d0]
          eigenvalues      = p_model%compute_x_eigenvalues(nodes_test)
          
          
          test_validated = .true.

          do k=1, ne
             test_loc = is_test_validated(eigenvalues(k),eigenvalues_test(k),detailled)
             test_validated = test_validated.and.test_loc
          end do

        end function test_compute_x_eigenvalues


        function test_compute_y_eigenvalues(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne) :: nodes_test
          real(rkind), dimension(ne) :: eigenvalues_test
          real(rkind), dimension(ne) :: eigenvalues
          type(pmodel_eq)            :: p_model
          integer                    :: k
          logical                    :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvalues_test = [0.739795918d0, 0.739795918d0, 0.107346331d0, 1.372245505d0]
          eigenvalues      = p_model%compute_y_eigenvalues(nodes_test)
          
          
          test_validated = .true.

          do k=1, ne
             test_loc = is_test_validated(eigenvalues(k),eigenvalues_test(k),detailled)
             test_validated = test_validated.and.test_loc
          end do

        end function test_compute_y_eigenvalues


        function test_compute_x_lefteigenvector(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne)    :: nodes_test
          real(rkind), dimension(ne,ne) :: eigenvector_test
          real(rkind), dimension(ne,ne) :: eigenvector
          type(pmodel_eq)               :: p_model
          integer                       :: k,l
          logical                       :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvector_test = reshape((/
     $         -0.075489379d0, 0.0d0        ,  0.102040816d0,  0.0d0,
     $          0.460522795d0, 0.527220796d0,  1.233016378d0, -1.666698001d0,
     $          0.207923704d0,-0.42166697d0 , -0.246598639d0,  0.333333333d0,
     $          0.007863121d0, 0.210782617d0, -0.246598639d0,  0.333333333d0
     $         /),
     $         (/ne,ne/))
          eigenvector      = p_model%compute_x_lefteigenvector(nodes_test)
          
          
          test_validated = .true.

          do k=1, ne
             do l=1,ne
                test_loc = is_test_validated(eigenvector(k,l),eigenvector_test(k,l),detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_compute_x_lefteigenvector


        function test_compute_x_righteigenvector(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne)    :: nodes_test
          real(rkind), dimension(ne,ne) :: eigenvector_test
          real(rkind), dimension(ne,ne) :: eigenvector
          type(pmodel_eq)               :: p_model
          integer                       :: k,l
          logical                       :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvector_test = reshape((/
     $         0.0d0, 1.0d0        ,  2.500047001d0, 2.500047001d0,
     $         0.0d0, 0.316326531d0, -0.790322499d0, 2.371984887d0,
     $         9.8d0, 0.739795918d0,  1.849524567d0, 1.849524567d0,
     $        7.25d0, 0.323680237d0,  1.809054945d0, 2.809376669d0
     $         /),
     $         (/ne,ne/))
          eigenvector      = p_model%compute_x_righteigenvector(nodes_test)
          
          
          test_validated = .true.

          do k=1, ne
             do l=1,ne
                test_loc = is_test_validated(eigenvector(k,l),eigenvector_test(k,l),detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_compute_x_righteigenvector


        function test_compute_y_lefteigenvector(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne)    :: nodes_test
          real(rkind), dimension(ne,ne) :: eigenvector_test
          real(rkind), dimension(ne,ne) :: eigenvector
          type(pmodel_eq)               :: p_model
          integer                       :: k,l
          logical                       :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvector_test = reshape((/
     $         -0.032278217d0,  0.102040816d0,   0.0d0        ,  0.0d0        ,
     $          0.460522795d0,  0.527220796d0,   1.233016378d0, -1.666698001d0,
     $          0.341835224d0, -0.105442177d0,  -0.562823433d0,  0.333333333d0,
     $         -0.126048399d0, -0.105442177d0,   0.069626154d0,  0.333333333d0
     $         /),
     $         (/ne,ne/))
          eigenvector      = p_model%compute_y_lefteigenvector(nodes_test)
          
          
          test_validated = .true.

          do k=1, ne
             do l=1,ne
                test_loc = is_test_validated(eigenvector(k,l),eigenvector_test(k,l),detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_compute_y_lefteigenvector


        function test_compute_y_righteigenvector(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne)    :: nodes_test
          real(rkind), dimension(ne,ne) :: eigenvector_test
          real(rkind), dimension(ne,ne) :: eigenvector
          type(pmodel_eq)               :: p_model
          integer                       :: k,l
          logical                       :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvector_test = reshape((/
     $         0.0d0, 1.0d0        , 2.500047001d0, 2.500047001d0,
     $         9.8d0, 0.316326531d0, 0.790831194d0, 0.790831194d0,
     $         0.0d0, 0.739795918d0, 0.268370874d0, 3.430678260d0,
     $         3.1d0, 0.323680237d0, 1.139484758d0, 3.478946855d0
     $         /),
     $         (/ne,ne/))
          eigenvector      = p_model%compute_y_righteigenvector(nodes_test)
          
          
          test_validated = .true.

          do k=1, ne
             do l=1,ne
                test_loc = is_test_validated(eigenvector(k,l),eigenvector_test(k,l),detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_compute_y_righteigenvector


        function test_compute_x_transM(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne)    :: nodes_test
          real(rkind), dimension(ne,ne) :: eigenvector_test
          real(rkind), dimension(ne,ne) :: eigenvector
          type(pmodel_eq)               :: p_model
          integer                       :: k,l
          logical                       :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvector_test = reshape((/
     $         0.0d0	    ,  0.0d0        , 1.0d0	   , 0.0d0        ,
     $        -0.234017076d0,  0.739795918d0, 0.316326531d0, 0.0d0        ,
     $        -0.331511176d0, -0.210884354d0, 0.986394558d0, 0.666666667d0,
     $        -0.523688312d0, -0.156011384d0, 0.558803623d0, 1.232993197d0
     $         /),
     $         (/ne,ne/))
          eigenvector      = p_model%compute_x_transM(nodes_test)
          
          
          test_validated = .true.

          do k=1, ne
             do l=1,ne
                test_loc = is_test_validated(eigenvector(k,l),eigenvector_test(k,l),detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_compute_x_transM


        function test_compute_y_transM(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne)    :: nodes_test
          real(rkind), dimension(ne,ne) :: eigenvector_test
          real(rkind), dimension(ne,ne) :: eigenvector
          type(pmodel_eq)               :: p_model
          integer                       :: k,l
          logical                       :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvector_test = reshape((/
     $         0.0d0        , 1.0d0        ,  0.0d0	   , 0.0d0        ,
     $         0.115724351d0, 0.421768707d0, -0.493197279d0, 0.666666667d0,
     $        -0.234017076d0, 0.739795918d0,  0.316326531d0, 0.0d0        ,
     $        -0.223921899d0, 0.856960641d0, -0.156011384d0, 0.527210884d0
     $         /),
     $         (/ne,ne/))
          eigenvector      = p_model%compute_y_transM(nodes_test)
          
          
          test_validated = .true.

          do k=1, ne
             do l=1,ne
                test_loc = is_test_validated(eigenvector(k,l),eigenvector_test(k,l),detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_compute_y_transM


        function test_compute_n1_lefteigenvector(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne)    :: nodes_test
          real(rkind), dimension(ne,ne) :: eigenvector_test
          real(rkind), dimension(ne,ne) :: eigenvector
          type(pmodel_eq)               :: p_model
          integer                       :: k,l
          logical                       :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvector_test = reshape((/
     $         -0.053883798d0,  0.051020408d0,  0.051020408d0,  0.0d0,
     $          0.460522795d0,  0.527220796d0,  1.233016378d0, -1.666698001d0,
     $          0.013203669d0, -0.329046873d0, -0.022993944d0,  0.333333333d0,
     $          0.202583156d0,  0.118162519d0, -0.470203335d0,  0.333333333d0
     $         /),
     $         (/ne,ne/))
          eigenvector      = p_model%compute_n1_lefteigenvector(nodes_test)
          
          
          test_validated = .true.

          do k=1, ne
             do l=1,ne
                test_loc = is_test_validated(eigenvector(k,l),eigenvector_test(k,l),detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_compute_n1_lefteigenvector


        function test_compute_n1_righteigenvector(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne)    :: nodes_test
          real(rkind), dimension(ne,ne) :: eigenvector_test
          real(rkind), dimension(ne,ne) :: eigenvector
          type(pmodel_eq)               :: p_model
          integer                       :: k,l
          logical                       :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvector_test = reshape((/
     $         0.0d0, 1.0d0        ,  2.500047001d0, 2.500047001d0,
     $         9.8d0, 0.316326531d0, -0.327213304d0, 1.908875693d0,
     $         9.8d0, 0.739795918d0,  2.967569065d0, 0.731480069d0,
     $       10.35d0, 0.323680237d0,  2.782673426d0, 1.835758188d0
     $         /),
     $         (/ne,ne/))
          eigenvector      = p_model%compute_n1_righteigenvector(nodes_test)
          
          
          test_validated = .true.

          do k=1, ne
             do l=1,ne
                test_loc = is_test_validated(eigenvector(k,l),eigenvector_test(k,l),detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_compute_n1_righteigenvector


        function test_compute_n2_lefteigenvector(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne)    :: nodes_test
          real(rkind), dimension(ne,ne) :: eigenvector_test
          real(rkind), dimension(ne,ne) :: eigenvector
          type(pmodel_eq)               :: p_model
          integer                       :: k,l
          logical                       :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvector_test = reshape((/
     $        -0.021605581d0, -0.051020408d0,  0.051020408d0,  0.0d0        ,
     $         0.460522795d0,  0.527220796d0,  1.233016378d0, -1.666698001d0,
     $         0.344047351d0, -0.329046873d0, -0.470203335d0,  0.333333333d0,
     $        -0.128260526d0,  0.118162519d0, -0.022993944d0,  0.333333333d0
     $         /),
     $         (/ne,ne/))
          eigenvector      = p_model%compute_n2_lefteigenvector(nodes_test)
          
          
          test_validated = .true.

          do k=1, ne
             do l=1,ne
                test_loc = is_test_validated(eigenvector(k,l),eigenvector_test(k,l),detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_compute_n2_lefteigenvector


        function test_compute_n2_righteigenvector(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne)    :: nodes_test
          real(rkind), dimension(ne,ne) :: eigenvector_test
          real(rkind), dimension(ne,ne) :: eigenvector
          type(pmodel_eq)               :: p_model
          integer                       :: k,l
          logical                       :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvector_test = reshape((/
     $          0.0d0, 1.0d0        ,  2.500047001d0, 2.500047001d0,
     $         -9.8d0, 0.316326531d0, -0.327213304d0, 1.908875693d0,
     $          9.8d0, 0.739795918d0,  0.731480069d0, 2.967569065d0,
     $         4.15d0, 0.323680237d0,  1.128423913d0, 3.490007701d0
     $         /),
     $         (/ne,ne/))
          eigenvector      = p_model%compute_n2_righteigenvector(nodes_test)
          
          
          test_validated = .true.

          do k=1, ne
             do l=1,ne
                test_loc = is_test_validated(eigenvector(k,l),eigenvector_test(k,l),detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_compute_n2_righteigenvector


        function test_compute_n1_transM(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne)    :: nodes_test
          real(rkind), dimension(ne,ne) :: eigenvector_test
          real(rkind), dimension(ne,ne) :: eigenvector
          type(pmodel_eq)               :: p_model
          integer                       :: k,l
          logical                       :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvector_test = reshape((/
     $         0.0d0	    , 0.707106781d0,  0.707106781d0, 0.0d0        ,
     $        -0.083645588d0, 0.821350224d0, -0.125066506d0, 0.471404521d0,
     $        -0.399888862d0, 0.373996954d0,  0.921162916d0, 0.471404521d0,
     $        -0.528640250d0, 0.495645973d0,  0.284817124d0, 1.244652242d0
     $         /),
     $         (/ne,ne/))
          eigenvector      = p_model%compute_n1_transM(nodes_test)
          
          test_validated = .true.

          do k=1, ne
             do l=1,ne
                test_loc = is_test_validated(eigenvector(k,l),eigenvector_test(k,l),detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_compute_n1_transM


        function test_compute_n2_transM(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          real(rkind), dimension(ne)    :: nodes_test
          real(rkind), dimension(ne,ne) :: eigenvector_test
          real(rkind), dimension(ne,ne) :: eigenvector
          type(pmodel_eq)               :: p_model
          integer                       :: k,l
          logical                       :: test_loc


          nodes_test       = [9.8d0,3.1d0,7.25d0,6.7d0]
          eigenvector_test = reshape((/
     $         0.0d0        , 0.707106781d0, -0.707106781d0, 0.0d0        ,
     $         0.247304535d0,-0.224879197d0, -0.572419775d0, 0.471404521d0,
     $         0.068938739d0, 0.672232467d0, -0.473809646d0,-0.471404521d0,
     $         0.211966864d0, 0.716279388d0, -0.505450539d0,-0.499063460d0
     $         /),
     $         (/ne,ne/))
          eigenvector      = p_model%compute_n2_transM(nodes_test)
          
          test_validated = .true.

          do k=1, ne
             do l=1,ne
                test_loc = is_test_validated(eigenvector(k,l),eigenvector_test(k,l),detailled)
                test_validated = test_validated.and.test_loc
             end do
          end do

        end function test_compute_n2_transM


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
     $         int(var*1e5)-
     $         int(cst*1e5)).le.1
          
        end function is_test_validated

      end program test_ns2d_eq
