      program test_ns_vdw2d_prim

        use check_data_module, only :
     $      is_real_validated,
     $      is_real_vector_validated,
     $      is_real_matrix_validated

        use dim2d_parameters, only :
     $       cv_r

        use ns_vdw2d_prim_module, only :
     $       compute_prim_var_ns_vdw2d,
     $       compute_cons_var_ns_vdw2d,
     $       compute_jacobian_prim_to_cons_ns_vdw2d,
     $       compute_jacobian_cons_to_prim_ns_vdw2d,
     $       compute_x_transM_ns_vdw2d,
     $       compute_y_transM_ns_vdw2d,
     $       compute_x_eigenvalues_ns_vdw2d,
     $       compute_y_eigenvalues_ns_vdw2d,
     $       compute_x_lefteigenvector_ns_vdw2d,
     $       compute_x_righteigenvector_ns_vdw2d,
     $       compute_y_lefteigenvector_ns_vdw2d,
     $       compute_y_righteigenvector_ns_vdw2d

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind


        implicit none

        logical :: test_loc
        logical :: test_validated
        logical :: detailled        

        
        detailled = .false.


        !test the input parameters
        call test_inputs()

        !test_compute_prim_var
        test_loc = test_compute_prim_var(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_prim_var: '',L1)', test_loc
        print '()'

        !test_compute_var_prim
        test_loc = test_compute_prim_var(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_prim_var: '',L1)', test_loc
        print '()'

        !test_compute_jacobian_prim_to_cons
        test_loc = test_compute_jacobian_prim_to_cons(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_jacobian_prim_to_cons: '',L1)', test_loc
        print '()'
        
        !test_compute_jacobian_cons_to_prim
        test_loc = test_compute_jacobian_cons_to_prim(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_jacobian_cons_to_prim: '',L1)', test_loc
        print '()'

        !test_x_transM
        test_loc = test_compute_x_transM(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_x_transM: '',L1)', test_loc
        print '()'
        
        !test_y_transM
        test_loc = test_compute_y_transM(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_y_transM: '',L1)', test_loc
        print '()'

        !test_compute_x_eigenvalues
        test_loc = test_compute_x_eigenvalues(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_x_eigenvalues: '',L1)', test_loc
        print '()'
        
        !test_compute_y_eigenvalues
        test_loc = test_compute_y_eigenvalues(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_y_eigenvalues: '',L1)', test_loc
        print '()'

        !test_compute_x_lefteigenvector
        test_loc = test_compute_x_lefteigenvector(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_x_lefteigenvector: '',L1)', test_loc
        print '()'

        !test_compute_x_righteigenvector
        test_loc = test_compute_x_righteigenvector(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_x_righteigenvector: '',L1)', test_loc
        print '()'

        !test_compute_y_lefteigenvector
        test_loc = test_compute_y_lefteigenvector(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_y_lefteigenvector: '',L1)', test_loc
        print '()'

        !test_compute_y_righteigenvector
        test_loc = test_compute_y_righteigenvector(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_y_righteigenvector: '',L1)', test_loc
        print '()'


        contains

        !test the inputs
        subroutine test_inputs()
        
          implicit none

          logical :: test_parameters

          test_parameters=.true.

          test_parameters=test_parameters.and.(ne.eq.4)
          test_parameters=test_parameters.and.is_real_validated(cv_r,2.5d0,.false.)

          if(.not.test_parameters) then

             !< print the dim2d parameters used for the test
             print '(''WARNING: this test is designed for:'')'
             print '(''ne            : '', I2)', 4
             print '()'
             print '(''cv_r          : '', F16.6)', 2.5
             print '()'
             print '(''it allows to see errors easily'')'
             print '('''')'

             stop 'dim2d_parameters not adapted for test'
          end if

        end subroutine test_inputs


        function test_compute_prim_var(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          test_validated = is_real_vector_validated(
     $         compute_prim_var_ns_vdw2d([4.2d0,5.8d0,2.6d0,2.2d0]),
     $         [4.2d0,1.380952381d0,0.619047619d0,-103.2304762d0],
     $         detailled)


        end function test_compute_prim_var


        function test_compute_cons_var(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          test_validated = is_real_vector_validated(
     $         compute_cons_var_ns_vdw2d([4.2d0,1.380952381d0,0.619047619d0,-103.2304762d0]),
     $         [4.2d0,5.8d0,2.6d0,2.2d0],
     $         detailled)


        end function test_compute_cons_var


        function test_compute_jacobian_prim_to_cons(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          test_validated = is_real_matrix_validated(
     $         compute_jacobian_prim_to_cons_ns_vdw2d(
     $             4.2d0,
     $             1.380952381d0,
     $             0.619047619d0,
     $            -103.2304762d0),
     $         reshape((/
     $         1.0d0        ,  0.0d0        ,  0.0d0        ,  0.0d0,
     $        -0.328798186d0,  0.238095238d0,  0.0d0        ,  0.0d0,
     $        -0.147392290d0,  0.0d0        ,  0.238095238d0,  0.0d0,
     $        -9.619727891d0,  1.380952381d0,  0.619047610d0, -1.0d0/),
     $         (/4,4/)),
     $         detailled)

        end function test_compute_jacobian_prim_to_cons


        function test_compute_jacobian_cons_to_prim(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          test_validated = is_real_matrix_validated(
     $         compute_jacobian_cons_to_prim_ns_vdw2d(
     $             4.2d0,
     $             1.380952381d0,
     $             0.619047619d0,
     $            -103.2304762d0),
     $         reshape((/
     $         1.0d0        ,  0.0d0,   0.0d0,   0.0d0,
     $         1.380952381d0,  4.2d0,   0.0d0,   0.0d0,
     $         0.619047619d0,  0.0d0,   4.2d0,   0.0d0,
     $        -7.329478458d0,  5.8d0,   2.6d0,  -1.0d0/),
     $         (/4,4/)),
     $         detailled)


        end function test_compute_jacobian_cons_to_prim


        function test_compute_x_transM(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          test_validated = is_real_matrix_validated(
     $         compute_x_transM_ns_vdw2d(
     $             4.2d0,
     $             0.619047619d0,
     $             4.089669525d0),
     $         reshape((/
     $         0.619047619d0, 0.0d0        , 4.2d0        , 0.0d0,
     $         0.0d0        , 0.619047619d0, 0.0d0        , 0.0d0,
     $         0.0d0	    , 0.0d0        , 0.619047619d0, 0.238095238d0,
     $         0.0d0        , 0.0d0        , 70.24666667d0, 0.619047619d0/),
     $         (/4,4/)),
     $         detailled)

        end function test_compute_x_transM


        function test_compute_y_transM(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          test_validated = is_real_matrix_validated(
     $         compute_y_transM_ns_vdw2d(
     $             4.2d0,
     $             1.380952381d0,
     $             4.089669525d0),
     $         reshape((/
     $         1.380952381d0, 4.2d0	   , 0.0d0, 	    0.0d0,
     $         0.0d0        , 1.380952381d0, 0.0d0,	    0.238095238d0,
     $         0.0d0        , 0.0d0        , 1.380952381d0, 0.0d0,
     $         0.0d0        , 70.24666667d0, 0.0d0        , 1.380952381d0/),
     $         (/4,4/)),
     $         detailled)

        end function test_compute_y_transM


        function test_compute_x_eigenvalues(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          test_validated = is_real_vector_validated(
     $         compute_x_eigenvalues_ns_vdw2d(1.380952381d0,4.089669525d0),
     $         [1.380952381d0, 1.380952381d0, -2.708717144d0, 5.470621906d0],
     $         detailled)

        end function test_compute_x_eigenvalues


        function test_compute_y_eigenvalues(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          test_validated = is_real_vector_validated(
     $         compute_y_eigenvalues_ns_vdw2d(0.619047619d0,4.089669525d0),
     $         [0.619047619d0, 0.619047619d0, -3.470621906d0, 4.708717144d0],
     $         detailled)

        end function test_compute_y_eigenvalues


        function test_compute_x_lefteigenvector(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          test_validated = is_real_matrix_validated(
     $         compute_x_lefteigenvector_ns_vdw2d(
     $             4.2d0,
     $             4.089669525d0),
     $         reshape((/
     $         0.0d0,    0.0d0,	        1.0d0,  0.0d0        ,
     $         1.0d0,    0.0d0,         0.0d0, -0.059789314d0,
     $         0.0d0,   -8.588306003d0,	0.0d0,  0.5d0        ,
     $         0.0d0,    8.588306003d0,	0.0d0, 	0.5d0/),
     $         (/4,4/)),
     $         detailled)

        end function test_compute_x_lefteigenvector


        function test_compute_x_righteigenvector(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          test_validated = is_real_matrix_validated(
     $         compute_x_righteigenvector_ns_vdw2d(
     $             4.2d0,
     $             4.089669525d0),
     $         reshape((/
     $         0.0d0, 1.0d0,  0.059789314d0, 0.059789314d0,
     $         0.0d0, 0.0d0, -0.058218699d0, 0.058218699d0,
     $         1.0d0, 0.0d0,  0.0d0, 0.0d0,
     $         0.0d0, 0.0d0,  1.0d0, 1.0d0/),
     $         (/4,4/)),
     $         detailled)

        end function test_compute_x_righteigenvector


        function test_compute_y_lefteigenvector(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          test_validated = is_real_matrix_validated(
     $         compute_y_lefteigenvector_ns_vdw2d(
     $             4.2d0,
     $             4.089669525d0),
     $         reshape((/
     $         0.0d0, 1.0d0,  0.0d0        ,  0.0d0        ,
     $         1.0d0, 0.0d0,  0.0d0        , -0.059789314d0,
     $         0.0d0, 0.0d0, -8.588306003d0,  0.5d0        ,
     $         0.0d0, 0.0d0,  8.588306003d0,  0.5d0/),
     $         (/4,4/)),
     $         detailled)

        end function test_compute_y_lefteigenvector


        function test_compute_y_righteigenvector(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          test_validated = is_real_matrix_validated(
     $         compute_y_righteigenvector_ns_vdw2d(
     $             4.2d0,
     $             4.089669525d0),
     $         reshape((/
     $         0.0d0, 1.0d0,  0.059789314d0, 0.059789314d0,
     $         1.0d0, 0.0d0,  0.0d0        , 0.0d0        ,
     $         0.0d0, 0.0d0, -0.058218699d0, 0.058218699d0,
     $         0.0d0, 0.0d0,  1.0d0        , 1.0d0        /),
     $         (/4,4/)),
     $         detailled)

        end function test_compute_y_righteigenvector

      end program test_ns_vdw2d_prim
