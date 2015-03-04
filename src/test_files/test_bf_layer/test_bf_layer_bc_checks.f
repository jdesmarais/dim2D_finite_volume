      program test_bf_layer_bc_checks

        use bf_layer_bc_checks_module, only :
     $     compute_edge_N,
     $     compute_edge_S,
     $     compute_edge_E,
     $     compute_edge_W

        use parameters_constant, only :
     $       bc_nodes_choice,
     $       bc_timedev_choice

        use parameters_input, only :
     $       nx,ny,ne,
     $       x_min, x_max,
     $       y_min, y_max,
     $       bc_N_type_choice,
     $       bc_S_type_choice,
     $       bc_E_type_choice,
     $       bc_W_type_choice

        use parameters_kind, only :
     $       rkind

        
        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_compute_edge_N(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_N: '',L1)', test_loc
        print '()'

        test_loc = test_compute_edge_S(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_S: '',L1)', test_loc
        print '()'

        test_loc = test_compute_edge_E(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_E: '',L1)', test_loc
        print '()'

        test_loc = test_compute_edge_W(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_compute_W: '',L1)', test_loc
        print '()'


        contains


        function test_compute_edge_N(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(6) :: y_test
          integer    , dimension(6) :: bc_type_test
          logical    , dimension(6) :: compute_edge_test
          logical                   :: compute_edge
          integer                   :: k

          test_validated = .true.

          y_test = [y_max-1, y_max, y_max+1,
     $              y_max-1, y_max, y_max+1]

          bc_type_test = [bc_timedev_choice, bc_timedev_choice, bc_timedev_choice,
     $                    bc_nodes_choice  , bc_nodes_choice  , bc_nodes_choice]
          
          compute_edge_test = [.true.,.true.,.true.,.true.,.false.,.false.]


          do k=1,size(y_test,1)

             !output
             compute_edge = compute_edge_N(y_test(k),bc_type_test(k))

             !validation
             test_loc = compute_edge.eqv.compute_edge_test(k)
             test_validated = test_validated.and.test_loc
             
             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed '')', k
             end if

          end do

        end function test_compute_edge_N


        function test_compute_edge_S(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(6) :: y_test
          integer    , dimension(6) :: bc_type_test
          logical    , dimension(6) :: compute_edge_test
          logical                   :: compute_edge
          integer                   :: k

          test_validated = .true.

          y_test = [y_min-1, y_min, y_min+1,
     $              y_min-1, y_min, y_min+1]

          bc_type_test = [bc_timedev_choice, bc_timedev_choice, bc_timedev_choice,
     $                    bc_nodes_choice  , bc_nodes_choice  , bc_nodes_choice]
          
          compute_edge_test = [.true.,.true.,.true.,.false.,.false.,.true.]


          do k=1,size(y_test,1)

             !output
             compute_edge = compute_edge_S(y_test(k),bc_type_test(k))

             !validation
             test_loc = compute_edge.eqv.compute_edge_test(k)
             test_validated = test_validated.and.test_loc
             
             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed '')', k
             end if

          end do

        end function test_compute_edge_S


        function test_compute_edge_E(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(6) :: x_test
          integer    , dimension(6) :: bc_type_test
          logical    , dimension(6) :: compute_edge_test
          logical                   :: compute_edge
          integer                   :: k

          test_validated = .true.

          x_test = [x_max-1, x_max, x_max+1,
     $              x_max-1, x_max, x_max+1]

          bc_type_test = [bc_timedev_choice, bc_timedev_choice, bc_timedev_choice,
     $                    bc_nodes_choice  , bc_nodes_choice  , bc_nodes_choice]
          
          compute_edge_test = [.true.,.true.,.true.,.true.,.false.,.false.]


          do k=1,size(x_test,1)

             !output
             compute_edge = compute_edge_E(x_test(k),bc_type_test(k))

             !validation
             test_loc = compute_edge.eqv.compute_edge_test(k)
             test_validated = test_validated.and.test_loc
             
             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed '')', k
             end if

          end do

        end function test_compute_edge_E


        function test_compute_edge_W(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(6) :: x_test
          integer    , dimension(6) :: bc_type_test
          logical    , dimension(6) :: compute_edge_test
          logical                   :: compute_edge
          integer                   :: k

          test_validated = .true.

          x_test = [x_min-1, x_min, x_min+1,
     $              x_min-1, x_min, x_min+1]

          bc_type_test = [bc_timedev_choice, bc_timedev_choice, bc_timedev_choice,
     $                    bc_nodes_choice  , bc_nodes_choice  , bc_nodes_choice]
          
          compute_edge_test = [.true.,.true.,.true.,.false.,.false.,.true.]


          do k=1,size(x_test,1)

             !output
             compute_edge = compute_edge_W(x_test(k),bc_type_test(k))

             !validation
             test_loc = compute_edge.eqv.compute_edge_test(k)
             test_validated = test_validated.and.test_loc
             
             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed '')', k
             end if

          end do

        end function test_compute_edge_W


        subroutine check_inputs()

          implicit none

          if(bc_N_type_choice.ne.bc_timedev_choice) then
             print '(''bc_N_type_choice.eq.bc_timedev_choice: '',L1)',
     $            bc_N_type_choice.eq.bc_timedev_choice
          end if
          
          if(bc_S_type_choice.ne.bc_timedev_choice) then
             print '(''bc_S_type_choice.eq.bc_timedev_choice: '',L1)',
     $            bc_S_type_choice.eq.bc_timedev_choice
          end if

          if(bc_E_type_choice.ne.bc_timedev_choice) then
             print '(''bc_E_type_choice.eq.bc_timedev_choice: '',L1)',
     $            bc_E_type_choice.eq.bc_timedev_choice
          end if

          if(bc_W_type_choice.ne.bc_timedev_choice) then
             print '(''bc_W_type_choice.eq.bc_timedev_choice: '',L1)',
     $            bc_W_type_choice.eq.bc_timedev_choice
          end if

        end subroutine check_inputs

      end program test_bf_layer_bc_checks
