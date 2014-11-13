      program test_bf_layer_newgrdpt_procedure

        use bf_layer_bc_procedure_module, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type

        use bf_layer_newgrdpt_procedure_module, only :
     $       no_gradient_type,
     $       gradient_I_type,
     $       gradient_L0_type,
     $       gradient_R0_type,
     $       get_newgrdpt_procedure

        use parameters_bf_layer, only :
     $       no_pt,
     $       bc_pt,
     $       bc_interior_pt

        use parameters_kind, only :
     $       ikind

        implicit none

        
        integer :: k
        integer :: nb_tests

        integer, dimension(:,:), allocatable :: grdpts_id_tested
        integer(ikind)                       :: i_tested
        integer(ikind)                       :: j_tested
        integer                              :: procedure_type_data
        integer                              :: gradient_type_data

        integer                              :: procedure_type
        integer                              :: gradient_type

        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        nb_tests = 22
        test_validated = .true.
        detailled = .true.


        do k=1, nb_tests

           !get the data to conduct the test
           call get_test_data(
     $          k,
     $          grdpts_id_tested,
     $          i_tested,
     $          j_tested,
     $          procedure_type_data,
     $          gradient_type_data)

           !carry the test
           call get_newgrdpt_procedure(
     $          i_tested,
     $          j_tested,
     $          grdpts_id_tested,
     $          procedure_type,
     $          gradient_type)
           
           !compare the results
           test_loc = (procedure_type.eq.procedure_type_data).and.
     $          (gradient_type.eq.gradient_type_data)
           test_validated = test_validated.and.test_loc

           !print message if problem
           if(detailled) then
              if(.not.test_loc) then

                 print '(''***test '',I2,'' failed***: '')', k
                 print '(''   - procedure_type: '',I2,''  -> '', I2)', procedure_type, procedure_type_data
                 print '(''   - gradient_type : '',I2,''  -> '', I2)', gradient_type , gradient_type_data
                            
              else

                 print '(''test '',I2,'' : '',L1)', k, test_loc

              end if
           end if

        end do

        print '()'
        print '(''test_validated: '',L1)', test_validated


        contains


        subroutine get_test_data(
     $       k,
     $       grdpts_id_tested,
     $       i_tested,j_tested,
     $       procedure_type_data,
     $       gradient_type_data)


          implicit none

          integer                             , intent(in)  :: k
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id_tested
          integer                             , intent(out) :: i_tested
          integer                             , intent(out) :: j_tested
          integer                             , intent(out) :: procedure_type_data
          integer                             , intent(out) :: gradient_type_data


          if(allocated(grdpts_id_tested)) then
             deallocate(grdpts_id_tested)
          end if

          
          select case(k)

            !  -------
            ! | 0 0*0 |
            ! | 3 3 3 |
            !  -------
            case(1)

               allocate(grdpts_id_tested(3,2))
               grdpts_id_tested = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,2/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = N_edge_type
               gradient_type_data  = gradient_I_type

            !  -----
            ! | 3 0*|
            ! | 3 3 |
            !  ----- 
            case(2)

               allocate(grdpts_id_tested(2,2))
               grdpts_id_tested = reshape( (/
     $              bc_pt,bc_pt,
     $              bc_pt,no_pt
     $              /),
     $              (/2,2/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = NE_edge_type
               gradient_type_data  = no_gradient_type

            !  -------
            ! | 0 0*0 |
            ! | 3 3 0 |
            !  -------
            case(3)

               allocate(grdpts_id_tested(3,2))
               grdpts_id_tested = reshape( (/
     $              bc_pt,bc_pt,no_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,2/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = N_edge_type
               gradient_type_data  = gradient_R0_type

            !  -----
            ! | 3 0 |
            ! | 3 0*|
            ! | 3 0 |
            !  -----
            case(4)

               allocate(grdpts_id_tested(2,3))
               grdpts_id_tested = reshape( (/
     $              bc_pt,no_pt,
     $              bc_pt,no_pt,
     $              bc_pt,no_pt
     $              /),
     $              (/2,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = E_edge_type
               gradient_type_data  = gradient_I_type

            !  -----
            ! | 3 0*|
            ! | 3 0 |
            !  -----
            case(5)

               allocate(grdpts_id_tested(2,2))
               grdpts_id_tested = reshape( (/
     $              bc_pt,no_pt,
     $              bc_pt,no_pt
     $              /),
     $              (/2,2/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = E_edge_type
               gradient_type_data  = gradient_R0_type

            !  -------
            ! ! 0 0 0*|
            ! | 3 3 0 |
            ! | 2 3 0 |
            !  -------
            case(6)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              bc_interior_pt,bc_pt,no_pt,
     $              bc_pt,bc_pt,no_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 3
               j_tested = 3

               procedure_type_data = NE_corner_type
               gradient_type_data  = no_gradient_type

            !  -------
            ! ! 0 0*3 |
            ! | 0 3 3 |
            !  -------
            case(7)

               allocate(grdpts_id_tested(3,2))
               grdpts_id_tested = reshape( (/
     $              no_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,bc_pt
     $              /),
     $              (/3,2/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = NW_edge_type
               gradient_type_data  = no_gradient_type

            !  -----
            ! ! 0*0 |
            ! | 3 3 |
            !  -----
            case(8)

               allocate(grdpts_id_tested(2,2))
               grdpts_id_tested = reshape( (/
     $              bc_pt,bc_pt,
     $              no_pt,no_pt
     $              /),
     $              (/2,2/))

               i_tested = 1
               j_tested = 2

               procedure_type_data = N_edge_type
               gradient_type_data  = gradient_L0_type


            !  -------
            ! | 0 3 3 |
            ! | 0 0*3 |
            ! | 0 0 3 |
            !  -------
            case(9)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,bc_pt,
     $              no_pt,no_pt,bc_pt,
     $              no_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = SW_edge_type
               gradient_type_data  = no_gradient_type

            !  -------
            ! | 0 0 3 |
            ! | 0 0*3 |
            ! | 0 0 3 |
            !  -------
            case(10)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,bc_pt,
     $              no_pt,no_pt,bc_pt,
     $              no_pt,no_pt,bc_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = W_edge_type
               gradient_type_data  = gradient_I_type

            !  -----
            ! | 0*3 |
            ! | 0 3 |
            !  -----
            case(11)

               allocate(grdpts_id_tested(2,2))
               grdpts_id_tested = reshape( (/
     $              no_pt,bc_pt,
     $              no_pt,bc_pt
     $              /),
     $              (/2,2/))

               i_tested = 1
               j_tested = 2

               procedure_type_data = W_edge_type
               gradient_type_data  = gradient_R0_type


            !  -----
            ! | 0*0 |
            ! | 0 3 |
            !  -----
            case(12)

               allocate(grdpts_id_tested(2,2))
               grdpts_id_tested = reshape( (/
     $              no_pt,bc_pt,
     $              no_pt,no_pt
     $              /),
     $              (/2,2/))

               i_tested = 1
               j_tested = 2

               procedure_type_data = NW_corner_type
               gradient_type_data  = no_gradient_type

            !  -----
            ! | 3 3 |
            ! | 3 0*|
            !  -----
            case(13)

               allocate(grdpts_id_tested(2,2))
               grdpts_id_tested = reshape( (/
     $              bc_pt,no_pt,
     $              bc_pt,bc_pt
     $              /),
     $              (/2,2/))

               i_tested = 2
               j_tested = 1

               procedure_type_data = SE_edge_type
               gradient_type_data  = no_gradient_type

            !  -----
            ! | 3 0 |
            ! | 3 0*|
            !  -----
            case(14)

               allocate(grdpts_id_tested(2,2))
               grdpts_id_tested = reshape( (/
     $              bc_pt,no_pt,
     $              bc_pt,no_pt
     $              /),
     $              (/2,2/))

               i_tested = 2
               j_tested = 1

               procedure_type_data = E_edge_type
               gradient_type_data  = gradient_L0_type

            !  -----
            ! | 3 3 |
            ! | 0*3 |
            !  -----
            case(15)

               allocate(grdpts_id_tested(2,2))
               grdpts_id_tested = reshape( (/
     $              no_pt,bc_pt,
     $              bc_pt,bc_pt
     $              /),
     $              (/2,2/))

               i_tested = 1
               j_tested = 1

               procedure_type_data = SW_edge_type
               gradient_type_data  = no_gradient_type

            !  -----
            ! | 0 3 |
            ! | 0*3 |
            !  -----
            case(16)

               allocate(grdpts_id_tested(2,2))
               grdpts_id_tested = reshape( (/
     $              no_pt,bc_pt,
     $              no_pt,bc_pt
     $              /),
     $              (/2,2/))

               i_tested = 1
               j_tested = 1

               procedure_type_data = W_edge_type
               gradient_type_data  = gradient_L0_type

            !  -------
            ! | 3 3 3 |
            ! | 0 0*0 |
            !  -------
            case(17)

               allocate(grdpts_id_tested(3,2))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,2/))

               i_tested = 2
               j_tested = 1

               procedure_type_data = S_edge_type
               gradient_type_data  = gradient_I_type

            !  -------
            ! | 3 3 0 |
            ! | 0 0*0 |
            !  -------
            case(18)

               allocate(grdpts_id_tested(3,2))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,no_pt
     $              /),
     $              (/3,2/))

               i_tested = 2
               j_tested = 1

               procedure_type_data = S_edge_type
               gradient_type_data  = gradient_R0_type

            !  -----
            ! | 3 0 |
            ! | 0 0*|
            !  -----
            case(19)

               allocate(grdpts_id_tested(2,2))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,
     $              bc_pt,no_pt
     $              /),
     $              (/2,2/))

               i_tested = 2
               j_tested = 1

               procedure_type_data = SE_corner_type
               gradient_type_data  = no_gradient_type

            !  -----
            ! | 0 3 |
            ! | 0*0 |
            !  -----
            case(20)

               allocate(grdpts_id_tested(2,2))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,
     $              no_pt,bc_pt
     $              /),
     $              (/2,2/))

               i_tested = 1
               j_tested = 1

               procedure_type_data = SW_corner_type
               gradient_type_data  = no_gradient_type

            !  -----
            ! | 3 3 |
            ! | 0*0 |
            !  -----
            case(21)

               allocate(grdpts_id_tested(2,2))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,
     $              bc_pt,bc_pt
     $              /),
     $              (/2,2/))

               i_tested = 1
               j_tested = 1

               procedure_type_data = S_edge_type
               gradient_type_data  = gradient_L0_type

            !  -------
            ! | 0 0 3 |
            ! | 0 0*0 |
            ! | 0 0 0 |
            !  -------
            case(22)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,
     $              no_pt,no_pt,bc_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = SW_corner_type
               gradient_type_data  = no_gradient_type

            case default
               print '(''test_bf_layer_newgrdpt_procedure'')'
               print '(''get_test_data'')'
               stop 'case not recognized'

          end select

        end subroutine get_test_data

      end program test_bf_layer_newgrdpt_procedure
