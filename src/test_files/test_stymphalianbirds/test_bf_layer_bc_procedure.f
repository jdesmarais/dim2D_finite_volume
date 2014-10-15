      program test_bf_layer_bc_procedure

        use bf_layer_bc_procedure_module, only :
     $     SW_corner_type,
     $     SE_corner_type,
     $     NW_corner_type,
     $     NE_corner_type,
     $     S_edge_type,
     $     E_edge_type,
     $     W_edge_type,
     $     N_edge_type,
     $     get_bc_interior_pt_procedure

        use parameters_bf_layer, only :
     $     interior_pt,
     $     bc_interior_pt,
     $     bc_pt

        implicit none

        integer                 :: test_id
        integer, dimension(3,3) :: grdpts_id
        integer                 :: procedure_type
        integer                 :: test_procedure_type
        logical                 :: test_validated


        do test_id=1,16

           call make_test_bc_interior_pt_procedure(
     $          test_id,
     $          grdpts_id,
     $          test_procedure_type)

           procedure_type = get_bc_interior_pt_procedure(
     $          2,
     $          2,
     $          grdpts_id)

           test_validated = procedure_type.eq.test_procedure_type

           print '(''test '',I2,'': '',L1)', test_id, test_validated

        end do
        

        
        contains

        subroutine make_test_bc_interior_pt_procedure(
     $       test_id,
     $       grdpts_id,
     $       test_procedure_type)

          implicit none

          integer                , intent(in)  :: test_id
          integer, dimension(3,3), intent(out) :: grdpts_id
          integer                , intent(out) :: test_procedure_type


          select case(test_id)

            !  -------
            ! | 1 1 1 |
            ! | 0 0 0 |
            ! |       |
            !  -------
            case(1)
               grdpts_id(1,1) = interior_pt
               grdpts_id(2,1) = interior_pt
               grdpts_id(3,1) = interior_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = bc_pt
               grdpts_id(2,3) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=N_edge_type
               
            !  -------
            ! |       |
            ! | 0 0 0 |
            ! | 1 1 1 |
            !  -------
            case(2)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = interior_pt
               grdpts_id(2,3) = interior_pt
               grdpts_id(3,3) = interior_pt

               test_procedure_type=S_edge_type

            !  -------
            ! |   0 1 |
            ! |   0 1 |
            ! |   0 1 |
            !  -------
            case(3)
               grdpts_id(1,1) = interior_pt
               grdpts_id(1,2) = interior_pt
               grdpts_id(1,3) = interior_pt
                                                     
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(2,3) = bc_interior_pt
                                                     
               grdpts_id(3,1) = bc_pt
               grdpts_id(3,2) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=E_edge_type

            !  -------
            ! | 1 0   |
            ! | 1 0   |
            ! | 1 0   |
            !  -------
            case(4)
               grdpts_id(1,1) = bc_pt
               grdpts_id(1,2) = bc_pt
               grdpts_id(1,3) = bc_pt
                                                     
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(2,3) = bc_interior_pt
                                                     
               grdpts_id(3,1) = interior_pt
               grdpts_id(3,2) = interior_pt
               grdpts_id(3,3) = interior_pt

               test_procedure_type=W_edge_type


            !  -------
            ! | 1 1 1 |
            ! | 0 0 0 |
            ! | 0     |
            !  -------
            case(5)
               grdpts_id(1,1) = bc_interior_pt
               grdpts_id(2,1) = interior_pt
               grdpts_id(3,1) = interior_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = bc_pt
               grdpts_id(2,3) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=N_edge_type
               
            !  -------
            ! |     0 |
            ! | 0 0 0 |
            ! | 1 1 1 |
            !  -------
            case(6)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = interior_pt
               grdpts_id(2,3) = interior_pt
               grdpts_id(3,3) = bc_interior_pt

               test_procedure_type=S_edge_type

            !  -------
            ! |   0 1 |
            ! |   0 1 |
            ! | 0 0 1 |
            !  -------
            case(7)
               grdpts_id(1,1) = bc_interior_pt
               grdpts_id(1,2) = interior_pt
               grdpts_id(1,3) = interior_pt
                                                     
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(2,3) = bc_interior_pt
                                                     
               grdpts_id(3,1) = bc_pt
               grdpts_id(3,2) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=E_edge_type

            !  -------
            ! | 1 0 0 |
            ! | 1 0   |
            ! | 1 0   |
            !  -------
            case(8)
               grdpts_id(1,1) = bc_pt
               grdpts_id(1,2) = bc_pt
               grdpts_id(1,3) = bc_pt
                                                     
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(2,3) = bc_interior_pt
                                                     
               grdpts_id(3,1) = interior_pt
               grdpts_id(3,2) = interior_pt
               grdpts_id(3,3) = bc_interior_pt

               test_procedure_type=W_edge_type

            !  -------
            ! | 1 1 1 |
            ! | 1 0 0 |
            ! | 1 0   |
            !  -------
            case(9)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(3,1) = interior_pt

               grdpts_id(1,2) = bc_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = bc_pt
               grdpts_id(2,3) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=NW_corner_type

            !  -------
            ! | 1 1 1 |
            ! | 0 0 1 |
            ! |   0 1 |
            !  -------
            case(10)
               grdpts_id(1,1) = interior_pt
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_pt

               grdpts_id(1,3) = bc_pt
               grdpts_id(2,3) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=NE_corner_type

            !  -------
            ! | 1 0   |
            ! | 1 0 0 |
            ! | 1 1 1 |
            !  -------
            case(11)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = bc_pt
               grdpts_id(2,3) = bc_interior_pt
               grdpts_id(3,3) = interior_pt

               test_procedure_type=SW_corner_type

            !  -------
            ! |   0 1 |
            ! | 0 0 1 |
            ! | 1 1 1 |
            !  -------
            case(12)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_pt

               grdpts_id(1,3) = interior_pt
               grdpts_id(2,3) = bc_interior_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=SE_corner_type

            !  -------
            ! | 1 1 1 |
            ! | 1 0 0 |
            ! | 0 0   |
            !  -------
            case(13)
               grdpts_id(1,1) = bc_interior_pt
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(3,1) = interior_pt

               grdpts_id(1,2) = bc_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = bc_pt
               grdpts_id(2,3) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=NW_corner_type

            !  -------
            ! | 0 1 1 |
            ! | 0 0 1 |
            ! |   0 1 |
            !  -------
            case(14)
               grdpts_id(1,1) = interior_pt
               grdpts_id(2,1) = bc_interior_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_pt

               grdpts_id(1,3) = bc_interior_pt
               grdpts_id(2,3) = bc_pt
               grdpts_id(3,3) = bc_pt

               test_procedure_type=NE_corner_type

            !  -------
            ! | 0 0   |
            ! | 1 0 0 |
            ! | 1 1 1 |
            !  -------
            case(15)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_interior_pt

               grdpts_id(1,3) = bc_interior_pt
               grdpts_id(2,3) = bc_interior_pt
               grdpts_id(3,3) = interior_pt

               test_procedure_type=SW_corner_type

            !  -------
            ! |   0 0 |
            ! | 0 0 1 |
            ! | 1 1 1 |
            !  -------
            case(16)
               grdpts_id(1,1) = bc_pt
               grdpts_id(2,1) = bc_pt
               grdpts_id(3,1) = bc_pt

               grdpts_id(1,2) = bc_interior_pt
               grdpts_id(2,2) = bc_interior_pt
               grdpts_id(3,2) = bc_pt

               grdpts_id(1,3) = interior_pt
               grdpts_id(2,3) = bc_interior_pt
               grdpts_id(3,3) = bc_interior_pt

               test_procedure_type=SE_corner_type

            case default
               print '(''test_bf_layer_bc_procedure'')'
               print '(''get_test_bc_interior_pt_procedure'')'
               print '(''test case not implemented'')'
               stop 'choose [1,16]'

          end select

        end subroutine make_test_bc_interior_pt_procedure

      end program test_bf_layer_bc_procedure
