      program test_bf_layer_bc_procedure

        use parameters_bf_layer, only :
     $     bc_pt,
     $     bc_interior_pt,
     $     interior_pt

        use bf_layer_bc_procedure_module, only :
     $     SW_corner_type,
     $     SE_corner_type,
     $     NW_corner_type,
     $     NE_corner_type,
     $     S_edge_type,
     $     E_edge_type,
     $     W_edge_type,
     $     N_edge_type,
     $     SE_edge_type,
     $     SW_edge_type,
     $     NE_edge_type,
     $     NW_edge_type,
     $     get_bc_interior_pt_procedure

        implicit none

        integer                              :: k
        integer, dimension(:,:), allocatable :: grdpts_id
        integer                              :: test_i
        integer                              :: test_j
        integer                              :: test_procedure_type
        integer                              :: test_i_proc
        integer                              :: test_j_proc
        integer                              :: procedure_type
        integer                              :: i_proc
        integer                              :: j_proc

        logical                              :: test_validated

        
        !test get_bc_interior_pt_procedure
        print '(''test get_bc_interior_pt_procedure()'')'
        do k=1, 28

           call make_test_bf_layer_bc_procedure(
     $          k,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_procedure_type,
     $          test_i_proc,
     $          test_j_proc)

           call get_bc_interior_pt_procedure(
     $          test_i,
     $          test_j,
     $          grdpts_id,
     $          procedure_type,
     $          i_proc,
     $          j_proc)

           test_validated = test_procedure_type.eq.procedure_type
           test_validated = test_validated.and.(test_i_proc.eq.i_proc)
           test_validated = test_validated.and.(test_j_proc.eq.j_proc)

           print '(''test '',I2,'': '', L1)', k, test_validated

           if(.not.test_validated) then
              print '(''    procedure_type: '', I2, 2X, I2)',
     $             test_procedure_type, procedure_type
              print '(''    i_proc        : '', I2, 2X, I2)',
     $             test_i_proc, i_proc
              print '(''    j_proc        : '', I2, 2X, I2)',
     $             test_j_proc, j_proc
              
           end if

           deallocate(grdpts_id) 

        end do 
        print '()'

        
        contains

        subroutine make_test_bf_layer_bc_procedure(
     $       test_id,
     $       grdpts_id,
     $       test_i,
     $       test_j,
     $       test_procedure,
     $       test_i_proc,
     $       test_j_proc)
        
          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          integer                             , intent(out) :: test_procedure
          integer                             , intent(out) :: test_i_proc
          integer                             , intent(out) :: test_j_proc

          
          allocate(grdpts_id(3,3))
          test_i = 2
          test_j = 2

          select case(test_id)
            case(1)
               grdpts_id = reshape(
     $              (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = SW_edge_type
               test_i_proc    = 1
               test_j_proc    = 1

            case(2)
               grdpts_id =  reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_procedure = NE_corner_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(3)
               grdpts_id = reshape(
     $              (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_procedure = NW_corner_type
               test_i_proc    = 1
               test_j_proc    = 2

            case(4)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = SE_edge_type
               test_i_proc    = 2
               test_j_proc    = 1

            case(5)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_procedure = E_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(6)
               grdpts_id = reshape( (/
     $              bc_interior_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = NW_edge_type
               test_i_proc    = 1
               test_j_proc    = 1

            case(7)
               grdpts_id = reshape( (/
     $              interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_procedure = NE_edge_type
               test_i_proc    = 2
               test_j_proc    = 1

            case(8)
               grdpts_id = reshape( (/
     $              bc_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = W_edge_type
               test_i_proc    = 2
               test_j_proc    = 2
               
            case(9)
               grdpts_id = reshape( (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_procedure = E_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(10)
               grdpts_id = reshape( (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = SW_edge_type
               test_i_proc    = 1
               test_j_proc    = 2

            case(11)
               grdpts_id = reshape( (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = SE_edge_type
               test_i_proc    = 2
               test_j_proc    = 2               

            case(12)
               grdpts_id = reshape( (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = W_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(13)
               grdpts_id = reshape( (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_procedure = E_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(14)
               grdpts_id = reshape( (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = W_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(15)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_procedure = SE_corner_type
               test_i_proc    = 2
               test_j_proc    = 1

            case(16)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = NW_edge_type
               test_i_proc    = 1
               test_j_proc    = 2

            case(17)
               grdpts_id = reshape( (/
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = SE_edge_type
               test_i_proc    = 1
               test_j_proc    = 1

            case(18)
               grdpts_id = reshape( (/
     $              bc_interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_procedure = N_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(19)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_procedure = N_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(20)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = SW_edge_type
               test_i_proc    = 2
               test_j_proc    = 1

            case(21)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = S_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(22)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_procedure = NE_edge_type
               test_i_proc    = 1
               test_j_proc    = 2

            case(23)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = NW_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(24)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = S_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(25)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = S_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(26)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_procedure = N_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case(27)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_procedure = SW_corner_type
               test_i_proc    = 1
               test_j_proc    = 1

            case(28)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_procedure = NE_edge_type
               test_i_proc    = 2
               test_j_proc    = 2

            case default
               print '(''test_bf_layer_bc_procedure'')'
               print '(''make_test_bf_layer_bc_procedure'')'
               print '(''test case not implemented: '', I2)', test_id
               stop ''

          end select

        end subroutine make_test_bf_layer_bc_procedure

      end program test_bf_layer_bc_procedure
