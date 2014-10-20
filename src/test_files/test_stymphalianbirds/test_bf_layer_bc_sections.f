      program test_bf_layer_bc_sections

        use bf_layer_bc_procedure_module, only : 
     $     N_edge_type,
     $     S_edge_type,
     $     E_edge_type,
     $     W_edge_type,
     $     NE_corner_type,
     $     NW_corner_type,
     $     SE_corner_type,
     $     SW_corner_type,
     $     NE_edge_type,
     $     NW_edge_type,
     $     SE_edge_type,
     $     SW_edge_type

        use bf_layer_bc_sections_class, only :
     $     bf_layer_bc_sections

        use parameters_bf_layer, only :
     $     interior_pt,
     $     bc_interior_pt,
     $     bc_pt,
     $     no_pt


        implicit none

        type(bf_layer_bc_sections)           :: bc_sections

        integer, dimension(5)                :: new_bc_section
        integer                              :: k,l
        logical                              :: compatible
        logical                              :: remove_ele

        integer                              :: test_i
        integer                              :: test_j
        integer, dimension(:,:), allocatable :: grdpts_id
        integer, dimension(5)                :: test_bc_section
        logical                              :: test_validated
        logical                              :: test_compatible
        logical                              :: test_remove_ele
        
        
        
        !test ini()
        print '(''test ini()'')'
        call bc_sections%ini()
        print '(''nb_ele_temp : '',I2)', bc_sections%get_nb_ele_temp()
        print '(''nb_ele_final: '',I2)', bc_sections%get_nb_ele_final()
        print '()'


        !test add_to_temporary_bc_sections() and
        !add_to_final_bc_sections()
        print '(''test add_to_temporary_bc_sections()'')'
        
        do k=1,7
           new_bc_section = get_test_bc_section(k)
           call bc_sections%add_to_temporary_bc_sections(new_bc_section)
           call bc_sections%add_to_final_bc_sections(new_bc_section)
        end do
        call bc_sections%print_bc_sections()


        !test deallocate_tables
        call bc_sections%deallocate_tables()

        !test add_to_bc_sections()
        print '(''test add_to_bc_sections()'')'

        call bc_sections%ini()
        do k=1,9
           new_bc_section = get_test_bc_section(k)
           call bc_sections%add_to_bc_sections(new_bc_section)
        end do
        call bc_sections%print_bc_sections()


        !test get_bc_section()
        print '(''test get_bc_section()'')'
        do k=1,26

           call make_test_get_bc_section(
     $          k,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_bc_section)

           new_bc_section = bc_sections%get_bc_section(
     $          test_i,
     $          test_j,
     $          grdpts_id)


           test_validated = test_bc_section(1).eq.new_bc_section(1)

           select case(test_bc_section(1))

             !if this is an edge procedure, only the 2:4
             !elements should be tested
             case(N_edge_type,S_edge_type,E_edge_type,W_edge_type)

                test_validated = test_validated.and.(
     $               test_bc_section(2).eq.new_bc_section(2))
                test_validated = test_validated.and.(
     $               test_bc_section(3).eq.new_bc_section(3))
                test_validated = test_validated.and.(
     $               test_bc_section(4).eq.new_bc_section(4))

                print '(''test '',I2,'':'',L1)', k, test_validated

                if(.not.test_validated) then
                   print '(''procedure_type: '',I2, 2X,I2)', test_bc_section(1), new_bc_section(1)
                   print '(''edge_min      : '',I2, 2X,I2)', test_bc_section(2), new_bc_section(2)
                   print '(''edge_max      : '',I2, 2X,I2)', test_bc_section(3), new_bc_section(3)
                   print '(''coord         : '',I2, 2X,I2)', test_bc_section(4), new_bc_section(4)
                end if

             !if this is a corner procedure, only the 2:3
             !elements should be tested
             case(NE_corner_type, NW_corner_type, SE_corner_type, SW_corner_type)

                test_validated = test_validated.and.(
     $               test_bc_section(2).eq.new_bc_section(2))
                test_validated = test_validated.and.(
     $               test_bc_section(3).eq.new_bc_section(3))

                print '(''test '',I2,'':'',L1)', k, test_validated

                if(.not.test_validated) then
                   print '(''procedure_type: '',I2, 2X,I2)', test_bc_section(1), new_bc_section(1)
                   print '(''i_min         : '',I2, 2X,I2)', test_bc_section(2), new_bc_section(2)
                   print '(''j_min         : '',I2, 2X,I2)', test_bc_section(3), new_bc_section(3)
                end if

             !if this is a special edge procedure, only the 2:3 and 5
             !elements should be tested
             case(NE_edge_type,NW_edge_type,SE_edge_type,SW_edge_type)

                test_validated = test_validated.and.(
     $               test_bc_section(2).eq.new_bc_section(2))
                test_validated = test_validated.and.(
     $               test_bc_section(3).eq.new_bc_section(3))
                test_validated = test_validated.and.(
     $               test_bc_section(5).eq.new_bc_section(5))

                print '(''test '',I2,'':'',L1)', k, test_validated

                if(.not.test_validated) then
                   print '(''procedure_type: '',I2, 2X,I2)', test_bc_section(1), new_bc_section(1)
                   print '(''i_min         : '',I2, 2X,I2)', test_bc_section(2), new_bc_section(2)
                   print '(''j_min         : '',I2, 2X,I2)', test_bc_section(3), new_bc_section(3)
                   print '(''match_nb      : '',I2, 2X,I2)', test_bc_section(5), new_bc_section(5)
                end if

             case default
                
                print '(''test_bf_layer_bc_sections'')'
                print '(''test_get_bc-sections()'')'
                print '(''case '', I2, ''not recognized'')', k
                stop ''

           end select

           deallocate(grdpts_id)

        end do
        print '()'


        !test analyse_grdpt_with_bc_section()
        print '(''test analyse_grdpt_with_bc_section()'')'
        do k=1,12

           call make_test_analyse_grdpt_with_bc_section(
     $          k,
     $          grdpts_id,
     $          test_i,
     $          test_j,
     $          test_bc_section,
     $          test_compatible,
     $          test_remove_ele)

           compatible = bc_sections%analyse_grdpt_with_bc_section(
     $          test_i,test_j,grdpts_id,test_bc_section,remove_ele)

           test_validated = compatible.eqv.test_compatible
           test_validated = test_validated.and.(remove_ele.eqv.test_remove_ele)

           print '(''test '',I2,'':'',L1)', k, test_validated

           if(.not.test_validated) then
              print '(''  compatible    : '',L1, 2X,L1)', test_compatible, compatible
              print '(''  remove_ele    : '',L1, 2X,L1)', test_remove_ele, remove_ele
           end if

           deallocate(grdpts_id)

        end do
        print '()'


        contains

        function get_test_bc_section(k)
     $       result(bc_section)

          implicit none

          integer, intent(in)   :: k
          integer, dimension(5) :: bc_section

          integer, dimension(12) :: proc_type

          proc_type = [
     $         N_edge_type,
     $         S_edge_type,
     $         E_edge_type,
     $         W_edge_type,
     $         NE_edge_type,
     $         NW_edge_type,
     $         SE_edge_type,
     $         SW_edge_type,
     $         NE_corner_type,
     $         NW_corner_type,
     $         SE_corner_type,
     $         SW_corner_type]

          bc_section = [
     $         proc_type(k),
     $         k,
     $         k,
     $         k+1,
     $         mod(k,3)]

        end function get_test_bc_section


      subroutine make_test_get_bc_section(
     $       test_id,
     $       grdpts_id,
     $       test_i,
     $       test_j,
     $       test_bc_section)
        
          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          integer, dimension(5)               , intent(out) :: test_bc_section

          
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
               test_bc_section(1) = SW_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(5) = 1

            case(2)
               grdpts_id =  reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NE_corner_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2

            case(3)
               grdpts_id = reshape(
     $              (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NW_corner_type
               test_bc_section(2) = 1
               test_bc_section(3) = 2

            case(4)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SE_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 1
               test_bc_section(5) = 1

            case(5)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = E_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(6)
               grdpts_id = reshape( (/
     $              bc_interior_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NW_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(5) = 1

            case(7)
               grdpts_id = reshape( (/
     $              interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NE_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 1
               test_bc_section(5) = 1

            case(8)
               grdpts_id = reshape( (/
     $              bc_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = W_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2
               
            case(9)
               grdpts_id = reshape( (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = E_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(10)
               grdpts_id = reshape( (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SW_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 2
               test_bc_section(5) = 1

            case(11)
               grdpts_id = reshape( (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SE_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(5) = 1

            case(12)
               grdpts_id = reshape( (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = W_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(13)
               grdpts_id = reshape( (/
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = E_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(14)
               grdpts_id = reshape( (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = W_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(15)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SE_corner_type
               test_bc_section(2) = 2
               test_bc_section(3) = 1

            case(16)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NW_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 2
               test_bc_section(5) = 1

            case(17)
               grdpts_id = reshape( (/
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SE_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(5) = 1

            case(18)
               grdpts_id = reshape( (/
     $              bc_interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = N_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(19)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = N_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(20)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SW_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 1
               test_bc_section(5) = 1

            case(21)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = S_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(22)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_interior_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NE_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 2
               test_bc_section(5) = 1

            case(23)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NW_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(5) = 1

            case(24)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = S_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(25)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = S_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(26)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = N_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 2

            case(27)
               grdpts_id = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              bc_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = SW_corner_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1

            case(28)
               grdpts_id = reshape( (/
     $              interior_pt,interior_pt,interior_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))
               test_bc_section(1) = NE_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(5) = 1

            case default
               print '(''test_bf_layer_bc_procedure'')'
               print '(''make_test_bf_layer_bc_procedure'')'
               print '(''test case not implemented: '', I2)', test_id
               stop ''

          end select

        end subroutine make_test_get_bc_section


        subroutine make_test_analyse_grdpt_with_bc_section(
     $       test_id,
     $       grdpts_id,
     $       test_i,
     $       test_j,
     $       test_bc_section,
     $       test_compatible,
     $       test_remove_ele)
        
          implicit none

          integer                             , intent(in)  :: test_id
          integer, dimension(:,:), allocatable, intent(out) :: grdpts_id
          integer                             , intent(out) :: test_i
          integer                             , intent(out) :: test_j
          integer, dimension(5)               , intent(out) :: test_bc_section
          logical                             , intent(out) :: test_compatible
          logical                             , intent(out) :: test_remove_ele

          
          allocate(grdpts_id(3,3))

          select case(test_id)

            !test W_edge: compatible

            !  -------    
            ! | 2 1 0 |
            ! | 2 1*0 |
            ! | 2 1 0 |
            !  -------
            case(1)
               grdpts_id = reshape(
     $              (/
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt,
     $              bc_pt,bc_interior_pt,interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 2

               test_bc_section(1) = W_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 2
               test_bc_section(5) = 0
               test_compatible    = .true.
               test_remove_ele    = .false.

            !test E_edge: compatible

            !  -------    
            ! | 1 2   |
            ! | 1*2   |
            ! | 1 2   |
            !  -------
            case(2)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_pt,no_pt,
     $              bc_interior_pt,bc_pt,no_pt,
     $              bc_interior_pt,bc_pt,no_pt
     $              /),
     $              (/3,3/))

               test_i             = 1
               test_j             = 2

               test_bc_section(1) = E_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 1
               test_bc_section(5) = 0
               test_compatible    = .true.
               test_remove_ele    = .false.

            !test N_edge: compatible

            !  -------    
            ! |       |
            ! | 2 2 2 |
            ! | 1 1*1 |
            !  -------
            case(3)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 1

               test_bc_section(1) = N_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 1
               test_bc_section(5) = 0
               test_compatible    = .true.
               test_remove_ele    = .false.

            !test S_edge: compatible

            !  -------    
            ! | 1 1*1 |
            ! | 2 2 2 |
            ! |       |
            !  -------
            case(4)
               grdpts_id = reshape(
     $              (/
     $              no_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 3

               test_bc_section(1) = S_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 3
               test_bc_section(5) = 0
               test_compatible    = .true.
               test_remove_ele    = .false.

            !test W_edge: not compatible

            !  -------
            ! | 1 1 1 |
            ! | 2 2 1*|
            ! |   2 1 |
            !  -------
            case(5)
               grdpts_id = reshape(
     $              (/
     $              no_pt,bc_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 3
               test_j             = 2

               test_bc_section(1) = W_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 3
               test_bc_section(5) = 0
               test_compatible    = .false.
               test_remove_ele    = .false.

            !test E_edge: not compatible

            !  -------    
            ! | 1 1 1 |
            ! | 1*2 2 |
            ! | 1 2   |
            !  -------
            case(6)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_pt,no_pt,
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 1
               test_j             = 2

               test_bc_section(1) = E_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 1
               test_bc_section(5) = 0
               test_compatible    = .false.
               test_remove_ele    = .false.

            !test N_edge: not compatible

            !  -------    
            ! |   2 1 |
            ! | 2 2 1 |
            ! | 1 1*1 |
            !  -------
            case(7)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              no_pt,bc_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 1

               test_bc_section(1) = N_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 1
               test_bc_section(5) = 0
               test_compatible    = .false.
               test_remove_ele    = .false.

            !test S_edge: not compatible

            !  -------    
            ! | 1 1*1 |
            ! | 2 2 1 |
            ! |   2 1 |
            !  -------
            case(8)
               grdpts_id = reshape(
     $              (/
     $              no_pt,bc_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 3

               test_bc_section(1) = S_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 1
               test_bc_section(4) = 3
               test_bc_section(5) = 0
               test_compatible    = .false.
               test_remove_ele    = .false.

            !test SW_edge: compatible

            !  -------    
            ! | 1 1 1 |
            ! | 2 2 1*|
            ! |   2 1 |
            !  -------
            case(9)
               grdpts_id = reshape(
     $              (/
     $              no_pt,bc_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 3
               test_j             = 2

               test_bc_section(1) = SW_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 0
               test_bc_section(5) = 1
               test_compatible    = .true.
               test_remove_ele    = .false.

            !test SE_edge: compatible

            !  -------    
            ! | 1 1 1 |
            ! | 1*2 2 |
            ! | 1 2   |
            !  -------
            case(10)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_pt,no_pt,
     $              bc_interior_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 1
               test_j             = 2

               test_bc_section(1) = SE_edge_type
               test_bc_section(2) = 1
               test_bc_section(3) = 2
               test_bc_section(4) = 0
               test_bc_section(5) = 2
               test_compatible    = .true.
               test_remove_ele    = .true.

            !test NW_edge: compatible

            !  -------    
            ! |   2 1 |
            ! | 2 2 1 |
            ! | 1 1*1 |
            !  -------
            case(11)
               grdpts_id = reshape(
     $              (/
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              bc_pt,bc_pt,bc_interior_pt,
     $              no_pt,bc_pt,bc_interior_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 1

               test_bc_section(1) = NW_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 1
               test_bc_section(4) = 0
               test_bc_section(5) = 2
               test_compatible    = .true.
               test_remove_ele    = .true.

            !test NE_edge: compatible

            !  -------    
            ! | 0 1 2 |
            ! | 0 1*1 |
            ! | 0 0 0 |
            !  -------
            case(12)
               grdpts_id = reshape(
     $              (/
     $              interior_pt,interior_pt,interior_pt,
     $              interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,bc_interior_pt,bc_pt
     $              /),
     $              (/3,3/))

               test_i             = 2
               test_j             = 2

               test_bc_section(1) = NE_edge_type
               test_bc_section(2) = 2
               test_bc_section(3) = 2
               test_bc_section(4) = 0
               test_bc_section(5) = 2
               test_compatible    = .true.
               test_remove_ele    = .true.

            end select

        end subroutine make_test_analyse_grdpt_with_bc_section

      end program test_bf_layer_bc_sections
