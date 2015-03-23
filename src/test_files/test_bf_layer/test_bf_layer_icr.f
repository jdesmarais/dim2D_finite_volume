      program test_bf_layer_icr
      
        use bf_layer_icr_class, only :
     $     bf_layer_icr

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       no_pt,
     $       interior_pt,
     $       BF_SUCCESS

        use parameters_input, only :
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        implicit none

        
        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        call check_inputs()


        test_loc = test_is_node_activated(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_is_node_activated: '',L1)', test_loc
        print '()'


        test_loc = test_get_bc_sections(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_bc_sections: '',L1)', test_loc
        print '()'


        contains


        function test_get_bc_sections(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_layer_icr)                   :: bf_layer_used
          integer, dimension(:,:), allocatable :: bc_sections
          integer, dimension(:,:), allocatable :: test_bc_sections

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,2

             !input
             call get_param_test_get_bc_section(
     $            k,
     $            bf_layer_used,
     $            test_bc_sections)

             !output
             call bf_layer_used%get_bc_sections(bc_sections)

             !validation
             test_loc = allocated(bc_sections).eqv.allocated(test_bc_sections)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test allocation('',I2,'') failed'')', k
             end if

             if(allocated(bc_sections)) then

                test_loc = is_int_matrix_validated(
     $               bc_sections,
     $               test_bc_sections,
     $               detailled)

                test_loc = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test bc_sections('',I2,'') failed'')', k
                end if

             end if

          end do

        end function test_get_bc_sections


        subroutine get_param_test_get_bc_section(
     $     test_id,
     $     bf_layer_used,
     $     test_bc_sections)

          implicit none

          integer                             , intent(in)    :: test_id
          type(bf_layer_icr)                  , intent(inout) :: bf_layer_used
          integer, dimension(:,:), allocatable, intent(inout) :: test_bc_sections

          integer :: i,j


          if(test_id.eq.2) then
             allocate(bf_layer_used%bc_sections(5,4))
             bf_layer_used%bc_sections = reshape((/
     $            ((i+5*j,i=1,5),j=1,4)/),
     $            (/5,4/))

             allocate(test_bc_sections(5,4))
             test_bc_sections = reshape((/
     $            ((i+5*j,i=1,5),j=1,4)/),
     $            (/5,4/))
             
          end if

        end subroutine get_param_test_get_bc_section


        function test_is_node_activated(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(bf_layer_icr)           :: bf_layer_used
          type(pmodel_eq)              :: p_model
          integer(ikind), dimension(2) :: loc_coords
          logical                      :: ierror
          logical                      :: node_activated
          logical                      :: test_ierror
          logical                      :: test_node_activated

          integer :: k
          logical :: test_loc


          test_validated = .true.


          do k=1,5

             !input
             call get_param_test_node_activated(
     $            k,
     $            bf_layer_used,
     $            loc_coords,
     $            test_ierror,
     $            test_node_activated)

             !output
             node_activated = bf_layer_used%is_node_activated(
     $            loc_coords,
     $            p_model,
     $            ierror)

             !validation
             test_loc = ierror.eqv.test_ierror
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test ierror('',I2,'') failed'')',k
             end if

             if(test_loc.and.(test_ierror.eqv.BF_SUCCESS)) then
                
                test_loc = node_activated.eqv.test_node_activated
                test_validated = test_validated.and.test_loc
                if(detailled.and.(.not.test_loc)) then
                   print '(''test node_activated('',I2,'') failed'')',k
                end if

             end if

          end do

        end function test_is_node_activated


        subroutine get_param_test_node_activated(
     $     test_id,
     $     bf_layer_used,
     $     loc_coords,
     $     test_ierror,
     $     test_node_activated)

          implicit none

          integer                     , intent(in)    :: test_id
          type(bf_layer_icr)          , intent(inout) :: bf_layer_used
          integer(ikind), dimension(2), intent(out)   :: loc_coords
          logical                     , intent(out)   :: test_ierror
          logical                     , intent(out)   :: test_node_activated

          integer :: i,j


          loc_coords = [1,1]
          test_node_activated = .false.

          select case(test_id)

            !nodes not allocated
            case(1)
               test_ierror = .not.BF_SUCCESS

            !nodes allocated but not enough gridpoints
            case(2)
               allocate(bf_layer_used%x_map(5))
               allocate(bf_layer_used%y_map(5))
               allocate(bf_layer_used%nodes(5,5,ne))
               allocate(bf_layer_used%grdpts_id(5,5))

               test_ierror = .not.BF_SUCCESS

            !nodes allocated and enough grid points but no_pt
            case(3)
               bf_layer_used%grdpts_id = reshape((/
     $              ((no_pt,i=1,5),j=1,5)/),(/5,5/))
               
               loc_coords = [2,2]

               test_ierror = BF_SUCCESS
               test_node_activated = .false.

            !nodes allocated, enough grdpts but desactivated
            case(4)
               bf_layer_used%grdpts_id = reshape((/
     $              ((interior_pt,i=1,5),j=1,5)/),(/5,5/))

               bf_layer_used%nodes(2,2,1) = 1.0d0

               loc_coords = [2,2]

               test_ierror = BF_SUCCESS
               test_node_activated = .false.

            !nodes allocated, enough grdpts and activated
            case(5)
               bf_layer_used%nodes(2,2,1) = -1.0d0

               loc_coords = [2,2]

               test_ierror = BF_SUCCESS
               test_node_activated = .true.

          end select

        end subroutine get_param_test_node_activated

        subroutine check_inputs()

          implicit none

          type(pmodel_eq)                :: p_model
          real(rkind), dimension(3)      :: x_map
          real(rkind), dimension(3)      :: y_map
          real(rkind), dimension(3,3,ne) :: nodes
          
          logical :: test_loc


          nodes(2,2,1) = -1.0d0
          test_loc = p_model%are_openbc_undermined(x_map,y_map,nodes)
          if(.not.test_loc) then
             print '(''the test requires: '')'
             print '(''openbc_undermined if nodes(2,2,1)<0'')'
             stop ''
          end if

          
          nodes(2,2,1) = 1.0d0
          test_loc = .not.p_model%are_openbc_undermined(x_map,y_map,nodes)
          if(.not.test_loc) then
             print '(''the test requires: '')'
             print '(''.not.openbc_undermined if nodes(2,2,1)>0'')'
             stop ''
          end if

        end subroutine check_inputs


      end program test_bf_layer_icr
