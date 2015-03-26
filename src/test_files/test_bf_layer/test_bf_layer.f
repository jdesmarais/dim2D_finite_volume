      program test_bf_layer

        use bf_layer_class, only :
     $     bf_layer

        use parameters_bf_layer, only :
     $       interior_pt,
     $       align_N, align_S,
     $       align_E, align_W,
     $       search_dcr

        use parameters_constant, only :
     $       N

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size,
     $       x_min,y_min

        use parameters_kind, only :
     $       ikind, rkind

        use pmodel_eq_class, only :
     $       pmodel_eq
        

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        detailled = .true.
        test_validated = .true.


        call check_inputs()
        

        test_loc = test_should_remain(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_should_remain: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_should_remain(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_layer)                   :: bf_layer_used
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          type(pmodel_eq)                  :: p_model

          integer :: i,j
          logical :: should_remain
          logical :: test_loc


          test_validated = .true.


          !input
          interior_nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,nx),j=1,ny)/),(/nx,ny/))

          call bf_layer_used%ini(N)

          call bf_layer_used%allocate_bf_layer(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         reshape((/
     $           align_W+5,align_N,align_W+6,align_N+1/),
     $           (/2,2/))
     $         )

          bf_layer_used%grdpts_id = reshape((/
     $         ((interior_pt,i=1,6),j=1,6)/),(/6,6/))

          bf_layer_used%nodes(3,5,1) =-1.0d0


          !output
          should_remain = bf_layer_used%should_remain(
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model)


          !validation
          test_loc = should_remain
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test should_remain(1) failed'')'
          end if

          test_loc = bf_layer_used%get_remain_status()
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test should_remain(2) failed'')'
          end if 

        end function test_should_remain


        subroutine check_inputs()

          implicit none


          type(pmodel_eq)                :: p_model
          real(rkind), dimension(3)      :: x_map
          real(rkind), dimension(3)      :: y_map
          real(rkind), dimension(3,3,ne) :: nodes
          
          logical :: test_loc

          if(search_dcr.ne.4) then

             print '(''the test requires: '')'
             print '(''search_dcr=4'')'
             stop ''
             
          end if

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

      end program test_bf_layer
