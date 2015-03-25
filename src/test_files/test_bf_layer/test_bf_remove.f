      program test_bf_remove

        use bf_remove_module, only :
     $     check_if_bf_layer_remains,
     $     grdpts_common_with_check_layer,     
     $     get_check_line_param,
     $     check_line_neighbors,
     $     check_layer_interior,
     $     check_layer_bf,
     $     no_nogrdpt

        use bf_layer_extract_module, only :
     $       get_bf_layer_match_table

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       align_N,align_S,
     $       align_E,align_W,
     $       no_pt,
     $       search_dcr,
     $       interior_pt,
     $       search_dcr

        use parameters_constant, only :
     $       N,S,E,W

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


        test_loc = test_no_nogrdpt(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_no_nogrdpt: '',L1)', test_loc
        print '()'


        test_loc = test_check_layer_bf(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_check_layer_bf: '',L1)', test_loc
        print '()'


        test_loc = test_check_layer_interior(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_check_layer_interior: '',L1)', test_loc
        print '()'


        test_loc = test_check_line_neighbors(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_check_line_neighbors: '',L1)', test_loc
        print '()'


        test_loc = test_get_check_line_param(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_check_line_param: '',L1)', test_loc
        print '()'


        test_loc = test_grdpts_common_with_check_layer(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_grdpts_common_with_check_layer: '',L1)', test_loc
        print '()'


        test_loc = test_check_if_bf_layer_remains(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_check_if_bf_layer_remains: '',L1)', test_loc
        print '()'


        contains


        function test_check_if_bf_layer_remains(detailled)
     $     result(test_validated)

          implicit none
          
          logical, intent(in) :: detailled
          logical             :: test_validated

          integer                          :: bf_localization
          integer    , dimension(2,2,2)    :: bf_alignment
          integer    , dimension(2,2)      :: bf_match_table
          integer    , dimension(6,6)      :: bf_grdpts_id
          real(rkind), dimension(6)        :: bf_x_map
          real(rkind), dimension(6)        :: bf_y_map
          real(rkind), dimension(6,6,ne)   :: bf_nodes
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          type(pmodel_eq)                  :: p_model
          logical, dimension(2)            :: test_bf_remains
          logical                          :: bf_remains
          
          integer(ikind) :: i,j

          integer :: k
          logical :: test_loc


          test_validated = .true.


          !input
          bf_localization = N

          bf_alignment = reshape((/
     $         align_W+5, align_N, align_W+6, align_N+1,
     $         align_W-7, align_N, align_W-6, align_N+1/),
     $         (/2,2,2/))

          bf_match_table = reshape((/
     $         align_W+5-3,align_N-3,
     $         align_W-7-3,align_N-3/),
     $         (/2,2/))

          bf_grdpts_id = reshape((/
     $         ((interior_pt,i=1,6),j=1,6)/),
     $         (/6,6/))

          bf_nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,6),j=1,6)/),
     $         (/6,6/))

          bf_nodes(2,5,1) = -1.0d0

          interior_nodes(:,:,1) = reshape((/
     $         ((1.0d0,i=1,nx),j=1,ny)/),
     $         (/nx,ny/))

          test_bf_remains = [.true.,.false.]


          do k=1,size(test_bf_remains,1)

             !output
             bf_remains = check_if_bf_layer_remains(
     $            bf_localization,
     $            bf_alignment(:,:,k),
     $            bf_match_table(:,k),
     $            bf_grdpts_id,
     $            bf_x_map,
     $            bf_y_map,
     $            bf_nodes,
     $            interior_x_map,
     $            interior_y_map,
     $            interior_nodes,
     $            p_model)


             !validation
             test_loc = bf_remains.eqv.test_bf_remains(k)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')',k
             end if

          end do

        end function test_check_if_bf_layer_remains


        function test_grdpts_common_with_check_layer(detailled)
     $     result(test_validated)

          implicit none
          
          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(2,2,2) :: test_bf_alignment
          logical, dimension(2)     :: test_grdpts_common
          logical                   :: grdpts_common
          
          integer :: k
          logical :: test_loc
          
          
          test_validated = .true.


          !input
          test_bf_alignment = reshape((/
     $         align_W+5, align_N, align_W+5, align_N,
     $         align_W-6, align_N, align_W-6, align_N/),
     $         (/2,2,2/))

          test_grdpts_common = [.true.,.false.]


          do k=1,size(test_grdpts_common,1)

             !output
             grdpts_common = grdpts_common_with_check_layer(
     $            test_bf_alignment(:,:,k))


             !validation
             test_loc = grdpts_common.eqv.test_grdpts_common(k)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')',k
             end if

          end do

        end function test_grdpts_common_with_check_layer


        function test_get_check_line_param(detailled)
     $     result(test_validated)

          implicit none
          
          logical, intent(in) :: detailled
          logical             :: test_validated


          integer, dimension(8)     :: test_bf_localization
          integer, dimension(2,2,8) :: test_bf_alignment
          integer, dimension(2,2,8) :: test_bf_coords
          integer, dimension(2,2,8) :: test_in_coords

          integer, dimension(2,2)   :: bf_coords
          integer, dimension(2,2)   :: in_coords

          integer :: k
          logical :: test_loc

          test_validated = .true.

          
          test_bf_localization = [N,N,N,S,S,S,W,E]

          test_bf_alignment = reshape((/
     $         align_W-3, align_N  , align_W+2, align_N+5,
     $         align_W+1, align_N  , align_E-1, align_N+5,
     $         align_E-2, align_N  , align_E+3, align_N+5,
     $                                          
     $         align_W-3, align_S-5, align_W+2, align_S  ,
     $         align_W+1, align_S-5, align_E-1, align_S  ,
     $         align_E-2, align_S-5, align_E+3, align_S  ,
     $                                          
     $         align_W-3, align_S+1, align_W  , align_S+2,
     $         align_E  , align_S+1, align_E+3, align_S+2/),
     $         (/2,2,8/))

          test_bf_coords = reshape((/
     $         3,3,   9,6,
     $         2,3,nx-1,6,
     $         2,3,   8,6,
     $         
     $         3,5,   9,8,
     $         2,5,nx-1,8,
     $         2,5,   8,8,
     $         
     $         3,2,6,5,
     $         3,2,6,5
     $         /),
     $         (/2,2,8/))

          test_in_coords = reshape((/
     $                 2, align_N-4,    6, align_N-1,
     $                 2, align_N-4, nx-1, align_N-1,
     $         align_E-4, align_N-4, nx-1, align_N-1,
     $         
     $                 2, align_S+1,    6, align_S+4,
     $                 2, align_S+1, nx-1, align_S+4,
     $         align_E-4, align_S+1, nx-1, align_S+4,
     $         
     $         align_W+1,2,align_W+4,6,
     $         align_E-4,2,align_E-1,6
     $         /),
     $         (/2,2,8/))


          do k=1, size(test_bf_alignment,3)

             !output
             call get_check_line_param(
     $            test_bf_localization(k),
     $            test_bf_alignment(:,:,k),
     $            get_bf_layer_match_table(
     $            test_bf_alignment(:,:,k)),
     $            test_bf_alignment(1,2,k)-test_bf_alignment(1,1,k)+2*bc_size+1,
     $            test_bf_alignment(2,2,k)-test_bf_alignment(2,1,k)+2*bc_size+1,
     $            bf_coords, in_coords)

             !validation
             test_loc = is_int_matrix_validated(
     $            bf_coords,
     $            test_bf_coords(:,:,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test bf('',I2,'') failed'')',k
             end if

             test_loc = is_int_matrix_validated(
     $            in_coords,
     $            test_in_coords(:,:,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test in('',I2,'') failed'')',k
             end if

             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')',k
             end if

          end do

        end function test_get_check_line_param

        
        function test_check_line_neighbors(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer    , dimension(2,2)     :: in_coords
          integer    , dimension(2,2)     :: bf_coords

          integer    , dimension(11,6)    :: bf_grdpts_id
          real(rkind), dimension(11)      :: bf_x_map
          real(rkind), dimension(6)       :: bf_y_map
          real(rkind), dimension(11,6,ne) :: bf_nodes

          real(rkind), dimension(nx)        :: interior_x_map
          real(rkind), dimension(ny)        :: interior_y_map
          real(rkind), dimension(nx,ny,ne)  :: interior_nodes
                                          
          real(rkind), dimension(27*ne)   :: bf_nodes_replaced
          real(rkind), dimension(15*ne)   :: interior_nodes_replaced
          type(pmodel_eq)                 :: p_model

          integer :: i,j,k
          logical :: bf_remains
          logical :: test_loc


          test_validated = .true.


          !test 1: none are activated
          !------------------------------------------------------------
          !input
          in_coords       = reshape((/2,2,6,4/),(/2,2/))
          interior_x_map  = (/ (x_min+(i-1)*0.1d0,i=1,nx) /)
          interior_y_map  = (/ (y_min+(j-1)*0.2d0,j=1,ny) /)
          interior_nodes  = reshape((/ (((1.0d0,i=1,nx),j=1,ny),k=1,ne) /), (/nx,ny,ne/))

          bf_coords       = reshape((/2,2,10,5/),(/2,2/))
          bf_grdpts_id    = reshape((/ ((interior_pt,i=1,11),j=1,6) /), (/11,6/))
          bf_x_map        = (/ (x_min+(i-1)*0.1d0,i=1,11) /)
          bf_y_map        = (/ (y_min+(j-1)*0.2d0,j=1,6) /)
          bf_nodes        = reshape((/ (((1.0d0,i=1,11),j=1,6),k=1,ne) /), (/11,6,ne/))

          !output
          call check_line_neighbors(
     $         bf_coords, in_coords,
     $         bf_grdpts_id,
     $         bf_x_map,
     $         bf_y_map,
     $         bf_nodes,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         bf_remains)

          !validation
          test_loc = .not.bf_remains
          test_validated = test_validated.and.test_loc

          !detailled
          if(detailled.and.(.not.test_loc)) then
             print '(''test(1) failed'')'
          end if


          !test 2: only interior is activated
          !------------------------------------------------------------
          !input
          interior_nodes_replaced    = (/((1.0d0,i=1,15),k=1,ne)/)
          interior_nodes_replaced(1) = -1.0d0
             
          interior_nodes(2:6,2:4,:) = reshape(interior_nodes_replaced, (/5,3,ne/))

          !output
          call check_line_neighbors(
     $         bf_coords, in_coords,
     $         bf_grdpts_id,
     $         bf_x_map,
     $         bf_y_map,
     $         bf_nodes,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         bf_remains)

          !validation
          test_loc = bf_remains
          test_validated = test_validated.and.test_loc

          !detailled
          if(detailled.and.(.not.test_loc)) then
             print '(''test(2) failed'')'
          end if

          
          !test 3: interior and buffer are activated
          !------------------------------------------------------------
          !input
          bf_nodes_replaced = (/((1.0d0,i=1,27),k=1,ne)/)
          bf_nodes_replaced(1) = -1.0d0

          bf_nodes(2:10,3:5,:) = reshape(bf_nodes_replaced, (/9,3,ne/))

          !output
          call check_line_neighbors(
     $         bf_coords, in_coords,
     $         bf_grdpts_id,
     $         bf_x_map,
     $         bf_y_map,
     $         bf_nodes,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         bf_remains)

          !validation
          test_loc = bf_remains
          test_validated = test_validated.and.test_loc

          !detailled
          if(detailled.and.(.not.test_loc)) then
             print '(''test(3) failed'')'
          end if


          !test 4: only buffer is activated
          !------------------------------------------------------------
          !input
          interior_nodes = reshape((/ (((1.0d0,i=1,nx),j=1,ny),k=1,ne) /), (/nx,ny,ne/))

          !output
          call check_line_neighbors(
     $         bf_coords, in_coords,
     $         bf_grdpts_id,
     $         bf_x_map,
     $         bf_y_map,
     $         bf_nodes,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         p_model,
     $         bf_remains)

          !validation
          test_loc = bf_remains
          test_validated = test_validated.and.test_loc

          !detailled
          if(detailled.and.(.not.test_loc)) then
             print '(''test(4) failed'')'
          end if


        end function test_check_line_neighbors


        function test_check_layer_interior(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer    , dimension(2,2)      :: pt_coords
          real(rkind), dimension(nx)       :: x_map
          real(rkind), dimension(ny)       :: y_map
          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(15*ne)    :: nodes_replaced
          type(pmodel_eq)                  :: p_model

          integer :: i,j,k
          logical :: bf_remains
          logical :: test_loc

          test_validated = .true.


          !test 1: no node activated
          !------------------------------------------------------------
          !input
          pt_coords = reshape((/2,2,6,4/),(/2,2/))
          x_map     = (/ (x_min+(i-1)*0.1d0,i=1,nx) /)
          y_map     = (/ (y_min+(j-1)*0.2d0,j=1,ny) /)
          nodes     = reshape((/ (((1.0d0,i=1,nx),j=1,ny),k=1,ne) /), (/nx,ny,ne/))
          
          call check_layer_interior(
     $         pt_coords,
     $         x_map,
     $         y_map,
     $         nodes,
     $         p_model,
     $         bf_remains)
          test_loc = .not.bf_remains
          test_validated = test_validated.and.test_loc

          if(detailled.and.(.not.test_loc)) then
             print '(''test(1) failed'')'
          end if


          !test 2-16: only one node activated
          !------------------------------------------------------------
          do k=1,15
             
             !input
             nodes_replaced    = (/((1.0d0,i=1,15),k=1,ne)/)
             nodes_replaced(k) = -1.0d0
             
             nodes(2:6,2:4,:) = reshape(nodes_replaced, (/5,3,ne/))

             !output
             call check_layer_interior(
     $            pt_coords,
     $            x_map,
     $            y_map,
     $            nodes,
     $            p_model,
     $            bf_remains)
             test_loc = bf_remains
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k+1
             end if

          end do


        end function test_check_layer_interior


        function test_check_layer_bf(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer    , dimension(2,2)    :: pt_coords
          integer    , dimension(7,8)    :: grdpts_id
          real(rkind), dimension(7)      :: x_map
          real(rkind), dimension(8)      :: y_map
          real(rkind), dimension(7,8,ne) :: nodes
          real(rkind), dimension(15*ne)  :: nodes_replaced
          type(pmodel_eq)                :: p_model

          integer :: i,j,k
          logical :: bf_remains
          logical :: test_loc

          test_validated = .true.

          !test 1: all nodes are desactivated
          !------------------------------------------------------------
          !input
          pt_coords = reshape((/2,2,6,4/),(/2,2/))
          grdpts_id = reshape((/ ((interior_pt,i=1,7),j=1,8) /), (/7,8/))
          x_map     = (/ (x_min+(i-1)*0.1d0,i=1,7) /)
          y_map     = (/ (y_min+(j-1)*0.2d0,j=1,8) /)
          nodes     = reshape((/ (((1.0d0,i=1,7),j=1,8),k=1,ne) /), (/7,8,ne/))
          
          call check_layer_bf(
     $         pt_coords,
     $         grdpts_id,
     $         x_map,
     $         y_map,
     $         nodes,
     $         p_model,
     $         bf_remains)
          test_loc = .not.bf_remains
          test_validated = test_validated.and.test_loc

          if(detailled.and.(.not.test_loc)) then
             print '(''test(1) failed'')'
          end if


          !test 2-16: only one node in the nodes checked is activated
          !------------------------------------------------------------
          do k=1,15
             
             !input
             nodes_replaced    = (/((1.0d0,i=1,15),k=1,ne)/)
             nodes_replaced(k) = -1.0d0
             
             nodes(2:6,2:4,:) = reshape(nodes_replaced, (/5,3,ne/))

             !output
             call check_layer_bf(
     $            pt_coords,
     $            grdpts_id,
     $            x_map,
     $            y_map,
     $            nodes,
     $            p_model,
     $            bf_remains)
             test_loc = bf_remains
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k+1
             end if

          end do

          !test 17: no node can be checked
          !------------------------------------------------------------
          !input
          grdpts_id = reshape((/ ((no_pt,i=1,7),j=1,8) /), (/7,8/))

          !output
          call check_layer_bf(
     $         pt_coords,
     $         grdpts_id,
     $         x_map,
     $         y_map,
     $         nodes,
     $         p_model,
     $         bf_remains)
          test_loc = .not.bf_remains
          test_validated = test_validated.and.test_loc
          
          !detailled
          if(detailled.and.(.not.test_loc)) then
             print '(''test(17) failed'')'
          end if

        end function test_check_layer_bf


        function test_no_nogrdpt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer                 :: i,k
          integer, dimension(9)   :: grdpts_line
          integer, dimension(3,3) :: grdpts

          logical :: test_loc

          
          test_validated = .true.


          do k=1,9
             
             grdpts_line    = (/(interior_pt,i=1,9)/)
             grdpts_line(k) = no_pt

             grdpts = reshape(grdpts_line, (/3,3/))

             test_loc = .not.(no_nogrdpt(grdpts,2,2))
             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')',k
             end if

          end do


          grdpts_line = (/(interior_pt,i=1,9)/)
          grdpts = reshape(grdpts_line,(/3,3/))
          test_loc = no_nogrdpt(grdpts,2,2)
          test_validated = test_validated.and.test_loc
          

          if(detailled.and.(.not.test_loc)) then
             print '(''test(10) failed'')'
          end if

        end function test_no_nogrdpt


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

      end program test_bf_remove
