      program test_bf_newgrdpt_verification

        use bf_newgrdpt_verification_module, only :
     $       are_grdpts_available,
     $       get_newgrdpt_verification_bounds

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
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
     $       SW_corner_type,
     $       no_gradient_type,
     $       gradient_I_type,
     $       gradient_L0_type,
     $       gradient_R0_type,
     $       gradient_xLR0_yI_type,
     $       gradient_xI_yLR0_type,
     $       gradient_xLR0_yLR0_type,
     $       
     $       interior_pt,
     $       no_pt

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled = .true.
        test_validated = .true.

        test_loc = test_are_grdpts_available(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_are_grdpts_available: '',L1)', test_loc
        print '()'

        test_loc = test_get_newgrdpt_verification_bounds(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_get_newgrdpt_verification_bounds: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated

        contains


        function test_are_grdpts_available(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer, dimension(4,6) :: grdpts_id
          integer, dimension(2,2) :: gen_coords
          integer, dimension(6)   :: grdpts_id_replaced

          integer :: i,j
          integer :: k


          test_validated = .true.


          gen_coords = reshape((/2,4,3,6/),(/2,2/))

          grdpts_id = reshape((/((no_pt,i=1,4),j=1,6)/),(/4,6/))
          grdpts_id_replaced = (/(interior_pt,i=1,6)/)
          grdpts_id(2:3,4:6) = reshape(grdpts_id_replaced,(/2,3/))

          test_loc = are_grdpts_available(grdpts_id,gen_coords)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test( 0) failed'')'
          end if

          do k=1,6
             grdpts_id_replaced    = (/(interior_pt,i=1,6)/)
             grdpts_id_replaced(k) = no_pt

             grdpts_id(2:3,4:6) = reshape(grdpts_id_replaced,(/2,3/))

             test_loc = .not.are_grdpts_available(grdpts_id,gen_coords)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')',k
             end if
          end do

          grdpts_id = reshape((/((interior_pt,i=1,4),j=1,6)/),(/4,6/))
          grdpts_id_replaced = (/(interior_pt,i=1,6)/)
          grdpts_id(2:3,4:6) = reshape(grdpts_id_replaced,(/2,3/))

          test_loc = are_grdpts_available(grdpts_id,gen_coords)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test( 7) failed'')'
          end if

          do k=1,6
             grdpts_id_replaced    = (/(interior_pt,i=1,6)/)
             grdpts_id_replaced(k) = no_pt

             grdpts_id(2:3,4:6) = reshape(grdpts_id_replaced,(/2,3/))

             test_loc = .not.are_grdpts_available(grdpts_id,gen_coords)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')',7+k
             end if
          end do

        end function test_are_grdpts_available


        function test_get_newgrdpt_verification_bounds(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(7,7,32) :: test_grdpts
          integer, dimension(32)     :: test_procedure_type
          integer, dimension(32)     :: test_gradient_type

          integer                   :: k,l
          integer                   :: nb_bounds
          integer, dimension(2,2,2) :: bounds

          integer, dimension(7,7)   :: grdpts
          integer                   :: i
          integer                   :: j
          
          test_validated = .true.

          !corners
          test_procedure_type(1) = SW_corner_type
          test_grdpts(:,:,1) = reshape((/
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,1,1,1,
     $         0,0,0,0,1,1,1,
     $         0,0,0,0,1,1,1/),
     $         (/7,7/))

          test_procedure_type(2) = SE_corner_type
          test_grdpts(:,:,2)     = reflect_x(test_grdpts(:,:,1))

         
          test_procedure_type(3) = NW_corner_type
          test_grdpts(:,:,3)     = reflect_y(test_grdpts(:,:,1))
          
          test_procedure_type(4) = NE_corner_type
          test_grdpts(:,:,4)     = reflect_y(test_grdpts(:,:,2))

          !E_edge
          test_procedure_type(5) = E_edge_type
          test_gradient_type(5) = gradient_I_type
          
          test_grdpts(:,:,5) = reshape((/
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,1,1,0,0,0,0,
     $         0,1,1,0,0,0,0,
     $         0,1,1,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0/),
     $         (/7,7/))

          test_procedure_type(6) = E_edge_type
          test_gradient_type(6)  = gradient_L0_type
          
          test_grdpts(:,:,6) = reshape((/
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,1,1,0,0,0,0,
     $         0,1,1,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0/),
     $         (/7,7/))

          test_procedure_type(7) = E_edge_type
          test_gradient_type(7) = gradient_R0_type
          
          test_grdpts(:,:,7) = reshape((/
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,1,1,0,0,0,0,
     $         0,1,1,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0/),
     $         (/7,7/))

          !W_edge
          test_procedure_type(8) = W_edge_type
          test_gradient_type(8) = gradient_I_type
          test_grdpts(:,:,8) = reflect_x(test_grdpts(:,:,5))

          test_procedure_type(9) = W_edge_type
          test_gradient_type(9) = gradient_L0_type
          test_grdpts(:,:,9) = reflect_x(test_grdpts(:,:,6))

          test_procedure_type(10) = W_edge_type
          test_gradient_type(10) = gradient_R0_type
          test_grdpts(:,:,10) = reflect_x(test_grdpts(:,:,7))

          !N_edge
          test_procedure_type(11) =N_edge_type
          test_gradient_type(11) = gradient_I_type
          test_grdpts(:,:,11) = transpose(test_grdpts(:,:,5))

          test_procedure_type(12) =N_edge_type
          test_gradient_type(12) = gradient_L0_type
          test_grdpts(:,:,12) = transpose(test_grdpts(:,:,6))

          test_procedure_type(13) =N_edge_type
          test_gradient_type(13) = gradient_R0_type
          test_grdpts(:,:,13) = transpose(test_grdpts(:,:,7))

          !S_edge
          test_procedure_type(14) =S_edge_type
          test_gradient_type(14) = gradient_I_type
          test_grdpts(:,:,14) = transpose(test_grdpts(:,:,8))

          test_procedure_type(15) =S_edge_type
          test_gradient_type(15) = gradient_L0_type
          test_grdpts(:,:,15) = transpose(test_grdpts(:,:,9))

          test_procedure_type(16) =S_edge_type
          test_gradient_type(16) = gradient_R0_type
          test_grdpts(:,:,16) = transpose(test_grdpts(:,:,10))

          !anti_corners
          !SW_anti_corner
          test_procedure_type(17) = SW_edge_type
          test_gradient_type(17) = gradient_I_type
          
          test_grdpts(:,:,17) = reshape((/
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,1,1,0,
     $         0,0,0,0,1,1,0,
     $         0,0,1,1,1,1,0,
     $         0,0,1,1,1,1,0,
     $         0,0,0,0,0,0,0/),
     $         (/7,7/))

          test_procedure_type(18) = SW_edge_type
          test_gradient_type(18)  = gradient_xI_yLR0_type
          
          test_grdpts(:,:,18) = reshape((/
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,1,1,0,
     $         0,0,1,1,1,1,0,
     $         0,0,1,1,1,1,0,
     $         0,0,0,0,0,0,0/),
     $         (/7,7/))

          test_procedure_type(19) = SW_edge_type
          test_gradient_type(19)  = gradient_xLR0_yI_type
          
          test_grdpts(:,:,19) = reshape((/
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,1,1,0,
     $         0,0,0,0,1,1,0,
     $         0,0,0,1,1,1,0,
     $         0,0,0,1,1,1,0,
     $         0,0,0,0,0,0,0/),
     $         (/7,7/))

          test_procedure_type(20) = SW_edge_type
          test_gradient_type(20)  = gradient_xLR0_yLR0_type
          
          test_grdpts(:,:,20) = reshape((/
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,0,0,0,
     $         0,0,0,0,1,1,0,
     $         0,0,0,1,1,1,0,
     $         0,0,0,1,1,1,0,
     $         0,0,0,0,0,0,0/),
     $         (/7,7/))
          
          !SE_anti_corner
          test_procedure_type(21) =SE_edge_type
          test_gradient_type(21) = gradient_I_type
          test_grdpts(:,:,21) = reflect_x(test_grdpts(:,:,17))

          test_procedure_type(22) =SE_edge_type
          test_gradient_type(22) = gradient_xI_yLR0_type
          test_grdpts(:,:,22) = reflect_x(test_grdpts(:,:,18))

          test_procedure_type(23) =SE_edge_type
          test_gradient_type(23) = gradient_xLR0_yI_type
          test_grdpts(:,:,23) = reflect_x(test_grdpts(:,:,19))

          test_procedure_type(24) =SE_edge_type
          test_gradient_type(24) = gradient_xLR0_yLR0_type
          test_grdpts(:,:,24) = reflect_x(test_grdpts(:,:,20))

          !NE_anti_corner
          test_procedure_type(25) =NE_edge_type
          test_gradient_type(25) = gradient_I_type
          test_grdpts(:,:,25) = reflect_y(test_grdpts(:,:,21))

          test_procedure_type(26) =NE_edge_type
          test_gradient_type(26) = gradient_xI_yLR0_type
          test_grdpts(:,:,26) = reflect_y(test_grdpts(:,:,22))

          test_procedure_type(27) =NE_edge_type
          test_gradient_type(27) = gradient_xLR0_yI_type
          test_grdpts(:,:,27) = reflect_y(test_grdpts(:,:,23))

          test_procedure_type(28) =NE_edge_type
          test_gradient_type(28) = gradient_xLR0_yLR0_type
          test_grdpts(:,:,28) = reflect_y(test_grdpts(:,:,24))

          !NW_anti_corner
          test_procedure_type(29) =NW_edge_type
          test_gradient_type(29) = gradient_I_type
          test_grdpts(:,:,29) = reflect_x(test_grdpts(:,:,25))

          test_procedure_type(30) =NW_edge_type
          test_gradient_type(30) = gradient_xI_yLR0_type
          test_grdpts(:,:,30) = reflect_x(test_grdpts(:,:,26))

          test_procedure_type(31) =NW_edge_type
          test_gradient_type(31) = gradient_xLR0_yI_type
          test_grdpts(:,:,31) = reflect_x(test_grdpts(:,:,27))

          test_procedure_type(32) =NW_edge_type
          test_gradient_type(32) = gradient_xLR0_yLR0_type
          test_grdpts(:,:,32) = reflect_x(test_grdpts(:,:,28))


          do k=1,32

             !output
             call get_newgrdpt_verification_bounds(
     $            test_procedure_type(k),
     $            test_gradient_type(k),
     $            nb_bounds,
     $            bounds)

             !initialize grdpts
             grdpts = reshape((/
     $            0,0,0,0,0,0,0,
     $            0,0,0,0,0,0,0,
     $            0,0,0,0,0,0,0,
     $            0,0,0,0,0,0,0,
     $            0,0,0,0,0,0,0,
     $            0,0,0,0,0,0,0,
     $            0,0,0,0,0,0,0/),
     $            (/7,7/))

             !set as 1 the grdpts tested
             do l=1,nb_bounds
                
                do j=4+bounds(2,1,l),4+bounds(2,2,l)
                   do i=4+bounds(1,1,l),4+bounds(1,2,l)

                      grdpts(i,j) = 1

                   end do
                end do

             end do

             !compare grdpts with test_grdpts(:,:,k)
             test_loc = is_int_matrix_validated(
     $            grdpts,
     $            test_grdpts(:,:,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             
             !detailled
             if((.not.test_loc).and.detailled) then
                print '(''test('',I2,''): '',L1)', k, test_loc
             end if

          end do

        end function test_get_newgrdpt_verification_bounds


        function reflect_x(grdpts)
     $     result(grdpts_r)

          implicit none

          integer, dimension(7,7), intent(in) :: grdpts
          integer, dimension(7,7)             :: grdpts_r

          integer :: i,j
          integer :: i_r

          do j=1,7
             do i=1,7
                i_r = 7-i+1
                grdpts_r(i,j) = grdpts(i_r,j)
             end do
          end do

        end function reflect_x


        function reflect_y(grdpts)
     $     result(grdpts_r)

          implicit none

          integer, dimension(7,7), intent(in) :: grdpts
          integer, dimension(7,7)             :: grdpts_r

          integer :: i,j
          integer :: j_r

          do j=1,7
             j_r = 7-j+1
             do i=1,7
                grdpts_r(i,j) = grdpts(i,j_r)
             end do
          end do

        end function reflect_y

      end program test_bf_newgrdpt_verification
