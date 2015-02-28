      program test_bf_bc_sections_overlap
      
        use bf_bc_sections_overlap_module, only :
     $     determine_corner_or_anti_corner_grdpts_computed,
     $     overlap_N,
     $     overlap_S,
     $     overlap_E,
     $     overlap_W,
     $     overlap_bc_section_by_integration_borders

        use check_data_module, only :
     $       is_int_matrix_validated

        use parameters_bf_layer, only :
     $       no_overlap,
     $       N_overlap,
     $       S_overlap,
     $       E_overlap,
     $       W_overlap,
     $       NE_overlap,
     $       NW_overlap,
     $       SE_overlap,
     $       SW_overlap,
     $       NS_overlap,
     $       EW_overlap,
     $       
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NE_corner_type,
     $       NW_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       NE_edge_type,
     $       NW_edge_type,
     $       SE_edge_type,
     $       SW_edge_type

        use parameters_kind, only :
     $       ikind

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        
        detailled      = .true.
        test_validated = .true.

        
        test_loc = test_determine_corner_or_anti_corner_grdpts_computed(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_bf_bc_sections_overlap: '',L1)', test_loc
        print '()'


        test_loc = test_overlap_N(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_overlap_N: '',L1)', test_loc
        print '()'


        test_loc = test_overlap_S(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_overlap_S: '',L1)', test_loc
        print '()'


        test_loc = test_overlap_E(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_overlap_E: '',L1)', test_loc
        print '()'


        test_loc = test_overlap_W(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_overlap_W: '',L1)', test_loc
        print '()'


        test_loc = test_overlap_bc_section_by_integration_borders(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_overlap_bc_section_by_integration_borders: '',L1)', test_loc
        print '()' 


        contains


        function test_determine_corner_or_anti_corner_grdpts_computed(
     $       detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer                  :: k
          integer, dimension(11)   :: overlap_test
          integer, dimension(4,11) :: compute_pt_test
          logical                  :: compute_point1
          logical                  :: compute_point2
          logical                  :: compute_point3
          logical                  :: compute_point4
          logical                  :: test_loc

          
          test_validated = .true.


          !input
          overlap_test = [
     $         no_overlap,
     $         N_overlap,
     $         S_overlap,
     $         E_overlap,
     $         W_overlap,
     $         NE_overlap,
     $         NW_overlap,
     $         SE_overlap,
     $         SW_overlap,
     $         NS_overlap,
     $         EW_overlap]

          compute_pt_test(:,1) = [.true.,.true.,.true.,.true.]     !no_overlap
          compute_pt_test(:,2) = [.true.,.true.,.false.,.false.]   !N
          compute_pt_test(:,3) = [.false.,.false.,.true.,.true.]   !S
          compute_pt_test(:,4) = [.true.,.false.,.true.,.false.]   !E
          compute_pt_test(:,5) = [.false.,.true.,.false.,.true.]   !W
          compute_pt_test(:,6) = [.true.,.false.,.false.,.false.]  !NE
          compute_pt_test(:,7) = [.false.,.true.,.false.,.false.]  !NW
          compute_pt_test(:,8) = [.false.,.false.,.true.,.false.]  !SE
          compute_pt_test(:,9) = [.false.,.false.,.false.,.true.]  !SW
          compute_pt_test(:,10)= [.false.,.false.,.false.,.false.] !NS
          compute_pt_test(:,11)= [.false.,.false.,.false.,.false.] !EW


          do k=1, size(overlap_test,1)


             !output
             call determine_corner_or_anti_corner_grdpts_computed(
     $            overlap_test(k),
     $            compute_point1,
     $            compute_point2,
     $            compute_point3,
     $            compute_point4)

             
             !validation
             test_loc =
     $            (compute_point1.eqv.compute_pt_test(1,k)).and.
     $            (compute_point2.eqv.compute_pt_test(2,k)).and.
     $            (compute_point3.eqv.compute_pt_test(3,k)).and.
     $            (compute_point4.eqv.compute_pt_test(4,k))
             test_validated = test_validated.and.test_loc
             

             !detailled
             if(detailled.and.(.not.test_loc)) then

                print '(''test '',I2,'' failed'')', k
                print '(''    compute_point1: '',L1,'' -> '',L1)', compute_point1, compute_pt_test(1,k)
                print '(''    compute_point2: '',L1,'' -> '',L1)', compute_point2, compute_pt_test(2,k)
                print '(''    compute_point3: '',L1,'' -> '',L1)', compute_point3, compute_pt_test(3,k)
                print '(''    compute_point4: '',L1,'' -> '',L1)', compute_point4, compute_pt_test(4,k)
                print '()'

             end if

          end do

        end function test_determine_corner_or_anti_corner_grdpts_computed


        function test_overlap_N(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(11) :: overlap_test
          integer, dimension(11) :: overlap_after_test
          integer                :: k
          logical                :: test_loc

          test_validated = .true.

          !input
          overlap_test = [
     $         no_overlap,
     $         N_overlap,
     $         S_overlap,
     $         E_overlap,
     $         W_overlap,
     $         NE_overlap,
     $         NW_overlap,
     $         SE_overlap,
     $         SW_overlap,
     $         NS_overlap,
     $         EW_overlap]

          overlap_after_test = [
     $         N_overlap,
     $         N_overlap,
     $         NS_overlap,
     $         NE_overlap,
     $         NW_overlap,
     $         NE_overlap,
     $         NW_overlap,
     $         NS_overlap,
     $         NS_overlap,
     $         NS_overlap,
     $         EW_overlap]

          do k=1,size(overlap_test,1)

             !output
             call overlap_N(overlap_test(k))

             !validation
             test_loc = overlap_test(k).eqv.overlap_after_test(k)
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                
                print '(''test('',I2,'') failed: '',L1,'' -> '',L1)',
     $               k, overlap_test(k), overlap_after_test(k)
                
             end if

          end do

        end function test_overlap_N


        function test_overlap_S(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(11) :: overlap_test
          integer, dimension(11) :: overlap_after_test
          integer                :: k
          logical                :: test_loc

          test_validated = .true.

          !input
          overlap_test = [
     $         no_overlap,
     $         N_overlap,
     $         S_overlap,
     $         E_overlap,
     $         W_overlap,
     $         NE_overlap,
     $         NW_overlap,
     $         SE_overlap,
     $         SW_overlap,
     $         NS_overlap,
     $         EW_overlap]

          overlap_after_test = [
     $         S_overlap,
     $         NS_overlap,
     $         S_overlap,
     $         SE_overlap,
     $         SW_overlap,
     $         NS_overlap,
     $         NS_overlap,
     $         SE_overlap,
     $         SW_overlap,
     $         NS_overlap,
     $         EW_overlap]

          do k=1,size(overlap_test,1)

             !output
             call overlap_S(overlap_test(k))

             !validation
             test_loc = overlap_test(k).eqv.overlap_after_test(k)
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                
                print '(''test('',I2,'') failed: '',L1,'' -> '',L1)',
     $               k, overlap_test(k), overlap_after_test(k)
                
             end if

          end do

        end function test_overlap_S


        function test_overlap_E(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(11) :: overlap_test
          integer, dimension(11) :: overlap_after_test
          integer                :: k
          logical                :: test_loc

          test_validated = .true.

          !input
          overlap_test = [
     $         no_overlap,
     $         N_overlap,
     $         S_overlap,
     $         E_overlap,
     $         W_overlap,
     $         NE_overlap,
     $         NW_overlap,
     $         SE_overlap,
     $         SW_overlap,
     $         NS_overlap,
     $         EW_overlap]

          overlap_after_test = [
     $         E_overlap,
     $         NE_overlap,
     $         SE_overlap,
     $         E_overlap,
     $         EW_overlap,
     $         NE_overlap,
     $         EW_overlap,
     $         SE_overlap,
     $         EW_overlap,
     $         NS_overlap,
     $         EW_overlap]

          do k=1,size(overlap_test,1)

             !output
             call overlap_E(overlap_test(k))

             !validation
             test_loc = overlap_test(k).eqv.overlap_after_test(k)
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                
                print '(''test('',I2,'') failed: '',L1,'' -> '',L1)',
     $               k, overlap_test(k), overlap_after_test(k)
                
             end if

          end do

        end function test_overlap_E


        function test_overlap_W(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(11) :: overlap_test
          integer, dimension(11) :: overlap_after_test
          integer                :: k
          logical                :: test_loc

          test_validated = .true.

          !input
          overlap_test = [
     $         no_overlap,
     $         N_overlap,
     $         S_overlap,
     $         E_overlap,
     $         W_overlap,
     $         NE_overlap,
     $         NW_overlap,
     $         SE_overlap,
     $         SW_overlap,
     $         NS_overlap,
     $         EW_overlap]

          overlap_after_test = [
     $         W_overlap,
     $         NW_overlap,
     $         SW_overlap,
     $         EW_overlap,
     $         W_overlap,
     $         EW_overlap,
     $         NW_overlap,
     $         EW_overlap,
     $         SW_overlap,
     $         NS_overlap,
     $         EW_overlap]

          do k=1,size(overlap_test,1)

             !output
             call overlap_W(overlap_test(k))

             !validation
             test_loc = overlap_test(k).eqv.overlap_after_test(k)
             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                
                print '(''test('',I2,'') failed: '',L1,'' -> '',L1)',
     $               k, overlap_test(k), overlap_after_test(k)
                
             end if

          end do

        end function test_overlap_W


        function test_overlap_bc_section_by_integration_borders(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer, dimension(15)     :: test_edge_overlap
          integer, dimension(15)     :: test_edge_overlap_after
          integer, dimension(2,2,15) :: test_edge_borders
          integer, dimension(2,2,15) :: test_edge_borders_after

          integer, dimension(25)     :: test_square_overlap
          integer, dimension(25)     :: test_square_overlap_after
          integer, dimension(2,2,25) :: test_square_borders
          integer, dimension(2,2,25) :: test_square_borders_after

          integer                        :: i
          integer(ikind), dimension(2,2) :: gen_borders
          logical                        :: test_loc

          test_validated = .true.

          ![1,10]x[0,20]
          gen_borders(1,1) = 2
          gen_borders(1,2) = 9
          gen_borders(2,1) = 1
          gen_borders(2,2) = 19          

          !test of N_edge and S_edge bc_sections
          test_edge_overlap = (/(no_overlap, i=1,15)/)

          test_edge_overlap_after = [
     $         NS_overlap,
     $         NS_overlap,
     $         NS_overlap,
     $         N_overlap,
     $         N_overlap,
     $         N_overlap,
     $         no_overlap,
     $         no_overlap,
     $         no_overlap,
     $         S_overlap,
     $         S_overlap,
     $         S_overlap,
     $         NS_overlap,
     $         NS_overlap,
     $         NS_overlap]

          test_edge_borders = reshape((/
     $         1, 20, 4, 21,
     $         7, 20, 9, 21,
     $         9, 20,12, 21,
     $         1, 19, 4, 20,
     $         7, 19, 9, 20,
     $         9, 19,12, 20,
     $         1, 10, 4, 11,
     $         7, 10, 9, 11,
     $         9, 10,12, 11,
     $         1,  0, 4,  1,
     $         7,  0, 9,  1,
     $         9,  0,12,  1,
     $         1, -1, 4,  0,
     $         7, -1, 9,  0,
     $         9, -1,12,  0/),
     $         (/2,2,15/))

          test_edge_borders_after = reshape((/
     $         2, 20, 4, 21,
     $         7, 20, 9, 21,
     $         9, 20, 9, 21,
     $         2, 19, 4, 20,
     $         7, 19, 9, 20,
     $         9, 19, 9, 20,
     $         2, 10, 4, 11,
     $         7, 10, 9, 11,
     $         9, 10, 9, 11,
     $         2,  0, 4,  1,
     $         7,  0, 9,  1,
     $         9,  0, 9,  1,
     $         2, -1, 4,  0,
     $         7, -1, 9,  0,
     $         9, -1, 9,  0/),
     $         (/2,2,15/))

          !N_edge test
          test_loc = test_overlap_bc_section(
     $         N_edge_type,
     $         test_edge_borders,
     $         test_edge_borders_after,
     $         test_edge_overlap,
     $         test_edge_overlap_after,
     $         gen_borders,
     $         detailled)
          test_validated = test_validated.and.test_loc

          !S_edge test
          test_loc = test_overlap_bc_section(
     $         S_edge_type,
     $         test_edge_borders,
     $         test_edge_borders_after,
     $         test_edge_overlap,
     $         test_edge_overlap_after,
     $         gen_borders,
     $         detailled)
          test_validated = test_validated.and.test_loc


          !test of E_edge and W_edge bc_sections
          test_edge_overlap = (/(no_overlap, i=1,15)/)

          test_edge_overlap_after = [
     $         EW_overlap,
     $         W_overlap,
     $         no_overlap,
     $         E_overlap,
     $         EW_overlap,
     $         EW_overlap,
     $         W_overlap,
     $         no_overlap,
     $         E_overlap,
     $         EW_overlap,
     $         EW_overlap,
     $         W_overlap,
     $         no_overlap,
     $         E_overlap,
     $         EW_overlap]

          test_edge_borders = reshape((/
     $         0, 18,  1, 21,
     $         1, 18,  2, 21,
     $         6, 18,  7, 21,
     $         9, 18, 10, 21,
     $        10, 18, 11, 21,
     $         0, 10,  1, 15,
     $         1, 10,  2, 15,
     $         6, 10,  7, 15,
     $         9, 10, 10, 15,
     $        10, 10, 11, 15,
     $         0, -1,  1,  2,
     $         1, -1,  2,  2,
     $         6, -1,  7,  2,
     $         9, -1, 10,  2,
     $        10, -1, 11,  2/),
     $         (/2,2,15/))

          test_edge_borders_after = reshape((/
     $         0, 18,  1, 19,
     $         1, 18,  2, 19,
     $         6, 18,  7, 19,
     $         9, 18, 10, 19,
     $        10, 18, 11, 19,
     $         0, 10,  1, 15,
     $         1, 10,  2, 15,
     $         6, 10,  7, 15,
     $         9, 10, 10, 15,
     $        10, 10, 11, 15,
     $         0,  1,  1,  2,
     $         1,  1,  2,  2,
     $         6,  1,  7,  2,
     $         9,  1, 10,  2,
     $        10,  1, 11,  2/),
     $         (/2,2,15/))

          !E_edge test
          test_loc = test_overlap_bc_section(
     $         E_edge_type,
     $         test_edge_borders,
     $         test_edge_borders_after,
     $         test_edge_overlap,
     $         test_edge_overlap_after,
     $         gen_borders,
     $         detailled)
          test_validated = test_validated.and.test_loc

          !W_edge test
          test_loc = test_overlap_bc_section(
     $         W_edge_type,
     $         test_edge_borders,
     $         test_edge_borders_after,
     $         test_edge_overlap,
     $         test_edge_overlap_after,
     $         gen_borders,
     $         detailled)
          test_validated = test_validated.and.test_loc

          !test of corners and anti_corners
          test_square_overlap = (/(no_overlap, i=1,25)/)

          test_square_overlap_after = [
     $         NS_overlap,
     $         NS_overlap,
     $         NS_overlap,
     $         NS_overlap,
     $         NS_overlap,
     $         
     $         EW_overlap,
     $         NW_overlap,
     $         N_overlap,
     $         NE_overlap,
     $         EW_overlap,
     $         
     $         EW_overlap,
     $         W_overlap,
     $         no_overlap,
     $         E_overlap,
     $         EW_overlap,
     $         
     $         EW_overlap,
     $         SW_overlap,
     $         S_overlap,
     $         SE_overlap,
     $         EW_overlap,
     $         
     $         NS_overlap,
     $         NS_overlap,
     $         NS_overlap,
     $         NS_overlap,
     $         NS_overlap]

          test_square_borders = reshape((/
     $         0, 20,  1, 21,
     $         1, 20,  2, 21,
     $         5, 20,  6, 21,
     $         9, 20, 10, 21,
     $        10, 20, 11, 21,
     $         
     $         0, 19,  1, 20,
     $         1, 19,  2, 20,
     $         5, 19,  6, 20,
     $         9, 19, 10, 20,
     $        10, 19, 11, 20,
     $         
     $         0, 10,  1, 15,
     $         1, 10,  2, 15,
     $         5, 10,  6, 15,
     $         9, 10, 10, 15,
     $        10, 10, 11, 15,
     $         
     $         0,  0,  1,  1,
     $         1,  0,  2,  1,
     $         5,  0,  6,  1,
     $         9,  0, 10,  1,
     $        10,  0, 11,  1,
     $         
     $         0, -1,  1,  0,
     $         1, -1,  2,  0,
     $         5, -1,  6,  0,
     $         9, -1, 10,  0,
     $        10, -1, 11,  0/),
     $         (/2,2,25/))

          test_square_borders_after = reshape((/
     $         0, 20,  1, 21,
     $         1, 20,  2, 21,
     $         5, 20,  6, 21,
     $         9, 20, 10, 21,
     $        10, 20, 11, 21,
     $         
     $         0, 19,  1, 20,
     $         1, 19,  2, 20,
     $         5, 19,  6, 20,
     $         9, 19, 10, 20,
     $        10, 19, 11, 20,
     $         
     $         0, 10,  1, 15,
     $         1, 10,  2, 15,
     $         5, 10,  6, 15,
     $         9, 10, 10, 15,
     $        10, 10, 11, 15,
     $         
     $         0,  0,  1,  1,
     $         1,  0,  2,  1,
     $         5,  0,  6,  1,
     $         9,  0, 10,  1,
     $        10,  0, 11,  1,
     $         
     $         0, -1,  1,  0,
     $         1, -1,  2,  0,
     $         5, -1,  6,  0,
     $         9, -1, 10,  0,
     $        10, -1, 11,  0/),
     $         (/2,2,25/))

          !NE_corner test
          test_loc = test_overlap_bc_section(
     $         NE_corner_type,
     $         test_square_borders,
     $         test_square_borders_after,
     $         test_square_overlap,
     $         test_square_overlap_after,
     $         gen_borders,
     $         detailled)
          test_validated = test_validated.and.test_loc

          !NW_corner test
          test_loc = test_overlap_bc_section(
     $         NW_corner_type,
     $         test_square_borders,
     $         test_square_borders_after,
     $         test_square_overlap,
     $         test_square_overlap_after,
     $         gen_borders,
     $         detailled)
          test_validated = test_validated.and.test_loc

          !SE_edge test
          test_loc = test_overlap_bc_section(
     $         SE_edge_type,
     $         test_square_borders,
     $         test_square_borders_after,
     $         test_square_overlap,
     $         test_square_overlap_after,
     $         gen_borders,
     $         detailled)
          test_validated = test_validated.and.test_loc

          !SW_edge test
          test_loc = test_overlap_bc_section(
     $         SW_edge_type,
     $         test_square_borders,
     $         test_square_borders_after,
     $         test_square_overlap,
     $         test_square_overlap_after,
     $         gen_borders,
     $         detailled)
          test_validated = test_validated.and.test_loc

          !NE_corner test
          test_loc = test_overlap_bc_section(
     $         NE_corner_type,
     $         test_square_borders,
     $         test_square_borders_after,
     $         test_square_overlap,
     $         test_square_overlap_after,
     $         gen_borders,
     $         detailled)
          test_validated = test_validated.and.test_loc

          !NW_corner test
          test_loc = test_overlap_bc_section(
     $         NW_corner_type,
     $         test_square_borders,
     $         test_square_borders_after,
     $         test_square_overlap,
     $         test_square_overlap_after,
     $         gen_borders,
     $         detailled)
          test_validated = test_validated.and.test_loc

          !SE_edge test
          test_loc = test_overlap_bc_section(
     $         SE_edge_type,
     $         test_square_borders,
     $         test_square_borders_after,
     $         test_square_overlap,
     $         test_square_overlap_after,
     $         gen_borders,
     $         detailled)
          test_validated = test_validated.and.test_loc

          !SW_edge test
          test_loc = test_overlap_bc_section(
     $         SW_edge_type,
     $         test_square_borders,
     $         test_square_borders_after,
     $         test_square_overlap,
     $         test_square_overlap_after,
     $         gen_borders,
     $         detailled)
          test_validated = test_validated.and.test_loc
          
        end function test_overlap_bc_section_by_integration_borders


        function test_overlap_bc_section(
     $     bc_section_type,
     $     test_bc_section_borders,
     $     test_bc_section_borders_after,
     $     test_bc_section_overlap,
     $     test_bc_section_overlap_after,
     $     gen_borders,
     $     detailled)
     $     result(test_validated)

          implicit none

          integer                              , intent(in)    :: bc_section_type
          integer       , dimension(:,:,:)     , intent(inout) :: test_bc_section_borders
          integer       , dimension(:,:,:)     , intent(in)    :: test_bc_section_borders_after
          integer       , dimension(:)         , intent(inout) :: test_bc_section_overlap
          integer       , dimension(:)         , intent(in)    :: test_bc_section_overlap_after
          integer(ikind), dimension(2,2)       , intent(in)    :: gen_borders
          logical                              , intent(in)    :: detailled
          logical                                              :: test_validated

          
          integer :: k

          test_validated = .true.

          
          do k=1,size(test_bc_section_borders,3)

             !output
             call overlap_bc_section_by_integration_borders(
     $            bc_section_type,
     $            test_bc_section_borders(:,:,k),
     $            test_bc_section_overlap(k),
     $            gen_borders)

             !validation
             test_loc = is_int_matrix_validated(
     $            test_bc_section_borders(:,:,k),
     $            test_bc_section_borders_after(:,:,k),
     $            detailled)

             test_loc = test_loc.and.
     $            (test_bc_section_overlap(k).eq.test_bc_section_overlap_after(k))

             test_validated = test_validated.and.test_loc

             !detailled
             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
                print '(''   - overlap: '',I2,'' -> '',I2)',
     $               test_bc_section_overlap(k),
     $               test_bc_section_overlap_after(k)
                print '()'
             end if

          end do

        end function test_overlap_bc_section

      end program test_bf_bc_sections_overlap
