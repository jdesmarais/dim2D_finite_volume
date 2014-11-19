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

        use bf_layer_class, only :
     $       bf_layer

        use bf_layer_newgrdpt_procedure_module, only :
     $       no_gradient_type,
     $       gradient_I_type,
     $       gradient_L0_type,
     $       gradient_R0_type,
     $       get_newgrdpt_procedure,
     $       get_interior_data_for_newgrdpt,
     $       are_intermediate_newgrdpt_data_needed,
     $       get_x_map_for_newgrdpt,
     $       get_y_map_for_newgrdpt

        use parameters_bf_layer, only :
     $       no_pt,
     $       bc_pt,
     $       bc_interior_pt,
     $       interior_pt

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne,bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        
        logical :: detailled
        logical :: test_validated

        
        detailled = .true.

        test_validated = test_get_newgrdpt_procedure(detailled)
        print '(''test_get_newgrdpt_procedure: '',L1)', test_validated
        print '()'

        test_validated = test_get_interior_data_for_newgrdpt(detailled)
        print '(''test_get_interior_data_for_newgrdpt: '', L1)', test_validated
        print '()'

        test_validated = test_are_intermediate_newgrdpt_data_needed(detailled)
        print '(''test_are_intermediate_newgrdpt_data_needed: '', L1)', test_validated
        print '()'
        
        test_validated = test_get_buffer_data_for_newgrdpt(detailled)
        print '(''test_get_buffer_data_for_newgrdpt: '', L1)', test_validated
        print '()'

        test_validated = test_get_x_map_for_newgrdpt(detailled)
        print '(''test_get_x_map_for_newgrdpt: '', L1)', test_validated
        print '()'

        test_validated = test_get_y_map_for_newgrdpt(detailled)
        print '(''test_get_y_map_for_newgrdpt: '', L1)', test_validated
        print '()'

        contains

        function test_get_buffer_data_for_newgrdpt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer       , dimension(:,:)     , allocatable      :: alignment0
          integer       , dimension(:,:)     , allocatable      :: alignment1
          integer       , dimension(:,:)     , allocatable      :: grdpts_id0
          integer       , dimension(:,:)     , allocatable      :: grdpts_id1
          real(rkind)   , dimension(:,:,:)   , allocatable      :: interior_nodes0
          real(rkind)   , dimension(:,:,:)   , allocatable      :: interior_nodes1
          integer       , dimension(2*bc_size+1,2*bc_size+1)    :: test_tmp_grdpts_id0
          real(rkind)   , dimension(2*bc_size+1,2*bc_size+1,ne) :: test_tmp_nodes0
          real(rkind)   , dimension(2*bc_size+1,2*bc_size+1,ne) :: test_tmp_nodes1
          integer       , dimension(2*bc_size+1,2*bc_size+1)    :: tmp_grdpts_id0
          real(rkind)   , dimension(2*bc_size+1,2*bc_size+1,ne) :: tmp_nodes0
          real(rkind)   , dimension(2*bc_size+1,2*bc_size+1,ne) :: tmp_nodes1
          integer(ikind), dimension(2,2)                        :: test_gen_coords
          integer(ikind), dimension(2,2)                        :: test_indices


          type(bf_layer) :: bf_layer_used
          integer        :: k
          integer        :: nb_tests
          logical        :: test_loc
          
          nb_tests = 4

          !allocations
          allocate(alignment0(2,2))
          allocate(alignment1(2,2))
          allocate(grdpts_id0(nx,ny))
          allocate(grdpts_id1(nx,ny))
          allocate(interior_nodes0(nx,ny,ne))
          allocate(interior_nodes1(nx,ny,ne))
          

          !initialize the nodes for the buffer layer
          alignment0(1,1) = bc_size+1
          alignment0(1,2) = nx-bc_size
          alignment0(2,1) = bc_size+1
          alignment0(2,2) = ny-bc_size

          alignment1(1,1) = bc_size+1
          alignment1(1,2) = nx-bc_size
          alignment1(2,1) = bc_size+1
          alignment1(2,2) = ny-bc_size

          call initialize_grdpts_id(grdpts_id0)
          call initialize_nodes(interior_nodes0,interior_nodes1)

          call bf_layer_used%set_alignment_tab(alignment1)
          call bf_layer_used%set_grdpts_id(grdpts_id1)
          call bf_layer_used%set_nodes(interior_nodes1)

          call bf_layer_used%bf_compute_used%set_alignment(alignment0)
          call bf_layer_used%bf_compute_used%set_grdpts_id(grdpts_id0)
          call bf_layer_used%bf_compute_used%set_nodes(interior_nodes0)


          test_validated = .true.

          !tests different gen_coords
          do k=1, nb_tests

             !get the data for the test
             call get_data_get_interior_data_for_newgrdpt(
     $            k,
     $            test_tmp_grdpts_id0,
     $            test_tmp_nodes0,
     $            test_tmp_nodes1,
     $            test_gen_coords,
     $            test_indices)
             
             !initialize the tmp_grdpts_id0 with no_pt
             tmp_grdpts_id0 = reshape((/
     $            no_pt,no_pt,no_pt,no_pt,no_pt,
     $            no_pt,no_pt,no_pt,no_pt,no_pt,
     $            no_pt,no_pt,no_pt,no_pt,no_pt,
     $            no_pt,no_pt,no_pt,no_pt,no_pt,
     $            no_pt,no_pt,no_pt,no_pt,no_pt/),
     $            (/5,5/))

             !get the data according to the fct tested
             call bf_layer_used%get_data_for_newgrdpt(
     $            tmp_grdpts_id0,
     $            tmp_nodes0,
     $            tmp_nodes1,
     $            test_gen_coords)

             !compare the data
             test_loc = compare_get_interior_data_for_newgrdpt(
     $            test_indices,
     $            tmp_grdpts_id0,
     $            tmp_nodes0,
     $            tmp_nodes1,
     $            test_tmp_grdpts_id0,
     $            test_tmp_nodes0,
     $            test_tmp_nodes1,
     $            detailled)

             test_validated = test_validated.and.test_loc

             !print the test
             if(detailled) then
                if(test_loc) then

                   print '(''test '',I2, '' validated'')', k

                else

                   print '(''*** test '',I2, '' failed ***'')', k

                end if
             end if

          end do

        end function test_get_buffer_data_for_newgrdpt



        function test_are_intermediate_newgrdpt_data_needed(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer :: k
          logical :: test_loc
          integer :: nb_tests

          integer, dimension(2,2) :: test_gen_coords
          integer                 :: test_mainlayer_id
          logical                 :: test_needed
          logical                 :: needed

          nb_tests = 12

          do k=1,nb_tests

             call get_data_are_intermediate_newgrdpt_data_needed(
     $            k,
     $            test_gen_coords,
     $            test_mainlayer_id,
     $            test_needed)

             needed = are_intermediate_newgrdpt_data_needed(
     $            test_mainlayer_id, test_gen_coords)

             test_loc = needed.eqv.test_needed
             test_validated=test_loc.and.test_validated

             if(detailled) then
                if(test_loc) then
                   print '(''test '', I2, '' validated'')', k
                else
                   print '(''**test '', I2, ''failed**'')', k
                   print '(L1,'' -> '',L1)', needed, test_needed
                end if
             end if

          end do

        end function test_are_intermediate_newgrdpt_data_needed


        subroutine get_data_are_intermediate_newgrdpt_data_needed(
     $     config,
     $     test_gen_coords,
     $     test_mainlayer_id,
     $     test_needed)

          implicit none

          integer                       , intent(in)  :: config
          integer(ikind), dimension(2,2), intent(out) :: test_gen_coords
          integer                       , intent(out) :: test_mainlayer_id
          logical                       , intent(out) :: test_needed


          select case(config)
            case(1)
               test_mainlayer_id    = N
               test_gen_coords(2,1) = ny+1
               test_needed          = .false.

            case(2)
               test_mainlayer_id    = N
               test_gen_coords(2,1) = ny
               test_needed          = .true.

            case(3)
               test_mainlayer_id    = S
               test_gen_coords(2,2) = 0
               test_needed          = .false.

            case(4)
               test_mainlayer_id    = S
               test_gen_coords(2,2) = 1
               test_needed          = .true.

            case(5)
               test_mainlayer_id    = E
               test_gen_coords(1,1) = nx+1
               test_gen_coords(2,1) = bc_size+1
               test_gen_coords(2,2) = ny-bc_size
               test_needed          = .false.

            case(6)
               test_mainlayer_id    = E
               test_gen_coords(1,1) = nx+1
               test_gen_coords(2,1) = bc_size
               test_gen_coords(2,2) = ny-bc_size
               test_needed          = .true.

            case(7)
               test_mainlayer_id    = E
               test_gen_coords(1,1) = nx+1
               test_gen_coords(2,1) = bc_size+1
               test_gen_coords(2,2) = ny
               test_needed          = .true.

            case(8)
               test_mainlayer_id    = E
               test_gen_coords(1,1) = nx
               test_gen_coords(2,1) = bc_size+1
               test_gen_coords(2,2) = ny-bc_size
               test_needed          = .true.

            case(9)
               test_mainlayer_id    = W
               test_gen_coords(1,2) = 0
               test_gen_coords(2,1) = bc_size+1
               test_gen_coords(2,2) = ny-bc_size
               test_needed          = .false.

            case(10)
               test_mainlayer_id    = W
               test_gen_coords(1,2) = 0
               test_gen_coords(2,1) = bc_size
               test_gen_coords(2,2) = ny-bc_size
               test_needed          = .true.

            case(11)
               test_mainlayer_id    = W
               test_gen_coords(1,2) = 0
               test_gen_coords(2,1) = bc_size+1
               test_gen_coords(2,2) = ny
               test_needed          = .true.

            case(12)
               test_mainlayer_id    = W
               test_gen_coords(1,2) = 1
               test_gen_coords(2,1) = bc_size+1
               test_gen_coords(2,2) = ny-bc_size
               test_needed          = .true.

            case default
              print '(''test_bf_layer_newgrdpt_procedure'')'
           print '(''get_data_are_intermediate_newgrdpt_data_needed'')'
              print '(''test not implemented: '',I2)', config
              stop ''

          end select

        end subroutine get_data_are_intermediate_newgrdpt_data_needed



        function test_get_interior_data_for_newgrdpt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind)    , dimension(nx,ny,ne)                   :: interior_nodes0
          real(rkind)    , dimension(nx,ny,ne)                   :: interior_nodes1
          integer        , dimension(2*bc_size+1,2*bc_size+1)    :: test_tmp_grdpts_id0
          real(rkind)    , dimension(2*bc_size+1,2*bc_size+1,ne) :: test_tmp_nodes0
          real(rkind)    , dimension(2*bc_size+1,2*bc_size+1,ne) :: test_tmp_nodes1
          integer        , dimension(2*bc_size+1,2*bc_size+1)    :: tmp_grdpts_id0
          real(rkind)    , dimension(2*bc_size+1,2*bc_size+1,ne) :: tmp_nodes0
          real(rkind)    , dimension(2*bc_size+1,2*bc_size+1,ne) :: tmp_nodes1
          integer(ikind) , dimension(2,2)                        :: test_gen_coords
          integer(ikind) , dimension(2,2)                        :: test_indices

          integer :: k
          integer :: nb_tests
          logical :: test_loc
          
          nb_tests = 4

          !initialize the nodes for interior_nodes0 and interior_nodes1
          call initialize_nodes(interior_nodes0,interior_nodes1)

          test_validated = .true.

          !tests different gen_coords
          do k=1, nb_tests

             !get the data for the test
             call get_data_get_interior_data_for_newgrdpt(
     $            k,
     $            test_tmp_grdpts_id0,
     $            test_tmp_nodes0,
     $            test_tmp_nodes1,
     $            test_gen_coords,
     $            test_indices)
             
             !get the data according to the fct tested
             call get_interior_data_for_newgrdpt(
     $            interior_nodes0,
     $            interior_nodes1,
     $            tmp_grdpts_id0,
     $            tmp_nodes0,
     $            tmp_nodes1,
     $            test_gen_coords)

             !compare the data
             test_loc = compare_get_interior_data_for_newgrdpt(
     $            test_indices,
     $            tmp_grdpts_id0,
     $            tmp_nodes0,
     $            tmp_nodes1,
     $            test_tmp_grdpts_id0,
     $            test_tmp_nodes0,
     $            test_tmp_nodes1,
     $            detailled)

             test_validated = test_validated.and.test_loc

             !print the test
             if(detailled) then
                if(test_loc) then

                   print '(''test '',I2, '' validated'')', k

                else

                   print '(''*** test '',I2, '' failed ***'')', k

                end if
             end if

          end do

        end function test_get_interior_data_for_newgrdpt

      
        function compare_get_interior_data_for_newgrdpt(
     $     test_indices,
     $     tmp_grdpts_id0,
     $     tmp_nodes0,
     $     tmp_nodes1,
     $     test_tmp_grdpts_id0,
     $     test_tmp_nodes0,
     $     test_tmp_nodes1,
     $     detailled)
     $     result(test_validated)

          integer(ikind), dimension(2,2)   , intent(in) :: test_indices
          integer       , dimension(5,5)   , intent(in) :: tmp_grdpts_id0
          real(rkind)   , dimension(5,5,ne), intent(in) :: tmp_nodes0
          real(rkind)   , dimension(5,5,ne), intent(in) :: tmp_nodes1
          integer       , dimension(5,5)   , intent(in) :: test_tmp_grdpts_id0
          real(rkind)   , dimension(5,5,ne), intent(in) :: test_tmp_nodes0
          real(rkind)   , dimension(5,5,ne), intent(in) :: test_tmp_nodes1
          logical                          , intent(in) :: detailled
          logical                                       :: test_validated


          logical        :: test_loc
          logical        :: test_grdpts_id
          logical        :: test_nodes0
          logical        :: test_nodes1
          integer(ikind) :: i,j
          integer        :: k
          
          
          test_validated = .true.

          !test grdpts_id
          do j=1,5
             do i=1,5
                test_loc = tmp_grdpts_id0(i,j).eq.test_tmp_grdpts_id0(i,j)
                test_validated = test_loc.and.test_validated
             end do
          end do

          test_grdpts_id = test_validated
          if(detailled) then
             print '(''-- test grdpts_id: '',L1)', test_grdpts_id

             if(.not.test_validated) then

                print '(''  - test_tmp_grdpts_id0:'')'
                do j=1,5
                   print '(''    '',5I2)', test_tmp_grdpts_id0(:,5-j+1)
                end do
                print '(''  - tmp_grdpts_id0:'')'
                do j=1,5
                   print '(''    '',5I2)', tmp_grdpts_id0(:,5-j+1)
                end do
                print ''

             end if
          end if


          test_validated = .true.

          !test nodes0
          do k=1,ne
             do j=test_indices(2,1),test_indices(2,2)
                do i=test_indices(1,1),test_indices(1,2)
                   test_loc = is_test_validated(
     $                  tmp_nodes0(i,j,k),
     $                  test_tmp_nodes0(i,j,k),
     $                  .false.)
                   test_validated = test_loc.and.test_validated
                   if(.not.test_loc) then
                      print '(F5.2,F5.2,3I2)',
     $                     tmp_nodes0(i,j,k),
     $                     test_tmp_nodes0(i,j,k),
     $                     i,j,k
                   end if
                end do
             end do
          end do

          test_nodes0 = test_validated
          if(detailled) then
             print '(''-- test nodes0: '',L1)', test_nodes0

             if(.not.test_validated) then

                print '(''  - test_tmp_nodes0:'')'
                do j=1,5
                   print '(''    '',5F6.1)', test_tmp_nodes0(:,5-j+1,1)
                end do
                print '(''  - tmp_nodes0:'')'
                do j=1,5
                   print '(''    '',5F6.1)', tmp_nodes0(:,5-j+1,1)
                end do
                print ''

             end if

          end if          

             
          test_validated = .true.

          !test nodes1
          do k=1,ne
             do j=test_indices(2,1),test_indices(2,2)
                do i=test_indices(1,1),test_indices(1,2)
                   test_loc = is_test_validated(
     $                  tmp_nodes1(i,j,k),
     $                  test_tmp_nodes1(i,j,k),
     $                  .false.)
                   test_validated = test_loc.and.test_validated
                   if(.not.test_loc) then
                      print '(F5.2,F5.2,3I2)',
     $                     tmp_nodes1(i,j,k),
     $                     test_tmp_nodes1(i,j,k),
     $                     i,j,k
                   end if
                end do
             end do
          end do

          test_nodes1 = test_validated
          if(detailled) then
             print '(''-- test nodes1: '',L1)', test_nodes1

             if(.not.test_validated) then

                print '(''  - test_tmp_nodes1:'')'
                do j=1,5
                   print '(''    '',5F6.1)', test_tmp_nodes1(:,5-j+1,1)
                end do
                print '(''  - tmp_nodes1:'')'
                do j=1,5
                   print '(''    '',5F6.1)', tmp_nodes1(:,5-j+1,1)
                end do
                print ''

             end if

          end if


          test_validated = test_grdpts_id.and.test_nodes0.and.test_nodes1

        end function compare_get_interior_data_for_newgrdpt


        subroutine get_data_get_interior_data_for_newgrdpt(
     $     config,
     $     test_tmp_grdpts_id0,
     $     test_tmp_nodes0,
     $     test_tmp_nodes1,
     $     test_gen_coords,
     $     test_indices)

          implicit none

          integer                                              , intent(in)  :: config
          integer       , dimension(2*bc_size+1,2*bc_size+1)   , intent(out) :: test_tmp_grdpts_id0
          real(rkind)   , dimension(2*bc_size+1,2*bc_size+1,ne), intent(out) :: test_tmp_nodes0
          real(rkind)   , dimension(2*bc_size+1,2*bc_size+1,ne), intent(out) :: test_tmp_nodes1
          integer(ikind), dimension(2,2)                       , intent(out) :: test_gen_coords
          integer(ikind), dimension(2,2)                       , intent(out) :: test_indices
            
          
          integer(ikind) :: i,j
          integer        :: k

          
          !---------------------------------------------
          !                       ___2____
          !      ___1____        |        |
          !     |        |       |        |
          !     |     ___|_______|___     |
          !     |    |   |       |___|____|
          !     |____|___|           |     
          !          |               |
          !          |               |   ___3____
          !          |               |  |        |
          !          |    ____4___   |  |        |
          !          |___|________|__|  |        |
          !              |        |     |________|
          !              |        |
          !              |________|
          !
          !---------------------------------------------          
          select case(config)

            case(1)

               test_gen_coords(1,1) = -1
               test_gen_coords(1,2) =  3
               test_gen_coords(2,1) =  8
               test_gen_coords(2,2) = 12

               test_indices(1,1) = 3
               test_indices(1,2) = 5
               test_indices(2,1) = 1
               test_indices(2,2) = 3

               test_tmp_grdpts_id0 = reshape((/
     $              no_pt,no_pt,bc_pt,bc_interior_pt,interior_pt,
     $              no_pt,no_pt,bc_pt,bc_interior_pt,bc_interior_pt,
     $              no_pt,no_pt,bc_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt/),
     $              (/5,5/))

            case(2)

               test_gen_coords(1,1) =  8
               test_gen_coords(1,2) =  12
               test_gen_coords(2,1) =  9
               test_gen_coords(2,2) = 13

               test_indices(1,1) = 1
               test_indices(1,2) = 3
               test_indices(2,1) = 1
               test_indices(2,2) = 2

               test_tmp_grdpts_id0 = reshape((/
     $              bc_interior_pt,bc_interior_pt,bc_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,bc_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt/),
     $              (/5,5/))

            case(3)
               
               test_gen_coords(1,1) =  12
               test_gen_coords(1,2) =  16
               test_gen_coords(2,1) =  3
               test_gen_coords(2,2) =  7
               
               test_indices(1,1) = 0
               test_indices(1,2) = -1
               test_indices(2,1) = 0
               test_indices(2,2) = -1
               
               test_tmp_grdpts_id0 = reshape((/
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt/),
     $              (/5,5/))

            case(4)
               
               test_gen_coords(1,1) =   4
               test_gen_coords(1,2) =   8
               test_gen_coords(2,1) =  -1
               test_gen_coords(2,2) =   3
               
               test_indices(1,1) = 1
               test_indices(1,2) = 5
               test_indices(2,1) = 3
               test_indices(2,2) = 5
               
               test_tmp_grdpts_id0 = reshape((/
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,bc_pt,bc_pt,bc_pt,
     $              bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,bc_interior_pt,
     $              interior_pt,interior_pt,interior_pt,interior_pt,interior_pt/),
     $              (/5,5/))

            case default

               print '(''test not implemented: '', I2)', k
               stop ''

          end select

          do k=1,ne
             do j=test_gen_coords(2,1), test_gen_coords(2,2)
                do i=test_gen_coords(1,1), test_gen_coords(1,2)

                   test_tmp_nodes0(
     $                  i-test_gen_coords(1,1)+1,
     $                  j-test_gen_coords(2,1)+1,
     $                  k) =
     $                  + (i-1) + 10*(j-1) + 100*(k-1)

                   test_tmp_nodes1(
     $                  i-test_gen_coords(1,1)+1,
     $                  j-test_gen_coords(2,1)+1,
     $                  k) =
     $                  -(i-1) - 10*(j-1) - 100*(k-1)

                end do
             end do
          end do

        end subroutine get_data_get_interior_data_for_newgrdpt


        subroutine initialize_grdpts_id(grdpts_id)

          implicit none

          integer, dimension(nx,ny), intent(out) :: grdpts_id

          integer(ikind) :: i,j

          j=1
          do i=1,nx
             grdpts_id(i,j) = bc_pt
          end do

          j=2
          i=1
          grdpts_id(i,j) = bc_pt
          do i=2, nx-1
             grdpts_id(i,j) = bc_interior_pt
          end do
          i=nx
          grdpts_id(i,j) = bc_pt

          do j=3,ny-2

             i=1
             grdpts_id(i,j) = bc_pt

             i=2
             grdpts_id(i,j) = bc_interior_pt

             do i=3,nx-2
                grdpts_id(i,j) = interior_pt
             end do

             i=nx-1
             grdpts_id(i,j) = bc_interior_pt

             i=nx
             grdpts_id(i,j) = bc_pt

          end do

          j=ny-1
          i=1
          grdpts_id(i,j) = bc_pt
          do i=2, nx-1
             grdpts_id(i,j) = bc_interior_pt
          end do
          i=nx
          grdpts_id(i,j) = bc_pt

          j=ny
          do i=1,nx
             grdpts_id(i,j) = bc_pt
          end do

        end subroutine initialize_grdpts_id


        subroutine initialize_nodes(interior_nodes0,interior_nodes1)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes0
          real(rkind), dimension(nx,ny,ne), intent(out) :: interior_nodes1

          integer(ikind) :: i,j
          integer        :: k


          if((nx.ne.10).and.(ny.ne.10)) then
             print '(''test_bf_layer_newgrdpt_procedure'')'
             print '(''initialize_nodes'')'
             print '(''the test requires (nx,ny)=(10,10)'')'
             stop ''
          end if

          do k=1,ne
             do j=1, ny
                do i=1, nx
                   interior_nodes0(i,j,k) = (i-1) + 10*(j-1) + 100*(k-1)
                   interior_nodes1(i,j,k) =-(i-1) - 10*(j-1) - 100*(k-1)
                end do
             end do
          end do

        end subroutine initialize_nodes


        function test_get_newgrdpt_procedure(detailled)
     $       result(test_validated)

          logical, intent(in) :: detailled
          logical             :: test_validated


          integer                              :: k
          integer                              :: nb_tests

          integer, dimension(:,:), allocatable :: grdpts_id_tested
          integer(ikind)                       :: i_tested
          integer(ikind)                       :: j_tested
          integer                              :: procedure_type_data
          integer                              :: gradient_type_data

          integer                              :: procedure_type
          integer                              :: gradient_type

          logical                              :: test_loc


          nb_tests = 22
          test_validated = .true.

          do k=1, nb_tests

             !get the data to conduct the test
             call get_test_data(
     $            k,
     $            grdpts_id_tested,
     $            i_tested,
     $            j_tested,
     $            procedure_type_data,
     $            gradient_type_data)

             !carry the test
             call get_newgrdpt_procedure(
     $            i_tested,
     $            j_tested,
     $            grdpts_id_tested,
     $            procedure_type,
     $            gradient_type)
             
             !compare the results
             test_loc = (procedure_type.eq.procedure_type_data).and.
     $            (gradient_type.eq.gradient_type_data)
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

        end function test_get_newgrdpt_procedure


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
            ! | 0 0 0 |
            ! | 0 0*0 |
            ! | 3 3 3 |
            !  -------
            case(1)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              bc_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = N_edge_type
               gradient_type_data  = gradient_I_type

            !  -------
            ! | 0 0 0 | 
            ! | 3 0*0 |
            ! | 3 3 0 |
            !  ------- 
            case(2)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              bc_pt,bc_pt,no_pt,
     $              bc_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = NE_edge_type
               gradient_type_data  = no_gradient_type

            !  -------
            ! | 0 0 0 |
            ! | 0 0*0 |
            ! | 3 3 0 |
            !  -------
            case(3)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              bc_pt,bc_pt,no_pt,
     $              no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = N_edge_type
               gradient_type_data  = gradient_R0_type

            !  -------
            ! | 3 0 0 |
            ! | 3 0*0 |
            ! | 3 0 0 |
            !  -------
            case(4)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              bc_pt,no_pt,no_pt,
     $              bc_pt,no_pt,no_pt,
     $              bc_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = E_edge_type
               gradient_type_data  = gradient_I_type

            !  -------
            ! | 0 0 0 |
            ! | 3 0*0 |
            ! | 3 0 0 |
            !  -------
            case(5)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              bc_pt,no_pt,no_pt,
     $              bc_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = E_edge_type
               gradient_type_data  = gradient_R0_type

            !  ---------
            ! | 0 0 0 0 |
            ! | 0 0 0*0 |
            ! | 3 3 0 0 |
            ! | 2 3 0 0 |
            !  ---------
            case(6)

               allocate(grdpts_id_tested(4,4))
               grdpts_id_tested = reshape( (/
     $              bc_interior_pt,bc_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,no_pt
     $              /),
     $              (/4,4/))

               i_tested = 3
               j_tested = 3

               procedure_type_data = NE_corner_type
               gradient_type_data  = no_gradient_type

            !  -------
            ! | 0 0 0 |
            ! | 0 0*3 |
            ! | 0 3 3 |
            !  -------
            case(7)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,bc_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = NW_edge_type
               gradient_type_data  = no_gradient_type

            !  -------
            ! | 0 0 0 |
            ! ! 0 0*0 |
            ! | 0 3 3 |
            !  -------
            case(8)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,bc_pt,bc_pt,
     $              no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
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

            !  -------
            ! | 0 0 0 |
            ! | 0 0*3 |
            ! | 0 0 3 |
            !  -------
            case(11)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,bc_pt,
     $              no_pt,no_pt,bc_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = W_edge_type
               gradient_type_data  = gradient_R0_type


            !  -------
            ! | 0 0 0 |
            ! | 0 0*0 |
            ! | 0 0 3 |
            !  -------
            case(12)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,bc_pt,
     $              no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = NW_corner_type
               gradient_type_data  = no_gradient_type

            !  -------
            ! | 3 3 0 |
            ! | 3 0*0 |
            ! | 0 0 0 |
            !  -------
            case(13)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,no_pt,
     $              bc_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = SE_edge_type
               gradient_type_data  = no_gradient_type

            !  -------
            ! | 3 0 0 |
            ! | 3 0*0 |
            ! | 0 0 0 |
            !  -------
            case(14)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,no_pt,
     $              bc_pt,no_pt,no_pt,
     $              bc_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = E_edge_type
               gradient_type_data  = gradient_L0_type

            !  -------
            ! | 0 3 3 |
            ! | 0 0*3 |
            ! | 0 0 0 |
            !  -------
            case(15)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,no_pt,
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
            ! | 0 0 0 |
            !  -------
            case(16)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,no_pt,
     $              no_pt,no_pt,bc_pt,
     $              no_pt,no_pt,bc_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = W_edge_type
               gradient_type_data  = gradient_L0_type

            !  -------
            ! | 3 3 3 |
            ! | 0 0*0 |
            ! | 0 0 0 |
            !  -------
            case(17)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = S_edge_type
               gradient_type_data  = gradient_I_type

            !  -------
            ! | 3 3 0 |
            ! | 0 0*0 |
            ! | 0 0 0 |
            !  -------
            case(18)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,
     $              bc_pt,bc_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = S_edge_type
               gradient_type_data  = gradient_R0_type

            !  -------
            ! | 3 0 0 |
            ! | 0 0*0 |
            ! | 0 0 0 |
            !  -------
            case(19)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,
     $              bc_pt,no_pt,no_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

               procedure_type_data = SE_corner_type
               gradient_type_data  = no_gradient_type

            !  -------
            ! | 0 0 3 |
            ! | 0 0*0 |
            ! | 0 0 0 |
            !  -------
            case(20)

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

            !  -------
            ! | 0 3 3 |
            ! | 0 0*0 |
            ! | 0 0 0 |
            !  -------
            case(21)

               allocate(grdpts_id_tested(3,3))
               grdpts_id_tested = reshape( (/
     $              no_pt,no_pt,no_pt,
     $              no_pt,no_pt,no_pt,
     $              no_pt,bc_pt,bc_pt
     $              /),
     $              (/3,3/))

               i_tested = 2
               j_tested = 2

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


        function is_test_validated(var,cst,detailled) result(test_validated)

          implicit none

          real(rkind), intent(in) :: var
          real(rkind), intent(in) :: cst
          logical    , intent(in) :: detailled
          logical                 :: test_validated

          if(detailled) then
             print *, int(var*1e5)
             print *, int(cst*1e5)
          end if
          
          test_validated=abs(
     $         int(var*1e5)-
     $         int(cst*1e5)).le.1
          
        end function is_test_validated


        function test_get_x_map_for_newgrdpt(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind)   , dimension(nx)          :: interior_x_map
          real(rkind)                            :: dx
          integer                                :: nb_tests
          integer                                :: k
          real(rkind)   , dimension(2*bc_size+1) :: tmp_x_map
          real(rkind)   , dimension(2*bc_size+1) :: tmp_x_map_data
          integer(ikind), dimension(2,2)         :: gen_coords

          logical :: test_loc

          
          dx = 0.1d0
          call initialize_map(interior_x_map,dx)

          nb_tests = 5

          do k=1, nb_tests

             !get the data for the test
             call get_data_test_get_map_for_newgrdpt(
     $            k,dx,
     $            tmp_x_map_data,
     $            gen_coords)

             !compute the temporary x_map
             tmp_x_map = get_x_map_for_newgrdpt(
     $            interior_x_map,
     $            gen_coords)

             !compare the results
             test_loc = compare_tmp_maps(
     $            tmp_x_map,
     $            tmp_x_map_data,
     $            .false.)
             test_validated = test_validated.and.test_loc

             !display the results of the comparison
             if(detailled) then
                if(test_loc) then
                   print '(''test '',I2, '' validated'')', k
                else
                   print '(''** test '',I2, '' failed **'')', k
                end if
             end if

          end do

        end function test_get_x_map_for_newgrdpt


        function test_get_y_map_for_newgrdpt(detailled)
     $     result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(ny)          :: interior_y_map
          real(rkind)                         :: dy
          integer                             :: nb_tests
          integer                             :: k
          real(rkind), dimension(2*bc_size+1) :: tmp_y_map
          real(rkind), dimension(2*bc_size+1) :: tmp_y_map_data
          integer(ikind), dimension(2,2)      :: gen_coords

          logical :: test_loc

          
          dy = 0.2d0
          call initialize_map(interior_y_map,dy)
          
          nb_tests = 5

          do k=1, nb_tests

             !get the data for the test
             call get_data_test_get_map_for_newgrdpt(
     $            k,dy,
     $            tmp_y_map_data,
     $            gen_coords)

             !compute the temporary x_map
             tmp_y_map = get_y_map_for_newgrdpt(
     $            interior_y_map,
     $            gen_coords)

             !compare the results
             test_loc = compare_tmp_maps(
     $            tmp_y_map,
     $            tmp_y_map_data,
     $            .false.)
             test_validated = test_validated.and.test_loc

             !display the results of the comparison
             if(detailled) then
                if(test_loc) then
                   print '(''test '',I2, '' validated'')', k
                else
                   print '(''** test '',I2, '' failed **'')', k
                end if
             end if

          end do

        end function test_get_y_map_for_newgrdpt


        subroutine get_data_test_get_map_for_newgrdpt(
     $     config,
     $     space_step,
     $     tmp_map_data,
     $     gen_coords)

          implicit none

          integer                            , intent(in)  :: config
          real(rkind)                        , intent(in)  :: space_step
          real(rkind), dimension(2*bc_size+1), intent(out) :: tmp_map_data
          integer(ikind), dimension(2,2)     , intent(out) :: gen_coords

          integer :: k

          select case(config)
            case(1)

               gen_coords(1,1) = -4
               gen_coords(1,2) =  0
               gen_coords(2,1) = -4
               gen_coords(2,2) =  0

               do k=1,2*bc_size+1
                  tmp_map_data(k) = (k-6)*space_step
               end do

            case(2)

               gen_coords(1,1) = -3
               gen_coords(1,2) =  1
               gen_coords(2,1) = -3
               gen_coords(2,2) =  1

               do k=1,2*bc_size+1
                  tmp_map_data(k) = (k-5)*space_step
               end do

            case(3)

               gen_coords(1,1) = 1
               gen_coords(1,2) = 5
               gen_coords(2,1) = 1
               gen_coords(2,2) = 5

               do k=1,2*bc_size+1
                  tmp_map_data(k) = (k-1)*space_step
               end do

            case(4)

               gen_coords(1,1) = 8
               gen_coords(1,2) = 13
               gen_coords(2,1) = 8
               gen_coords(2,2) = 13

               do k=1,2*bc_size+1
                  tmp_map_data(k) = (k+6)*space_step
               end do

            case(5)

               gen_coords(1,1) = 11
               gen_coords(1,2) = 15
               gen_coords(2,1) = 11
               gen_coords(2,2) = 15

               do k=1,2*bc_size+1
                  tmp_map_data(k) = (k+9)*space_step
               end do

            case default
               print '(''test_bf_layer_newgrdpt_procedure'')'
               print '(''get_data_test_get_map_for_newgrdpt'')'
               print '(''test not yet implemented: '',I2)',k
               stop ''

          end select

        end subroutine get_data_test_get_map_for_newgrdpt


        function compare_tmp_maps(tmp_map,tmp_map_data,detailled)
     $     result(test_validated)

          implicit none

          real(rkind), dimension(2*bc_size), intent(in) :: tmp_map
          real(rkind), dimension(2*bc_size), intent(in) :: tmp_map_data
          logical                          , intent(in) :: detailled
          logical                                       :: test_validated

          integer :: k
          logical :: test_loc

          test_validated = .true.

          do k=1, 2*bc_size+1

             test_loc = is_test_validated(
     $            tmp_map(k),
     $            tmp_map_data(k),
     $            detailled)

             test_validated = test_validated.and.test_loc

          end do

        end function compare_tmp_maps

      
        subroutine initialize_map(interior_map,space_step)

          implicit none

          real(rkind), dimension(:), intent(out) :: interior_map
          real(rkind)              , intent(in)  :: space_step


          integer     :: i

          do i=1, size(interior_map,1)
             interior_map(i) = (i-1)*space_step
          end do

        end subroutine initialize_map

      end program test_bf_layer_newgrdpt_procedure
