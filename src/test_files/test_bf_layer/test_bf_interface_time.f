      program test_bf_interface_time

        use bc_operators_class, only :
     $       bc_operators

        use bf_interface_time_class, only :
     $       bf_interface_time

        use bf_sublayer_class, only :
     $       bf_sublayer

        use check_data_module, only :
     $       is_int_vector_validated,
     $       is_int_matrix_validated,
     $       is_real_matrix3D_validated

        use parameters_bf_layer, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NW_corner_type,
     $       NE_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       
     $       align_N, align_S,
     $       align_E, align_W,
     $       no_overlap,
     $       
     $       interior_pt,
     $       bc_interior_pt,
     $       bc_pt

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use pmodel_eq_class, only :
     $       pmodel_eq

        use sd_operators_class, only :
     $       sd_operators

        use td_operators_class ,only :
     $       td_operators

        use rk3tvd_steps_module, only :
     $       compute_1st_step,
     $       compute_1st_step_nopt

        implicit none

        real(rkind), parameter :: x_min = -10.0d0
        real(rkind), parameter :: y_min = -10.0d0
        real(rkind), parameter :: dx    =  0.2d0
        real(rkind), parameter :: dy    =  0.2d0
        real(rkind), parameter :: t     =  0.0d0
        real(rkind), parameter :: dt    =  0.05d0
        

        logical, parameter :: generate_small_domain = .false.
        
        logical :: detailled
        logical :: test_loc
        logical :: test_validated


        if(generate_small_domain) then

           call generate_small_domain_results()

        else
           
           call check_inputs()
           
           detailled = .true.
           test_validated = .true.

           test_loc = test_initialize_before_timeInt(detailled)
           test_validated = test_validated.and.test_loc
           print '(''test_initialize_before_timeInt: '',L1)', test_loc
           print '()'

           test_loc = test_finalize_after_timeInt(detailled)
           test_validated = test_validated.and.test_loc
           print '(''test_finalize_after_timeInt: '',L1)', test_loc
           print '()'

           test_loc = test_compute_time_dev(detailled)
           test_validated = test_validated.and.test_loc
           print '(''test_compute_time_dev: '',L1)', test_loc
           print '()'

           test_loc = test_compute_integration_step(detailled)
           test_validated = test_validated.and.test_loc
           print '(''test_compute_integration_step: '',L1)', test_loc
           print '()'
           
        end if


        contains

        subroutine generate_small_domain_results()

          implicit none

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          real(rkind), dimension(nx,ny,ne) :: interior_nodes_tmp
          real(rkind), dimension(nx,ny,ne) :: time_dev

          type(sd_operators) :: s
          type(pmodel_eq)    :: p_model
          type(bc_operators) :: bc_used
          type(td_operators) :: td_operators_used

          integer(ikind) :: i,j

          integer :: ios

          call check_inputs_small_domain()


          !> initialization
          interior_x_map = (/(x_min+dx*(i-1),i=1,nx)/)
          interior_y_map = (/(y_min+dy*(j-1),j=1,ny)/)
          
          call apply_ic(
     $         interior_nodes,
     $         interior_x_map,
     $         interior_y_map)

          !> write the nodes0 on an output file
          open(unit=2,
     $         file='nodes0.out',
     $         action="write",
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             write(unit=2, iostat=ios) interior_nodes
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          print '(''nodes0 written in nodes0.out'')'


          !> compute the time derivatives
          time_dev = td_operators_used%compute_time_dev(
     $         t,
     $         interior_nodes,
     $         interior_x_map,
     $         interior_y_map,
     $         s,
     $         p_model,
     $         bc_used)

          !> write the time derivatives on an output file
          open(unit=2,
     $         file='timedev.out',
     $         action="write", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             write(unit=2, iostat=ios) time_dev
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          print '(''timedev written in timedev.out'')'


          !> compute the time integration step
          call compute_1st_step(
     $         interior_nodes,
     $         dt,
     $         interior_nodes_tmp,
     $         time_dev,
     $         x_borders=[1,nx],
     $         y_borders=[1,ny])

          !> write the nodes1st on an output file
          open(unit=2,
     $         file='nodes1st.out',
     $         action="write",
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             write(unit=2, iostat=ios) interior_nodes
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          print '(''nodes1st written in nodes1st.out'')'

        end subroutine generate_small_domain_results


        function test_initialize_before_timeInt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_interface_time)              :: bf_interface_used
          integer, dimension(:,:), allocatable :: interior_bc_sections

          integer, dimension(:,:), allocatable :: bc_sections
          integer, dimension(2,4)              :: x_borders_test
          integer, dimension(2,4)              :: y_borders_test

          type(bf_sublayer), pointer :: bf_sublayer_ptr

          logical :: test_loc
          integer :: k


          test_validated = .true.


          !input
          call ini_bf_interface_for_tests(bf_interface_used)

          !output
          call bf_interface_used%initialize_before_timeInt(interior_bc_sections)

          !validation
          !test the interior bc_sections
          test_loc = .not.allocated(interior_bc_sections)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test interior_bc_sections failed'')'
          end if

          !test the N bc_sections
          allocate(bc_sections(5,5))
          bc_sections(:,1) = [W_edge_type   , 1 , 3, 40        , no_overlap]
          bc_sections(:,2) = [E_edge_type   , 99, 3, 40        , no_overlap]
          bc_sections(:,3) = [NW_corner_type, 1 ,41, no_overlap, no_overlap]
          bc_sections(:,4) = [N_edge_type   , 3 ,41, 98        , no_overlap]
          bc_sections(:,5) = [NE_corner_type, 99,41, no_overlap, no_overlap]

          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%bc_sections,
     $         bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test N_bc_sections failed'')'
          end if

          !test the S bc_sections
          bc_sections(:,1) = [SW_corner_type,  1, 1, no_overlap, no_overlap]
          bc_sections(:,2) = [S_edge_type   ,  3, 1, 98        , no_overlap]
          bc_sections(:,3) = [SE_corner_type, 99, 1, no_overlap, no_overlap]
          bc_sections(:,4) = [W_edge_type   ,  1, 3, 20        , no_overlap]
          bc_sections(:,5) = [E_edge_type   , 99, 3, 20        , no_overlap]

          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%ptr%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%bc_sections,
     $         bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test S_bc_sections failed'')'
          end if
          deallocate(bc_sections)

          !test the E bc_sections
          allocate(bc_sections(5,1))
          bc_sections(:,1) = [E_edge_type, 11, 3, 52, no_overlap]

          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%bc_sections,
     $         bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test E_bc_sections failed'')'
          end if

          !test the W bc_sections
          bc_sections(:,1) = [W_edge_type, 1, 3, 52, no_overlap]

          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%ptr%get_head_sublayer()
          test_loc = is_int_matrix_validated(
     $         bf_sublayer_ptr%bc_sections,
     $         bc_sections,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test W_bc_sections failed'')'
          end if
          deallocate(bc_sections)

          !test the x_borders and y_borders: N
          x_borders_test(:,N) = [1,100]
          y_borders_test(:,N) = [3,42]
          x_borders_test(:,S) = [1,100]
          y_borders_test(:,S) = [1,20]
          x_borders_test(:,E) = [3,12]
          y_borders_test(:,E) = [3,52]
          x_borders_test(:,W) = [1,30]
          y_borders_test(:,W) = [3,52]

          do k=1,4
             
             bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(k)%ptr%get_head_sublayer()
             test_loc = is_int_vector_validated(
     $            bf_sublayer_ptr%x_borders,
     $            x_borders_test(:,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test x_borders('',I2,'') failed'')',k
             end if

             test_loc = is_int_vector_validated(
     $            bf_sublayer_ptr%y_borders,
     $            y_borders_test(:,k),
     $            detailled)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test y_borders('',I2,'') failed'')',k
             end if

             test_loc = allocated(bf_sublayer_ptr%bf_compute_used%nodes_tmp)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test nodes_tmp allocated('',I2,'') failed'')',k
             end if
             
             test_loc = allocated(bf_sublayer_ptr%bf_compute_used%time_dev)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test time_dev allocated('',I2,'') failed'')',k
             end if

          end do

        end function test_initialize_before_timeInt


        function test_finalize_after_timeInt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          type(bf_interface_time) :: bf_interface_used
          integer, dimension(:,:), allocatable :: interior_bc_sections

          type(bf_sublayer), pointer :: bf_sublayer_ptr

          logical :: test_loc
          integer :: k


          test_validated = .true.


          !input
          call ini_bf_interface_for_tests(bf_interface_used)
          call bf_interface_used%initialize_before_timeInt(interior_bc_sections)

          !output
          call bf_interface_used%finalize_after_timeInt()

          !validation
          do k=1,4
             
             bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(k)%ptr%get_head_sublayer()

             test_loc = .not.allocated(bf_sublayer_ptr%bf_compute_used%nodes_tmp)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test nodes_tmp deallocated('',I2,'') failed'')',k
             end if
             
             test_loc = .not.allocated(bf_sublayer_ptr%bf_compute_used%time_dev)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test time_dev deallocated('',I2,'') failed'')',k
             end if

          end do

        end function test_finalize_after_timeInt


        function test_compute_time_dev(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes          

          type(bf_interface_time)              :: bf_interface_used
          integer, dimension(:,:), allocatable :: interior_bc_sections

          type(td_operators) :: td_operators_used
          type(sd_operators) :: sd_operators_used
          type(pmodel_eq)    :: p_model
          type(bc_operators) :: bc_used

          type(bf_sublayer), pointer :: bf_sublayer_ptr

          real(rkind), dimension(100,110,ne) :: time_dev_test


          logical        :: test_loc
          integer(ikind) :: i,j
          integer        :: k

          integer :: ios


          test_validated = .true.


          !input
          interior_x_map = (/(x_min+(28+i-1)*dx,i=1,nx)/)
          interior_y_map = (/(y_min+(18+j-1)*dy,j=1,ny)/)

          call apply_ic(
     $            interior_nodes,
     $            interior_x_map,
     $            interior_y_map)

          call ini_bf_interface_for_tests(bf_interface_used)
          call bf_interface_used%initialize_before_timeInt(interior_bc_sections)

          do k=1,4

             bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(k)%ptr%get_head_sublayer()

             call apply_ic(
     $            bf_sublayer_ptr%nodes,
     $            bf_sublayer_ptr%x_map,
     $            bf_sublayer_ptr%y_map)

          end do


          !output
          call bf_interface_used%compute_time_dev(
     $         td_operators_used,
     $         t,sd_operators_used,p_model,bc_used,
     $         interior_nodes)


          !validation: read the results from the small domain
          open(unit=2,
     $         file='timedev.out',
     $         action="read", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             read(unit=2, iostat=ios) time_dev_test
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          
          !north buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%bf_compute_used%time_dev(:,3:42,:),
     $         time_dev_test(:,71:110,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''time_dev(N) failed'')'
          end if

          !south buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%bf_compute_used%time_dev(:,1:20,:),
     $         time_dev_test(:,1:20,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''time_dev(S) failed'')'
          end if

          !west buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%bf_compute_used%time_dev(1:30,3:52,:),
     $         time_dev_test(1:30,21:70,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''time_dev(W) failed'')'
          end if

          !east buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%bf_compute_used%time_dev(3:12,3:52,:),
     $         time_dev_test(91:100,21:70,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''time_dev(E) failed'')'
          end if


        end function test_compute_time_dev


        function test_compute_integration_step(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes
          real(rkind), dimension(nx,ny,ne) :: interior_nodes_tmp
          real(rkind), dimension(nx,ny,ne) :: interior_time_dev

          type(bf_interface_time)              :: bf_interface_used
          integer, dimension(:,:), allocatable :: interior_bc_sections

          type(td_operators) :: td_operators_used
          type(sd_operators) :: sd_operators_used
          type(pmodel_eq)    :: p_model
          type(bc_operators) :: bc_used

          type(bf_sublayer), pointer :: bf_sublayer_ptr

          real(rkind), dimension(100,110,ne) :: nodes0_test
          real(rkind), dimension(100,110,ne) :: nodes1st_test
          real(rkind), dimension(100,110,ne) :: time_dev_test


          logical        :: test_loc
          integer(ikind) :: i,j
          integer        :: k

          integer :: ios


          test_validated = .true.


          !input
          interior_x_map = (/(x_min+(28+i-1)*dx,i=1,nx)/)
          interior_y_map = (/(y_min+(18+j-1)*dy,j=1,ny)/)

          call apply_ic(
     $            interior_nodes,
     $            interior_x_map,
     $            interior_y_map)

          call ini_bf_interface_for_tests(bf_interface_used)
          call bf_interface_used%initialize_before_timeInt(interior_bc_sections)

          do k=1,4

             bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(k)%ptr%get_head_sublayer()

             call apply_ic(
     $            bf_sublayer_ptr%nodes,
     $            bf_sublayer_ptr%x_map,
     $            bf_sublayer_ptr%y_map)

          end do

          interior_time_dev = td_operators_used%compute_time_dev(
     $         t,
     $         interior_nodes,
     $         interior_x_map,
     $         interior_y_map,
     $         sd_operators_used,
     $         p_model,
     $         bc_used,
     $         bc_sections=interior_bc_sections)

          call bf_interface_used%compute_time_dev(
     $         td_operators_used,
     $         t,sd_operators_used,p_model,bc_used,
     $         interior_nodes)


          !output
          call compute_1st_step(
     $         interior_nodes,
     $         dt,
     $         interior_nodes_tmp,
     $         interior_time_dev)

          call bf_interface_used%compute_integration_step(
     $         dt, compute_1st_step_nopt,
     $         interior_nodes)


          !nodes 0 validation
          !------------------------------------------------------------
          !validation: read the results from the small domain: 0th step
          open(unit=2,
     $         file='nodes0.out',
     $         action="read", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             read(unit=2, iostat=ios) nodes0_test
             close(unit=2)
          else
             stop 'file opening pb'
          end if


          !interior nodes
          test_loc = is_real_matrix3D_validated(
     $         interior_nodes_tmp(3:62,3:52,:),
     $         nodes0_test(31:90,21:70,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes0 interior failed'')'
          end if
          

          !north buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%bf_compute_used%nodes_tmp,
     $         nodes0_test(:,69:110,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes0(N) failed'')'
          end if


          !south buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%bf_compute_used%nodes_tmp,
     $         nodes0_test(:,1:22,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes0(S) failed'')'
          end if


          !east buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%bf_compute_used%nodes_tmp,
     $         nodes0_test(89:100,19:72,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes0(E) failed'')'
          end if

          !west buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%bf_compute_used%nodes_tmp,
     $         nodes0_test(1:32,19:72,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes0(W) failed'')'
          end if


          !time dev validation
          !------------------------------------------------------------
          !validation: read the results from the small domain: time_dev
          open(unit=2,
     $         file='timedev.out',
     $         action="read", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             read(unit=2, iostat=ios) time_dev_test
             close(unit=2)
          else
             stop 'file opening pb'
          end if


          !interior nodes
          test_loc = is_real_matrix3D_validated(
     $         interior_time_dev(3:62,3:52,:),
     $         time_dev_test(31:90,21:70,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''time_dev interior failed'')'
          end if


          !north buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%bf_compute_used%time_dev(:,3:42,:),
     $         time_dev_test(:,71:110,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''time_dev(N) failed'')'
          end if

          !south buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%bf_compute_used%time_dev(:,1:20,:),
     $         time_dev_test(:,1:20,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''time_dev(S) failed'')'
          end if

          !east buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%bf_compute_used%time_dev(3:12,3:52,:),
     $         time_dev_test(91:100,21:70,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''time_dev(E) failed'')'
          end if

          !west buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%bf_compute_used%time_dev(1:30,3:52,:),
     $         time_dev_test(1:30,21:70,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''time_dev(W) failed'')'
          end if


          !nodes 1st step validation
          !------------------------------------------------------------
          !validation: read the results from the small domain: 1st step
          open(unit=2,
     $         file='nodes1st.out',
     $         action="read", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             read(unit=2, iostat=ios) nodes1st_test
             close(unit=2)
          else
             stop 'file opening pb'
          end if


          !interior nodes
          test_loc = is_real_matrix3D_validated(
     $         interior_nodes(:,:,:),
     $         nodes1st_test(29:92,19:72,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes1st interior failed'')'
          end if


          !north buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%nodes(:,1:42,:),
     $         nodes1st_test(:,69:110,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes1st(N) failed'')'
          end if

          !south buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%nodes(:,1:22,:),
     $         nodes1st_test(:,1:22,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes1st(S) failed'')'
          end if

          !west buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%nodes,
     $         nodes1st_test(1:32,19:72,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes1st(W) failed'')'
          end if

          !east buffer layer
          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
          test_loc = is_real_matrix3D_validated(
     $         bf_sublayer_ptr%nodes,
     $         nodes1st_test(89:100,19:72,:),
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes1st(E) failed'')'
          end if

        end function test_compute_integration_step          
          

        subroutine ini_bf_interface_for_tests(bf_interface_used)

          implicit none

          type(bf_interface_time), intent(inout) :: bf_interface_used

          type(bf_sublayer), pointer     :: added_sublayer
          
          real(rkind), dimension(nx)       :: interior_x_map
          real(rkind), dimension(ny)       :: interior_y_map
          real(rkind), dimension(nx,ny,ne) :: interior_nodes

          integer(ikind), dimension(2,2) :: bf_alignment_tmp
          integer, dimension(:,:), allocatable :: grdpts_id

          integer(ikind) :: i,j

          
          interior_x_map = (/(x_min+(28+i-1)*dx,i=1,nx)/)
          interior_y_map = (/(y_min+(18+j-1)*dy,j=1,ny)/)


          call bf_interface_used%ini(
     $         interior_x_map,
     $         interior_y_map)


          !north buffer layer
          bf_alignment_tmp = reshape((/
     $         align_W+1, align_N, align_E-1, align_N/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         N,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          bf_alignment_tmp = reshape((/
     $         align_W-27, align_N, align_E+7, align_N+37/),
     $         (/2,2/))

          call bf_interface_used%reallocate_sublayer(
     $         added_sublayer,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          allocate(grdpts_id(100,42))
          grdpts_id = reshape((/
     $         ((interior_pt,i=1,100),j=1,42)/),
     $         (/100,42/))
          grdpts_id(2,:)   = bc_interior_pt
          grdpts_id(99,:)  = bc_interior_pt
          grdpts_id(:,41)  = bc_interior_pt
          grdpts_id(1,:)   = bc_pt
          grdpts_id(100,:) = bc_pt
          grdpts_id(:,42)  = bc_pt
          added_sublayer%grdpts_id = grdpts_id
          deallocate(grdpts_id)


          !south buffer layer
          bf_alignment_tmp = reshape((/
     $         align_W+1, align_S, align_E-1, align_S/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         S,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          bf_alignment_tmp = reshape((/
     $         align_W-27, align_S-17, align_E+7, align_S/),
     $         (/2,2/))

          call bf_interface_used%reallocate_sublayer(
     $         added_sublayer,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          allocate(grdpts_id(100,22))
          grdpts_id = reshape((/
     $         ((interior_pt,i=1,100),j=1,22)/),
     $         (/100,22/))
          grdpts_id(2,:)   = bc_interior_pt
          grdpts_id(99,:)  = bc_interior_pt
          grdpts_id(:,2)   = bc_interior_pt
          grdpts_id(1,:)   = bc_pt
          grdpts_id(100,:) = bc_pt
          grdpts_id(:,1)   = bc_pt
          added_sublayer%grdpts_id = grdpts_id
          deallocate(grdpts_id)


          !west buffer layer
          bf_alignment_tmp = reshape((/
     $         align_W-27, align_S+1, align_W, align_N-1/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         W,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          allocate(grdpts_id(32,54))
          grdpts_id = reshape((/
     $         ((interior_pt,i=1,32),j=1,54)/),
     $         (/32,54/))
          grdpts_id(1,:) = bc_pt
          grdpts_id(2,:) = bc_interior_pt
          added_sublayer%grdpts_id = grdpts_id
          deallocate(grdpts_id)


          !east buffer layer
          bf_alignment_tmp = reshape((/
     $         align_E, align_S+1, align_E+7, align_N-1/),
     $         (/2,2/))

          added_sublayer => bf_interface_used%allocate_sublayer(
     $         E,
     $         interior_x_map,
     $         interior_y_map,
     $         interior_nodes,
     $         bf_alignment_tmp)

          allocate(grdpts_id(12,54))
          grdpts_id = reshape((/
     $         ((interior_pt,i=1,12),j=1,54)/),
     $         (/12,54/))
          grdpts_id(11,:) = bc_interior_pt
          grdpts_id(12,:) = bc_pt
          added_sublayer%grdpts_id = grdpts_id
          deallocate(grdpts_id)

        end subroutine ini_bf_interface_for_tests


        subroutine apply_ic(nodes,x_map,y_map)
        
          implicit none

          real(rkind), dimension(:,:,:), intent(out) :: nodes
          real(rkind), dimension(:)    , intent(in)  :: x_map
          real(rkind), dimension(:)    , intent(in)  :: y_map          

          integer(ikind) :: i,j
          real(rkind)    :: x,y

          
          do j=1, size(y_map,1)

             y = y_map(j)/10.0d0*(2.0d0*ACOS(-1.0d0))

             do i=1,size(x_map,1)

                x = x_map(i)/10.0d0*(2.0d0*ACOS(-1.0d0))

                nodes(i,j,1) = SIN(x)*SIN(y)
                nodes(i,j,2) = COS(x)*SIN(y)
                nodes(i,j,3) = SIN(x)*COS(y)

             end do

          end do

        end subroutine apply_ic



        subroutine check_inputs_small_domain()

          if(.not.(
     $       (nx.eq.100).and.
     $       (ny.eq.110).and.
     $       (ne.eq.3))) then

             print '(''the test requires:'')'
             print '(''   - nx=100'')'
             print '(''   - ny=110'')'
             print '(''   - ne=3'')'
             stop ''

          end if

        end subroutine check_inputs_small_domain


        subroutine check_inputs()

          if(.not.(
     $       (nx.eq.64).and.
     $       (ny.eq.54).and.
     $       (ne.eq.3))) then

             print '(''the test requires:'')'
             print '(''   - nx=64'')'
             print '(''   - ny=54'')'
             print '(''   - ne=3'')'
             stop ''

          end if

        end subroutine check_inputs

      end program test_bf_interface_time
