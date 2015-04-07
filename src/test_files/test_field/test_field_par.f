      !  _____________                    _______________
      ! |             |                  |       |       |
      ! |             |                  |       |       |
      ! |             |  compared with   |_______|_______|
      ! |             |                  |       |       |
      ! |             |                  |       |       |
      ! |_____________|                  |_______|_______|
      !
      !      field                           field_par
      !------------------------------------------------------------
      program test_field_par

        use bc_operators_class, only :
     $       bc_operators

        use check_data_module, only :
     $       is_int_vector_validated,
     $       is_int_matrix_validated,
     $       is_real_validated,
     $       is_real_matrix3D_validated

        use field_class, only :
     $       field

        use field_par_class, only :
     $       field_par

        use mpi

        use mpi_process_class, only :
     $       mpi_process

        use parameters_constant, only :
     $       N,S,E,W,
     $       sincos,
     $       reflection_xy_choice

        use parameters_input, only :
     $       nx,ny,ne,
     $       x_min,
     $       y_min,
     $       x_max,
     $       y_max,
     $       ic_choice,
     $       bc_choice,
     $       npx,npy,
     $       ntx,nty

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

        real(rkind), parameter :: t     =  0.0d0
        real(rkind), parameter :: dt    =  0.05d0
        

        logical, parameter :: generate_small_domain = .false.
        
        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        type(mpi_process) :: mpi_process_used
        integer :: ierror
        integer :: rank


        if(generate_small_domain) then

           call generate_small_domain_results()

        else
           
           call check_inputs()
           
           detailled = .true.
           test_validated = .true.

           test_loc = test_ini(detailled)
           test_validated = test_validated.and.test_loc
           call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
           if(rank.eq.0) then
              print '(''test_ini: '',L1)', test_loc
              print '()'
           end if

           test_loc = test_compute_integration_step(detailled)
           test_validated = test_validated.and.test_loc
           call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
           if(rank.eq.0) then
              print '(''test_compute_integration_step: '',L1)', test_loc
              print '()'
           end if
c$$$
c$$$           test_loc = test_integrate(detailled)
c$$$           test_validated = test_validated.and.test_loc
c$$$           call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
c$$$           if(rank.eq.0) then
c$$$              print '(''test_integrate: '',L1)', test_loc
c$$$              print '()'
c$$$           end if
           
           
           call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
           if(rank.eq.0) then
              print '(''test_validated: '',L1)', test_validated
              print '()'
           end if

           call mpi_process_used%finalize_mpi()


        end if


        contains


c$$$        function test_integrate(detailled) result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          logical, intent(in) :: detailled
c$$$          logical             :: test_validated
c$$$
c$$$
c$$$          type(field_extended)              :: field_used
c$$$          type(bf_sublayer), pointer        :: bf_sublayer_ptr
c$$$          real(rkind), dimension(100,110,3) :: nodesInt
c$$$
c$$$          logical :: test_loc
c$$$          integer :: ios
c$$$
c$$$
c$$$          test_validated = .true.
c$$$
c$$$
c$$$          !output
c$$$          !------------------------------------------------------------
c$$$          call field_used%ini()
c$$$          call ini_bf_interface_for_tests(field_used)
c$$$          call field_used%apply_initial_conditions()
c$$$          call field_used%integrate(dt)
c$$$
c$$$
c$$$          !validation
c$$$          !------------------------------------------------------------
c$$$
c$$$          !validation: read the results from the small domain
c$$$          open(unit=2,
c$$$     $         file='field_nodesInt.out',
c$$$     $         action="read", 
c$$$     $         status="unknown",
c$$$     $         form='unformatted',
c$$$     $         access='sequential',
c$$$     $         position='rewind',
c$$$     $         iostat=ios)
c$$$          
c$$$          if(ios.eq.0) then
c$$$             read(unit=2, iostat=ios) nodesInt
c$$$             close(unit=2)
c$$$          else
c$$$             stop 'file opening pb'
c$$$          end if
c$$$
c$$$          !interior
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         field_used%nodes(1:64,1:54,:),
c$$$     $         nodesInt(29:92,19:72,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodesInt(interior) failed'')'
c$$$          end if
c$$$
c$$$          !north buffer layer
c$$$          bf_sublayer_ptr => field_used%domain_extension%mainlayer_pointers(N)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%nodes(:,1:42,:),
c$$$     $         nodesInt(:,69:110,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodesInt(N) failed'')'
c$$$          end if
c$$$
c$$$          !south buffer layer
c$$$          bf_sublayer_ptr => field_used%domain_extension%mainlayer_pointers(S)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%nodes(:,1:22,:),
c$$$     $         nodesInt(:,1:22,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodesInt(S) failed'')'
c$$$          end if
c$$$
c$$$          !west buffer layer
c$$$          bf_sublayer_ptr => field_used%domain_extension%mainlayer_pointers(W)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%nodes(1:32,1:54,:),
c$$$     $         nodesInt(1:32,19:72,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodesInt(W) failed'')'
c$$$          end if
c$$$
c$$$          !east buffer layer
c$$$          bf_sublayer_ptr => field_used%domain_extension%mainlayer_pointers(E)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%nodes(1:12,1:54,:),
c$$$     $         nodesInt(89:100,19:72,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodesInt(E) failed'')'
c$$$          end if
c$$$
c$$$        end function test_integrate


        function test_compute_integration_step(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(field_par)                   :: field_used
          real(rkind), dimension(100,110,3) :: nodes1

          real(rkind), dimension(nx,ny,ne)  :: interior_nodes_tmp
          real(rkind), dimension(nx,ny,ne)  :: interior_timedev

          logical :: test_loc
          integer :: ios

          integer :: ierror
          integer :: i
          integer :: rank
          logical, dimension(4) :: test_loc_gathered


          test_validated = .true.


          !output
          !------------------------------------------------------------
          call field_used%ini()

          interior_timedev = field_used%compute_time_dev()

          call field_used%compute_integration_step(
     $         dt,
     $         interior_nodes_tmp,
     $         interior_timedev,
     $         compute_1st_step)

          call field_used%apply_bc_on_nodes()


          !validation
          !------------------------------------------------------------
          !validation: read the results from the small domain
          open(unit=2,
     $         file='field_nodes1.out',
     $         action="read", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             read(unit=2, iostat=ios) nodes1
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          select case(field_used%usr_rank)
      
             !SW tile
             case(0)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodes1(1:52,1:57,:),
     $               detailled)

             !NW tile
             case(1)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodes1(1:52,54:110,:),
     $               detailled)
                
             !SE tile
             case(2)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodes1(49:100,1:57,:),
     $               detailled)
                
             !NE tile
             case(3)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodes1(49:100,54:110,:),
     $               detailled)

          end select

          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes('',I2,'') failed'')', field_used%usr_rank
          end if


          ! the master processor gathers the test results
          call MPI_GATHER(
     $         test_loc,1,MPI_LOGICAL,
     $         test_loc_gathered,1,MPI_LOGICAL,
     $         0, MPI_COMM_WORLD, ierror)

          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

          if(rank.eq.0) then

             test_loc = test_loc_gathered(1)
             do i=2,4
                test_loc = test_loc.and.test_loc_gathered(i)
             end do
             test_validated = test_validated.and.test_loc
          end if

        end function test_compute_integration_step


c$$$        function test_compute_time_dev_ext(detailled) result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          logical, intent(in) :: detailled
c$$$          logical             :: test_validated
c$$$
c$$$
c$$$          type(field_extended)              :: field_used
c$$$          type(bf_sublayer), pointer        :: bf_sublayer_ptr
c$$$          real(rkind), dimension(100,110,3) :: timedev
c$$$          real(rkind), dimension(nx,ny,ne)  :: timedev_ext
c$$$
c$$$          logical :: test_loc
c$$$          integer :: ios
c$$$
c$$$
c$$$          test_validated = .true.
c$$$
c$$$
c$$$          !output
c$$$          !------------------------------------------------------------
c$$$          call field_used%ini()
c$$$          call ini_bf_interface_for_tests(field_used)
c$$$          call field_used%apply_initial_conditions()
c$$$
c$$$          call field_used%domain_extension%initialize_before_timeInt(
c$$$     $         field_used%interior_bc_sections)
c$$$
c$$$          timedev_ext = field_used%compute_time_dev_ext()
c$$$
c$$$
c$$$          !validation
c$$$          !------------------------------------------------------------
c$$$
c$$$          !validation: read the results from the small domain
c$$$          open(unit=2,
c$$$     $         file='field_timedev.out',
c$$$     $         action="read", 
c$$$     $         status="unknown",
c$$$     $         form='unformatted',
c$$$     $         access='sequential',
c$$$     $         position='rewind',
c$$$     $         iostat=ios)
c$$$          
c$$$          if(ios.eq.0) then
c$$$             read(unit=2, iostat=ios) timedev
c$$$             close(unit=2)
c$$$          else
c$$$             stop 'file opening pb'
c$$$          end if
c$$$
c$$$          !interior
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         timedev_ext(3:62,3:52,:),
c$$$     $         timedev(31:90,21:70,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(interior) failed'')'
c$$$          end if
c$$$
c$$$          !north buffer layer
c$$$          bf_sublayer_ptr => field_used%domain_extension%mainlayer_pointers(N)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%time_dev(:,3:42,:),
c$$$     $         timedev(:,71:110,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(N) failed'')'
c$$$          end if
c$$$
c$$$          !south buffer layer
c$$$          bf_sublayer_ptr => field_used%domain_extension%mainlayer_pointers(S)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%time_dev(:,1:20,:),
c$$$     $         timedev(:,1:20,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(S) failed'')'
c$$$          end if
c$$$
c$$$          !west buffer layer
c$$$          bf_sublayer_ptr => field_used%domain_extension%mainlayer_pointers(W)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%time_dev(1:30,3:52,:),
c$$$     $         timedev(1:30,21:70,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(W) failed'')'
c$$$          end if
c$$$
c$$$          !east buffer layer
c$$$          bf_sublayer_ptr => field_used%domain_extension%mainlayer_pointers(E)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%time_dev(3:12,3:52,:),
c$$$     $         timedev(91:100,21:70,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(E) failed'')'
c$$$          end if
c$$$
c$$$        end function test_compute_time_dev_ext


        function test_ini(detailled) result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mpi_process)                 :: mpi_process_used
          type(field_par)                   :: field_used
          real(rkind), dimension(100,110,3) :: nodes0

          logical :: test_loc
          integer :: ios

          integer :: ierror
          integer :: rank
          integer :: i

          logical, dimension(4) :: test_loc_gathered


          test_validated = .true.


          call mpi_process_used%ini_mpi()


          !output
          !------------------------------------------------------------
          call field_used%ini()


          !validation
          !------------------------------------------------------------

          !validation: read the results from the small domain
          open(unit=2,
     $         file='field_nodes0.out',
     $         action="read", 
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             read(unit=2, iostat=ios) nodes0
             close(unit=2)
          else
             stop 'file opening pb'
          end if


          select case(field_used%usr_rank)
      
             !SW tile
             case(0)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodes0(1:52,1:57,:),
     $               detailled)

             !NW tile
             case(1)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodes0(1:52,54:110,:),
     $               detailled)
                
             !SE tile
             case(2)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodes0(49:100,1:57,:),
     $               detailled)
                
             !NE tile
             case(3)
                test_loc = is_real_matrix3D_validated(
     $               field_used%nodes(:,:,:),
     $               nodes0(49:100,54:110,:),
     $               detailled)

          end select

          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''nodes('',I2,'') failed'')', field_used%usr_rank
          end if


          ! the master processor gathers the test results
          call MPI_GATHER(
     $         test_loc,1,MPI_LOGICAL,
     $         test_loc_gathered,1,MPI_LOGICAL,
     $         0, MPI_COMM_WORLD, ierror)

          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

          if(rank.eq.0) then

             test_loc = test_loc_gathered(1)
             do i=2,4
                test_loc = test_loc.and.test_loc_gathered(i)
             end do
             test_validated = test_validated.and.test_loc
          end if
          

        end function test_ini


        subroutine generate_small_domain_results()

          implicit none

          type(field) :: field_used

          real(rkind), dimension(nx,ny,ne) :: time_dev
          real(rkind), dimension(nx,ny,ne) :: interior_nodes_tmp

          integer :: ios


          call check_inputs_small_domain()


          !> initialization
          call field_used%ini()


          !> write the nodes0 on an output file
          open(unit=2,
     $         file='field_nodes0.out',
     $         action="write",
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             write(unit=2, iostat=ios) field_used%nodes
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          print '(''field%nodes written in field_nodes0.out'')'


          !> compute the time derivatives
          time_dev = field_used%compute_time_dev()

          !> write the time derivatives on an output file
          open(unit=2,
     $         file='field_timedev.out',
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

          print '(''time derivatives written in field_timedev.out'')'


          !> compute the time integration step
          call field_used%compute_integration_step(
     $         dt,
     $         interior_nodes_tmp,
     $         time_dev,
     $         compute_1st_step)

          !> write the nodes1st on an output file
          open(unit=2,
     $         file='field_nodes1st.out',
     $         action="write",
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             write(unit=2, iostat=ios) field_used%nodes
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          print '(''field%nodes written in field_nodes1st.out'')'


          !> re-initialization
          call field_used%ini()

          !> complete integration cycle
          call field_used%integrate(dt)

          !> write the nodes1st on an output file
          open(unit=2,
     $         file='field_nodesInt.out',
     $         action="write",
     $         status="unknown",
     $         form='unformatted',
     $         access='sequential',
     $         position='rewind',
     $         iostat=ios)
          
          if(ios.eq.0) then
             write(unit=2, iostat=ios) field_used%nodes
             close(unit=2)
          else
             stop 'file opening pb'
          end if

          print '(''field%nodes after integration written in field_nodesInt.out'')'

        end subroutine generate_small_domain_results


c$$$        function test_initialize_before_timeInt(detailled)
c$$$     $       result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          logical, intent(in) :: detailled
c$$$          logical             :: test_validated
c$$$
c$$$          type(bf_interface_time)              :: bf_interface_used
c$$$          integer, dimension(:,:), allocatable :: interior_bc_sections
c$$$
c$$$          integer, dimension(:,:), allocatable :: bc_sections
c$$$          integer, dimension(2,4)              :: x_borders_test
c$$$          integer, dimension(2,4)              :: y_borders_test
c$$$
c$$$          type(bf_sublayer), pointer :: bf_sublayer_ptr
c$$$
c$$$          logical :: test_loc
c$$$          integer :: k
c$$$
c$$$
c$$$          test_validated = .true.
c$$$
c$$$
c$$$          !input
c$$$          call ini_bf_interface_for_tests(bf_interface_used)
c$$$
c$$$          !output
c$$$          call bf_interface_used%initialize_before_timeInt(interior_bc_sections)
c$$$
c$$$          !validation
c$$$          !test the interior bc_sections
c$$$          test_loc = .not.allocated(interior_bc_sections)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''test interior_bc_sections failed'')'
c$$$          end if
c$$$
c$$$          !test the N bc_sections
c$$$          allocate(bc_sections(5,5))
c$$$          bc_sections(:,1) = [W_edge_type   , 1 , 3, 40        , no_overlap]
c$$$          bc_sections(:,2) = [E_edge_type   , 99, 3, 40        , no_overlap]
c$$$          bc_sections(:,3) = [NW_corner_type, 1 ,41, no_overlap, no_overlap]
c$$$          bc_sections(:,4) = [N_edge_type   , 3 ,41, 98        , no_overlap]
c$$$          bc_sections(:,5) = [NE_corner_type, 99,41, no_overlap, no_overlap]
c$$$
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
c$$$          test_loc = is_int_matrix_validated(
c$$$     $         bf_sublayer_ptr%bc_sections,
c$$$     $         bc_sections,
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''test N_bc_sections failed'')'
c$$$          end if
c$$$
c$$$          !test the S bc_sections
c$$$          bc_sections(:,1) = [SW_corner_type,  1, 1, no_overlap, no_overlap]
c$$$          bc_sections(:,2) = [S_edge_type   ,  3, 1, 98        , no_overlap]
c$$$          bc_sections(:,3) = [SE_corner_type, 99, 1, no_overlap, no_overlap]
c$$$          bc_sections(:,4) = [W_edge_type   ,  1, 3, 20        , no_overlap]
c$$$          bc_sections(:,5) = [E_edge_type   , 99, 3, 20        , no_overlap]
c$$$
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%ptr%get_head_sublayer()
c$$$          test_loc = is_int_matrix_validated(
c$$$     $         bf_sublayer_ptr%bc_sections,
c$$$     $         bc_sections,
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''test S_bc_sections failed'')'
c$$$          end if
c$$$          deallocate(bc_sections)
c$$$
c$$$          !test the E bc_sections
c$$$          allocate(bc_sections(5,1))
c$$$          bc_sections(:,1) = [E_edge_type, 11, 3, 52, no_overlap]
c$$$
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
c$$$          test_loc = is_int_matrix_validated(
c$$$     $         bf_sublayer_ptr%bc_sections,
c$$$     $         bc_sections,
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''test E_bc_sections failed'')'
c$$$          end if
c$$$
c$$$          !test the W bc_sections
c$$$          bc_sections(:,1) = [W_edge_type, 1, 3, 52, no_overlap]
c$$$
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%ptr%get_head_sublayer()
c$$$          test_loc = is_int_matrix_validated(
c$$$     $         bf_sublayer_ptr%bc_sections,
c$$$     $         bc_sections,
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''test W_bc_sections failed'')'
c$$$          end if
c$$$          deallocate(bc_sections)
c$$$
c$$$          !test the x_borders and y_borders: N
c$$$          x_borders_test(:,N) = [1,100]
c$$$          y_borders_test(:,N) = [3,42]
c$$$          x_borders_test(:,S) = [1,100]
c$$$          y_borders_test(:,S) = [1,20]
c$$$          x_borders_test(:,E) = [3,12]
c$$$          y_borders_test(:,E) = [3,52]
c$$$          x_borders_test(:,W) = [1,30]
c$$$          y_borders_test(:,W) = [3,52]
c$$$
c$$$          do k=1,4
c$$$             
c$$$             bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(k)%ptr%get_head_sublayer()
c$$$             test_loc = is_int_vector_validated(
c$$$     $            bf_sublayer_ptr%x_borders,
c$$$     $            x_borders_test(:,k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$             if(detailled.and.(.not.test_loc)) then
c$$$                print '(''test x_borders('',I2,'') failed'')',k
c$$$             end if
c$$$
c$$$             test_loc = is_int_vector_validated(
c$$$     $            bf_sublayer_ptr%y_borders,
c$$$     $            y_borders_test(:,k),
c$$$     $            detailled)
c$$$             test_validated = test_validated.and.test_loc
c$$$             if(detailled.and.(.not.test_loc)) then
c$$$                print '(''test y_borders('',I2,'') failed'')',k
c$$$             end if
c$$$
c$$$             test_loc = allocated(bf_sublayer_ptr%bf_compute_used%nodes_tmp)
c$$$             test_validated = test_validated.and.test_loc
c$$$             if(detailled.and.(.not.test_loc)) then
c$$$                print '(''test nodes_tmp allocated('',I2,'') failed'')',k
c$$$             end if
c$$$             
c$$$             test_loc = allocated(bf_sublayer_ptr%bf_compute_used%time_dev)
c$$$             test_validated = test_validated.and.test_loc
c$$$             if(detailled.and.(.not.test_loc)) then
c$$$                print '(''test time_dev allocated('',I2,'') failed'')',k
c$$$             end if
c$$$
c$$$          end do
c$$$
c$$$        end function test_initialize_before_timeInt
c$$$
c$$$
c$$$        function test_finalize_after_timeInt(detailled)
c$$$     $       result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          logical, intent(in) :: detailled
c$$$          logical             :: test_validated
c$$$
c$$$          type(bf_interface_time) :: bf_interface_used
c$$$          integer, dimension(:,:), allocatable :: interior_bc_sections
c$$$
c$$$          type(bf_sublayer), pointer :: bf_sublayer_ptr
c$$$
c$$$          logical :: test_loc
c$$$          integer :: k
c$$$
c$$$
c$$$          test_validated = .true.
c$$$
c$$$
c$$$          !input
c$$$          call ini_bf_interface_for_tests(bf_interface_used)
c$$$          call bf_interface_used%initialize_before_timeInt(interior_bc_sections)
c$$$
c$$$          !output
c$$$          call bf_interface_used%finalize_after_timeInt()
c$$$
c$$$          !validation
c$$$          do k=1,4
c$$$             
c$$$             bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(k)%ptr%get_head_sublayer()
c$$$
c$$$             test_loc = .not.allocated(bf_sublayer_ptr%bf_compute_used%nodes_tmp)
c$$$             test_validated = test_validated.and.test_loc
c$$$             if(detailled.and.(.not.test_loc)) then
c$$$                print '(''test nodes_tmp deallocated('',I2,'') failed'')',k
c$$$             end if
c$$$             
c$$$             test_loc = .not.allocated(bf_sublayer_ptr%bf_compute_used%time_dev)
c$$$             test_validated = test_validated.and.test_loc
c$$$             if(detailled.and.(.not.test_loc)) then
c$$$                print '(''test time_dev deallocated('',I2,'') failed'')',k
c$$$             end if
c$$$
c$$$          end do
c$$$
c$$$        end function test_finalize_after_timeInt
c$$$
c$$$
c$$$        function test_compute_time_dev(detailled)
c$$$     $       result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          logical, intent(in) :: detailled
c$$$          logical             :: test_validated
c$$$
c$$$          real(rkind), dimension(nx)       :: interior_x_map
c$$$          real(rkind), dimension(ny)       :: interior_y_map
c$$$          real(rkind), dimension(nx,ny,ne) :: interior_nodes          
c$$$
c$$$          type(bf_interface_time)              :: bf_interface_used
c$$$          integer, dimension(:,:), allocatable :: interior_bc_sections
c$$$
c$$$          type(td_operators) :: td_operators_used
c$$$          type(sd_operators) :: sd_operators_used
c$$$          type(pmodel_eq)    :: p_model
c$$$          type(bc_operators) :: bc_used
c$$$
c$$$          type(bf_sublayer), pointer :: bf_sublayer_ptr
c$$$
c$$$          real(rkind), dimension(100,110,ne) :: time_dev_test
c$$$
c$$$
c$$$          logical        :: test_loc
c$$$          integer(ikind) :: i,j
c$$$          integer        :: k
c$$$
c$$$          integer :: ios
c$$$
c$$$
c$$$          test_validated = .true.
c$$$
c$$$
c$$$          !input
c$$$          interior_x_map = (/(x_min+(28+i-1)*dx,i=1,nx)/)
c$$$          interior_y_map = (/(y_min+(18+j-1)*dy,j=1,ny)/)
c$$$
c$$$          call apply_ic(
c$$$     $            interior_nodes,
c$$$     $            interior_x_map,
c$$$     $            interior_y_map)
c$$$
c$$$          call ini_bf_interface_for_tests(bf_interface_used)
c$$$          call bf_interface_used%initialize_before_timeInt(interior_bc_sections)
c$$$
c$$$          do k=1,4
c$$$
c$$$             bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(k)%ptr%get_head_sublayer()
c$$$
c$$$             call apply_ic(
c$$$     $            bf_sublayer_ptr%nodes,
c$$$     $            bf_sublayer_ptr%x_map,
c$$$     $            bf_sublayer_ptr%y_map)
c$$$
c$$$          end do
c$$$
c$$$
c$$$          !output
c$$$          call bf_interface_used%compute_time_dev(
c$$$     $         td_operators_used,
c$$$     $         t,sd_operators_used,p_model,bc_used,
c$$$     $         interior_nodes)
c$$$
c$$$
c$$$          !validation: read the results from the small domain
c$$$          open(unit=2,
c$$$     $         file='timedev.out',
c$$$     $         action="read", 
c$$$     $         status="unknown",
c$$$     $         form='unformatted',
c$$$     $         access='sequential',
c$$$     $         position='rewind',
c$$$     $         iostat=ios)
c$$$          
c$$$          if(ios.eq.0) then
c$$$             read(unit=2, iostat=ios) time_dev_test
c$$$             close(unit=2)
c$$$          else
c$$$             stop 'file opening pb'
c$$$          end if
c$$$
c$$$          
c$$$          !north buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%time_dev(:,3:42,:),
c$$$     $         time_dev_test(:,71:110,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(N) failed'')'
c$$$          end if
c$$$
c$$$          !south buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%time_dev(:,1:20,:),
c$$$     $         time_dev_test(:,1:20,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(S) failed'')'
c$$$          end if
c$$$
c$$$          !west buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%time_dev(1:30,3:52,:),
c$$$     $         time_dev_test(1:30,21:70,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(W) failed'')'
c$$$          end if
c$$$
c$$$          !east buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%time_dev(3:12,3:52,:),
c$$$     $         time_dev_test(91:100,21:70,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(E) failed'')'
c$$$          end if
c$$$
c$$$
c$$$        end function test_compute_time_dev
c$$$
c$$$
c$$$        function test_compute_integration_step(detailled)
c$$$     $       result(test_validated)
c$$$
c$$$          implicit none
c$$$
c$$$          logical, intent(in) :: detailled
c$$$          logical             :: test_validated
c$$$
c$$$          real(rkind), dimension(nx)       :: interior_x_map
c$$$          real(rkind), dimension(ny)       :: interior_y_map
c$$$          real(rkind), dimension(nx,ny,ne) :: interior_nodes
c$$$          real(rkind), dimension(nx,ny,ne) :: interior_nodes_tmp
c$$$          real(rkind), dimension(nx,ny,ne) :: interior_time_dev
c$$$
c$$$          type(bf_interface_time)              :: bf_interface_used
c$$$          integer, dimension(:,:), allocatable :: interior_bc_sections
c$$$
c$$$          type(td_operators) :: td_operators_used
c$$$          type(sd_operators) :: sd_operators_used
c$$$          type(pmodel_eq)    :: p_model
c$$$          type(bc_operators) :: bc_used
c$$$
c$$$          type(bf_sublayer), pointer :: bf_sublayer_ptr
c$$$
c$$$          real(rkind), dimension(100,110,ne) :: nodes0_test
c$$$          real(rkind), dimension(100,110,ne) :: nodes1st_test
c$$$          real(rkind), dimension(100,110,ne) :: time_dev_test
c$$$
c$$$
c$$$          logical        :: test_loc
c$$$          integer(ikind) :: i,j
c$$$          integer        :: k
c$$$
c$$$          integer :: ios
c$$$
c$$$
c$$$          test_validated = .true.
c$$$
c$$$
c$$$          !input
c$$$          interior_x_map = (/(x_min+(28+i-1)*dx,i=1,nx)/)
c$$$          interior_y_map = (/(y_min+(18+j-1)*dy,j=1,ny)/)
c$$$
c$$$          call apply_ic(
c$$$     $            interior_nodes,
c$$$     $            interior_x_map,
c$$$     $            interior_y_map)
c$$$
c$$$          call ini_bf_interface_for_tests(bf_interface_used)
c$$$          call bf_interface_used%initialize_before_timeInt(interior_bc_sections)
c$$$
c$$$          do k=1,4
c$$$
c$$$             bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(k)%ptr%get_head_sublayer()
c$$$
c$$$             call apply_ic(
c$$$     $            bf_sublayer_ptr%nodes,
c$$$     $            bf_sublayer_ptr%x_map,
c$$$     $            bf_sublayer_ptr%y_map)
c$$$
c$$$          end do
c$$$
c$$$          interior_time_dev = td_operators_used%compute_time_dev(
c$$$     $         t,
c$$$     $         interior_nodes,
c$$$     $         interior_x_map,
c$$$     $         interior_y_map,
c$$$     $         sd_operators_used,
c$$$     $         p_model,
c$$$     $         bc_used,
c$$$     $         bc_sections=interior_bc_sections)
c$$$
c$$$          call bf_interface_used%compute_time_dev(
c$$$     $         td_operators_used,
c$$$     $         t,sd_operators_used,p_model,bc_used,
c$$$     $         interior_nodes)
c$$$
c$$$
c$$$          !output
c$$$          call compute_1st_step(
c$$$     $         interior_nodes,
c$$$     $         dt,
c$$$     $         interior_nodes_tmp,
c$$$     $         interior_time_dev)
c$$$
c$$$          call bf_interface_used%compute_integration_step(
c$$$     $         dt, compute_1st_step_nopt,
c$$$     $         interior_nodes)
c$$$
c$$$
c$$$          !nodes 0 validation
c$$$          !------------------------------------------------------------
c$$$          !validation: read the results from the small domain: 0th step
c$$$          open(unit=2,
c$$$     $         file='nodes0.out',
c$$$     $         action="read", 
c$$$     $         status="unknown",
c$$$     $         form='unformatted',
c$$$     $         access='sequential',
c$$$     $         position='rewind',
c$$$     $         iostat=ios)
c$$$          
c$$$          if(ios.eq.0) then
c$$$             read(unit=2, iostat=ios) nodes0_test
c$$$             close(unit=2)
c$$$          else
c$$$             stop 'file opening pb'
c$$$          end if
c$$$
c$$$
c$$$          !interior nodes
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         interior_nodes_tmp(3:62,3:52,:),
c$$$     $         nodes0_test(31:90,21:70,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodes0 interior failed'')'
c$$$          end if
c$$$          
c$$$
c$$$          !north buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%nodes_tmp,
c$$$     $         nodes0_test(:,69:110,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodes0(N) failed'')'
c$$$          end if
c$$$
c$$$
c$$$          !south buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%nodes_tmp,
c$$$     $         nodes0_test(:,1:22,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodes0(S) failed'')'
c$$$          end if
c$$$
c$$$
c$$$          !east buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%nodes_tmp,
c$$$     $         nodes0_test(89:100,19:72,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodes0(E) failed'')'
c$$$          end if
c$$$
c$$$          !west buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%nodes_tmp,
c$$$     $         nodes0_test(1:32,19:72,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodes0(W) failed'')'
c$$$          end if
c$$$
c$$$
c$$$          !time dev validation
c$$$          !------------------------------------------------------------
c$$$          !validation: read the results from the small domain: time_dev
c$$$          open(unit=2,
c$$$     $         file='timedev.out',
c$$$     $         action="read", 
c$$$     $         status="unknown",
c$$$     $         form='unformatted',
c$$$     $         access='sequential',
c$$$     $         position='rewind',
c$$$     $         iostat=ios)
c$$$          
c$$$          if(ios.eq.0) then
c$$$             read(unit=2, iostat=ios) time_dev_test
c$$$             close(unit=2)
c$$$          else
c$$$             stop 'file opening pb'
c$$$          end if
c$$$
c$$$
c$$$          !interior nodes
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         interior_time_dev(3:62,3:52,:),
c$$$     $         time_dev_test(31:90,21:70,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev interior failed'')'
c$$$          end if
c$$$
c$$$
c$$$          !north buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%time_dev(:,3:42,:),
c$$$     $         time_dev_test(:,71:110,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(N) failed'')'
c$$$          end if
c$$$
c$$$          !south buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%time_dev(:,1:20,:),
c$$$     $         time_dev_test(:,1:20,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(S) failed'')'
c$$$          end if
c$$$
c$$$          !east buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%time_dev(3:12,3:52,:),
c$$$     $         time_dev_test(91:100,21:70,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(E) failed'')'
c$$$          end if
c$$$
c$$$          !west buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%bf_compute_used%time_dev(1:30,3:52,:),
c$$$     $         time_dev_test(1:30,21:70,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''time_dev(W) failed'')'
c$$$          end if
c$$$
c$$$
c$$$          !nodes 1st step validation
c$$$          !------------------------------------------------------------
c$$$          !validation: read the results from the small domain: 1st step
c$$$          open(unit=2,
c$$$     $         file='nodes1st.out',
c$$$     $         action="read", 
c$$$     $         status="unknown",
c$$$     $         form='unformatted',
c$$$     $         access='sequential',
c$$$     $         position='rewind',
c$$$     $         iostat=ios)
c$$$          
c$$$          if(ios.eq.0) then
c$$$             read(unit=2, iostat=ios) nodes1st_test
c$$$             close(unit=2)
c$$$          else
c$$$             stop 'file opening pb'
c$$$          end if
c$$$
c$$$
c$$$          !interior nodes
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         interior_nodes(:,:,:),
c$$$     $         nodes1st_test(29:92,19:72,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodes1st interior failed'')'
c$$$          end if
c$$$
c$$$
c$$$          !north buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(N)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%nodes(:,1:42,:),
c$$$     $         nodes1st_test(:,69:110,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodes1st(N) failed'')'
c$$$          end if
c$$$
c$$$          !south buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(S)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%nodes(:,1:22,:),
c$$$     $         nodes1st_test(:,1:22,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodes1st(S) failed'')'
c$$$          end if
c$$$
c$$$          !west buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(W)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%nodes,
c$$$     $         nodes1st_test(1:32,19:72,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodes1st(W) failed'')'
c$$$          end if
c$$$
c$$$          !east buffer layer
c$$$          bf_sublayer_ptr => bf_interface_used%mainlayer_pointers(E)%ptr%get_head_sublayer()
c$$$          test_loc = is_real_matrix3D_validated(
c$$$     $         bf_sublayer_ptr%nodes,
c$$$     $         nodes1st_test(89:100,19:72,:),
c$$$     $         detailled)
c$$$          test_validated = test_validated.and.test_loc
c$$$          if(detailled.and.(.not.test_loc)) then
c$$$             print '(''nodes1st(E) failed'')'
c$$$          end if
c$$$
c$$$        end function test_compute_integration_step          
          

        subroutine check_inputs_small_domain()

          if(.not.(
     $       (npx.eq.1).and.
     $       (npy.eq.1).and.
     $       (ntx.eq.100).and.
     $       (nty.eq.110).and.
     $       (ne.eq.3).and.
     $       is_real_validated(x_min,-10.0d0,.true.).and.
     $       is_real_validated(x_max,-0.50d0,.true.).and.
     $       is_real_validated(y_min,-10.0d0,.true.).and.
     $       is_real_validated(y_max, 11.0d0,.true.).and.
     $       (ic_choice.eq.sincos).and.
     $       (bc_choice.eq.reflection_xy_choice))) then

             print '(''the test requires:'')'
             print '(''   - npx=1'')'
             print '(''   - npy=1'')'
             print '(''   - nx=100'')'
             print '(''   - ny=110'')'
             print '(''   - ne=3'')'
             print '(''   - x_min=-10'')'
             print '(''   - x_max=-0.5'')'
             print '(''   - y_min=-10'')'
             print '(''   - y_max= 11'')'
             print '(''   - wave2d model'')'
             print '(''   - sincos i.c.'')'
             print '(''   - reflection i.c.'')'
             stop ''

          end if

        end subroutine check_inputs_small_domain


        subroutine check_inputs()

          if(.not.(
     $       (npx.eq.2).and.
     $       (npy.eq.2).and.
     $       (ntx.eq.104).and.
     $       (nty.eq.114).and.
     $       (ne.eq.3).and.
     $       is_real_validated(x_min,-10.0d0,.true.).and.
     $       is_real_validated(x_max,-0.50d0,.true.).and.
     $       is_real_validated(y_min,-10.0d0,.true.).and.
     $       is_real_validated(y_max, 11.0d0,.true.).and.
     $       (ic_choice.eq.sincos).and.
     $       (bc_choice.eq.reflection_xy_choice))) then

             print '(''the test requires:'')'
             print '(''   - npx=2'')'
             print '(''   - npy=2'')'
             print '(''   - ntx=104'')'
             print '(''   - nty=114'')'
             print '(''   - ne=3'')'
             print '(''   - x_min=-10.0'')'
             print '(''   - x_max=-0.50'')'
             print '(''   - y_min=-10.0'')'
             print '(''   - y_max= 11.0'')'
             print '(''   - wave2d model'')'
             print '(''   - sincos i.c.'')'
             print '(''   - reflection i.c.'')'
             stop ''

          end if

        end subroutine check_inputs

      end program test_field_par
