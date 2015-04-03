      program test_mpi_interface

        use check_data_module, only :
     $     is_int_vector_validated,
     $     is_int_matrix_validated,
     $     is_real_matrix3D_validated

        use mpi

        use mpi_interface_class, only :
     $       mpi_interface

        use mpi_process_class, only :
     $       mpi_process

        use parameters_bf_layer, only :
     $       N_edge_type,
     $       S_edge_type,
     $       E_edge_type,
     $       W_edge_type,
     $       NW_corner_type,
     $       NE_corner_type,
     $       SE_corner_type,
     $       SW_corner_type,
     $       no_overlap

        use parameters_constant, only :
     $       reflection_xy_choice

        use parameters_input, only :
     $       npx,npy,
     $       nx,ny,ne,
     $       bc_choice

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated
        
        integer :: rank
        integer :: ierror


        detailled = .true.
        test_validated = .true.


        ! test_ini
        test_loc = test_ini(detailled)
        test_validated = test_validated.and.test_loc

        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_ini: '',L1)', test_loc
           print '()'
        end if


        ! test ini_for_timeInt
        test_loc = test_ini_for_timeInt(detailled)
        test_validated = test_validated.and.test_loc

        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_ini_for_timeInt: '',L1)', test_loc
           print '()'
        end if


        ! test MPI_ISENDRECV_XDIR
        test_loc = test_mpi_isendrecv_xdir(detailled)
        test_validated = test_validated.and.test_loc

        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_MPI_ISENDRECV_XDIR: '',L1)', test_loc
           print '()'
        end if


        ! test MPI_ISENDRECV_YDIR
        test_loc = test_mpi_isendrecv_ydir(detailled)
        test_validated = test_validated.and.test_loc

        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_MPI_ISENDRECV_YDIR: '',L1)', test_loc
           print '()'
        end if


        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_validated: '',L1)', test_validated
        end if


        call MPI_FINALIZE(ierror)


        contains


        function test_ini(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mpi_process)     :: mpi_process_used
          integer               :: ierror
          integer               :: nb_procs
          integer               :: comm2d
          integer               :: rank
          type(mpi_interface)   :: mpi_interface_used
          integer, dimension(6) :: test_nb_mpi_requests_x
          integer, dimension(6) :: test_nb_mpi_requests_y

          logical, dimension(6) :: test_validated_gathered

          logical :: test_loc
          integer :: i


          test_validated = .true.

          
          ! in this test we will first create the 
          ! mpi interface to exchange data with the
          ! N,S,E,W neighbors.
          ! we will test the number of mpi requests
          ! to be send in each direction
          ! we use the reflection b.c.
          !
          !     ________ ______ _______ 
          !    |       |       |       |
          !    |       |       |       |
          !    |   1 /_|_\ 3 /_|_\ 5   |        
          !    |   . \ | / . \ | / .   |
          !    |  /|\  |  /|\  |  /|\  |
          !    |___|___|__ |___|___|___|
          !    |   |   |   |   |   |   |
          !    |  \|/  |  \|/  |  \|/  |
          !    |   '   |   '   |   '   |
          !    |   0 /_|_\ 2 /_|_\ 4   |
          !    |     \ | /   \ | /     |
          !    |_______|_______|_______|
          !                         
          !
          ! we test whether the number of mpi requests is correct:
          !
          ! | proc_id | nb_mpi_requests_x | nb_mpi_requests_y |
          ! ---------------------------------------------------
          ! |       0 |                 2 |                 2 |
          ! |       1 |                 2 |                 2 |
          ! |       2 | (only_exchange) 0 |                 2 |
          ! |       3 | (only_exchange) 0 |                 2 |
          ! |       4 |                 2 |                 2 |
          ! |       5 |                 2 |                 2 |
          ! ---------------------------------------------------
          !          
          !------------------------------------------------------------

          ! input
          !============================================================
          call mpi_process_used%ini_mpi()

          call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,ierror)
          if(.not.(
     $         (nb_procs.eq.6).and.
     $         (npx.eq.3).and.
     $         (npy.eq.2).and.
     $         (nx.eq.6).and.
     $         (ny.eq.8).and.
     $         (ne.eq.2).and.
     $         (bc_choice.eq.reflection_xy_choice)))then

             print '(''the test requires:'')'
             print '(''   - nb_procs=6'')'
             print '(''   - npx=3'')'
             print '(''   - npy=2'')'
             print '(''   - nx=6'')'
             print '(''   - ny=8'')'
             print '(''   - ne=2'')'
             print '(''   - bc_choice=reflection_xy_choice'')'
             call MPI_FINALIZE(ierror)
             stop ''

          end if

          call mpi_process_used%ini_cartesian_communicator(
     $         comm2d, rank)


          test_nb_mpi_requests_x = [2,2,0,0,2,2]
          test_nb_mpi_requests_y = [2,2,2,2,2,2]


          ! output
          !============================================================
          call mpi_interface_used%ini(comm2d)


          ! validation
          !============================================================
          test_loc = mpi_interface_used%nb_mpi_requests_x.eq.test_nb_mpi_requests_x(rank+1)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_mpi_requests_x('',I2,'') failed'')', rank
          end if

          test_loc = mpi_interface_used%nb_mpi_requests_y.eq.test_nb_mpi_requests_y(rank+1)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test nb_mpi_requests_y('',I2,'') failed'')', rank
          end if

          ! the master processor gathers the test results 
          !------------------------------------------------------------
          call MPI_GATHER(
     $         test_validated,1,MPI_LOGICAL,
     $         test_validated_gathered,1,MPI_LOGICAL,
     $         0, MPI_COMM_WORLD, ierror)

          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

          if(rank.eq.0) then

             test_validated = test_validated_gathered(1)
             do i=2,6
                test_validated = test_validated.and.test_validated_gathered(i)
             end do

          end if          

        end function test_ini


        function test_ini_for_timeInt(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mpi_process)                           :: mpi_process_used
          integer                                     :: ierror
          integer                                     :: nb_procs
          integer                                     :: comm2d
          integer                                     :: rank
          type(mpi_interface)                         :: mpi_interface_used

          integer(ikind), dimension(2)                :: x_borders
          integer(ikind), dimension(2)                :: y_borders
          integer(ikind), dimension(:,:), allocatable :: bc_sections_x
          integer(ikind), dimension(:,:), allocatable :: bc_sections_y

          integer(ikind), dimension(2)                :: test_x_borders
          integer(ikind), dimension(2)                :: test_y_borders
          integer(ikind), dimension(:,:), allocatable :: test_bc_sections_x
          integer(ikind), dimension(:,:), allocatable :: test_bc_sections_y

          logical, dimension(6) :: test_validated_gathered

          logical :: test_loc
          integer :: i


          test_validated = .true.

          
          ! in this test we will first create the 
          ! mpi interface to exchange data with the
          ! N,S,E,W neighbors.
          !
          ! we then initialize the borders of integration
          ! as well as the bc_sections
          !
          ! we use the reflection b.c. to test the
          ! initialization
          !
          ! each tile is 6x8
          !
          !     ________ ______ _______ 
          !    |       |       |       |
          !    |       |       |       |
          !    |   1 /_|_\ 3 /_|_\ 5   |        
          !    |   . \ | / . \ | / .   |
          !    |  /|\  |  /|\  |  /|\  |
          !    |___|___|__ |___|___|___|
          !    |   |   |   |   |   |   |
          !    |  \|/  |  \|/  |  \|/  |
          !    |   '   |   '   |   '   |
          !    |   0 /_|_\ 2 /_|_\ 4   |
          !    |     \ | /   \ | /     |
          !    |_______|_______|_______|
          !                         
          !
          ! we test whether the integration borders are corrects:
          !
          ! -----------------------------------
          ! | proc_id | x_borders | y_borders |
          ! -----------------------------------
          ! |       0 |     [1,4] |     [1,6] |
          ! |       1 |     [1,4] |     [3,8] |
          ! |       2 |     [3,4] |     [1,6] |
          ! |       3 |     [3,4] |     [3,8] |
          ! |       4 |     [3,6] |     [1,6] |
          ! |       5 |     [3,6] |     [3,8] |
          ! -----------------------------------
          !
          ! we test whether the bc_sections are correct
          !
          ! -----------------------------------------------
          ! | proc_id | bc_sections_x |     bc_sections_y |
          ! -----------------------------------------------
          ! |       0 |        W_edge | SW_corner, S_edge |
          ! |       1 |        W_edge | NW_corner, N_edge |
          ! |       2 |             - |            S_edge |
          ! |       3 |             - |            N_edge |
          ! |       4 |        E_edge | S_edge, SE_corner |
          ! |       5 |        E_edge | N_edge, NE_corner |
          ! -----------------------------------------------          
          !          
          !------------------------------------------------------------

          ! input
          !============================================================
          call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,ierror)
          if(.not.(
     $         (nb_procs.eq.6).and.
     $         (npx.eq.3).and.
     $         (npy.eq.2).and.
     $         (nx.eq.6).and.
     $         (ny.eq.8).and.
     $         (ne.eq.2).and.
     $         (bc_choice.eq.reflection_xy_choice)))then

             print '(''the test requires:'')'
             print '(''   - nb_procs=6'')'
             print '(''   - npx=3'')'
             print '(''   - npy=2'')'
             print '(''   - nx=6'')'
             print '(''   - ny=8'')'
             print '(''   - ne=2'')'
             print '(''   - bc_choice=reflection_xy_choice'')'
             call MPI_FINALIZE(ierror)
             stop ''

          end if

          call mpi_process_used%ini_cartesian_communicator(
     $         comm2d, rank)


          select case(rank)
            case(0)
               test_x_borders = [1,4]
               test_y_borders = [1,6]

               allocate(test_bc_sections_x(5,1))
               test_bc_sections_x(:,1) = [W_edge_type,1,3,6,no_overlap]

               allocate(test_bc_sections_y(5,2))
               test_bc_sections_y = reshape((/
     $              SW_corner_type,1,1,no_overlap,no_overlap,
     $              S_edge_type   ,3,1,6         , no_overlap/),
     $              (/5,2/))

            case(1)
               test_x_borders = [1,4]
               test_y_borders = [3,8]

               allocate(test_bc_sections_x(5,1))
               test_bc_sections_x(:,1) = [W_edge_type,1,3,6,no_overlap]

               allocate(test_bc_sections_y(5,2))
               test_bc_sections_y = reshape((/
     $              NW_corner_type,1,7,no_overlap,no_overlap,
     $              N_edge_type   ,3,7,6         , no_overlap/),
     $              (/5,2/))

            case(2)
               test_x_borders = [3,4]
               test_y_borders = [1,6]

               allocate(test_bc_sections_y(5,1))
               test_bc_sections_y(:,1) = [S_edge_type,1,1,6,no_overlap]

            case(3)
               test_x_borders = [3,4]
               test_y_borders = [3,8]

               allocate(test_bc_sections_y(5,1))
               test_bc_sections_y(:,1) = [N_edge_type,1,7,6,no_overlap]

            case(4)
               test_x_borders = [3,6]
               test_y_borders = [1,6]

               allocate(test_bc_sections_x(5,1))
               test_bc_sections_x(:,1) = [E_edge_type,5,3,6,no_overlap]

               allocate(test_bc_sections_y(5,2))
               test_bc_sections_y = reshape((/
     $              S_edge_type   ,1,1,4         ,no_overlap,
     $              SE_corner_type,5,1,no_overlap,no_overlap/),
     $              (/5,2/))

            case(5)
               test_x_borders = [3,6]
               test_y_borders = [3,8]

               allocate(test_bc_sections_x(5,1))
               test_bc_sections_x(:,1) = [E_edge_type,5,3,6,no_overlap]

               allocate(test_bc_sections_y(5,2))
               test_bc_sections_y = reshape((/
     $              N_edge_type   ,1,7,4         ,no_overlap,
     $              NE_corner_type,5,7,no_overlap,no_overlap/),
     $              (/5,2/))

          end select


          ! output
          !============================================================
          call mpi_interface_used%ini(comm2d)
          call mpi_interface_used%ini_for_timeInt(
     $         comm2d,
     $         x_borders,
     $         y_borders,
     $         bc_sections_x,
     $         bc_sections_y)


          ! validation
          !============================================================
          ! x_borders
          !------------------------------------------------------------
          test_loc = is_int_vector_validated(x_borders,test_x_borders,detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test x_borders('',I2,'') failed'')', rank
          end if

          ! y_borders
          !------------------------------------------------------------
          test_loc = is_int_vector_validated(y_borders,test_y_borders,detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test y_borders('',I2,'') failed'')', rank
          end if

          ! bc_sections_x
          !------------------------------------------------------------
          if(allocated(bc_sections_x)) then
             test_loc = is_int_matrix_validated(
     $            bc_sections_x,
     $            test_bc_sections_x,
     $            detailled)
          else
             test_loc = .not.allocated(test_bc_sections_x)
          end if
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test bc_sections_x('',I2,'') failed'')', rank
          end if

          ! bc_sections_y
          !------------------------------------------------------------
          if(allocated(bc_sections_y)) then
             test_loc = is_int_matrix_validated(
     $            bc_sections_y,
     $            test_bc_sections_y,
     $            detailled)
          else
             test_loc = .not.allocated(test_bc_sections_y)
          end if
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test bc_sections_y('',I2,'') failed'')', rank
          end if          


          ! the master processor gathers the test results 
          !------------------------------------------------------------
          call MPI_GATHER(
     $         test_validated,1,MPI_LOGICAL,
     $         test_validated_gathered,1,MPI_LOGICAL,
     $         0, MPI_COMM_WORLD, ierror)

          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

          if(rank.eq.0) then

             test_validated = test_validated_gathered(1)
             do i=2,6
                test_validated = test_validated.and.test_validated_gathered(i)
             end do

          end if          

        end function test_ini_for_timeInt


        function test_mpi_isendrecv_xdir(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mpi_process)                           :: mpi_process_used
          integer                                     :: ierror
          integer                                     :: nb_procs
          integer                                     :: comm2d
          integer                                     :: rank
          type(mpi_interface)                         :: mpi_interface_used

          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(nx,ny,ne) :: test_nodes

          logical, dimension(6) :: test_validated_gathered

          logical :: test_loc
          integer :: i,j,k


          test_validated = .true.

          
          ! in this test we will first create the 
          ! mpi interface to exchange data with the
          ! N,S,E,W neighbors.
          !
          ! we then initialize the borders of integration
          ! as well as the bc_sections
          !
          ! we use the reflection b.c. to test the
          ! exchnage of data in the x-direction
          !
          ! each tile is 6x8
          !
          !     ________ ______ _______ 
          !    |       |       |       |
          !    |       |       |       |
          !    |   1 /_|_\ 3 /_|_\ 5   |        
          !    |   . \ | / . \ | / .   |
          !    |  /|\  |  /|\  |  /|\  |
          !    |___|___|__ |___|___|___|
          !    |   |   |   |   |   |   |
          !    |  \|/  |  \|/  |  \|/  |
          !    |   '   |   '   |   '   |
          !    |   0 /_|_\ 2 /_|_\ 4   |
          !    |     \ | /   \ | /     |
          !    |_______|_______|_______|
          !                         
          !
          ! we test whether the exchange of data is correct in the 
          ! x-dir
          !
          !  -------------------------------------------
          !  | 1 1 1 1 1 1 | 3 3 3 3 3 3 | 5 5 5 5 5 5 |
          !  | 1 1 1 1 1 1 | 3 3 3 3 3 3 | 5 5 5 5 5 5 |
          !  | 1 1 1 1 3 3 | 1 1 3 3 5 5 | 3 3 5 5 5 5 |
          !  | 1 1 1 1 3 3 | 1 1 3 3 5 5 | 3 3 5 5 5 5 |
          !  | 1 1 1 1 3 3 | 1 1 3 3 5 5 | 3 3 5 5 5 5 |
          !  | 1 1 1 1 3 3 | 1 1 3 3 5 5 | 3 3 5 5 5 5 |
          !  | 1 1 1 1 1 1 | 3 3 3 3 3 3 | 5 5 5 5 5 5 |
          !  | 1 1 1 1 1 1 | 3 3 3 3 3 3 | 5 5 5 5 5 5 |
          !  -------------------------------------------
          !  | 0 0 0 0 0 0 | 2 2 2 2 2 2 | 4 4 4 4 4 4 |             
          !  | 0 0 0 0 0 0 | 2 2 2 2 2 2 | 4 4 4 4 4 4 |           
          !  | 0 0 0 0 2 2 | 0 0 2 2 4 4 | 2 2 4 4 4 4 |
          !  | 0 0 0 0 2 2 | 0 0 2 2 4 4 | 2 2 4 4 4 4 |
          !  | 0 0 0 0 2 2 | 0 0 2 2 4 4 | 2 2 4 4 4 4 |
          !  | 0 0 0 0 2 2 | 0 0 2 2 4 4 | 2 2 4 4 4 4 |
          !  | 0 0 0 0 0 0 | 2 2 2 2 2 2 | 4 4 4 4 4 4 |
          !  | 0 0 0 0 0 0 | 2 2 2 2 2 2 | 4 4 4 4 4 4 |
          !  -------------------------------------------
          !          
          !------------------------------------------------------------

          ! input
          !============================================================
          call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,ierror)
          if(.not.(
     $         (nb_procs.eq.6).and.
     $         (npx.eq.3).and.
     $         (npy.eq.2).and.
     $         (nx.eq.6).and.
     $         (ny.eq.8).and.
     $         (ne.eq.2).and.
     $         (bc_choice.eq.reflection_xy_choice)))then

             print '(''the test requires:'')'
             print '(''   - nb_procs=6'')'
             print '(''   - npx=3'')'
             print '(''   - npy=2'')'
             print '(''   - nx=6'')'
             print '(''   - ny=8'')'
             print '(''   - ne=2'')'
             print '(''   - bc_choice=reflection_xy_choice'')'
             call MPI_FINALIZE(ierror)
             stop ''

          end if

          call mpi_process_used%ini_cartesian_communicator(
     $         comm2d, rank)


          select case(rank)

            case(0)
               nodes = reshape((/(((0,i=1,6),j=1,8),k=1,2)/),(/6,8,2/))
               test_nodes = nodes
               test_nodes(5:6,3:6,:) = reshape((/
     $              (((2,i=1,2),j=1,4),k=1,2)/), (/2,4,2/))

            case(1)
               nodes = reshape((/(((1,i=1,6),j=1,8),k=1,2)/),(/6,8,2/))
               test_nodes = nodes
               test_nodes(5:6,3:6,:) = reshape((/
     $              (((3,i=1,2),j=1,4),k=1,2)/), (/2,4,2/))

            case(2)
               nodes = reshape((/(((2,i=1,6),j=1,8),k=1,2)/),(/6,8,2/))
               test_nodes = nodes
               test_nodes(1:2,3:6,:) = reshape((/
     $              (((0,i=1,2),j=1,4),k=1,2)/), (/2,4,2/))
               test_nodes(5:6,3:6,:) = reshape((/
     $              (((4,i=1,2),j=1,4),k=1,2)/), (/2,4,2/))

            case(3)
               nodes = reshape((/(((3,i=1,6),j=1,8),k=1,2)/),(/6,8,2/))
               test_nodes = nodes
               test_nodes(1:2,3:6,:) = reshape((/
     $              (((1,i=1,2),j=1,4),k=1,2)/), (/2,4,2/))
               test_nodes(5:6,3:6,:) = reshape((/
     $              (((5,i=1,2),j=1,4),k=1,2)/), (/2,4,2/))
               
            case(4)
               nodes = reshape((/(((4,i=1,6),j=1,8),k=1,2)/),(/6,8,2/))
               test_nodes = nodes
               test_nodes(1:2,3:6,:) = reshape((/
     $              (((2,i=1,2),j=1,4),k=1,2)/), (/2,4,2/))
               
            case(5)
               nodes = reshape((/(((5,i=1,6),j=1,8),k=1,2)/),(/6,8,2/))
               test_nodes = nodes
               test_nodes(1:2,3:6,:) = reshape((/
     $              (((3,i=1,2),j=1,4),k=1,2)/), (/2,4,2/))
               
          end select


          ! output
          !============================================================
          call mpi_interface_used%ini(comm2d)
          call mpi_interface_used%MPI_ISENDRECV_XDIR(comm2d,nodes)
          call mpi_interface_used%MPI_WAITALL_XDIR()

          ! validation
          !============================================================
          test_loc = is_real_matrix3D_validated(
     $         nodes,
     $         test_nodes,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test send_recv('',I2,'') failed'')', rank
          end if


          ! the master processor gathers the test results 
          !------------------------------------------------------------
          call MPI_GATHER(
     $         test_validated,1,MPI_LOGICAL,
     $         test_validated_gathered,1,MPI_LOGICAL,
     $         0, MPI_COMM_WORLD, ierror)

          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

          if(rank.eq.0) then

             test_validated = test_validated_gathered(1)
             do i=2,6
                test_validated = test_validated.and.test_validated_gathered(i)
             end do

          end if          

        end function test_mpi_isendrecv_xdir


        function test_mpi_isendrecv_ydir(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mpi_process)                           :: mpi_process_used
          integer                                     :: ierror
          integer                                     :: nb_procs
          integer                                     :: comm2d
          integer                                     :: rank
          type(mpi_interface)                         :: mpi_interface_used

          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(nx,ny,ne) :: test_nodes

          logical, dimension(6) :: test_validated_gathered

          logical :: test_loc
          integer :: i,j,k


          test_validated = .true.

          
          ! in this test we will first create the 
          ! mpi interface to exchange data with the
          ! N,S,E,W neighbors.
          !
          ! we then initialize the borders of integration
          ! as well as the bc_sections
          !
          ! we use the reflection b.c. to test the
          ! exchnage of data in the y-direction
          !
          ! each tile is 6x8
          !
          !     ________ ______ _______ 
          !    |       |       |       |
          !    |       |       |       |
          !    |   1 /_|_\ 3 /_|_\ 5   |        
          !    |   . \ | / . \ | / .   |
          !    |  /|\  |  /|\  |  /|\  |
          !    |___|___|__ |___|___|___|
          !    |   |   |   |   |   |   |
          !    |  \|/  |  \|/  |  \|/  |
          !    |   '   |   '   |   '   |
          !    |   0 /_|_\ 2 /_|_\ 4   |
          !    |     \ | /   \ | /     |
          !    |_______|_______|_______|
          !                         
          !
          ! we test whether the exchange of data is correct in the 
          ! y-dir
          !
          !  -------------------------------------------
          !  | 1 1 1 1 1 1 | 3 3 3 3 3 3 | 5 5 5 5 5 5 |
          !  | 1 1 1 1 1 1 | 3 3 3 3 3 3 | 5 5 5 5 5 5 |
          !  | 1 1 1 1 1 1 | 3 3 3 3 3 3 | 5 5 5 5 5 5 |
          !  | 1 1 1 1 1 1 | 3 3 3 3 3 3 | 5 5 5 5 5 5 |
          !  | 1 1 1 1 1 1 | 3 3 3 3 3 3 | 5 5 5 5 5 5 |
          !  | 1 1 1 1 1 1 | 3 3 3 3 3 3 | 5 5 5 5 5 5 |
          !  | 0 0 0 0 0 0 | 2 2 2 2 2 2 | 4 4 4 4 4 4 |
          !  | 0 0 0 0 0 0 | 2 2 2 2 2 2 | 4 4 4 4 4 4 |
          !  -------------------------------------------
          !  | 1 1 1 1 1 1 | 3 3 3 3 3 3 | 5 5 5 5 5 5 |
          !  | 1 1 1 1 1 1 | 3 3 3 3 3 3 | 5 5 5 5 5 5 |
          !  | 0 0 0 0 0 0 | 2 2 2 2 2 2 | 4 4 4 4 4 4 |
          !  | 0 0 0 0 0 0 | 2 2 2 2 2 2 | 4 4 4 4 4 4 |
          !  | 0 0 0 0 0 0 | 2 2 2 2 2 2 | 4 4 4 4 4 4 |
          !  | 0 0 0 0 0 0 | 2 2 2 2 2 2 | 4 4 4 4 4 4 |
          !  | 0 0 0 0 0 0 | 2 2 2 2 2 2 | 4 4 4 4 4 4 |
          !  | 0 0 0 0 0 0 | 2 2 2 2 2 2 | 4 4 4 4 4 4 |
          !  -------------------------------------------
          !          
          !------------------------------------------------------------

          ! input
          !============================================================
          call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,ierror)
          if(.not.(
     $         (nb_procs.eq.6).and.
     $         (npx.eq.3).and.
     $         (npy.eq.2).and.
     $         (nx.eq.6).and.
     $         (ny.eq.8).and.
     $         (ne.eq.2).and.
     $         (bc_choice.eq.reflection_xy_choice)))then

             print '(''the test requires:'')'
             print '(''   - nb_procs=6'')'
             print '(''   - npx=3'')'
             print '(''   - npy=2'')'
             print '(''   - nx=6'')'
             print '(''   - ny=8'')'
             print '(''   - ne=2'')'
             print '(''   - bc_choice=reflection_xy_choice'')'
             call MPI_FINALIZE(ierror)
             stop ''

          end if

          call mpi_process_used%ini_cartesian_communicator(
     $         comm2d, rank)


          select case(rank)

            case(0)
               nodes = reshape((/(((0,i=1,6),j=1,8),k=1,2)/),(/6,8,2/))
               test_nodes = nodes
               test_nodes(1:6,7:8,:) = reshape((/
     $              (((1,i=1,6),j=7,8),k=1,2)/), (/6,2,2/))

            case(1)
               nodes = reshape((/(((1,i=1,6),j=1,8),k=1,2)/),(/6,8,2/))
               test_nodes = nodes
               test_nodes(1:6,1:2,:) = reshape((/
     $              (((0,i=1,6),j=1,2),k=1,2)/), (/6,2,2/))

            case(2)
               nodes = reshape((/(((2,i=1,6),j=1,8),k=1,2)/),(/6,8,2/))
               test_nodes = nodes
               test_nodes(1:6,7:8,:) = reshape((/
     $              (((3,i=1,6),j=1,2),k=1,2)/), (/6,2,2/))

            case(3)
               nodes = reshape((/(((3,i=1,6),j=1,8),k=1,2)/),(/6,8,2/))
               test_nodes = nodes
               test_nodes(1:6,1:2,:) = reshape((/
     $              (((2,i=1,6),j=1,2),k=1,2)/), (/6,2,2/))
               
            case(4)
               nodes = reshape((/(((4,i=1,6),j=1,8),k=1,2)/),(/6,8,2/))
               test_nodes = nodes
               test_nodes(1:6,7:8,:) = reshape((/
     $              (((5,i=1,6),j=1,2),k=1,2)/), (/6,2,2/))
               
            case(5)
               nodes = reshape((/(((5,i=1,6),j=1,8),k=1,2)/),(/6,8,2/))
               test_nodes = nodes
               test_nodes(1:6,1:2,:) = reshape((/
     $              (((4,i=1,6),j=1,2),k=1,2)/), (/6,2,2/))
               
          end select


          ! output
          !============================================================
          call mpi_interface_used%ini(comm2d)
          call mpi_interface_used%MPI_ISENDRECV_YDIR(comm2d,nodes)
          call mpi_interface_used%MPI_WAITALL_YDIR()

          ! validation
          !============================================================
          test_loc = is_real_matrix3D_validated(
     $         nodes,
     $         test_nodes,
     $         detailled)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test send_recv('',I2,'') failed'')', rank
          end if


          ! the master processor gathers the test results 
          !------------------------------------------------------------
          call MPI_GATHER(
     $         test_validated,1,MPI_LOGICAL,
     $         test_validated_gathered,1,MPI_LOGICAL,
     $         0, MPI_COMM_WORLD, ierror)

          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

          if(rank.eq.0) then

             test_validated = test_validated_gathered(1)
             do i=2,6
                test_validated = test_validated.and.test_validated_gathered(i)
             end do

          end if          

        end function test_mpi_isendrecv_ydir


      end program test_mpi_interface
