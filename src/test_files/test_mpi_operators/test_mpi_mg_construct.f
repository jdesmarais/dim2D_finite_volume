      program test_mpi_mg_construct

        use check_data_module, only :
     $       is_real_matrix3D_validated

        use mpi

        use mpi_mg_bc_class, only :
     $       mpi_mg_bc

        use mpi_mg_construct, only :
     $       update_mpi_derived_types

        use mpi_mg_ini_bc_proc, only :
     $       ini_bc_procedures

        use mpi_process_class, only :
     $       mpi_process

        use mpi_requests_module, only :
     $       create_requests_for_one_direction,
     $       only_exchange_twice

        use parameters_constant, only :
     $       N,S,E,W,
     $       only_exchange_proc

        use parameters_input, only :
     $       nx,ny,ne,
     $       npx,npy

        use parameters_kind, only :
     $       rkind

        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated
        
        integer :: rank
        integer :: ierror


        detailled = .true.
        test_validated = .true.


        test_loc = test_update_mpi_derived_types(detailled)
        test_validated = test_validated.and.test_loc

        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_update_mpi_derived_types: '',L1)', test_loc
           print '()'
        end if


        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_validated: '',L1)', test_validated
        end if


        call MPI_FINALIZE(ierror)


        contains


        function test_update_mpi_derived_types(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mpi_process) :: mpi_process_used
          integer           :: ierror
          integer           :: nb_procs
          integer           :: comm2d
          integer           :: rank
          type(mpi_mg_bc)   :: mpi_mg_bc_used

          integer               :: proc_x_choice
          integer               :: proc_y_choice
          integer, dimension(2) :: exchange_id

          integer, dimension(2)                 :: mpi_requests
          integer, dimension(MPI_STATUS_SIZE,2) :: status

          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(nx,ny,ne) :: test_nodes

          logical, dimension(6) :: test_validated_gathered

          logical :: test_loc
          integer :: i

          test_validated = .true.


          ! in this test we will first create the 
          ! mpi derived to exchange data with the
          ! N,S,E,W neighbors.
          ! then, because of the use of periodic
          ! b.c. some mpi derived types can be combined
          ! into one derived type to send it to its
          ! neighbor
          !     ___|____ __|___ ___|___ 
          !    |  \|/  |  \|/  |  \|/  |
          !    |   '   |   '   |   '   |
          !   _|_\ 1 /_|_\ 3 /_|_\ 5 /_|_       
          !    | / . \ | / . \ | / . \ |
          !    |  /|\  |  /|\  |  /|\  |
          !    |___|___|__ |___|___|___|
          !    |   |   |   |   |   |   |
          !    |  \|/  |  \|/  |  \|/  |
          !    |   '   |   '   |   '   |
          !   _|_\ 0 /_|_\ 2 /_|_\ 4 /_|_
          !    | / . \ | / . \ | / . \ |
          !    |__/|\__|__/|\__|__/|\__|
          !        |       |       |
          !
          ! we test whether the merge of mpi_derived types
          ! work by using them to exchange data between the
          ! processors. If the data present in each processor
          ! after the exchange match the predicted data, the
          ! test is validated
          !------------------------------------------------------------

          ! inputs
          !============================================================
          call mpi_process_used%ini_mpi()

          call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,ierror)
          if(.not.(
     $         (nb_procs.eq.6).and.
     $         (npx.eq.3).and.
     $         (npy.eq.2).and.
     $         (nx.eq.6).and.
     $         (ny.eq.8).and.
     $         (ne.eq.2)))then

             print '(''the test requires:'')'
             print '(''   - nb_procs=6'')'
             print '(''   - npx=3'')'
             print '(''   - npy=2'')'
             print '(''   - nx=6'')'
             print '(''   - ny=8'')'
             print '(''   - ne=2'')'
             print '(''   - bc_choice=periodic_xy_choice'')'
             call MPI_FINALIZE(ierror)
             stop ''

          end if

          call mpi_process_used%ini_cartesian_communicator(
     $         comm2d, rank)

          call mpi_mg_bc_used%ini(comm2d)

          call ini_bc_procedures(
     $         comm2d, proc_x_choice, proc_y_choice, exchange_id)


          ! output
          !============================================================
          call update_mpi_derived_types(
     $         mpi_mg_bc_used%com_recv,
     $         mpi_mg_bc_used%com_send,
     $         mpi_mg_bc_used%com_rank,
     $         proc_x_choice,
     $         proc_y_choice)


          ! validation
          !============================================================

          ! initialize the nodes before the exchange
          !------------------------------------------------------------
          nodes = initialize_data_exchanged(rank)


          ! test the exchange along x
          !------------------------------------------------------------
          !exchange the data along the x-direction
          if(proc_x_choice.eq.only_exchange_proc) then
             call only_exchange_twice(
     $            mpi_mg_bc_used, comm2d, rank, nodes, nb_procs, [E,W])
          else
             print '(''test proc_x_choice('',I2,'') failed'')',rank
          end if


          ! test the nodes after the exchange along x
          test_nodes = get_test_data_after_x_exchange(
     $         rank,
     $         mpi_mg_bc_used%com_rank)

          ! test the nodes
          test_loc = is_real_matrix3D_validated(
     $         nodes,
     $         test_nodes,
     $         .false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test exchange_x('',I2,'') failed'')', rank
          end if


          ! test the exchange along y
          !------------------------------------------------------------
          !exchange the data along the y-direction
          if(proc_y_choice.eq.only_exchange_proc) then
             mpi_requests = create_requests_for_one_direction(
     $            mpi_mg_bc_used, comm2d, rank, nodes, N)
             call MPI_WAITALL(2,mpi_requests,status,ierror)
          else
             print '(''test proc_y_choice('',I2,'') failed'')',rank
          end if

          ! test the nodes after the exchange along y
          test_nodes = get_test_data_after_y_exchange(
     $         rank,
     $         mpi_mg_bc_used%com_rank)

          ! test the nodes
          test_loc = is_real_matrix3D_validated(
     $         nodes,
     $         test_nodes,
     $         .false.)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test exchange_y('',I2,'') failed'')', rank
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

        end function test_update_mpi_derived_types


        function initialize_data_exchanged(rank)
     $     result(data_ini)

          implicit none

          integer          , intent(in) :: rank
          real(rkind), dimension(6,8,2) :: data_ini

          integer :: i,j

          data_ini(:,:,1) = reshape((/
     $         ((rank*10,i=1,6),j=1,8)/),
     $         (/6,8/))

          data_ini(:,:,2) = reshape((/
     $         ((rank*10+1,i=1,6),j=1,8)/),
     $         (/6,8/))

        end function initialize_data_exchanged


        function get_test_data_after_x_exchange(rank,com_rank)
     $     result(data_ex)

          implicit none

          integer                      , intent(in) :: rank
          integer    , dimension(4)    , intent(in) :: com_rank
          real(rkind), dimension(6,8,2)             :: data_ex

          integer :: i,j

          data_ex = initialize_data_exchanged(rank)

          data_ex(1:2,3:6,1) = reshape((/
     $         ((com_rank(W)*10,i=1,2),j=3,6)/),
     $         (/2,4/))

          data_ex(5:6,3:6,1) = reshape((/
     $         ((com_rank(E)*10,i=5,6),j=3,6)/),
     $         (/2,4/))

          data_ex(1:2,3:6,2) = reshape((/
     $         ((com_rank(W)*10+1,i=1,2),j=3,6)/),
     $         (/2,4/))

          data_ex(5:6,3:6,2) = reshape((/
     $         ((com_rank(E)*10+1,i=5,6),j=3,6)/),
     $         (/2,4/))

        end function get_test_data_after_x_exchange


        function get_test_data_after_y_exchange(rank,com_rank)
     $     result(data_ex)

          implicit none

          integer                      , intent(in) :: rank
          integer    , dimension(4)    , intent(in) :: com_rank
          real(rkind), dimension(6,8,2)             :: data_ex

          integer :: i,j
          integer :: i_add

          data_ex = get_test_data_after_x_exchange(rank,com_rank)

          
          if(mod(rank,2).eq.0) then
             i_add =  1
          else
             i_add = -1
          end if


          ! data_ex(:,:,1)
          !============================================================
          !south
          !------------------------------------------------------------
          data_ex(1:2,1:2,1) = reshape((/
     $         (((com_rank(W)+i_add)*10,i=1,2),j=1,2)/),
     $         (/2,2/))

          data_ex(3:4,1:2,1) = reshape((/
     $         ((com_rank(S)*10,i=1,2),j=1,2)/),
     $         (/2,2/))

          data_ex(5:6,1:2,1) = reshape((/
     $         (((com_rank(E)+i_add)*10,i=1,2),j=1,2)/),
     $         (/2,2/))

          !north
          !------------------------------------------------------------
          data_ex(1:2,7:8,1) = reshape((/
     $         (((com_rank(W)+i_add)*10,i=1,2),j=1,2)/),
     $         (/2,2/))

          data_ex(3:4,7:8,1) = reshape((/
     $         ((com_rank(N)*10,i=1,2),j=1,2)/),
     $         (/2,2/))

          data_ex(5:6,7:8,1) = reshape((/
     $         (((com_rank(E)+i_add)*10,i=1,2),j=1,2)/),
     $         (/2,2/))

          
          ! data_ex(:,:,2)
          !============================================================
          !south
          !------------------------------------------------------------
          data_ex(1:2,1:2,2) = reshape((/
     $         (((com_rank(W)+i_add)*10+1,i=1,2),j=1,2)/),
     $         (/2,2/))

          data_ex(3:4,1:2,2) = reshape((/
     $         ((com_rank(S)*10+1,i=1,2),j=1,2)/),
     $         (/2,2/))

          data_ex(5:6,1:2,2) = reshape((/
     $         (((com_rank(E)+i_add)*10+1,i=1,2),j=1,2)/),
     $         (/2,2/))

          !north
          !------------------------------------------------------------
          data_ex(1:2,7:8,2) = reshape((/
     $         (((com_rank(W)+i_add)*10+1,i=1,2),j=1,2)/),
     $         (/2,2/))

          data_ex(3:4,7:8,2) = reshape((/
     $         ((com_rank(N)*10+1,i=1,2),j=1,2)/),
     $         (/2,2/))

          data_ex(5:6,7:8,2) = reshape((/
     $         (((com_rank(E)+i_add)*10+1,i=1,2),j=1,2)/),
     $         (/2,2/))

        end function get_test_data_after_y_exchange

      end program test_mpi_mg_construct

