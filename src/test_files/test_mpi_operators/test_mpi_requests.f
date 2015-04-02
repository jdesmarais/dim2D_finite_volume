      program test_mpi_requests

        use check_data_module, only :
     $     is_real_matrix3D_validated

        use mpi

        use mpi_mg_bc_class, only :
     $       mpi_mg_bc

        use mpi_process_class, only :
     $       mpi_process

        use mpi_requests_module, only :
     $       create_requests_for_one_direction,
     $       only_exchange_twice

        use parameters_constant, only :
     $       N,S

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


        test_loc = test_create_requests_for_one_direction(detailled)
        test_validated = test_validated.and.test_loc
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_create_requests_for_one_direction: '',L1)', test_loc
           print '()'
        end if


        test_loc = test_only_exchange_twice(detailled)
        test_validated = test_validated.and.test_loc
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_only_exchange_twice: '',L1)', test_loc
           print '()'
        end if


        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_validated: '',L1)', test_validated
        end if


        call MPI_FINALIZE(ierror)

        
        contains


        function test_create_requests_for_one_direction(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer :: nb_procs
          integer :: comm2d
          integer :: rank
          integer :: ngh_rank

          type(mpi_process)                :: mpi_process_used
          type(mpi_mg_bc)                  :: mpi_mg_bc_used
          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(nx,ny,ne) :: test_data_exchanged

          integer, dimension(2)                 :: mpi_requests
          integer, dimension(MPI_STATUS_SIZE,2) :: status

          logical, dimension(6) :: test_loc_gathered
          logical               :: test_loc

          integer :: i


          test_validated = .true.


          ! input
          call mpi_process_used%ini_mpi()
          call mpi_process_used%ini_cartesian_communicator(
     $         comm2d, rank)
          call mpi_mg_bc_used%ini(comm2d)

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
             call MPI_FINALIZE(ierror)
             stop ''

          end if


          ! data for validation of the requests
          if(mod(rank,2).eq.0) then
             ngh_rank = rank+1
          else
             ngh_rank = rank-1
          end if
          test_data_exchanged = get_data_N_exchanges(rank,ngh_rank)


          ! output
          nodes = initialize_data_exchanged(rank)
          mpi_requests = create_requests_for_one_direction(
     $         mpi_mg_bc_used, comm2d, rank, nodes, N)

          call MPI_WAITALL(2, mpi_requests, status, ierror)


          ! validation
          test_loc = is_real_matrix3D_validated(
     $         nodes,
     $         test_data_exchanged,
     $         .false.)
          if(detailled.and.(.not.test_loc)) then
             print '(''test rank '',I2,''failed'')',rank
          end if


          ! the master processor gathered the results of the
          ! local tests
          call MPI_GATHER(
     $         test_loc,1,MPI_LOGICAL,
     $         test_loc_gathered,1,MPI_LOGICAL,
     $         0, MPI_COMM_WORLD, ierror)

          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

          if(rank.eq.0) then

             test_loc = test_loc_gathered(1)
             do i=2,6
                test_loc = test_loc.and.test_loc_gathered(i)
             end do
             if(detailled.and.(.not.test_loc)) then
                print '(''test com_rank failed'')'
             end if
             test_validated = test_validated.and.test_loc

          else
             test_validated = test_validated.and.test_loc
          end if

        end function test_create_requests_for_one_direction


        function test_only_exchange_twice(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated

          integer :: nb_procs
          integer :: comm2d
          integer :: rank
          integer :: ngh_rank

          type(mpi_process)                :: mpi_process_used
          type(mpi_mg_bc)                  :: mpi_mg_bc_used
          real(rkind), dimension(nx,ny,ne) :: nodes
          real(rkind), dimension(nx,ny,ne) :: test_data_exchanged

          logical, dimension(6) :: test_loc_gathered
          logical               :: test_loc

          integer :: i


          test_validated = .true.


          ! input
          call mpi_process_used%ini_cartesian_communicator(
     $         comm2d, rank)
          call mpi_mg_bc_used%ini(comm2d)

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
             call MPI_FINALIZE(ierror)
             stop ''

          end if


          ! data for validation of the requests
          if(mod(rank,2).eq.0) then
             ngh_rank = rank+1
          else
             ngh_rank = rank-1
          end if
          test_data_exchanged = get_data_NS_exchanges(rank,ngh_rank)


          ! output
          nodes = initialize_data_exchanged(rank)
          call only_exchange_twice(
     $         mpi_mg_bc_used, comm2d, rank, nodes, nb_procs, [N,S])


          ! validation
          test_loc = is_real_matrix3D_validated(
     $         nodes,
     $         test_data_exchanged,
     $         .false.)
          if(detailled.and.(.not.test_loc)) then
             print '(''test rank '',I2,''failed'')',rank
          end if


          ! the master processor gathered the results of the
          ! local tests
          call MPI_GATHER(
     $         test_loc,1,MPI_LOGICAL,
     $         test_loc_gathered,1,MPI_LOGICAL,
     $         0, MPI_COMM_WORLD, ierror)

          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)

          if(rank.eq.0) then

             test_loc = test_loc_gathered(1)
             do i=2,6
                test_loc = test_loc.and.test_loc_gathered(i)
             end do
             if(detailled.and.(.not.test_loc)) then
                print '(''test com_rank failed'')'
             end if
             test_validated = test_validated.and.test_loc

          else
             test_validated = test_validated.and.test_loc
          end if

        end function test_only_exchange_twice


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


        function get_data_N_exchanges(rank,other_rank)
     $     result(data_ex)

          implicit none

          integer          , intent(in) :: rank
          integer          , intent(in) :: other_rank
          real(rkind), dimension(6,8,2) :: data_ex

          integer :: i,j

          data_ex = initialize_data_exchanged(rank)

          data_ex(1:6,7:8,1) = reshape((/
     $         ((other_rank*10,i=1,6),j=1,2)/),
     $         (/6,2/))

          data_ex(1:6,7:8,2) = reshape((/
     $         ((other_rank*10+1,i=1,6),j=1,2)/),
     $         (/6,2/))

        end function get_data_N_exchanges


        function get_data_NS_exchanges(rank,other_rank)
     $     result(data_ex)

          implicit none

          integer          , intent(in) :: rank
          integer          , intent(in) :: other_rank
          real(rkind), dimension(6,8,2) :: data_ex

          integer :: i,j

          data_ex = get_data_N_exchanges(rank,other_rank)

          data_ex(1:6,1:2,1) = reshape((/
     $         ((other_rank*10,i=1,6),j=1,2)/),
     $         (/6,2/))

          data_ex(1:6,1:2,2) = reshape((/
     $         ((other_rank*10+1,i=1,6),j=1,2)/),
     $         (/6,2/))

        end function get_data_NS_exchanges

      end program test_mpi_requests
