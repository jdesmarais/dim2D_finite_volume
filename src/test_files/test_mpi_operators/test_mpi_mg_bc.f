      program test_mpi_mg_bc

        use check_data_module, only :
     $     is_int_vector_validated,
     $     is_real_matrix3D_validated

        use mpi

        use mpi_mg_bc_class, only :
     $       mpi_mg_bc

        use mpi_process_class, only :
     $       mpi_process

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       npx,npy,
     $       nx,ny,ne

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


        test_loc = test_ini(detailled)
        test_validated = test_validated.and.test_loc
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_ini: '',L1)', test_loc
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


          type(mpi_process) :: mpi_process_used

          integer :: nb_procs
          integer :: ierror
          integer :: comm2d
          integer :: rank
          integer :: ngh_rank

          integer, dimension(4,6) :: test_com_rank
          logical, dimension(6)   :: test_loc_gathered

          integer, parameter :: tile_nx = 6
          integer, parameter :: tile_ny = 8
          integer, parameter :: tile_ne = 2
          
          real(rkind), dimension(tile_nx,tile_ny,tile_ne)   :: data_exchanged
          real(rkind), dimension(tile_nx,tile_ny,tile_ne,4) :: test_data_exchanged
          integer, parameter :: tag1=1111
          integer, parameter :: tag2=200

          type(mpi_mg_bc) :: mpi_mg_bc_used

          logical :: test_loc
          integer :: i
          integer :: k


          test_validated = .true.


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
             call MPI_FINALIZE(ierror)
             stop ''

          end if

          call mpi_process_used%ini_cartesian_communicator(
     $         comm2d, rank)

          ! data for validation of com_rank
          test_com_rank(:,1) = [1,1,2,4]
          test_com_rank(:,2) = [0,0,3,5]
          test_com_rank(:,3) = [3,3,4,0]
          test_com_rank(:,4) = [2,2,5,1]
          test_com_rank(:,5) = [5,5,0,2]
          test_com_rank(:,6) = [4,4,1,3]

          !data for validation of com_recv and com_send
          if(mod(rank,2).eq.0) then
             ngh_rank = rank+1
          else
             ngh_rank = rank-1
          end if
          test_data_exchanged(:,:,:,N) = get_data_N_exchanges(rank,ngh_rank)
          test_data_exchanged(:,:,:,S) = get_data_S_exchanges(rank,ngh_rank)
          test_data_exchanged(:,:,:,E) = get_data_E_exchanges(rank,ngh_rank)
          test_data_exchanged(:,:,:,W) = get_data_W_exchanges(rank,ngh_rank)


          ! output
          !============================================================
          call mpi_mg_bc_used%ini(comm2d)


          ! validation
          !============================================================

          ! validation of the com_rank
          !------------------------------------------------------------
          test_loc = is_int_vector_validated(
     $         mpi_mg_bc_used%com_rank,
     $         test_com_rank(:,rank+1),
     $         detailled)

          ! gather the test_loc of all the processors
          ! to the processor of rank 0
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

          ! validation of com_recv and com_send
          !------------------------------------------------------------
          do k=1,4

             data_exchanged = initialize_data_exchanged(rank)

             call MPI_SENDRECV (
     $            data_exchanged, 1, mpi_mg_bc_used%com_send(k),
     $            ngh_rank, tag1,
     $            data_exchanged, 1, mpi_mg_bc_used%com_recv(k),
     $            ngh_rank, tag1,
     $            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)

             test_loc = is_real_matrix3D_validated(
     $            data_exchanged,
     $            test_data_exchanged(:,:,:,k),
     $            detailled)

             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test com_recv/send('',I2,'') failed'')', k
             end if

          end do

        end function test_ini


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


        function get_data_S_exchanges(rank,other_rank)
     $     result(data_ex)

          implicit none

          integer          , intent(in) :: rank
          integer          , intent(in) :: other_rank
          real(rkind), dimension(6,8,2) :: data_ex

          integer :: i,j

          data_ex = initialize_data_exchanged(rank)

          data_ex(1:6,1:2,1) = reshape((/
     $         ((other_rank*10,i=1,6),j=1,2)/),
     $         (/6,2/))

          data_ex(1:6,1:2,2) = reshape((/
     $         ((other_rank*10+1,i=1,6),j=1,2)/),
     $         (/6,2/))

        end function get_data_S_exchanges


        function get_data_E_exchanges(rank,other_rank)
     $     result(data_ex)

          implicit none

          integer          , intent(in) :: rank
          integer          , intent(in) :: other_rank
          real(rkind), dimension(6,8,2) :: data_ex

          integer :: i,j

          data_ex = initialize_data_exchanged(rank)

          data_ex(5:6,3:6,1) = reshape((/
     $         ((other_rank*10,i=5,6),j=3,6)/),
     $         (/2,4/))

          data_ex(5:6,3:6,2) = reshape((/
     $         ((other_rank*10+1,i=5,6),j=3,6)/),
     $         (/2,4/))

        end function get_data_E_exchanges


        function get_data_W_exchanges(rank,other_rank)
     $     result(data_ex)

          implicit none

          integer          , intent(in) :: rank
          integer          , intent(in) :: other_rank
          real(rkind), dimension(6,8,2) :: data_ex

          integer :: i,j

          data_ex = initialize_data_exchanged(rank)

          data_ex(1:2,3:6,1) = reshape((/
     $         ((other_rank*10,i=1,2),j=3,6)/),
     $         (/2,4/))

          data_ex(1:2,3:6,2) = reshape((/
     $         ((other_rank*10+1,i=1,2),j=3,6)/),
     $         (/2,4/))

        end function get_data_W_exchanges

      end program test_mpi_mg_bc
