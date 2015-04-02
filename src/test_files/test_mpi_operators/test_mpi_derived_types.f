      program test_mpi_derived_types

        use check_data_module, only :
     $       is_real_matrix3D_validated

        use mpi
        
        use mpi_mg_derived_types, only :
     $       ini_mpi_derived_types

        use mpi_process_class, only :
     $       mpi_process

        use parameters_constant, only :
     $       N,S,E,W

        use parameters_input, only :
     $       bc_size

        use parameters_kind, only :
     $       ikind,
     $       rkind


        implicit none


        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        type(mpi_process) :: mpi_process_used
        integer           :: rank
        integer           :: ierror


        detailled = .true.
        test_validated  = .true.


        test_loc = test_ini_mpi_derived_types(detailled)
        test_validated = test_validated.and.test_loc
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_ini_mpi_derived_types: '',L1)', test_loc
           print '()'
        end if


        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_validated: '',L1)', test_validated
        end if

        call mpi_process_used%finalize_mpi()

        contains


        function test_ini_mpi_derived_types(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mpi_process) :: mpi_process_used

          integer, parameter :: tile_nx = 6
          integer, parameter :: tile_ny = 8
          integer, parameter :: tile_ne = 2
          
          real(rkind), dimension(tile_nx,tile_ny,tile_ne)   :: data_exchanged
          real(rkind), dimension(tile_nx,tile_ny,tile_ne,4) :: test_data_exchanged


          integer, dimension(4) :: com_send
          integer, dimension(4) :: com_recv

          integer, parameter                  :: tag1=1111
          integer, parameter                  :: tag2=200
          integer                             :: nb_procs
          integer                             :: usr_rank
          integer                             :: ngh_rank
          !integer, dimension(MPI_STATUS_SIZE) :: status
          integer                             :: ierror

          integer :: k
          logical :: test_loc


          test_validated = .true.


          ! inputs
          !============================================================
          ! initialization of the MPI processes
          call mpi_process_used%ini_mpi()

          ! verification of the number of processors
          call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,ierror)
          if(ierror.ne.MPI_SUCCESS) then
             print '(''MPI_COMM_SIZE FAILED'')'
          end if
          if(nb_procs.ne.2) then
             print '(''the test requires 2 processors'')'
             test_validated = .false.
             call mpi_process_used%finalize_mpi()
          end if

          ! determination of the rank of the current
          ! processor
          call MPI_COMM_RANK(MPI_COMM_WORLD,usr_rank,ierror)
          if(ierror.ne.MPI_SUCCESS) then
             print '(''MPI_COMM_RANK FAILED'')'
          end if
          ngh_rank = mod(usr_rank+1,2)

          !determination of the test data
          test_data_exchanged(:,:,:,N) = get_data_N_exchanges(usr_rank)
          test_data_exchanged(:,:,:,S) = get_data_S_exchanges(usr_rank)
          test_data_exchanged(:,:,:,E) = get_data_E_exchanges(usr_rank)
          test_data_exchanged(:,:,:,W) = get_data_W_exchanges(usr_rank)


          ! output
          !============================================================
          ! creation of the derived types to be tested
          call ini_mpi_derived_types(
     $         tile_nx,tile_ny,tile_ne,
     $         com_send, com_recv)

          
          ! validation
          !============================================================
          do k=1,4

             data_exchanged = initialize_data_exchanged(usr_rank)

             call MPI_SENDRECV (
     $            data_exchanged, 1, com_send(k),
     $            ngh_rank, tag1,
     $            data_exchanged, 1, com_recv(k),
     $            ngh_rank, tag1,
     $            MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierror)

             test_loc = is_real_matrix3D_validated(
     $            data_exchanged,
     $            test_data_exchanged(:,:,:,k),
     $            detailled)

             test_validated = test_validated.and.test_loc

             if(detailled.and.(.not.test_loc)) then
                print '(''test('',I2,'') failed'')', k
             end if

          end do

        end function test_ini_mpi_derived_types


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


        function get_data_N_exchanges(rank)
     $     result(data_ex)

          implicit none

          integer          , intent(in) :: rank
          real(rkind), dimension(6,8,2) :: data_ex

          integer :: other_rank
          integer :: i,j

          data_ex = initialize_data_exchanged(rank)

          other_rank = mod(rank+1,2)

          data_ex(1:6,7:8,1) = reshape((/
     $         ((other_rank*10,i=1,6),j=1,2)/),
     $         (/6,2/))

          data_ex(1:6,7:8,2) = reshape((/
     $         ((other_rank*10+1,i=1,6),j=1,2)/),
     $         (/6,2/))

        end function get_data_N_exchanges


        function get_data_S_exchanges(rank)
     $     result(data_ex)

          implicit none

          integer          , intent(in) :: rank
          real(rkind), dimension(6,8,2) :: data_ex

          integer :: other_rank
          integer :: i,j

          data_ex = initialize_data_exchanged(rank)

          other_rank = mod(rank+1,2)

          data_ex(1:6,1:2,1) = reshape((/
     $         ((other_rank*10,i=1,6),j=1,2)/),
     $         (/6,2/))

          data_ex(1:6,1:2,2) = reshape((/
     $         ((other_rank*10+1,i=1,6),j=1,2)/),
     $         (/6,2/))

        end function get_data_S_exchanges


        function get_data_E_exchanges(rank)
     $     result(data_ex)

          implicit none

          integer          , intent(in) :: rank
          real(rkind), dimension(6,8,2) :: data_ex

          integer :: other_rank
          integer :: i,j

          data_ex = initialize_data_exchanged(rank)

          other_rank = mod(rank+1,2)

          data_ex(5:6,3:6,1) = reshape((/
     $         ((other_rank*10,i=5,6),j=3,6)/),
     $         (/2,4/))

          data_ex(5:6,3:6,2) = reshape((/
     $         ((other_rank*10+1,i=5,6),j=3,6)/),
     $         (/2,4/))

        end function get_data_E_exchanges


        function get_data_W_exchanges(rank)
     $     result(data_ex)

          implicit none

          integer          , intent(in) :: rank
          real(rkind), dimension(6,8,2) :: data_ex

          integer :: other_rank
          integer :: i,j

          data_ex = initialize_data_exchanged(rank)

          other_rank = mod(rank+1,2)

          data_ex(1:2,3:6,1) = reshape((/
     $         ((other_rank*10,i=1,2),j=3,6)/),
     $         (/2,4/))

          data_ex(1:2,3:6,2) = reshape((/
     $         ((other_rank*10+1,i=1,2),j=3,6)/),
     $         (/2,4/))

        end function get_data_W_exchanges

      end program test_mpi_derived_types
