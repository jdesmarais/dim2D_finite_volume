      program test_mpi_process

        use mpi
        use mpi_process_class, only :
     $       mpi_process


        implicit none

        
        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        detailled = .true.
        test_validated = .true.


        test_loc = test_ini_mpi(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_ini_mpi: '',L1)', test_loc
        print '()'


        test_loc = test_ini_cartesian_communicator(detailled)
        test_validated = test_validated.and.test_loc
        print '(''test_cartesian_communicator: '',L1)', test_loc
        print '()'


        print '(''test_validated: '',L1)', test_validated
        print '()'


        contains


        function test_ini_mpi(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mpi_process) :: mpi_process_used
          integer           :: rank
          integer           :: ierror


          call mpi_process_used%ini_mpi()

          call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
          test_validated = ierror.eq.MPI_SUCCESS
          if(detailled.and.(.not.test_validated)) then
             print '(''test MPI_COMM_RANK failed'')'
          end if

        end function test_ini_mpi


        function test_ini_cartesian_communicator(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mpi_process) :: mpi_process_used
          integer           :: comm2d
          integer           :: rank
          integer           :: rank_test
          integer           :: ierror

          
          test_validated = .true.


          call mpi_process_used%ini_cartesian_communicator(
     $         comm2d, rank)

          call MPI_COMM_RANK(comm2d,rank_test,ierror)

          test_loc = ierror.eq.MPI_SUCCESS
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_validated)) then
             print '(''test MPI_COMM_RANK failed'')'
          end if

          test_loc = rank_test.eq.rank
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_validated)) then
             print '(''test rank failed'')'
          end if          

          call mpi_process_used%finalize_mpi()

        end function test_ini_cartesian_communicator

      end program test_mpi_process
