      program test_mpi_mg_neighbours

        use check_data_module, only :
     $     is_int_vector_validated

        use mpi

        use mpi_mg_neighbours, only :
     $       ini_neighbours_proc_id

        use mpi_process_class, only :
     $       mpi_process

        use parameters_input, only :
     $       npx, npy

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        integer           :: rank
        integer           :: ierror
        
        detailled = .true.
        test_validated = .true.


        test_loc = test_ini_neighbours_proc_id(detailled)
        test_validated = test_validated.and.test_loc
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_ini_neighbours_proc_id: '',L1)', test_loc
           print '()'
        end if        
        
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_validated: '',L1)', test_validated
        end if

        call MPI_FINALIZE(ierror)


        contains


        function test_ini_neighbours_proc_id(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mpi_process)       :: mpi_process_used
          integer                 :: nb_procs

          integer                 :: ierror
          integer                 :: comm2d
          integer                 :: rank
          integer, dimension(4)   :: com_rank
          integer, dimension(4,6) :: test_com_rank

          logical                 :: test_loc
          integer                 :: i
          logical, dimension(6)   :: test_loc_gathered


          ! input
          call MPI_INIT(ierror)

          call MPI_COMM_SIZE(MPI_COMM_WORLD,nb_procs,ierror)
          if(.not.(
     $         (nb_procs.eq.6).and.
     $         (npx.eq.3).and.
     $         (npy.eq.2)))then

             print '(''the test requires nb_procs=6'')'
             call MPI_FINALIZE(ierror)
             stop ''

          end if

          call mpi_process_used%ini_cartesian_communicator(
     $         comm2d, rank)

          test_com_rank(:,1) = [1,1,2,4]
          test_com_rank(:,2) = [0,0,3,5]
          test_com_rank(:,3) = [3,3,4,0]
          test_com_rank(:,4) = [2,2,5,1]
          test_com_rank(:,5) = [5,5,0,2]
          test_com_rank(:,6) = [4,4,1,3]


          ! output
          call ini_neighbours_proc_id(
     $         comm2d,
     $         [npx,npy],
     $         com_rank)


          ! validation
          test_loc = is_int_vector_validated(
     $         com_rank,
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
             test_validated = test_loc_gathered(1)
             do i=2,6
                test_validated = test_validated.and.test_loc_gathered(i)
             end do
          else
             test_validated = test_loc
          end if



        end function test_ini_neighbours_proc_id

      end program test_mpi_mg_neighbours
