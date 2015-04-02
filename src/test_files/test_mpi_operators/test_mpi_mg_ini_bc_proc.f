      program test_mpi_mg_ini_bc_proc

        use mpi

        use mpi_mg_ini_bc_proc, only :
     $     ini_bc_procedures

        use mpi_process_class, only :
     $       mpi_process

        use mpi_mg_ini_bc_proc, only :
     $       ini_bc_procedures

        use parameters_constant, only :
     $       N,S,E,W,
     $       reflection_xy_choice,
     $       only_exchange_proc,
     $       compute_and_exchange_proc

        use parameters_input, only :
     $       nx,ny,ne,
     $       npx,npy,
     $       bc_choice

        implicit none

        logical :: detailled
        logical :: test_loc
        logical :: test_validated

        integer :: rank
        integer :: ierror


        detailled = .true.
        test_validated = .true.


        test_loc = test_ini_bc_procedures(detailled)
        test_validated = test_validated.and.test_loc
        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_ini_bc_procedures: '',L1)', test_loc
           print '()'
        end if


        call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierror)
        if(rank.eq.0) then
           print '(''test_validated: '',L1)', test_validated
        end if


        call MPI_FINALIZE(ierror)


        contains

        
        function test_ini_bc_procedures(detailled)
     $       result(test_validated)

          implicit none

          logical, intent(in) :: detailled
          logical             :: test_validated


          type(mpi_process) :: mpi_process_used
          integer           :: ierror
          integer           :: nb_procs
          integer           :: comm2d
          integer           :: rank

          integer, dimension(6)   :: test_proc_x_choice
          integer, dimension(6)   :: test_proc_y_choice
          integer, dimension(2,6) :: test_exchange_id

          integer               :: proc_x_choice
          integer               :: proc_y_choice
          integer, dimension(2) :: exchange_id

          logical               :: test_loc
          logical, dimension(6) :: test_validated_gathered
          
          integer :: i


          test_validated = .true.


          ! we initialize the test with 6 processors with 3 processors
          ! along the x-direction and 2 processors along the y-direction
          !------------------------------------------------------------
          !   ___ ___ ___
          !  | 1 | 3 | 5 |
          !  |___|___|___|
          !  | 0 | 2 | 4 |
          !  |___|___|___|
          !
          ! we use reflection_xy b.c. such that the exchanges should
          ! look like
          !  ______ _______ _______ 
          ! |      |       |       |
          ! |  1 /_|_\ 3 /_|_\ 5   |       
          ! |  . \ | / . \ | / .   |
          ! | /|\  |  /|\  |  /|\  |
          ! |__|___|__ |___|___|___|
          ! |  |   |   |   |   |   |
          ! | \|/  |  \|/  |  \|/  |
          ! |  '   |   '   |   '   |
          ! |  0 /_|_\ 2 /_|_\ 4   |
          ! |    \ | /   \ | /     |
          ! |______|_______|_______|
          !
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
     $         (ne.eq.2).and.
     $         (bc_choice.eq.reflection_xy_choice))) then

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

          test_proc_x_choice = [compute_and_exchange_proc,
     $                          compute_and_exchange_proc,
     $                          only_exchange_proc,
     $                          only_exchange_proc,
     $                          compute_and_exchange_proc,
     $                          compute_and_exchange_proc]

          test_proc_y_choice = [compute_and_exchange_proc,
     $                          compute_and_exchange_proc,
     $                          compute_and_exchange_proc,
     $                          compute_and_exchange_proc,
     $                          compute_and_exchange_proc,
     $                          compute_and_exchange_proc]

          test_exchange_id   = reshape((/
     $                          E,N,
     $                          E,S,
     $                          0,N,
     $                          0,S,
     $                          W,N,
     $                          W,S/),
     $                          (/2,6/))

          ! output
          !============================================================
          call ini_bc_procedures(
     $         comm2d, proc_x_choice, proc_y_choice, exchange_id)


          ! validation
          !============================================================
          ! validation of proc_x_choice
          !------------------------------------------------------------
          test_loc = proc_x_choice.eq.test_proc_x_choice(rank+1)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test proc_x_choice('',I2,'') failed'')', rank
          end if

          
          ! validation of proc_y_choice
          !------------------------------------------------------------
          test_loc = proc_y_choice.eq.test_proc_y_choice(rank+1)
          test_validated = test_validated.and.test_loc
          if(detailled.and.(.not.test_loc)) then
             print '(''test proc_y_choice('',I2,'') failed'')', rank
          end if


          ! validation of exchange_id
          !------------------------------------------------------------
          if (proc_x_choice.ne.only_exchange_proc) then
             test_loc = exchange_id(1).eq.test_exchange_id(1,rank+1)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test exchange_id_x('',I2,'') failed'')', rank
             end if
          end if


          ! validation of exchange_id
          !------------------------------------------------------------
          if (proc_y_choice.ne.only_exchange_proc) then
             test_loc = exchange_id(2).eq.test_exchange_id(2,rank+1)
             test_validated = test_validated.and.test_loc
             if(detailled.and.(.not.test_loc)) then
                print '(''test exchange_id_y('',I2,'') failed'')', rank
             end if
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

        end function test_ini_bc_procedures

      end program test_mpi_mg_ini_bc_proc
