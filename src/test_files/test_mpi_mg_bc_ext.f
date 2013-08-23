      !> @file
      !> test file for the object 'mpi_mg_bc_ext'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the initialization of the mpi attributes of 'mpi_mg_bc_ext'
      !
      !> @date
      ! 23_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_mpi_mg_bc_ext

        use cg_operators_class , only : cg_operators
        use field_par_class    , only : field_par
        use mpi
        use mpi_mg_bc_ext_class, only : mpi_mg_bc_ext
        use mpi_process_class  , only : mpi_process
        use parameters_constant, only : periodic_xy_choice,
     $                                  reflection_xy_choice,
     $                                  N,S,E,W,
     $                                  compute_and_exchange_proc,
     $                                  only_exchange_proc
        use parameters_input   , only : npx,npy,bc_choice

        implicit none


        !< operators tested
        type(field_par)    :: f_tested
        type(mpi_process)  :: mpi_op
        type(mpi_mg_bc_ext):: mpi_mg
        type(cg_operators) :: s_op

        
        !< intermediate variables
        integer :: k
        integer :: proc_rank
        logical :: test_validated
        integer :: ierror
        

        !< test data
        integer               :: test_proc_x_choice
        integer               :: test_proc_y_choice
        integer, dimension(2) :: test_exchange_id


        !< the test is designed for (npx,npy)=(2,2)
        !> and periodic boundary conditions
        if((npx.ne.2).or.(npy.ne.2)) then
           stop '(''the test needs (npx,npy)=(2,2)'')'
        end if


        !< initialization of the mpi process
        call mpi_op%ini_mpi()


        !< initialization of the cartesian communicator
        call f_tested%ini_cartesian_communicator()


        !< get the rank of the processor computing this test
        call MPI_COMM_RANK(f_tested%comm_2d, proc_rank, ierror)
        if(ierror.ne.MPI_SUCCESS) then
           stop 'test_mpi_messenger_bc: MPI_COMM_RANK failed'
        end if


        !< initialize the test data
        select case(bc_choice)
     
          case(periodic_xy_choice)
             select case(proc_rank)
               case(0) 
                  test_proc_x_choice = only_exchange_proc
                  test_proc_y_choice = only_exchange_proc
                  test_exchange_id   = [E,N]
               case(1)
                  test_proc_x_choice = only_exchange_proc
                  test_proc_y_choice = only_exchange_proc
                  test_exchange_id   = [E,S]
               case(2)
                  test_proc_x_choice = only_exchange_proc
                  test_proc_y_choice = only_exchange_proc
                  test_exchange_id   = [W,N]
               case(3)
                  test_proc_x_choice = only_exchange_proc
                  test_proc_y_choice = only_exchange_proc
                  test_exchange_id   = [W,S]
               case default
                  stop 'proc rank not correct'
               end select

          case(reflection_xy_choice)
             select case(proc_rank)
               case(0) 
                  test_proc_x_choice = compute_and_exchange_proc
                  test_proc_y_choice = compute_and_exchange_proc
                  test_exchange_id   = [E,N]
               case(1)
                  test_proc_x_choice = compute_and_exchange_proc
                  test_proc_y_choice = compute_and_exchange_proc
                  test_exchange_id   = [E,S]
               case(2)
                  test_proc_x_choice = compute_and_exchange_proc
                  test_proc_y_choice = compute_and_exchange_proc
                  test_exchange_id   = [W,N]
               case(3)
                  test_proc_x_choice = compute_and_exchange_proc
                  test_proc_y_choice = compute_and_exchange_proc
                  test_exchange_id   = [W,S]
               case default
                  stop 'proc rank not correct'
             end select

          case default
             stop 'bc_choice not recognized for the test'
        end select


        !< test the initialization of the object 'mpi_mg_bc_ext'
        call mpi_mg%initialize(f_tested,s_op)


        !< compare the attributes of 'mpi_mg_bc_ext' with the test data
        test_validated=mpi_mg%proc_x_choice.eq.test_proc_x_choice
        print '(''proc, '',I2,'' test_proc_x_choice: '', L1)', 
     $       proc_rank, test_validated

        test_validated=mpi_mg%proc_y_choice.eq.test_proc_y_choice
        print '(''proc, '',I2,'' test_proc_y_choice: '', L1)', 
     $       proc_rank, test_validated

        if(bc_choice.ne.periodic_xy_choice) then
           test_validated=.true.
           k=1
           do while(test_validated.and.(k.le.2))
              test_validated=mpi_mg%exchange_id(k).eq.test_exchange_id(k)
              k=k+1
           end do
           print '(''proc, '',I2,'' test_exchange_id: '', L1)', 
     $          proc_rank, test_validated
        end if


        !< finalization of the mpi process
        call mpi_op%finalize_mpi()

      end program test_mpi_mg_bc_ext
