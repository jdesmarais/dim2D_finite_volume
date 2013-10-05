      !> @file
      !> test file for the object 'field_par'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the initialization of the mpi attributes of 'field_par'
      !
      !> @date
      ! 21_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_field_par

        use field_par_class  , only : field_par
        use mpi
        use mpi_process_class, only : mpi_process
        use parameters_kind  , only : rkind

        
        !< operators tested
        type(field_par)   :: f_tested
        type(mpi_process) :: mpi_op

        real(rkind) :: x_min, x_max, y_min, y_max

        !< initialize the mpi processes
        call mpi_op%ini_mpi()


        !< initialize the cartesian communicator
        call f_tested%ini_cartesian_communicator()
        

        !< initialize the coordinate table
        call f_tested%ini_coordinates()


        !< test the data
        print '(''I, proc '', I2, '' I belongs to the comm '', I2)',
     $       f_tested%usr_rank, f_tested%comm_2d


        !< finalize the mpi processes
        call mpi_op%finalize_mpi()

      end program test_field_par
