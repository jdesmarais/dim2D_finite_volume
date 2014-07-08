      !> @file
      !> test file for the object 'nf90_operators_wr_par'
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> test the subroutines writing the netcdf files
      !> in a parallel distributed memory system
      !
      !> @date
      ! 28_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program test_nf90_operators_par

        use cg_operators_class         , only : cg_operators
        use dim2d_eq_class             , only : dim2d_eq
        use field_par_class            , only : field_par
        use mpi_process_class          , only : mpi_process
        use nf90_operators_wr_par_class, only : nf90_operators_wr_par
        use parameters_input           , only : ne
        use parameters_kind            , only : rkind


        implicit none

        
        !<operators tested
        type(field_par)            :: field_tested
        type(cg_operators)         :: sd_op
        type(dim2d_eq)             :: p_model
        type(mpi_process)          :: mpi_op
        type(nf90_operators_wr_par):: nf90_writer
        real(rkind), parameter     :: time=3.0
        real(rkind) :: x_min, x_max, y_min, y_max

        !<CPU recorded times
        real    :: time1, time2

        
        !<if nx<4, ny<4 then the test cannot be done
        if(ne.ne.4) then
           stop 'the test needs: ne=4'
        end if
        

        !<get the initial CPU time
        call CPU_TIME(time1)


        !< initialization of the mpi process
        call mpi_op%ini_mpi()


        !< initialize the tables for the field
        x_min   = 0.
        x_max   = 1.
        y_min   = 0.
        y_max   = 1.

        call field_tested%ini_cartesian_communicator()
        call field_tested%ini_coordinates()
        call p_model%apply_ic(field_tested)


        !< write the data
        call nf90_writer%initialize(field_tested, sd_op, 3)
        call nf90_writer%write_data(field_tested, p_model, time)


        !< finalize the mpi process
        call mpi_op%finalize_mpi()


        !< compute the time elapsed
        call CPU_TIME(time2)
        print '(''time elapsed: '', F10.6)', time2-time1

      end program test_nf90_operators_par
