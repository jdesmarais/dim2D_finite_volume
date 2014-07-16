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

        use pmodel_eq_class       , only : pmodel_eq
        use mpi_process_class     , only : mpi_process
        use io_operators_par_class, only : io_operators_par
        use parameters_input      , only : nx,ny,ne
        use parameters_kind       , only : rkind


        implicit none

        
        !<operators tested
        integer :: comm_2d
        integer :: usr_rank
        real(rkind), dimension(nx,ny,ne) :: nodes
        real(rkind), dimension(nx)       :: x_map
        real(rkind), dimension(ny)       :: y_map
        real(rkind)                      :: dx
        real(rkind)                      :: dy
        integer                          :: i,j        
        type(pmodel_eq)                  :: p_model
        type(mpi_process)                :: mpi_op
        type(io_operators_par)           :: nf90_writer
        real(rkind), parameter           :: time=3.0
        real(rkind)                      :: x_min
        real(rkind)                      :: x_max
        real(rkind)                      :: y_min
        real(rkind)                      :: y_max

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

        dx=0.6
        dy=0.7

        do j=1, ny
           do i=1, nx
              nodes(i,j,1) = i + (j-1)*nx
           end do
        end do

        do i=1, nx
           x_map(i)=(i-1)*dx
        end do

        do j=1, ny
           y_map(j)=(j-1)*dy
        end do

        call mpi_op%ini_cartesian_communicator(comm_2d, usr_rank)


        !< write the data
        call nf90_writer%ini(comm_2d, usr_rank)
        call nf90_writer%write_data(
     $       comm_2d,
     $       nodes, x_map, y_map,
     $       p_model, time)


        !< finalize the mpi process
        call mpi_op%finalize_mpi()


        !< compute the time elapsed
        call CPU_TIME(time2)
        print '(''time elapsed: '', F10.6)', time2-time1

      end program test_nf90_operators_par
