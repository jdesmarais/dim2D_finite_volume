      !> @file
      !> program for the simulation of the Diffuse Interface
      !> Model on a fixed 2D-domain using a parallel memory
      !> distributed system: one domain = 2D-grid of tiles
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> program for the simulation of the Diffuse Interface
      !> Model on a fixed 2D-domain using a parallel memory
      !> distributed system: one domain = 2D-grid of tiles
      !
      !> @date
      !> 02_04_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program sim_dim2d_par

        use mpi_process_class, only :
     $       mpi_process

        use field_par_class, only :
     $       field_par

        use parameters_input, only :
     $       t_max,
     $       dt, 
     $       detail_print

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none


        ! operators needed for the simulation
        type(field_par)    :: f_simulated ! field simulated
        type(mpi_process)  :: mpi_op      ! mpi process

        ! intermediate variables for the simulation
        integer(ikind) :: nt,t,output_print
        real(rkind)    :: time

        ! CPU recorded times
        real :: time1, time2, time3


        ! get the initial CPU time
        !------------------------------------------------------
        ! this CPU time will be used as clock reference to
        ! compute the total duration of the simulation
        !------------------------------------------------------
        call CPU_TIME(time1)


        ! initialize the mpi process
        !------------------------------------------------------
        call mpi_op%ini_mpi()


        ! initialize intermediate variables
        !------------------------------------------------------
        ! initialize the variables determining the total number
        ! of timesteps for the simulation
        !------------------------------------------------------
        nt = int(t_max/dt)
        if(detail_print.eq.0) then
           output_print=0
        else
           output_print = int(1.0d0/detail_print)
        end if


        ! initialize the field
        !------------------------------------------------------
        ! initialize the communicator between the tiles
        ! initialize the coordinate tables for ech tile
        ! apply the initial conditions using the physical model
        !------------------------------------------------------
        time = 0
        call f_simulated%ini()
        call f_simulated%apply_bc_on_nodes()
        call f_simulated%write_data()


        ! initialization time
        !------------------------------------------------------
        ! compute the initialization time
        !------------------------------------------------------
        call CPU_TIME(time2)
        print *, 'time_elapsed: ', time2-time1


        ! integrate the field until t=t_max
        !------------------------------------------------------     
        do t=1, nt

           !< compute the simulation time
           time=(t-1)*dt

           !DEC$ FORCEINLINE RECURSIVE
           call f_simulated%integrate(dt)

           !< write the output data
           if((output_print.eq.1).or.
     $        ((output_print.ne.0).and.(mod(t,output_print).eq.0))) then
              call f_simulated%write_data()
           end if

        end do


        ! print the time needed for the simulation
        !------------------------------------------------------
        ! the total simulation time is evaluated by comparing
        ! the processor clock time at the begining and at the
        ! end of the simulation
        !------------------------------------------------------
        call CPU_TIME(time3)
        print *, 'time_elapsed: ', time3-time1


        ! write the last timestep
        !------------------------------------------------------
        ! even if the user asked no output, the initial state
        ! and the last state of the simulation are written
        !------------------------------------------------------        
        if((output_print.eq.0).or.(mod(nt,output_print).ne.0)) then
           call f_simulated%write_data()
        end if


        ! finalize the MPI processes
        !------------------------------------------------------
        call mpi_op%finalize_mpi()

      end program sim_dim2d_par
