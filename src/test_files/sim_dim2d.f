      !> @file
      !> program for the simulation of the Diffuse Interface
      !> Model in 2D
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> run the simulation using the Diffuse Interface Model,
      !> the cg_operators for the spatial discretization,
      !> and finite volume and the Runge-Kutta 3rd order TVD
      !> as time integration methods
      !
      !> @date
      ! 19_08_2013 - initial version              - J.L. Desmarais
      ! 15_07_2014 - composition over inheritance - J.L. Desmarais
      !-----------------------------------------------------------------
      program sim_dim2d

        use field_class            , only : field
        use parameters_input       , only : t_max,dt,detail_print
        use parameters_kind        , only : ikind

        implicit none


        !<operators needed for the simulation
        type(field) :: f_simulated !< field simulated
        
        !<intermediate variables for the simulation
        integer(ikind) :: nt, output_print
        integer(ikind) :: t

        !<CPU recorded times
        real :: time1, time2, time3


        !<get the initial CPU time
        call CPU_TIME(time1)


        !<allocate the field
        nt           = int(t_max/dt)
        output_print = int(1.0d0/detail_print)


        !<initialize the field
        call f_simulated%ini()
        call f_simulated%apply_bc_on_nodes()
        call f_simulated%write_data()


        !<initialization time
        call CPU_TIME(time2)
        print *, 'time_elapsed: ', time2-time1


        !<integrate the field until t=t_max
        do t=1, nt

           !DEC$ FORCEINLINE RECURSIVE
           call f_simulated%integrate(dt)

           !< write the output data
           if((output_print.eq.1).or.
     $        ((output_print.ne.0).and.(mod(t,output_print).eq.0))) then
              call f_simulated%write_data()
           end if

        end do


        !<print the time needed for the simulation
        call CPU_TIME(time3)
        print *, 'time_elapsed: ', time3-time1


        !<write the last timestep
        if((output_print.eq.0).or.(mod(nt,output_print).ne.0)) then
           call f_simulated%write_data()
        end if

      end program sim_dim2d
