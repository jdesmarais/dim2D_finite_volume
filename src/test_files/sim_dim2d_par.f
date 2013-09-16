      !> @file
      !> program for the simulation of the Diffuse Interface
      !> Model in 2D on a parallel memory distributed system
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> run the simulation using the Diffuse Interface Model,
      !> the cg_operators for the spatial discretization,
      !> and finite volume and the Runge-Kutta 3rd order TVD
      !> as time integration methods on a parallel memory
      !> distributed system
      !
      !> @date
      ! 28_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program sim_dim2d_par

        use bc_operators_par_class     , only : bc_operators_par
        use cg_operators_class         , only : cg_operators
        use dim2d_eq_class             , only : dim2d_eq
        use field_par_class            , only : field_par
        use fv_operators_class         , only : fv_operators
        use mpi_process_class          , only : mpi_process
        use nf90_operators_wr_par_class, only : nf90_operators_wr_par
        use parameters_input           , only : ne,t_max,dt,detail_print
        use parameters_kind            , only : ikind, rkind
        use rk3tvd_par_class           , only : rk3tvd_par

        implicit none


        !<operators needed for the simulation
        type(field_par)             :: f_simulated!< field simulated
        type(cg_operators)          :: sd         !< space discretisation
        type(dim2d_eq)              :: p_model    !< physical model
        type(bc_operators_par)      :: bc_used    !< boundary conditions 
        type(fv_operators)          :: td         !< time discretisation
        type(rk3tvd_par)            :: ti         !< time integration
        type(nf90_operators_wr_par) :: io_writer  !< output management
        type(mpi_process)           :: mpi_op     !< mpi process

        !<intermediate variables for the simulation
        integer(ikind) :: nt,t,output_print
        integer        :: bc_size
        real(rkind)    :: time

        !<CPU recorded times
        real :: time1, time2, time3


        !< get the initial CPU time
        !>------------------------------------------------------
        !> this CPU time will be used as clock reference to
        !> compute the total duration of the simulation
        !>------------------------------------------------------
        call CPU_TIME(time1)


        !< check the number of governing equations
        !>------------------------------------------------------
        !> ne which is the third dimension of the main table
        !> needs to correspond with the number of governing
        !> equations defined in the physical model
        !>------------------------------------------------------
        if(ne.ne.p_model%get_eq_nb()) then
           stop 'ne is not correct : check the physical model'
        end if


        !< initialize the mpi processes
        !>------------------------------------------------------
        !> start mpi on the processors computing the simulation
        !>------------------------------------------------------
        call mpi_op%ini_mpi()


        !< initialize intermediate variables
        !>------------------------------------------------------
        !> initialize the variables determining the total number
        !> of timesteps for the simulation
        !>------------------------------------------------------
        bc_size      = sd%get_bc_size()
        nt           = int(t_max/dt)
        if(detail_print.eq.0) then
           output_print=0
        else
           output_print = int(1.0d0/detail_print)
        end if


        !< initialize the field
        !>------------------------------------------------------
        !> initialize the communicator between the tiles
        !> initialize the coordinate tables for ech tile
        !> apply the initial conditions using the physical model
        !>------------------------------------------------------
        time = 0
        call f_simulated%ini_cartesian_communicator()
        call f_simulated%ini_coordinates(bc_size)
        call p_model%apply_ic(f_simulated)


        !< initialize the boundary conditions
        !>------------------------------------------------------
        !> initialize the procedures for the computation of the 
        !> boundary layers: ghost cells between tiles and boundary
        !> conditions at the edge of the computational domain
        !>------------------------------------------------------
        call bc_used%initialize(f_simulated,sd)
        call bc_used%apply_bc_on_nodes(
     $       f_simulated, f_simulated%nodes,sd,p_model)


        !< initialize the output writer
        !>------------------------------------------------------
        !> initialize the subarray limits of the nodes written by
        !> the current tile on the output file for the entire
        !> computational domain
        !> write the initial state on 'data0.nc'
        !>------------------------------------------------------
        call io_writer%initialize(f_simulated,sd)
        call io_writer%write_data(f_simulated,p_model,bc_size,time)


        !<initialization time
        !>------------------------------------------------------
        !> compute the initialization time
        !>------------------------------------------------------
        call CPU_TIME(time2)
        print *, 'time_elapsed: ', time2-time1


        !< integrate the field until t=t_max
        !>------------------------------------------------------
        !> use the runge-kutta 3rd order time integration scheme,
        !> finite volume method, the Diffuse Interface Model in 2D
        !> and the space discretisation method developed by Cockburn 
        !> and Gau
        !>------------------------------------------------------        
        do t=1, nt

           !< compute the simulation time
           time=(t-1)*dt

           !DEC$ FORCEINLINE RECURSIVE
           call ti%integrate(f_simulated,sd,p_model,td,bc_used,dt)

           !< write the output data
           if((output_print.eq.1).or.
     $        ((output_print.ne.0).and.(mod(t,output_print).eq.0))) then
              call io_writer%write_data(f_simulated,p_model,bc_size,time)
           end if

        end do


        !< write the last timestep
        !>------------------------------------------------------
        !> even if the user asked no output, the initial state
        !> and the last state of the simulation are written
        !>------------------------------------------------------        
        if((output_print.eq.0).or.(mod(nt,output_print).ne.0)) then
           call io_writer%write_data(f_simulated,p_model,bc_size,time)
        end if


        !< print the time needed for the simulation
        !>------------------------------------------------------
        !> the total simulation time is evaluated by comparing
        !> the processor clock time at the begining and at the
        !> end of the simulation
        !>------------------------------------------------------
        call CPU_TIME(time3)
        print *, 'time_elapsed: ', time3-time1


        !< finalize the MPI processes
        !>------------------------------------------------------
        call mpi_op%finalize_mpi()

      end program sim_dim2d_par
