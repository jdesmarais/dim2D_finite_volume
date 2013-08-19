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
      ! 19_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      program sim_dim2d

        use cg_operators_class     , only : cg_operators
        use dim2d_eq_class         , only : dim2d_eq
        use field_class            , only : field
        use fv_operators_class     , only : fv_operators
        use nf90_operators_wr_class, only : nf90_operators_wr
        use parameters_kind        , only : ikind, rkind
        use rk3tvd_class           , only : rk3tvd

        implicit none


        !<inputs for the simulation
        real(rkind) :: x_min, x_max, dx
        real(rkind) :: y_min, y_max, dy
        real(rkind) :: t_max, dt
        real(rkind) :: detail_print

        !<operators needed for the simulation
        type(field)             :: f_simulated !< field simulated
        type(cg_operators)      :: s           !< spatial discretisation
        type(dim2d_eq)          :: p_model     !< physical model
        type(fv_operators)      :: td          !< time discretisation
        type(rk3tvd)            :: ti          !< time integration
        type(nf90_operators_wr) :: io_writer   !< output management

        !<intermediate variables for the simulation
        integer(ikind) :: nx,ny,nt, output_print
        integer        :: ne, bc_size
        integer(ikind) :: i,j,t
        real(rkind)    :: time

        !<CPU recorded times
c$$$        real :: time1, time2



        !<get the initial CPU time
        !call CPU_TIME(time1)


        !<read the inputs
        dx           =  0.01
        x_min        = -0.4
        x_max        =  0.4
                       
        dy           =  0.01
        y_min        = -0.4
        y_max        =  0.4
                       
        t_max        =  0.15
        dt           =  0.00005
        detail_print =  0.0


        !<allocate the field
        bc_size      = s%get_bc_size()
        nx           = (x_max-x_min)/dx + 2*bc_size
        ny           = (y_max-y_min)/dy + 2*bc_size
        nt           = int(t_max/dt)
        output_print = int(1.0d0/detail_print)
        ne           = p_model%get_eq_nb()

        call f_simulated%allocate_tables(nx,ny,ne)


        !<initialize the field
        time = 0

        do i=1, nx
           f_simulated%x_map(i)=x_min + (i-(bc_size+1))*f_simulated%dx
        end do

        do j=1, ny
           f_simulated%y_map(j)=y_min + (j-(bc_size+1))*f_simulated%dy
        end do

        call p_model%apply_ic(f_simulated)


        !<write the initial state in an output file
        !call io_writer%initialize()
        !call io_writer%write_data(f_simulated,p_model,time)


        !<integrate the field until t=t_max
        do t=1, nt
           time=(t-1)*dt

           call ti%integrate(f_simulated,s,p_model,td,dt)

c$$$           if((output_print.ne.0).and.(mod(t,output_print).eq.0)) then
c$$$              call io_writer%write_data(f_simulated,p_model,time)
c$$$           end if

        end do


        !<write the last timestep
c$$$        if(mod(nt,output_print).ne.0) then
c$$$           call io_writer%write_data(f_simulated,p_model,time)
c$$$        end if


        !<print the time needed for the simulation
c$$$        call CPU_TIME(time2)
c$$$        print *, 'time_elapsed: ', time2-time1

      end program sim_dim2d
