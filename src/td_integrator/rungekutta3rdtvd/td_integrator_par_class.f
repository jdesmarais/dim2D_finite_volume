      !> @file
      !> class encapsulating subroutines to integrate
      !> the governing equations using Runge-Kutta 3rd
      !> order time integration scheme developed in
      !> “Efficient implementation of essentially non-
      !> oscillatory shock-capturing methods”, J. Comput.
      !> Phys., 77 (1988), pp. 439-471, C.-W. Shu and
      !> S. Osher, on a parallel memory distributed system
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to integrate
      !> the governing equations using Runge-Kutta 3rd
      !> order time integration scheme developed in
      !> “Efficient implementation of essentially non-
      !> oscillatory shock-capturing methods”, J. Comput.
      !> Phys., 77 (1988), pp. 439-471, C.-W. Shu and
      !> S. Osher, on a parallel memory distributed system
      !
      !> @date
      !> 27_08_2013 - initial version - J.L. Desmarais
      !
      !> \f{eqnarray*}{
      !> u_1     &=& u_n + \Delta t*\frac{d u_n}{dt} \\\
      !> u_2     &=& \frac{1}{4}u_n + \frac{3}{4} \left(
      !>             u_1 + \Delta t * \frac{d u_1}{dt} \right) \\\
      !> u_{n+1} &=& \frac{1}{3}u_n + \frac{2}{3} \left(
      !>             u_2 + \Delta t * \frac{d u_2}{dt} \right) \\\
      !> \f}
      !-----------------------------------------------------------------
      module td_integrator_par_class

        use field_abstract_par_class, only :
     $     field_abstract_par

        use parameters_input, only :
     $       nx,ny,ne

        use parameters_kind, only :
     $       rkind, ikind

        use td_integrator_abstract_par_class, only :
     $       td_integrator_abstract_par

        use rk3tvd_steps_module, only :
     $       compute_1st_step,
     $       compute_2nd_step,
     $       compute_3rd_step

        implicit none


        !> @class td_integrator_par
        !> class encapsulating subroutines to integrate
        !> the governing equations using Runge-Kutta 3rd
        !> order time integration scheme on a parallel 
        !> memory distributed system
        !>
        !> @param integrate
        !> integrate the computational field for dt
        !> \f{eqnarray*}{
        !> u_1     &=& u_n + \Delta t*\frac{d u_n}{dt} \\\
        !> u_2     &=& \frac{1}{4}u_n + \frac{3}{4} \left(
        !>             u_1 + \Delta t * \frac{d u_1}{dt} \right) \\\
        !> u_{n+1} &=& \frac{1}{3}u_n + \frac{2}{3} \left(
        !>             u_2 + \Delta t * \frac{d u_2}{dt}\right)\\\
        !> \f}
        !---------------------------------------------------------------
        type, extends(td_integrator_abstract_par) :: td_integrator_par

          contains
          procedure, nopass :: integrate

        end type td_integrator_par


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to integrate the governing equations using
        !> the numerical scheme developed by C.W.Shu and S.Osher
        !> on a parallel memory distributed system
        !
        !> @date
        !> 27_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !
        !>@param dt
        !> time step integrated
        !--------------------------------------------------------------
        subroutine integrate(field_used, dt)

          implicit none

          class(field_abstract_par), intent(inout) :: field_used
          real(rkind)              , intent(in)    :: dt

          real(rkind), dimension(nx,ny,ne) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne) :: time_dev


          !<runge-kutta first step
          !> u_1 = u_n + dt*d/dt(u_n)
          !> u_n is saved in nodes_tmp
          !> u_1 is saved in field_used%nodes
          !DEC$ FORCEINLINE RECURSIVE
          time_dev = field_used%compute_time_dev()

          !DEC$ FORCEINLINE RECURSIVE
          call field_used%compute_integration_step(
     $         dt, nodes_tmp, time_dev, compute_1st_step)
            
          !<apply the boundary conditions
          !DEC$ FORCEINLINE RECURSIVE
          call field_used%apply_bc_on_nodes()


          !<runge-kutta second step
          !> u_2 = 1/4*u_n + 3/4*(u_1 + dt * du_1/dt)
          !> u_n is saved in nodes_tmp
          !> u_2 is saved in field_used%nodes
          !DEC$ FORCEINLINE RECURSIVE
          time_dev = field_used%compute_time_dev()

          !DEC$ FORCEINLINE RECURSIVE
          call field_used%compute_integration_step(
     $         dt, nodes_tmp, time_dev, compute_2nd_step)
          
          !<apply the boundary conditions
          !DEC$ FORCEINLINE RECURSIVE
          call field_used%apply_bc_on_nodes()


          !<runge-kutta third step
          !> u_{n+1} = 1/3*u_n + 2/3*(u_2 + dt du_2/dt
          !> u_n is saved in nodes_tmp
          !> u_{n+1} is saved in field_used%nodes
          !DEC$ FORCEINLINE RECURSIVE
          time_dev = field_used%compute_time_dev()

          !DEC$ FORCEINLINE RECURSIVE
          call field_used%compute_integration_step(
     $         dt, nodes_tmp, time_dev, compute_3rd_step)

          !<apply the boundary conditions
          !DEC$ FORCEINLINE RECURSIVE
          call field_used%apply_bc_on_nodes()

        end subroutine integrate

      end module td_integrator_par_class


c$$$      print *,'first step'
c$$$
c$$$          do j=1, nodes_profile(2)
c$$$             do i=1, nodes_profile(1)
c$$$                print '(''('',I2,'','',I2,'')'', 3X,
c$$$     $                  ''u1='',F10.2)',
c$$$     $               i,j,
c$$$     $               field_used%nodes(i,j,1)
c$$$             end do
c$$$          end do

c$$$         print *,'second step'
c$$$
c$$$          do j=1, nodes_profile(2)
c$$$             do i=1, nodes_profile(1)
c$$$                print '(''('',I2,'','',I2,'')'', 3X,
c$$$     $                  ''u2='',F10.2)',
c$$$     $               i,j,
c$$$     $               field_used%nodes(i,j,1)
c$$$             end do
c$$$          end do

c$$$         do j=1, nodes_profile(2)-2*bc_size
c$$$             do i=1, nodes_profile(1)-2*bc_size
c$$$                print '(''('',I2,'','',I2,'')'', 3X,
c$$$     $                  ''time_dev='',F10.2)',
c$$$     $               i,j,
c$$$     $               time_dev(i,j,1)
c$$$             end do
c$$$          end do
