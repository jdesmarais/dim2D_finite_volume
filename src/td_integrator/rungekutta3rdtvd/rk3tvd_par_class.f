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
      module rk3tvd_par_class

        use bc_operators_par_class , only : bc_operators_par
        use cg_operators_class     , only : cg_operators
        use dim2d_eq_class         , only : dim2d_eq
        use field_par_class        , only : field_par
        use fv_operators_par_class , only : fv_operators_par
        use parameters_input       , only : nx,ny,ne
        use parameters_kind        , only : rkind, ikind
        use td_integrator_par_class, only : td_integrator_par
        use rk3tvd_steps_class     , only : rk3tvd_steps

        implicit none


        !> @class rk3tvd_par
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
        type, extends(td_integrator_par) :: rk3tvd_par

          contains
          procedure, nopass :: integrate

        end type rk3tvd_par


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
        !>@param sd
        !> space discretization operators
        !
        !>@param p_model
        !> physical model
        !
        !>@param td
        !> time discretisation operators
        !
        !>@param dt
        !> time step integrated
        !--------------------------------------------------------------
        subroutine integrate(
     $       field_used, sd, p_model, td_par, bc_par_used, dt)

          implicit none

          type(field_par)       , intent(inout) :: field_used
          type(cg_operators)    , intent(in)    :: sd
          type(dim2d_eq)        , intent(in)    :: p_model
          type(fv_operators_par), intent(in)    :: td_par
          type(bc_operators_par), intent(in)    :: bc_par_used
          real(rkind)           , intent(in)    :: dt

          integer                          :: bc_size
          real(rkind), dimension(nx,ny,ne) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne) :: time_dev
          type(rk3tvd_steps)               :: rk3tvd_int


          !<initialization of local variables for the allocation
          !>of nodes_tmp, a temporary table
          bc_size = sd%get_bc_size()

          
          !<runge-kutta first step
          !> u_1 = u_n + dt*d/dt(u_n)
          !> u_n is saved in nodes_tmp
          !> u_1 is saved in field_used%nodes
          !DEC$ FORCEINLINE RECURSIVE
          time_dev = td_par%compute_time_dev(
     $         field_used, sd, p_model, bc_par_used)

          !DEC$ FORCEINLINE RECURSIVE
          call rk3tvd_int%compute_1st_step(
     $         dt, time_dev, bc_size,
     $         field_used, nodes_tmp)
            
          !<apply the boundary conditions
          !DEC$ FORCEINLINE RECURSIVE
          call bc_par_used%apply_bc_on_nodes(
     $         field_used, field_used%nodes, sd, p_model)


          !<runge-kutta second step
          !> u_2 = 1/4*u_n + 3/4*(u_1 + dt * du_1/dt)
          !> u_n is saved in nodes_tmp
          !> u_2 is saved in field_used%nodes
          !DEC$ FORCEINLINE RECURSIVE
          time_dev = td_par%compute_time_dev(
     $         field_used, sd, p_model, bc_par_used)

          !DEC$ FORCEINLINE RECURSIVE
          call rk3tvd_int%compute_2nd_step(
     $         dt, time_dev, bc_size,
     $         field_used, nodes_tmp)
          
          !<apply the boundary conditions
          !DEC$ FORCEINLINE RECURSIVE
          call bc_par_used%apply_bc_on_nodes(
     $         field_used, field_used%nodes, sd, p_model)


          !<runge-kutta third step
          !> u_{n+1} = 1/3*u_n + 2/3*(u_2 + dt du_2/dt
          !> u_n is saved in nodes_tmp
          !> u_{n+1} is saved in field_used%nodes
          !DEC$ FORCEINLINE RECURSIVE
          time_dev = td_par%compute_time_dev(
     $         field_used, sd, p_model, bc_par_used)

          !DEC$ FORCEINLINE RECURSIVE
          call rk3tvd_int%compute_3rd_step(
     $         dt, time_dev, bc_size,
     $         field_used, nodes_tmp)

          !<apply the boundary conditions
          !DEC$ FORCEINLINE RECURSIVE
          call bc_par_used%apply_bc_on_nodes(
     $         field_used, field_used%nodes, sd, p_model)

        end subroutine integrate

      end module rk3tvd_par_class


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
