      !> @file
      !> class encapsulating subroutines to integrate
      !> the governing equations using Runge-Kutta 3rd
      !> order time integration scheme developed in
      !> “Efficient implementation of essentially non-
      !> oscillatory shock-capturing methods”, J. Comput.
      !> Phys., 77 (1988), pp. 439-471, C.-W. Shu and
      !> S. Osher
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
      !> S. Osher
      !
      !> @date
      !> 13_08_2013 - initial version                   - J.L. Desmarais
      !
      !> \f{eqnarray*}{
      !> u_1     &=& u_n + \Delta t*\frac{d u_n}{dt} \\\
      !> u_2     &=& \frac{1}{4}u_n + \frac{3}{4} \left(
      !>             u_1 + \Delta t * \frac{d u_1}{dt} \right) \\\
      !> u_{n+1} &=& \frac{1}{3}u_n + \frac{2}{3} \left(
      !>             u_2 + \Delta t * \frac{d u_2}{dt}\right)\\\
      !> \f}
      !-----------------------------------------------------------------
      module rk3tvd_class

        use bc_operators_class , only : bc_operators
        use cg_operators_class , only : cg_operators
        use field_class        , only : field
        use parameters_input   , only : nx,ny,ne
        use parameters_kind    , only : rkind, ikind
        use phy_model_eq_class , only : phy_model_eq
        use td_integrator_class, only : td_integrator
        use td_operators_class , only : td_operators

        implicit none


        !> @class rk3tvd
        !> class encapsulating subroutines to integrate
        !> the governing equations using Runge-Kutta 3rd
        !> order time integration scheme
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
        type, extends(td_integrator) :: rk3tvd

          contains
          procedure, nopass :: integrate

        end type rk3tvd


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to integrate the governing equations using
        !> the numerical scheme developed by C.W.Shu and S.Osher
        !
        !> @date
        !> 13_08_2013 - initial version - J.L. Desmarais
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
        subroutine integrate(field_used, sd, p_model, td, dt)

          implicit none

          class(field)       , intent(inout) :: field_used
          type(cg_operators) , intent(in)    :: sd
          class(phy_model_eq), intent(in)    :: p_model
          class(td_operators), intent(in)    :: td
          real(rkind)        , intent(in)    :: dt

          real(rkind) :: b2 !<Runge-Kutta scheme coeff
          real(rkind) :: b3 !<Runge-Kutta scheme coeff

          integer(ikind) :: i,j,k
          integer        :: bc_size

          real(rkind), dimension(nx,ny,ne) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne) :: time_dev

          type(bc_operators) :: bc_used !<boundary conditions
          

          !<initialization of the coeff for the Runge-Kutta scheme
          if(rkind.eq.8) then
             b2 = 0.75d0
             b3 = 1.0d0/3.0d0
          else
             b2 = 0.75
             b3 = 1.0/3.0
          end if


          !<initialization of local variables for the allocation
          !>of nodes_tmp, a temporary table
          bc_size = sd%get_bc_size()

          
          !<runge-kutta first step
          !> u_1 = u_n + dt*d/dt(u_n)
          !> u_n is saved in nodes_tmp
          !> u_1 is saved in field_used%nodes
          !DEC$ FORCEINLINE RECURSIVE
          time_dev = td%compute_time_dev(field_used, sd, p_model)

          do k=1, ne
             do j=bc_size+1, ny-bc_size
                do i=bc_size+1, nx-bc_size
                   nodes_tmp(i,j,k)        = field_used%nodes(i,j,k)
                   field_used%nodes(i,j,k) = field_used%nodes(i,j,k) +
     $                                       dt*time_dev(i,j,k)
                end do
             end do
          end do
            
          !<apply the boundary conditions
          !DEC$ FORCEINLINE RECURSIVE
          call bc_used%apply_bc_on_nodes(field_used,sd)


          !<runge-kutta second step
          !> u_2 = 1/4*u_n + 3/4*(u_1 + dt * du_1/dt)
          !> u_n is saved in nodes_tmp
          !> u_2 is saved in field_used%nodes
          !DEC$ FORCEINLINE RECURSIVE
          time_dev = td%compute_time_dev(field_used, sd, p_model)

          if(rkind.eq.8) then
             do k=1, ne
                do j=bc_size+1, ny-bc_size
                   !DEC$ VECTOR ALIGNED
                   do i=bc_size+1, nx-bc_size
                      field_used%nodes(i,j,k) = b2*nodes_tmp(i,j,k) +
     $                     (1.0d0-b2)*(field_used%nodes(i,j,k)+
     $                     dt*time_dev(i,j,k))
                   end do
                end do
             end do
          else
             do k=1, ne
                do j=bc_size+1, ny-bc_size
                   !DEC$ VECTOR ALIGNED
                   do i=bc_size+1, nx-bc_size
                      field_used%nodes(i,j,k) = b2*nodes_tmp(i,j,k) +
     $                     (1.0-b2)*(field_used%nodes(i,j,k)+
     $                     dt*time_dev(i,j,k))
                   end do
                end do
             end do
          end if
          
          !<apply the boundary conditions
          !DEC$ FORCEINLINE RECURSIVE
          call bc_used%apply_bc_on_nodes(field_used,sd)


          !<runge-kutta third step
          !> u_{n+1} = 1/3*u_n + 2/3*(u_2 + dt du_2/dt
          !> u_n is saved in nodes_tmp
          !> u_{n+1} is saved in field_used%nodes
          !DEC$ FORCEINLINE RECURSIVE
          time_dev = td%compute_time_dev(field_used, sd, p_model)

          if(rkind.eq.8) then
             do k=1 ,ne
                do j=bc_size+1, ny-bc_size
                   !DEC$ VECTOR ALIGNED
                   do i=bc_size+1, nx-bc_size
                      field_used%nodes(i,j,k) = b3*nodes_tmp(i,j,k) +
     $                     (1.0d0-b3)*(field_used%nodes(i,j,k)+
     $                     dt*time_dev(i,j,k))
                   end do
                end do
             end do
          else
             do k=1 ,ne
                do j=bc_size+1, ny-bc_size
                   !DEC$ VECTOR ALIGNED
                   do i=bc_size+1, nx-bc_size
                      field_used%nodes(i,j,k) = b3*nodes_tmp(i,j,k) +
     $                     (1.0-b3)*(field_used%nodes(i,j,k)+
     $                     dt*time_dev(i,j,k))
                   end do
                end do
             end do
          end if

          !<apply the boundary conditions
          !DEC$ FORCEINLINE RECURSIVE
          call bc_used%apply_bc_on_nodes(field_used,sd)

        end subroutine integrate

      end module rk3tvd_class


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
