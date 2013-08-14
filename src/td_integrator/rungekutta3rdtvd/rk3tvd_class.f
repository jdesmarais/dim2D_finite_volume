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
      !> order time integration scheme
      !
      !> @date
      !> 13_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module rk3tvd_class

        use field_class        , only : field
        use parameters_kind    , only : rkind, ikind
        use td_integrator_class, only : td_integrator

        implicit none


        !> @class rk3tvd
        !> class encapsulating subroutines to integrate
        !> the governing equations using Runge-Kutta 3rd
        !> order time integration scheme
        !>
        !> @param integrate
        !> integrate the computational field for dt
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
        !> interface to integrate the governing equations using
        !> space discretisation operators, physical model, 
        !> time discretisation operators and boundary conditions
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
        !>@param p
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

          real(rkind), parameter :: b2 = 0.75d0 !<coeff for the Runge-Kutta scheme
          real(rkind), parameter :: b3 = 1./3.  !<coeff for the Runge-Kutta scheme

          integer(ikind) :: i,j,k
          integer(ikind) :: nx,ny,ne,bc_size

          real(rkind), dimension(:,:,:), allocatable :: nodes_tmp
          real(rkind), dimension(:,:,:), allocatable :: time_dev

          type(bc_operators) :: bc_used !<boundary conditions
          

          !<initialization of local variables for the allocation
          !>of nodes_tmp, a temporary table
          nx      = size(field_used%nodes,1)
          ny      = size(field_used%nodes,2)
          ne      = size(field_used%nodes,3)
          bc_size = sd%get_bc_size()

          
          !<allocate the temporary tables
          allocate(nodes_tmp(nx,ny,ne))
          allocate(time_dev(nx,ny,ne))


          !<runge-kutta first step
          !>u_1 = u_n + dt*d/dt(u_n)
          call td%compute_time_dev(field_used, sd, p_model, time_dev)

          do k=1, nodes_profile(3)
             do j=bc_size+1, ny-bc_size
                do i=bc_size+1, nx-bc_size
                   nodes_tmp(i,j,k) = field_used%nodes(i,j,k) + dt*time_dev(i,j,k)
                end do
             end do
          end do
            
          call bc_used%apply_bc_on_nodes(nodes_tmp,sd)


          !runge-kutta second step
          !u_2 = b2*u_n +(1-b2)*(u_1 + *dt*d/dt(u_1))
          call field_bc_used%t(nodes_tmp, time_dev)

          do j=bc_size+1, nodes_profile(2)-bc_size
             do i=bc_size+1, nodes_profile(1)-bc_size
                do k=1, nodes_profile(3)
                   nodes_tmp(i,j,k) =
     $                  b2*field_bc_used%nodes(i,j,k) +
     $                  (1-b2)*(nodes_tmp(i,j,k) +
     $                  dt*time_dev(i,j,k))
                end do
             end do
          end do
          
          deallocate(time_dev)
          call field_bc_used%apply_bc(nodes_tmp)


          !runge-kutta third step
          !u_{n+1} = b3*u_n +(1-b3)*(u_1 + *dt*d/dt(u_1))
          call field_bc_used%t(nodes_tmp, time_dev)

          do j=bc_size+1, nodes_profile(2)-bc_size
             do i=bc_size+1, nodes_profile(1)-bc_size
                do k=1, nodes_profile(3)
                   field_bc_used%nodes(i,j,k) =
     $                  b3*field_bc_used%nodes(i,j,k) +
     $                  (1-b3)*(nodes_tmp(i,j,k) +
     $                  dt*time_dev(i,j,k))
                end do
             end do
          end do

          deallocate(nodes_tmp)
          deallocate(time_dev)
          call field_bc_used%apply_bc(field_bc_used%nodes)


        end subroutine integrate


      end module rk3tvd_class


c$$$      print *,'first step'
c$$$
c$$$          do j=1, nodes_profile(2)
c$$$             do i=1, nodes_profile(1)
c$$$                print '(''('',I2,'','',I2,'')'', 3X,
c$$$     $                  ''u1='',F10.2)',
c$$$     $               i,j,
c$$$     $               nodes_tmp(i,j,1)
c$$$             end do
c$$$          end do

c$$$         print *,'second step'
c$$$
c$$$          do j=1, nodes_profile(2)
c$$$             do i=1, nodes_profile(1)
c$$$                print '(''('',I2,'','',I2,'')'', 3X,
c$$$     $                  ''u2='',F10.2)',
c$$$     $               i,j,
c$$$     $               nodes_tmp(i,j,1)
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
