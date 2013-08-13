      module rungekutta3rdtvd_intg_class

        use field_bc_class     , only : field_bc
        use parameters_kind    , only : rkind, ikind
        use td_integrator_class, only : td_integrator


        implicit none


        type, extends(td_integrator) :: rungekutta3rdtvd_intg

          contains
          procedure, nopass :: integrate

        end type rungekutta3rdtvd_intg


        contains


        subroutine integrate(field_bc_used, dt)

          implicit none


          class(field_bc), intent(inout) :: field_bc_used
          real(rkind)    , intent(in)    :: dt


          !local variables
          real(rkind)                                   :: b2 = 3./4.
          real(rkind)                                   :: b3 = 1./3.

          integer(ikind)                                :: i,j,k
          integer(ikind), dimension(3)                  :: nodes_profile
          integer(ikind)                                :: bc_size

          real(rkind)   , dimension(:,:,:), allocatable :: nodes_tmp
          real(rkind)   , dimension(:,:,:), allocatable :: time_dev
          

          !initialization of local variables for the allocation
          !of nodes_tmp, a temporary table
          nodes_profile = field_bc_used%get_nodes_profile()
          bc_size       = field_bc_used%get_bc_size()


          !runge-kutta first step
          !u_1 = u_n + dt*d/dt(u_n)
          allocate(nodes_tmp(nodes_profile(1), nodes_profile(2), nodes_profile(3)))
          call field_bc_used%t(field_bc_used%nodes, time_dev)

          do j=bc_size+1, nodes_profile(2)-bc_size
             do i=bc_size+1, nodes_profile(1)-bc_size
                do k=1, nodes_profile(3)
                   nodes_tmp(i,j,k) =
     $                  field_bc_used%nodes(i,j,k) +
     $                  dt*time_dev(i,j,k)
                end do
             end do
          end do
            
          deallocate(time_dev)
          call field_bc_used%apply_bc(nodes_tmp)


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


      end module rungekutta3rdtvd_intg_class


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
