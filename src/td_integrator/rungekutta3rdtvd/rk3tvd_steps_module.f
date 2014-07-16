      module rk3tvd_steps_module

        use parameters_input, only : nx,ny,ne, bc_size
        use parameters_kind , only : ikind, rkind

        implicit none

        private
        public :: compute_1st_step,
     $            compute_2nd_step,
     $            compute_3rd_step

        contains


                !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to compute the 1st runge-kutta step using
        !> the numerical scheme developed by C.W.Shu and S.Osher
        !> \f$ u_1 = u_n + \Delta t*\frac{d u_n}{dt}\f$
        !
        !> @date
        !> 23_09_2013 - initial version - J.L. Desmarais
        !
        !>@param dt
        !> time step
        !
        !>@param time_dev
        !> table containing the time derivative at \f$ t=n \f$
        !
        !>@param field_used
        !> object encapsulating the main variables at \f$ t=\Delta t \f$
        !
        !>@param nodes_tmp
        !> nodes at \f$ t=n \f$
        !--------------------------------------------------------------
        subroutine compute_1st_step(
     $     nodes,
     $     dt,
     $     nodes_tmp,
     $     time_dev)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          real(rkind)                     , intent(in)    :: dt 
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          
          integer        :: k
          integer(ikind) :: i,j

          
          do k=1, ne
             do j=bc_size+1, ny-bc_size
                do i=bc_size+1, nx-bc_size
                   nodes_tmp(i,j,k) = nodes(i,j,k)
                   nodes(i,j,k)     = nodes(i,j,k) + dt*time_dev(i,j,k)
                end do
             end do
          end do

        end subroutine compute_1st_step


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to compute the 2nd runge-kutta step using
        !> the numerical scheme developed by C.W.Shu and S.Osher
        !> \f$u_2 = \frac{1}{4}u_n + \frac{3}{4} \left(
        !> u_1 + \Delta t * \frac{d u_1}{dt} \right) \f$
        !
        !> @date
        !> 23_09_2013 - initial version - J.L. Desmarais
        !
        !>@param dt
        !> time step
        !
        !>@param time_dev
        !> table containing the time derivative at \f$ t=n \f$
        !
        !>@param field_used
        !> object encapsulating the main variables at \f$ t=\Delta t \f$
        !
        !>@param nodes_tmp
        !> nodes at \f$ t=n \f$
        !--------------------------------------------------------------
        subroutine compute_2nd_step(
     $     nodes,
     $     dt,
     $     nodes_tmp,
     $     time_dev)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          
          integer        :: k
          integer(ikind) :: i,j

          real(rkind), parameter :: b2 = 0.75d0
          
          if(rkind.eq.8) then
             do k=1, ne
                do j=bc_size+1, ny-bc_size
                   do i=bc_size+1, nx-bc_size
                      nodes(i,j,k) =
     $                     b2*nodes_tmp(i,j,k)+
     $                     (1.0d0-b2)*(nodes(i,j,k)+dt*time_dev(i,j,k))
                   end do
                end do
             end do
          else
             do k=1, ne
                do j=bc_size+1, ny-bc_size
                   do i=bc_size+1, nx-bc_size
                      nodes(i,j,k) = 
     $                     b2*nodes_tmp(i,j,k)+
     $                     (1.0-b2)*(nodes(i,j,k)+dt*time_dev(i,j,k))
                   end do
                end do
             end do
          end if

        end subroutine compute_2nd_step


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to compute the 3rd runge-kutta step using
        !> the numerical scheme developed by C.W.Shu and S.Osher
        !> \f$ u_{n+1} = \frac{1}{3}u_n + \frac{2}{3} \left(
        !>             u_2 + \Delta t * \frac{d u_2}{dt}\right) \f$
        !> @date
        !> 23_09_2013 - initial version - J.L. Desmarais
        !
        !>@param dt
        !> time step
        !
        !>@param time_dev
        !> table containing the time derivative at \f$ t=n \f$
        !
        !>@param field_used
        !> object encapsulating the main variables at \f$ t=\Delta t \f$
        !
        !>@param nodes_tmp
        !> nodes at \f$ t=n \f$
        !--------------------------------------------------------------
        subroutine compute_3rd_step(
     $     nodes,
     $     dt,
     $     nodes_tmp,
     $     time_dev)

          implicit none

          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes
          real(rkind)                     , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev


          integer        :: k
          integer(ikind) :: i,j

          real(rkind), parameter :: b3 = 1.0d0/3.0d0

          if(rkind.eq.8) then

             do k=1 ,ne
                do j=bc_size+1, ny-bc_size
                   do i=bc_size+1, nx-bc_size
                      nodes(i,j,k) =
     $                     b3*nodes_tmp(i,j,k)+
     $                     (1.0d0-b3)*(nodes(i,j,k)+dt*time_dev(i,j,k))
                   end do
                end do
             end do

          else

             do k=1 ,ne
                do j=bc_size+1, ny-bc_size
                   do i=bc_size+1, nx-bc_size
                      nodes(i,j,k) =
     $                     b3*nodes_tmp(i,j,k)+
     $                     (1.0-b3)*(nodes(i,j,k)+dt*time_dev(i,j,k))
                   end do
                end do
             end do

          end if

        end subroutine compute_3rd_step

      end module rk3tvd_steps_module
