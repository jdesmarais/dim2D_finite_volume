      !> @file
      !> class encapsulating subroutines for the
      !> step computations of the Runge-Kutta 3rd
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
      !> class encapsulating subroutines for the
      !> step computations of the Runge-Kutta 3rd
      !> order time integration scheme developed in
      !> “Efficient implementation of essentially non-
      !> oscillatory shock-capturing methods”, J. Comput.
      !> Phys., 77 (1988), pp. 439-471, C.-W. Shu and
      !> S. Osher, on a parallel memory distributed system
      !
      !> @date
      !> 23_09_2013 - initial version - J.L. Desmarais
      !
      !> \f{eqnarray*}{
      !> u_1     &=& u_n + \Delta t*\frac{d u_n}{dt} \\\
      !> u_2     &=& \frac{1}{4}u_n + \frac{3}{4} \left(
      !>             u_1 + \Delta t * \frac{d u_1}{dt} \right) \\\
      !> u_{n+1} &=& \frac{1}{3}u_n + \frac{2}{3} \left(
      !>             u_2 + \Delta t * \frac{d u_2}{dt}\right) \\\
      !> \f}
      !-----------------------------------------------------------------
      module rk3tvd_steps_class

        use field_class     , only : field
        use parameters_input, only : nx,ny,ne
        use parameters_kind , only : ikind, rkind

        implicit none


        private
        public :: rk3tvd_steps


        !> @class rk3tvd_steps
        !> class encapsulating subroutines to integrate
        !> the governing equations using Runge-Kutta 3rd
        !> order time integration scheme
        !>
        !> @param b2
        !> 1st coeff for the Runge-Kutta scheme
        !> @param b3
        !> 2nd coeff for the Runge-Kutta scheme
        !---------------------------------------------------------------
        type :: rk3tvd_steps

          real(rkind) :: b2 = 0.75d0
          real(rkind) :: b3 = 1.0d0/3.0d0

          contains

          procedure, nopass :: compute_1st_step
          procedure,   pass :: compute_2nd_step
          procedure,   pass :: compute_3rd_step

        end type rk3tvd_steps


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
        !>@param bc_size
        !> size of the boundary layer
        !
        !>@param field_used
        !> object encapsulating the main variables at \f$ t=\Delta t \f$
        !
        !>@param nodes_tmp
        !> nodes at \f$ t=n \f$
        !--------------------------------------------------------------
        subroutine compute_1st_step(
     $       dt,
     $       time_dev,
     $       bc_size,
     $       field_used,
     $       nodes_tmp)

          implicit none

          real(rkind)                     , intent(in)    :: dt 
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          integer                         , intent(in)    :: bc_size
          class(field)                    , intent(inout) :: field_used
          real(rkind), dimension(nx,ny,ne), intent(inout) :: nodes_tmp

          
          integer        :: k
          integer(ikind) :: i,j

          
          do k=1, ne
             do j=bc_size+1, ny-bc_size
                do i=bc_size+1, nx-bc_size
                   nodes_tmp(i,j,k)        = field_used%nodes(i,j,k)
                   field_used%nodes(i,j,k) = field_used%nodes(i,j,k) +
     $                                       dt*time_dev(i,j,k)
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
        !>@param bc_size
        !> size of the boundary layer
        !
        !>@param field_used
        !> object encapsulating the main variables at \f$ t=\Delta t \f$
        !
        !>@param nodes_tmp
        !> nodes at \f$ t=n \f$
        !--------------------------------------------------------------
        subroutine compute_2nd_step(
     $     this,
     $     dt,
     $     time_dev,
     $     bc_size,
     $     field_used,
     $     nodes_tmp)

          implicit none

          class(rk3tvd_steps)        , intent(in)    :: this
          real(rkind)                     , intent(in)    :: dt 
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          integer                         , intent(in)    :: bc_size
          class(field)                    , intent(inout) :: field_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes_tmp

          
          integer        :: k
          integer(ikind) :: i,j

          
          if(rkind.eq.8) then
             do k=1, ne
                do j=bc_size+1, ny-bc_size
                   do i=bc_size+1, nx-bc_size
                      field_used%nodes(i,j,k) = this%b2*nodes_tmp(i,j,k)+
     $                     (1.0d0-this%b2)*(field_used%nodes(i,j,k)+
     $                     dt*time_dev(i,j,k))
                   end do
                end do
             end do
          else
             do k=1, ne
                do j=bc_size+1, ny-bc_size
                   do i=bc_size+1, nx-bc_size
                      field_used%nodes(i,j,k) = this%b2*nodes_tmp(i,j,k)+
     $                     (1.0-this%b2)*(field_used%nodes(i,j,k)+
     $                     dt*time_dev(i,j,k))
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
        !>@param bc_size
        !> size of the boundary layer
        !
        !>@param field_used
        !> object encapsulating the main variables at \f$ t=\Delta t \f$
        !
        !>@param nodes_tmp
        !> nodes at \f$ t=n \f$
        !--------------------------------------------------------------
        subroutine compute_3rd_step(
     $     this,
     $     dt,
     $     time_dev,
     $     bc_size,
     $     field_used,
     $     nodes_tmp)

          implicit none

          class(rk3tvd_steps)        , intent(in)    :: this
          real(rkind)                     , intent(in)    :: dt 
          real(rkind), dimension(nx,ny,ne), intent(in)    :: time_dev
          integer                         , intent(in)    :: bc_size
          class(field)                    , intent(inout) :: field_used
          real(rkind), dimension(nx,ny,ne), intent(in)    :: nodes_tmp

          
          integer        :: k
          integer(ikind) :: i,j

          
          if(rkind.eq.8) then

             do k=1 ,ne
                do j=bc_size+1, ny-bc_size
                   do i=bc_size+1, nx-bc_size
                      field_used%nodes(i,j,k) = this%b3*nodes_tmp(i,j,k)+
     $                     (1.0d0-this%b3)*(field_used%nodes(i,j,k)+
     $                     dt*time_dev(i,j,k))
                   end do
                end do
             end do

          else

             do k=1 ,ne
                do j=bc_size+1, ny-bc_size
                   do i=bc_size+1, nx-bc_size
                      field_used%nodes(i,j,k) = this%b3*nodes_tmp(i,j,k)+
     $                     (1.0-this%b3)*(field_used%nodes(i,j,k)+
     $                     dt*time_dev(i,j,k))
                   end do
                end do
             end do

          end if

        end subroutine compute_3rd_step

      end module rk3tvd_steps_class
