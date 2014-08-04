      module rk3tvd_steps_module

        use parameters_bf_layer, only : no_pt
        use parameters_input   , only : nx,ny,ne, bc_size
        use parameters_kind    , only : ikind, rkind

        implicit none

        private
        public :: compute_1st_step, compute_1st_step_nopt,
     $            compute_2nd_step, compute_2nd_step_nopt,
     $            compute_3rd_step, compute_3rd_step_nopt

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
        !>@param nodes
        !> array for the grid points \f$ u_1 \f$
        !
        !>@param dt
        !> time step
        !
        !>@param nodes_tmp
        !> array for the grid points \f$ u_n \f$
        !
        !>@param time_dev
        !> table containing the time derivative \f$ \frac{d u_n}{dt} \f$
        !--------------------------------------------------------------
        subroutine compute_1st_step(
     $       nodes,
     $       dt,
     $       nodes_tmp,
     $       time_dev,
     $       x_borders,
     $       y_borders)

          implicit none

          real(rkind)   , dimension(nx,ny,ne)   , intent(inout) :: nodes
          real(rkind)                           , intent(in)    :: dt 
          real(rkind)   , dimension(nx,ny,ne)   , intent(inout) :: nodes_tmp
          real(rkind)   , dimension(nx,ny,ne)   , intent(in)    :: time_dev
          integer(ikind), dimension(2), optional, intent(in)    :: x_borders
          integer(ikind), dimension(2), optional, intent(in)    :: y_borders
          
          integer        :: k
          integer(ikind) :: i,j
          integer(ikind) :: i_min, j_min, i_max, j_max


          if(present(x_borders)) then
             i_min = x_borders(1)
             i_max = x_borders(2)
          else
             i_min = bc_size+1
             i_max = nx-bc_size
          end if

          if(present(y_borders)) then
             j_min = y_borders(1)
             j_max = y_borders(2)
          else
             j_min = bc_size+1
             j_max = ny-bc_size
          end if

          
          do k=1, ne
             do j=j_min, j_max
                do i=i_min, i_max
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
        !> subroutine to compute the 1st runge-kutta step using
        !> the numerical scheme developed by C.W.Shu and S.Osher
        !> \f$ u_1 = u_n + \Delta t*\frac{d u_n}{dt}\f$
        !
        !> @date
        !> 23_09_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array for the grid points \f$ u_n \f$ (in) and
        !> \f$ u_1 \f$ (out)
        !
        !>@param dt
        !> time step
        !
        !>@param nodes_tmp
        !> array for the grid points \f$ u_n \f$ (out)
        !
        !>@param time_dev
        !> table containing the time derivative \f$ \frac{d u_n}{dt} \f$
        !        
        !>@param grdpts_id
        !> mask array identifying the role of the grid points
        !--------------------------------------------------------------
        subroutine compute_1st_step_nopt(
     $     nodes,
     $     dt,
     $     nodes_tmp,
     $     time_dev,
     $     grdpts_id,
     $     x_borders,
     $     y_borders)

          implicit none

          real(rkind), dimension(:,:,:)         , intent(inout) :: nodes
          real(rkind)                           , intent(in)    :: dt 
          real(rkind), dimension(:,:,:)         , intent(inout) :: nodes_tmp
          real(rkind), dimension(:,:,:)         , intent(in)    :: time_dev
          integer    , dimension(:,:)           , intent(in)    :: grdpts_id
          integer(ikind), dimension(2), optional, intent(in)    :: x_borders
          integer(ikind), dimension(2), optional, intent(in)    :: y_borders

          
          integer        :: k
          integer(ikind) :: i,j
          integer(ikind) :: i_min, j_min, i_max, j_max


          if(present(x_borders)) then
             i_min = x_borders(1)
             i_max = x_borders(2)
          else
             i_min = bc_size+1
             i_max = size(nodes,1)-bc_size
          end if

          if(present(y_borders)) then
             j_min = y_borders(1)
             j_max = y_borders(2)
          else
             j_min = bc_size+1
             j_max = size(nodes,2)-bc_size
          end if

          
          do k=1, ne
             do j=j_min, j_max
                do i=i_min, i_max
                   if(grdpts_id(i,j).ne.no_pt) then
                      nodes_tmp(i,j,k) = nodes(i,j,k)
                      nodes(i,j,k)     = nodes(i,j,k) + dt*time_dev(i,j,k)
                   end if
                end do
             end do
          end do

        end subroutine compute_1st_step_nopt


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to compute the 2nd runge-kutta step using
        !> the numerical scheme developed by C.W.Shu and S.Osher
        !> \f$u_2 = \frac{3}{4}u_n + \frac{1}{4} \left(
        !> u_1 + \Delta t * \frac{d u_1}{dt} \right) \f$
        !
        !> @date
        !> 23_09_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array for the grid points \f$ u_1 \f$ (in) and
        !> \f$ u_2 \f$ (out)
        !
        !>@param dt
        !> time step
        !
        !>@param nodes_tmp
        !> array for the grid points \f$ u_1 \f$ (in)
        !
        !>@param time_dev
        !> table containing the time derivative \f$ \frac{d u_1}{dt} \f$
        !--------------------------------------------------------------
        subroutine compute_2nd_step(
     $     nodes,
     $     dt,
     $     nodes_tmp,
     $     time_dev,
     $     x_borders, y_borders)

          implicit none

          real(rkind), dimension(nx,ny,ne)      , intent(inout) :: nodes
          real(rkind)                           , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne)      , intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne)      , intent(in)    :: time_dev
          integer(ikind), dimension(2), optional, intent(in)    :: x_borders
          integer(ikind), dimension(2), optional, intent(in)    :: y_borders

          
          integer        :: k
          integer(ikind) :: i,j

          real(rkind), parameter :: b2 = 0.75d0

          integer(ikind) :: i_min, j_min, i_max, j_max


          if(present(x_borders)) then
             i_min = x_borders(1)
             i_max = x_borders(2)
          else
             i_min = bc_size+1
             i_max = nx-bc_size
          end if

          if(present(y_borders)) then
             j_min = y_borders(1)
             j_max = y_borders(2)
          else
             j_min = bc_size+1
             j_max = ny-bc_size
          end if

          
          if(rkind.eq.8) then
             do k=1, ne
                do j=j_min, j_max
                   do i=i_min, i_max
                      nodes(i,j,k) =
     $                     b2*nodes_tmp(i,j,k)+
     $                     (1.0d0-b2)*(nodes(i,j,k)+dt*time_dev(i,j,k))
                   end do
                end do
             end do
          else
             do k=1, ne
                do j=j_min, j_max
                   do i=i_min, i_max
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
        !> subroutine to compute the 2nd runge-kutta step using
        !> the numerical scheme developed by C.W.Shu and S.Osher
        !> \f$u_2 = \frac{3}{4}u_n + \frac{1}{4} \left(
        !> u_1 + \Delta t * \frac{d u_1}{dt} \right) \f$
        !
        !> @date
        !> 23_09_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array for the grid points \f$ u_1 \f$ (in) and
        !> \f$ u_2 \f$ (out)
        !
        !>@param dt
        !> time step
        !
        !>@param nodes_tmp
        !> array for the grid points \f$ u_1 \f$ (in)
        !
        !>@param time_dev
        !> table containing the time derivative \f$ \frac{d u_1}{dt} \f$
        !        
        !>@param grdpts_id
        !> mask array identifying the role of the grid points
        !--------------------------------------------------------------
        subroutine compute_2nd_step_nopt(
     $     nodes,
     $     dt,
     $     nodes_tmp,
     $     time_dev,
     $     grdpts_id,
     $     x_borders,
     $     y_borders)

          implicit none

          real(rkind), dimension(:,:,:)         , intent(inout) :: nodes
          real(rkind)                           , intent(in)    :: dt
          real(rkind), dimension(:,:,:)         , intent(inout) :: nodes_tmp
          real(rkind), dimension(:,:,:)         , intent(in)    :: time_dev
          integer    , dimension(:,:)           , intent(in)    :: grdpts_id
          integer(ikind), dimension(2), optional, intent(in)    :: x_borders
          integer(ikind), dimension(2), optional, intent(in)    :: y_borders

          
          integer        :: k
          integer(ikind) :: i,j

          real(rkind), parameter :: b2 = 0.75d0

          integer(ikind) :: i_min, j_min, i_max, j_max

          if(present(x_borders)) then
             i_min = x_borders(1)
             i_max = x_borders(2)
          else
             i_min = bc_size+1
             i_max = size(nodes,1)-bc_size
          end if

          if(present(y_borders)) then
             j_min = y_borders(1)
             j_max = y_borders(2)
          else
             j_min = bc_size+1
             j_max = size(nodes,2)-bc_size
          end if

          
          if(rkind.eq.8) then
             do k=1, ne
                do j=j_min, j_max
                   do i=i_min, i_max
                      if(grdpts_id(i,j).ne.no_pt) then
                         nodes(i,j,k) =
     $                        b2*nodes_tmp(i,j,k)+
     $                        (1.0d0-b2)*(nodes(i,j,k)+dt*time_dev(i,j,k))
                      end if
                   end do
                end do
             end do
          else
             do k=1, ne
                do j=j_min, j_max
                   do i=i_min, i_max
                      if(grdpts_id(i,j).ne.no_pt) then
                         nodes(i,j,k) = 
     $                        b2*nodes_tmp(i,j,k)+
     $                        (1.0-b2)*(nodes(i,j,k)+dt*time_dev(i,j,k))
                      end if
                   end do
                end do
             end do
          end if

        end subroutine compute_2nd_step_nopt
        


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
        !>@param nodes
        !> array for the grid points \f$ u_2 \f$ (in) and
        !> \f$ u_{n+1} \f$ (out)
        !
        !>@param dt
        !> time step
        !
        !>@param nodes_tmp
        !> array for the grid points \f$ u_n \f$ (in)
        !
        !>@param time_dev
        !> table containing the time derivative \f$ \frac{d u_2}{dt} \f$
        !        
        !>@param grdpts_id
        !> mask array identifying the role of the grid points
        !--------------------------------------------------------------
        subroutine compute_3rd_step(
     $     nodes,
     $     dt,
     $     nodes_tmp,
     $     time_dev,
     $     x_borders,
     $     y_borders)

          implicit none

          real(rkind), dimension(nx,ny,ne)      , intent(inout) :: nodes
          real(rkind)                           , intent(in)    :: dt
          real(rkind), dimension(nx,ny,ne)      , intent(inout) :: nodes_tmp
          real(rkind), dimension(nx,ny,ne)      , intent(in)    :: time_dev
          integer(ikind), dimension(2), optional, intent(in)    :: x_borders
          integer(ikind), dimension(2), optional, intent(in)    :: y_borders

          integer        :: k
          integer(ikind) :: i,j

          real(rkind), parameter :: b3 = 1.0d0/3.0d0

          integer(ikind) :: i_min, j_min, i_max, j_max

          if(present(x_borders)) then
             i_min = x_borders(1)
             i_max = x_borders(2)
          else
             i_min = bc_size+1
             i_max = nx-bc_size
          end if

          if(present(y_borders)) then
             j_min = y_borders(1)
             j_max = y_borders(2)
          else
             j_min = bc_size+1
             j_max = ny-bc_size
          end if

          if(rkind.eq.8) then

             do k=1 ,ne
                do j=j_min, j_max
                   do i=i_min, i_max
                      nodes(i,j,k) =
     $                     b3*nodes_tmp(i,j,k)+
     $                     (1.0d0-b3)*(nodes(i,j,k)+dt*time_dev(i,j,k))
                   end do
                end do
             end do

          else

             do k=1 ,ne
                do j=j_min, j_max
                   do i=i_min, i_max
                      nodes(i,j,k) =
     $                     b3*nodes_tmp(i,j,k)+
     $                     (1.0-b3)*(nodes(i,j,k)+dt*time_dev(i,j,k))
                   end do
                end do
             end do

          end if

        end subroutine compute_3rd_step


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
        !>@param nodes
        !> array for the grid points \f$ u_2 \f$ (in) and
        !> \f$ u_{n+1} \f$ (out)
        !
        !>@param dt
        !> time step
        !
        !>@param nodes_tmp
        !> array for the grid points \f$ u_n \f$ (in)
        !
        !>@param time_dev
        !> table containing the time derivative \f$ \frac{d u_2}{dt} \f$
        !        
        !>@param grdpts_id
        !> mask array identifying the role of the grid points
        !--------------------------------------------------------------
        subroutine compute_3rd_step_nopt(
     $     nodes,
     $     dt,
     $     nodes_tmp,
     $     time_dev,
     $     grdpts_id,
     $     x_borders,
     $     y_borders)

          implicit none

          real(rkind), dimension(:,:,:)         , intent(inout) :: nodes
          real(rkind)                           , intent(in)    :: dt
          real(rkind), dimension(:,:,:)         , intent(inout) :: nodes_tmp
          real(rkind), dimension(:,:,:)         , intent(in)    :: time_dev
          integer    , dimension(:,:)           , intent(in)    :: grdpts_id
          integer(ikind), dimension(2), optional, intent(in)    :: x_borders
          integer(ikind), dimension(2), optional, intent(in)    :: y_borders

          integer        :: k
          integer(ikind) :: i,j

          real(rkind), parameter :: b3 = 1.0d0/3.0d0

          integer(ikind) :: i_min, j_min, i_max, j_max

          if(present(x_borders)) then
             i_min = x_borders(1)
             i_max = x_borders(2)
          else
             i_min = bc_size+1
             i_max = size(nodes,1)-bc_size
          end if

          if(present(y_borders)) then
             j_min = y_borders(1)
             j_max = y_borders(2)
          else
             j_min = bc_size+1
             j_max = size(nodes,2)-bc_size
          end if


          if(rkind.eq.8) then

             do k=1 ,ne
                do j=j_min,j_max
                   do i=i_min,i_max
                      if(grdpts_id(i,j).ne.no_pt) then
                         nodes(i,j,k) =
     $                        b3*nodes_tmp(i,j,k)+
     $                        (1.0d0-b3)*(nodes(i,j,k)+dt*time_dev(i,j,k))
                      end if
                   end do
                end do
             end do

          else

             do k=1 ,ne
                do j=j_min,j_max
                   do i=i_min,i_max
                      if(grdpts_id(i,j).ne.no_pt) then
                         nodes(i,j,k) =
     $                        b3*nodes_tmp(i,j,k)+
     $                        (1.0-b3)*(nodes(i,j,k)+dt*time_dev(i,j,k))
                      end if
                   end do
                end do
             end do

          end if

        end subroutine compute_3rd_step_nopt

      end module rk3tvd_steps_module
