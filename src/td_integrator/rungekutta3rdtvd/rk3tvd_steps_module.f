      module rk3tvd_steps_module
      
        use check_data_module, only :
     $       is_real_validated

        use parameters_bf_layer, only :
     $       no_pt

        use parameters_input, only :
     $       nx,ny,ne,
     $       bc_size,
     $       debug_initialize_timedev,
     $       debug_real

        use parameters_kind, only :
     $       ikind,
     $       rkind

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
             i_min = 1
             i_max = nx
          else
             i_min = 1
             i_max = nx
          end if

          if(present(y_borders)) then
             j_min = 1
             j_max = ny
          else
             j_min = 1
             j_max = ny
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

          real(rkind), dimension(:,:,:)          , intent(inout) :: nodes
          real(rkind)                            , intent(in)    :: dt 
          real(rkind), dimension(:,:,:)          , intent(inout) :: nodes_tmp
          real(rkind), dimension(:,:,:)          , intent(in)    :: time_dev
          integer    , dimension(:,:)            , intent(in)    :: grdpts_id
          integer(ikind), dimension(2) , optional, intent(in)    :: x_borders
          integer(ikind), dimension(2) , optional, intent(in)    :: y_borders

          
          integer        :: k
          integer(ikind) :: i,j
          integer(ikind) :: x_s
          integer(ikind) :: y_s


          if(present(x_borders)) then
             x_s = x_borders(1)
          end if

          if(present(y_borders)) then
             y_s = y_borders(1)
          end if
          
          do k=1, ne
             do j=1, size(nodes,2)
                do i=1, size(nodes,1)
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

          real(rkind), dimension(:,:,:)                        , intent(inout) :: nodes
          real(rkind)                                          , intent(in)    :: dt
          real(rkind), dimension(:,:,:)                        , intent(inout) :: nodes_tmp
          real(rkind), dimension(:,:,:)                        , intent(in)    :: time_dev
          integer    , dimension(:,:)                          , intent(in)    :: grdpts_id
          integer(ikind), dimension(2)               , optional, intent(in)    :: x_borders
          integer(ikind), dimension(2)               , optional, intent(in)    :: y_borders

          
          integer        :: k
          integer(ikind) :: i,j

          real(rkind)    :: b2
          real(rkind)    :: b2_m

          integer(ikind) :: i_min
          integer(ikind) :: j_min
          integer(ikind) :: i_max
          integer(ikind) :: j_max


          !coefficients
          if(rkind.eq.8) then
             b2   = 0.75d0
             b2_m = 0.25d0
          else
             b2   = 0.75
             b2_m = 0.25
          end if

          
          !borders
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


          !integration
          do k=1, ne
             do j=j_min, j_max
                do i=i_min, i_max

                   if(grdpts_id(i,j).ne.no_pt) then

                      if(debug_initialize_timedev) then
                         if(is_real_validated(time_dev(i,j,k),debug_real,.false.)) then
                            print '(''rk3tvd_steps_module'')'
                            print '(''compute_2nd_step_nopt'')'
                            print '(''this time derivative is not computed:'')'
                            print '(''timedev('',3I4,'')'')',i,j,k
                            print '()'
                            print *, time_dev(i,j,k)
                            print '()'
                            print '(''size(nodes,1): '',I4)', size(nodes,1)
                            print '(''size(nodes,2): '',I4)', size(nodes,2)
                            print '(''x_borders: '',2I4)', x_borders
                            print '(''y_borders: '',2I4)', y_borders
                            print '()'
                            stop ''
                         end if
                      end if

                      nodes(i,j,k) =
     $                     b2*nodes_tmp(i,j,k)+
     $                     b2_m*(nodes(i,j,k)+dt*time_dev(i,j,k))
                   end if

                end do
             end do
          end do

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

          real(rkind), dimension(:,:,:)          , intent(inout) :: nodes
          real(rkind)                            , intent(in)    :: dt
          real(rkind), dimension(:,:,:)          , intent(inout) :: nodes_tmp
          real(rkind), dimension(:,:,:)          , intent(in)    :: time_dev
          integer    , dimension(:,:)            , intent(in)    :: grdpts_id
          integer(ikind), dimension(2) , optional, intent(in)    :: x_borders
          integer(ikind), dimension(2) , optional, intent(in)    :: y_borders

          integer        :: k
          integer(ikind) :: i,j

          real(rkind) :: b3
          real(rkind) :: b3_m

          integer(ikind) :: i_min
          integer(ikind) :: i_max
          integer(ikind) :: j_min
          integer(ikind) :: j_max


          !coefficients
          if(rkind.eq.8) then
             b3   = 1.0d0/3.0d0
             b3_m = 2.0d0/3.0d0
          else
             b3   = 1.0/3.0
             b3_m = 2.0/3.0
          end if


          !borders
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
          

          !integration                
          do k=1, ne
             do j=j_min,j_max
                do i=i_min,i_max

                   if(grdpts_id(i,j).ne.no_pt) then

                      if(debug_initialize_timedev) then
                         if(is_real_validated(time_dev(i,j,k),debug_real,.false.)) then
                            print '(''rk3tvd_steps_module'')'
                            print '(''compute_3rd_step_nopt'')'
                            print '(''this time derivative is not computed:'')'
                            print '(''timedev('',3I4,'')'')',i,j,k
                            print '()'
                            print *, time_dev(i,j,k)
                            print '()'
                            print '(''size(nodes,1): '',I4)', size(nodes,1)
                            print '(''size(nodes,2): '',I4)', size(nodes,2)
                            print '(''x_borders: '',2I4)', x_borders
                            print '(''y_borders: '',2I4)', y_borders
                            print '()'
                            stop ''
                         end if
                      end if
                      
                      nodes(i,j,k) =
     $                     b3*nodes_tmp(i,j,k)+
     $                     b3_m*(nodes(i,j,k)+dt*time_dev(i,j,k))

                   end if

                end do
             end do
          end do

        end subroutine compute_3rd_step_nopt

      end module rk3tvd_steps_module
