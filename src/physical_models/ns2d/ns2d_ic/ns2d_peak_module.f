      !> @file
      !> module encapsulating subroutines to compute
      !> initial conditions with a peak for
      !> the Navier-Stokes equations
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute
      !> initial conditions with a peak for
      !> the Navier-Stokes equations
      !
      !> @date
      !> 14_08_2014 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module ns2d_peak_module

        use ns2d_parameters, only :
     $       gamma, mach_infty

        use parameters_kind, only :
     $       ikind, rkind

        implicit none

        private
        public :: apply_peak_ic


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> for a steady state
        !
        !> @date
        !> 12_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !---------------------------------------------------------------
        subroutine apply_peak_ic(nodes,x_map,y_map,velocity_meanflow)

          implicit none

          real(rkind), dimension(:,:,:)      , intent(inout) :: nodes
          real(rkind), dimension(:)          , intent(in)    :: x_map
          real(rkind), dimension(:)          , intent(in)    :: y_map
          real(rkind), dimension(2), optional, intent(in)    :: velocity_meanflow
          

          integer(ikind) :: i,j

          real(rkind) :: x_center, y_center, amplitude, period
          real(rkind) :: u0_mean_flow, v0_mean_flow
          real(rkind) :: p_infty


          x_center  = 0.0d0
          y_center  = 0.0d0
          amplitude = 0.1d0
          period    = 0.5d0
          

          if(present(velocity_meanflow)) then
             u0_mean_flow = velocity_meanflow(1)
             v0_mean_flow = velocity_meanflow(2)
          else
             if(rkind.eq.8) then
                u0_mean_flow = 0.0d0
                v0_mean_flow = 1.0d0
             else
                u0_mean_flow = 0.0
                v0_mean_flow = 1.0
             end if
          end if

          p_infty = 1.0d0/(gamma*mach_infty**2)


          do j=1, size(y_map,1)
             do i=1, size(x_map,1)
                
                nodes(i,j,1) = 1.0d0 + 
     $               peak(
     $               amplitude,
     $               period,
     $               x_map(i)-x_center,
     $               y_map(j)-y_center)
                nodes(i,j,2) = nodes(i,j,1)*u0_mean_flow
                nodes(i,j,3) = nodes(i,j,1)*v0_mean_flow
                nodes(i,j,4) = 0.5d0/nodes(i,j,1)*(
     $               nodes(i,j,2)**2+nodes(i,j,3)**2)
     $               + p_infty/(gamma-1.0d0)
                   
             end do
          end do

        end subroutine apply_peak_ic


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the peak needed for the initial conditions
        !> of the NS equations
        !
        !> @date
        !> 14_08_2014 - initial version - J.L. Desmarais
        !
        !>@param amplitude
        !> amplitude of the peak
        !
        !>@param period
        !> characteristic size for the compact support of the peak
        !
        !>@param x
        !> x-coordinate identifying where the initial condition is
        !> evaluated
        !
        !>@param y
        !> y-coordinate identifying where the initial condition is
        !> evaluated
        !---------------------------------------------------------------
        function peak(amplitude, period, x, y)

          implicit none

          real(rkind), intent(in) :: amplitude
          real(rkind), intent(in) :: period
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: peak


          real(rkind) :: radius
          real(rkind) :: radius_max
          real(rkind) :: omega
          

          radius = Sqrt(x**2+y**2)

          if(rkind.eq.8) then

             radius_max = period/2.0d0
             omega      = 2.0d0*ACOS(-1.0d0)/period

             if(radius.le.radius_max) then
                peak = amplitude*(1.0d0 + cos(omega*radius))
             else
                peak = 0.0d0                
             end if

          else

             radius_max = period/2.0
             omega      = 2.0*ACOS(-1.0)/period

             if(radius.le.radius_max) then
                peak = amplitude*(1.0 + cos(omega*radius))
             else
                peak = 0.0
             end if
             
          end if

        end function peak

      end module ns2d_peak_module
