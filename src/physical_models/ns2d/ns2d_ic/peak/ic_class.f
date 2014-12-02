      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a convected peak
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a convected peak
      !
      !> @date
      !> 05_09_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_class
       
        use ns2d_parameters, only :
     $     gamma, mach_infty

        use ic_abstract_class, only :
     $     ic_abstract

        use parameters_constant, only :
     $       x_direction,
     $       y_direction,
     $       xy_direction

        use parameters_input, only :
     $       flow_direction,
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: ic


        real(rkind) :: u0_x_flow  = 1.0d0
        real(rkind) :: u0_y_flow  = 0.0d0
        real(rkind) :: u0_xy_flow = 0.5d0*SQRT(2.0d0)
        real(rkind) :: v0_x_flow  = 0.0d0
        real(rkind) :: v0_y_flow  = 1.0d0
        real(rkind) :: v0_xy_flow = 0.5d0*SQRT(2.0d0)


        !> @class ic
        !> class encapsulating operators to set the initial
        !> conditions and the conditions enforced at the edge of the
        !> computational domain for a convected peak
        !
        !> @param apply_initial_conditions
        !> set the initial conditions
        !
        !> @param get_mach_ux_infty
        !> get the mach number along the x-direction in the far field
        !
        !> @param get_mach_uy_infty
        !> get the mach number along the y-direction in the far field
        !
        !> @param get_u_in
        !> get the x-component of the velocity at the edge of the
        !> computational domain
        !
        !> @param get_v_in
        !> get the y-component of the velocity at the edge of the
        !> computational domain
        !
        !> @param get_T_in
        !> get the temperature at the edge of the computational
        !> domain
        !
        !> @param get_P_out
        !> get the pressure at the edge of the computational domain
        !---------------------------------------------------------------
        type, extends(ic_abstract) :: ic

          contains

          procedure, nopass :: apply_ic
          procedure, nopass :: get_mach_ux_infty
          procedure, nopass :: get_mach_uy_infty
          procedure, nopass :: get_u_in
          procedure, nopass :: get_v_in
          procedure, nopass :: get_T_in
          procedure, nopass :: get_P_out
          procedure, nopass :: get_far_field

        end type ic


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
        !
        !>@param x_map
        !> map of x-coordinates
        !
        !>@param y_map
        !> map of y-coordinates
        !---------------------------------------------------------------
        subroutine apply_ic(nodes,x_map,y_map)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map          

          integer(ikind) :: i,j

          real(rkind) :: x_center, y_center, amplitude, period
          real(rkind) :: u0_mean_flow, v0_mean_flow
          real(rkind) :: mass_infty, p_infty

          x_center  = 0.0d0
          y_center  = 0.0d0
          amplitude = 0.1d0
          period    = 0.5d0
          
          select case(flow_direction)

            case(x_direction)
               u0_mean_flow = u0_x_flow
               v0_mean_flow = v0_x_flow

            case(y_direction)
               u0_mean_flow = u0_y_flow
               v0_mean_flow = v0_y_flow

            case(xy_direction)
               u0_mean_flow = u0_xy_flow
               v0_mean_flow = v0_xy_flow

            case default
               print '(''ns2d_ic'')'
               print '(''ic_class/peak/apply_ic'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop

          end select

          mass_infty = 1.0d0
          p_infty = 1.0d0/(gamma*mach_infty**2)

          do j=1, size(y_map,1)
             do i=1, size(x_map,1)
                
                nodes(i,j,1) = mass_infty + 
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

        end subroutine apply_ic


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_ux_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical :: side_s

          side_s = side

          select case(flow_direction)

            case(x_direction)
               var = u0_x_flow*mach_infty/Sqrt(u0_x_flow**2+v0_x_flow**2)
               
            case(y_direction)
               var = u0_y_flow*mach_infty/Sqrt(u0_y_flow**2+v0_y_flow**2)
               
            case(xy_direction)
               var = u0_xy_flow*mach_infty/Sqrt(u0_xy_flow**2+v0_xy_flow**2)

            case default
               print '(''ns2d_ic'')'
               print '(''ic_class/peak/get_mach_ux_infty'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop

          end select

        end function get_mach_ux_infty


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_uy_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical :: side_s

          side_s = side

          select case(flow_direction)

            case(x_direction)
               var = v0_x_flow*mach_infty/Sqrt(u0_x_flow**2+v0_x_flow**2)
               
            case(y_direction)
               var = v0_y_flow*mach_infty/Sqrt(u0_y_flow**2+v0_y_flow**2)

            case(xy_direction)
               var = v0_xy_flow*mach_infty/Sqrt(u0_xy_flow**2+v0_xy_flow**2)

            case default
               print '(''ns2d_ic'')'
               print '(''ic_class/peak/get_mach_uy_infty'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop

          end select

        end function get_mach_uy_infty


        !get the x-component of the velocity enforced
        !at the edge of the computational domain
        function get_u_in(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          
          real(rkind) :: t_s,x_s,y_s

          t_s = t
          x_s = x
          y_s = y

          select case(flow_direction)

            case(x_direction)
               var = u0_x_flow
               
            case(y_direction)
               var = u0_y_flow

            case(xy_direction)
               var = u0_xy_flow

            case default
               print '(''ns2d_ic'')'
               print '(''ic_class/peak/get_u_in'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop

          end select

        end function get_u_in


        !get the y-component of the velocity enforced
        !at the edge of the computational domain
        function get_v_in(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          
          real(rkind) :: t_s,x_s,y_s

          t_s = t
          x_s = x
          y_s = y

          select case(flow_direction)

            case(x_direction)
               var = v0_x_flow

            case(y_direction)
               var = v0_y_flow

            case(xy_direction)
               var = v0_xy_flow

            case default
               print '(''ns2d_ic'')'
               print '(''ic_class/peak/get_v_in'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop

          end select

        end function get_v_in

      
        !get the temperature enforced at the edge of the
        !computational domain
        function get_T_in(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          
          real(rkind) :: t_s,x_s,y_s

          t_s = t
          x_s = x
          y_s = y

          if(rkind.eq.8) then
             var = 1.0d0
          else
             var = 1.0
          end if

        end function get_T_in


        !get the pressure enforced at the edge of the
        !computational domain
        function get_P_out(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          
          real(rkind) :: t_s,x_s,y_s

          t_s = t
          x_s = x
          y_s = y

          if(rkind.eq.8) then
             var = 1.0d0/(gamma*mach_infty**2)
          else
             var = 1.0/(gamma*mach_infty**2)
          end if

        end function get_P_out


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the governing variables imposed in the far field
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !>@param t
        !> time
        !
        !>@param x
        !> x-coordinate
        !
        !>@param y
        !> y-coordinate
        !
        !>@return var
        !> governing variables in the far-field
        !--------------------------------------------------------------
        function get_far_field(t,x,y) result(var)

            implicit none

            real(rkind)   , intent(in) :: t
            real(rkind)   , intent(in) :: x
            real(rkind)   , intent(in) :: y
            real(rkind), dimension(ne) :: var

            if(rkind.eq.8) then

               var(1) = 1.0d0
               var(2) = var(1)*get_u_in(t,x,y)
               var(3) = var(1)*get_v_in(t,x,y)
               var(4) = 0.5d0*(var(2)**2+var(3)**2)/var(1) + get_P_out(t,x,y)/(gamma-1.0d0)

            else

               var(1) = 1.0
               var(2) = var(1)*get_u_in(t,x,y)
               var(3) = var(1)*get_v_in(t,x,y)
               var(4) = 0.5*(var(2)**2+var(3)**2)/var(1) + get_P_out(t,x,y)/(gamma-1.0)

            end if

          end function get_far_field


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

      end module ic_class
