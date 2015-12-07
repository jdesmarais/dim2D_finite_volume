      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a bubble nucleating at
      !> a wall
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a bubble nucleating at
      !> a wall
      !
      !> @date
      !> 05_06_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_class

        use dim2d_dropbubble_module, only :
     $       mass_density_ellipsoid,
     $       total_energy_ellipsoid

        use dim2d_parameters, only :
     $       cv_r,
     $       we

        use dim2d_prim_module, only :
     $       speed_of_sound

        use dim2d_state_eq_module, only :
     $       get_mass_density_liquid,
     $       get_mass_density_vapor,
     $       get_interface_length

        use ic_abstract_class, only :
     $       ic_abstract

        use parameters_constant, only :
     $       liquid,
     $       vapor,
     $       x_direction,
     $       y_direction,
     $       parabolic_profile,
     $       linear_profile

        use parameters_input, only :
     $       nx,ny,ne,
     $       T0,
     $       flow_velocity,
     $       flow_direction,
     $       flow_x_side,
     $       flow_y_side,
     $       flow_profile,
     $       phase_at_center,
     $       wall_micro_contact_angle,
     $       x_min,
     $       x_max,
     $       y_min,
     $       y_max,
     $       inflow_bubble_ac,
     $       inflow_bubble_x_center,
     $       inflow_bubble_y_center,
     $       inflow_bubble_radius

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: ic

        !flow velocities for the different flow configurations
        real(rkind), parameter :: u0_x_flow  = flow_velocity*flow_x_side                   !<@brief mean flow x-velocity if the flow is in the x-direction       
        real(rkind), parameter :: u0_y_flow  = 0.0d0                                       !<@brief mean flow x-velocity if the flow is in the y-direction       
        real(rkind), parameter :: u0_xy_flow = 0.5d0*SQRT(2.0d0)*flow_velocity*flow_x_side !<@brief mean flow x-velocity if the flow is in the diagonal direction
                                                                                                                                                                 
        real(rkind), parameter :: v0_x_flow  = 0.0d0                                       !<@brief mean flow y-velocity if the flow is in the x-direction       
        real(rkind), parameter :: v0_y_flow  = flow_velocity*flow_y_side                   !<@brief mean flow y-velocity if the flow is in the y-direction       
        real(rkind), parameter :: v0_xy_flow = 0.5d0*SQRT(2.0d0)*flow_velocity*flow_y_side !<@brief mean flow y-velocity if the flow is in the diagonal direction


        !> @class ic
        !> class encapsulating operators to set the initial
        !> conditions and the conditions enforced at the edge of the
        !> computational domain for a domain extension test
        !---------------------------------------------------------------
        type, extends(ic_abstract) :: ic

          character(18) :: name = 'bubble_nucleation' !<@brief name of the initial conditions

          contains

          procedure, nopass :: apply_ic          !<@brief set the initial conditions                                                 
          procedure, nopass :: get_mach_ux_infty !<@brief get the Mach number along the x-direction in the far-field                 
          procedure, nopass :: get_mach_uy_infty !<@brief get the Mach number along the y-direction in the far-field                 
          procedure, nopass :: get_u_in          !<@brief get the x-component of the velocity at the edge of the computational domain
          procedure, nopass :: get_v_in          !<@brief get the y-component of the velocity at the edge of the computational domain
          procedure, nopass :: get_T_in          !<@brief get the temperature at the edge of the computational domain                
          procedure, nopass :: get_P_out         !<@brief get the pressure at the edge of the computational domain                   
          procedure,   pass :: get_far_field     !<@brief get the governing variables imposed at the edge of the computational domain

        end type ic

        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> apply the initial conditions
        !> with saturated liquid everywhere or
        !> with an initial vapor bubble away from the wall
        !> in saturated liquid to see its influence
        !> on the nucleation of the second bubble (vapor-jet)
        !> \f[
        !> \begin{pmatrix} 
        !> \rho \\\ \rho u \\\ \rho v \\\ \rho E
        !> \end{pmatrix}(x,y) =
        !> \begin{pmatrix}
        !> \rho_\textrm{nucl}(x,y) \\\
        !> \rho(x,y) u_0(x,y) \\\
        !> \rho(x,y) v_0(x,y) \\\
        !> \rho(x,y) \left[ \frac{8}{3} c_v T_0 - 3 \rho(x,y) \right]
        !> + \frac{1}{2} \rho(x,y) (u_0(x,y)^2 + v_0(x,y)^2)
        !> + \frac{1}{2 \textrm{We}} | \nabla \rho(x,y) |^2
        !> \end{pmatrix}
        !> \f]
        !> 
        !> where the initial mass density field is:
        !> \f[
        !> \rho(x,y) = 
        !> \begin{cases}
        !> \rho_\textrm{liq} & \mbox{(for nucleation in homogeneous saturated liquid)} \\\
        !> \rho_\textrm{bubble}(r,r_c) & \mbox{for nucleation at the wall next to another bubble}
        !> \end{cases}
        !> \f]
        !> 
        !> For the initial conditions with the initial bubble away
        !> from the wall, we have:
        !> \f[ \rho_\textrm{bubble}(r,r_c) = \frac{\rho_\textrm{liq} + \rho_\textrm{vap}}{2}
        !> + \frac{\rho_\textrm{liq} - \rho_\textrm{vap}}{2} \tanh \left( \frac{2 (r-r_c)}{L_i} \right)\f]
        !> and \f$ L_i \f$ is the width of the interface
        !> The radius of the bubble, \f$ r_c\f$, is:
        !> \f[ r_c = \textrm{inflow\_bubble\_radius} \f]
        !> The radius coordinate is:
        !> \f[ r(x,y) = \sqrt{(x-x_c)^2 + (y-y_c)^2}\f]
        !> where the coordinates of the bubble center are:
        !> \f[ \{x_c,y_c\} = \{\textrm{inflow\_bubble\_x\_center},\textrm{inflow\_bubble\_y\_center} \} \f]
        !>
        !> The initial velocity field is given by:
        !> \f[ u_0(x,y) =
        !> \begin{cases}
        !> \displaystyle{u_{0x} \left( \frac{y - y_\textrm{min}}{y_\textrm{max} - y_\textrm{min}} \right)^\alpha} & \mbox{(for flow in the x-direction)} \\\
        !> \displaystyle{u_{0y} \left( \frac{x - x_\textrm{min}}{x_\textrm{max} - x_\textrm{min}} \right)^\alpha} & \mbox{(for flow in the y-direction)}
        !> \end{cases} \f]
        !> \f[ v_0(x,y) =
        !> \begin{cases}
        !> \displaystyle{v_{0x} \left( \frac{x - x_\textrm{min}}{x_\textrm{max} - x_\textrm{min}} \right)^\alpha} & \mbox{(for flow in the x-direction)} \\\
        !> \displaystyle{v_{0y} \left( \frac{y - y_\textrm{min}}{y_\textrm{max} - y_\textrm{min}} \right)^\alpha} & \mbox{(for flow in the y-direction)}
        !> \end{cases} \f]
        !> where
        !> \f[ \alpha =
        !> \begin{cases}
        !> 1 & \mbox{(linear profile)} \\\
        !> 2 & \mbox{(parabolic profile)}
        !> \end{cases} \f]
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data, \f$ (\rho, \rho u, \rho v ,\rho E) \f$
        !
        !>@param x_map
        !> array with the x-coordinates
        !
        !>@param y_map
        !> array with the y-coordinates                
        !--------------------------------------------------------------
        subroutine apply_ic(nodes,x_map,y_map)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          real(rkind), dimension(:)    , intent(in)    :: x_map
          real(rkind), dimension(:)    , intent(in)    :: y_map          

          
          !local variables for the droplet/bubble
          real(rkind)    :: xc,yc,a,b
          real(rkind)    :: dliq,dvap,li
          real(rkind)    :: dout

          !local variables for the vortices
          real(rkind)    :: velocity_x, velocity_y

          !local variables for the initialization
          integer(ikind) :: i,j
          real(rkind)    :: x,y

          real(rkind) :: s


          ! get the mass densities corresponding to the
          ! liquid and vapor phases for the initial
          ! temperature field
          dliq = get_mass_density_liquid(T0)
          dvap = get_mass_density_vapor(T0)

          ! get the interface length corresponding
          ! to the initial temperature field
          li = get_interface_length(T0)

          ! set the major and minor axes of the bubble ellipse
          a=inflow_bubble_radius
          b=a

          ! set the center of the droplet
          xc=inflow_bubble_x_center
          yc=inflow_bubble_y_center

          ! determine the flow velocities
          if(phase_at_center.eq.liquid) then
             dout = dvap
          else
             dout = dliq
          end if

          ! angle = wall_micro_contact_angle*ACOS(-1.0d0)/180.0d0
          ! x_pinned = 0.040

          ! initialize the mass, momentum and total energy fields
          if(inflow_bubble_ac) then

             do j=1, size(y_map,1)
                do i=1, size(x_map,1)

                   ! coordinates
                   x = x_map(i)
                   y = y_map(j)

                   ! velocities
                   velocity_x = get_velocity_x(x,y)
                   velocity_y = get_velocity_y(x,y)

                   !inflow bubble surrounding by saturated liquid
                   s = smoothing_fct(x,[xc-1.5d0*a,xc+1.5d0*a],li)*
     $                 smoothing_fct(y,[yc-1.5d0*a,yc+1.5d0*a],li)

                   nodes(i,j,1) =
     $                  mass_density_ellipsoid(
     $                  x,y,xc,yc,a,b,li,dliq,dvap,phase_at_center)*s+
     $                  dout*(1.0d0-s)

                   nodes(i,j,2) =
     $                  nodes(i,j,1)*velocity_x*s+
     $                  dout*velocity_x*(1.0d0-s)
                   
                   nodes(i,j,3) =
     $                  nodes(i,j,1)*velocity_y*s+
     $                  dout*velocity_y*(1.0d0-s)
                   
                   nodes(i,j,4) =
     $                  (total_energy_ellipsoid(
     $                  x,y,xc,yc,a,b,li,dliq,dvap,
     $                  nodes(i,j,1),T0)
     $                  +
     $                  0.5d0*nodes(i,j,1)*(velocity_x**2+velocity_y**2))*s+
     $                  (dout*(8.0d0/3.0d0*cv_r*T0-3.0d0*dout)+
     $                  0.5d0*dout*(velocity_x**2+velocity_y**2))*(1.0d0-s)

                end do
             end do

          else

             do j=1, size(y_map,1)
                do i=1, size(x_map,1)

                   ! coordinates
                   x = x_map(i)
                   y = y_map(j)

                   ! velocities
                   velocity_x = get_velocity_x(x,y)
                   velocity_y = get_velocity_y(x,y)

                   !constant field
                   nodes(i,j,1) = dliq
                   nodes(i,j,2) = nodes(i,j,1)*velocity_x
                   nodes(i,j,3) = nodes(i,j,1)*velocity_y
                   nodes(i,j,4) = 0.5d0*nodes(i,j,1)*(velocity_x**2+velocity_y**2)
     $                  + nodes(i,j,1)*(8.0d0/3.0d0*cv_r*T0-3.0d0*nodes(i,j,1))

                end do
             end do

          end if         

c$$$                !inclined bubble
c$$$                x1 = (x-x_pinned)*Sin(angle) - y*Cos(angle)
c$$$
c$$$                nodes(i,j,1) = 0.5d0*(dliq+dvap)
c$$$     $                       + 0.5d0*(dliq-dvap)*Tanh((2.0d0*x1)/li)
c$$$
c$$$                nodes(i,j,2) = nodes(i,j,1)*velocity_x
c$$$
c$$$                nodes(i,j,3) = nodes(i,j,1)*velocity_y
c$$$
c$$$                s = 0.5d0*(dliq-dvap)*(1.0d0-(Tanh((2.0d0*x1)/li))**2)*2.0d0/li
c$$$
c$$$                nodes(i,j,4) = 0.5d0*nodes(i,j,1)*(velocity_x**2+velocity_y**2)
c$$$     $                       + nodes(i,j,1)*(8.0d0/3.0d0*cv_r*T0-3.0d0*nodes(i,j,1))
c$$$     $                       + 1.0d0/We*s**2

        end subroutine apply_ic


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the Mach number imposed in the far-field
        !> for the velocity in the x-direction
        !> \f[ \textrm{Ma}_x = \frac{u_0(x_\textrm{max},y_\textrm{max})}{c}\f]
        !> where \f$u_0\f$ is the velocity of the mean flow
        !> in the x-direction and \f$c\f$ is the speed of sound
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param side
        !> left or right side
        !
        !>@return
        !> Mach number for the velocity in the x-direction,
        !> \f$ \textrm{Ma}_x \f$
        !--------------------------------------------------------------
        function get_mach_ux_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical     :: side_s
          real(rkind) :: velocity_x
          real(rkind) :: c

          side_s = side

          velocity_x = get_velocity_x(x_max,y_max)
          c          = get_speed_of_sound()

          var = velocity_x/c

        end function get_mach_ux_infty


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the Mach number imposed in the far-field
        !> for the velocity in the y-direction
        !> \f[ \textrm{Ma}_y = \frac{v_0(x_\textrm{max},y_\textrm{max})}{c}\f]
        !> where \f$v_0\f$ is the velocity of the mean flow
        !> in the y-direction and \f$c\f$ is the speed of sound
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param side
        !> left or right side
        !
        !>@return
        !> Mach number for the velocity in the y-direction,
        !> \f$ \textrm{Ma}_y \f$
        !--------------------------------------------------------------
        function get_mach_uy_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical     :: side_s
          real(rkind) :: velocity_y
          real(rkind) :: c

          side_s = side

          velocity_y = get_velocity_y(x_max,y_max)
          c          = get_speed_of_sound()

          var = velocity_y/c

        end function get_mach_uy_infty


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the value of the velocity
        !> in the x-direction imposed in the far-field
        !> \f[ u_\infty(t,x,y) = u_0(x,y) \f]
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        !>@return
        !> velocity along the x-direction imposed in the far-field,
        !> \f$ u_\infty(t,x,y) \f$
        !--------------------------------------------------------------
        function get_u_in(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          real(rkind) :: t_s,x_s,y_s

          t_s    = t
          x_s    = x
          y_s    = y

          var = get_velocity_x(x,y)

        end function get_u_in


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the value of the velocity
        !> in the y-direction imposed in the far-field
        !> \f[ v_\infty(t,x,y) = v_0(x,y) \f]
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        !>@return
        !> velocity along the y-direction imposed in the far-field,
        !> \f$ v_\infty(t,x,y) \f$
        !--------------------------------------------------------------
        function get_v_in(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          real(rkind) :: t_s,x_s,y_s

          t_s    = t
          x_s    = x
          y_s    = y

          var = get_velocity_y(x,y)

        end function get_v_in

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the value of the
        !> temperature imposed in the far-field
        !> \f[ T_\infty(t,x,y) = T_0 \f]
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        !>@return
        !> temperature imposed in the far-field,
        !> \f$ T_\infty(t,x,y) \f$
        !--------------------------------------------------------------
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

          var = T0

        end function get_T_in


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the value of the
        !> pressure imposed in the far-field
        !> \f[ P_\infty(t,x,y) =
        !> \frac{8 \rho_\textrm{liq}(T_0) T_0}{3 - \rho_\textrm{liq}(T_0)}
        !> - 3 \rho_\textrm{liq}^2(T_0) \f]
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        !>@return
        !> pressure imposed in the far-field,
        !> \f$ P_\infty(t,x,y) \f$
        !--------------------------------------------------------------
        function get_P_out(t,x,y) result(var)

          implicit none

          real(rkind), intent(in) :: t
          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: var
          
          
          real(rkind) :: t_s,x_s,y_s
          real(rkind) :: mass

          t_s = t
          x_s = x
          y_s = y

          mass = get_mass_far_field(T0)

          if(rkind.eq.8) then
             var = 8.0d0*mass*T0/(3.0d0-mass) - 3.0d0*mass**2
          else
             var = 8.0*mass*T0/(3.0-mass) - 3.0*mass**2
          end if

        end function get_P_out


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the value of the variables
        !> imposed at the edge of the computational domain
        !> depending on time and coordinates as well as the
        !> state of the object
        !> \f[ 
        !> \begin{pmatrix}
        !> \rho_\infty \\\
        !> {\rho u}_\infty \\\
        !> {\rho v}_\infty \\\
        !> {\rho E}_\infty \\\
        !> \end{pmatrix} =
        !> \begin{pmatrix}
        !> \rho_\textrm{liq}(T_0) \\\
        !> \rho_\textrm{liq}(T_0) u_0(x,y) \\\
        !> \rho_\textrm{liq}(T_0) v_0(x,y) \\\
        !> \rho_\textrm{liq}(T_0) \left[ \frac{8}{3} c_v T_0 - 3 \rho_\textrm{liq}(T_0) \right]
        !> + \frac{1}{2} \rho_\textrm{liq}(T_0) \left( u_0(x,y)^2 + v_0(x,y)^2 \right)
        !> \end{pmatrix}
        !> \f]
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !
        !>@param this
        !> object encapsulating the initial conditions and
        !> the state of the conditions imposed in the far-field
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
        !>@return
        !> variable imposed at the edge of the computational
        !> domain
        !--------------------------------------------------------------
        function get_far_field(this,t,x,y) result(var)

          implicit none

          class(ic)     , intent(in) :: this
          real(rkind)   , intent(in) :: t
          real(rkind)   , intent(in) :: x
          real(rkind)   , intent(in) :: y
          real(rkind), dimension(ne) :: var


          real(rkind) :: t_s,x_s,y_s
          real(rkind) :: mass
          real(rkind) :: velocity_x
          real(rkind) :: velocity_y
          real(rkind) :: temperature

          character(19) :: name_s
          
          t_s = t
          x_s = x
          y_s = y
          name_s = this%name

          velocity_x  = get_velocity_x(x,y)
          velocity_y  = get_velocity_y(x,y)
          temperature = T0

          mass = get_mass_far_field(temperature)

          if(rkind.eq.8) then

             var(1) = mass
             var(2) = mass*velocity_x
             var(3) = mass*velocity_y
             var(4) = 0.5d0*mass*(velocity_x**2+velocity_y**2) +
     $                mass*(8.0d0/3.0d0*cv_r*temperature-3.0d0*mass)

          else

             var(1) = mass
             var(2) = mass*velocity_x
             var(3) = mass*velocity_y
             var(4) = 0.5*mass*(velocity_x**2+velocity_y**2) +
     $                mass*(8.0/3.0*cv_r*temperature-3.0*mass)

          end if

        end function get_far_field


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the initial velocity of the mean flow
        !> in the x-direction
        !> \f[ u(x,y) =
        !> \begin{cases}
        !> \displaystyle{u_{0x} \left( \frac{y - y_\textrm{min}}{y_\textrm{max} - y_\textrm{min}} \right)^\alpha} & \mbox{(for flow in the x-direction)} \\\
        !> \displaystyle{u_{0y} \left( \frac{x - x_\textrm{min}}{x_\textrm{max} - x_\textrm{min}} \right)^\alpha} & \mbox{(for flow in the y-direction)}
        !> \end{cases} \f]
        !> where
        !> \f[ \alpha =
        !> \begin{cases}
        !> 1 & \mbox{(linear profile)} \\\
        !> 2 & \mbox{(parabolic profile)}
        !> \end{cases} \f]
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !>@param x
        !> x-coordinate
        !
        !>@param y
        !> y-coordinate
        !
        !>@return
        !> flow velocity in the x-direction at \f$\{x,y\}\f$
        !--------------------------------------------------------------
        function get_velocity_x(x,y) result(velocity_x)

          implicit none

          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: velocity_x


          select case(flow_profile)

            case(parabolic_profile)

               select case(flow_direction)

                 case(x_direction)
                    velocity_x = u0_x_flow*((y-y_min)/(y_max-y_min))**2

                 case(y_direction)
                    velocity_x = u0_y_flow*((x-x_min)/(x_max-x_min))**2

                 case default
                    print '(''bubble_spherical_cap/ic_class.f'')'
                    print '(''get_velocity_x'')'
                    print '(''flow_direction not recognized'')'
                    print '(''flow_direction: '',I2)', flow_direction
                    stop ''
                 end select

            case(linear_profile)

               select case(flow_direction)

                 case(x_direction)
                    velocity_x = u0_x_flow*((y-y_min)/(y_max-y_min))

                 case(y_direction)
                    velocity_x = u0_y_flow*((x-x_min)/(x_max-x_min))

                 case default
                    print '(''bubble_spherical_cap/ic_class.f'')'
                    print '(''get_velocity_x'')'
                    print '(''flow_direction not recognized'')'
                    print '(''flow_direction: '',I2)', flow_direction
                    stop ''
               end select
               
            case default

               print '(''bubble_spherical_cap/ic_class.f'')'
               print '(''get_velocity_x'')'
               print '(''flow_profile not recognized'')'
               print '(''flow_profile: '',I2)', flow_profile
               stop ''

          end select
          
        end function get_velocity_x

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the initial velocity of the mean flow
        !> in the y-direction
        !> \f[ v(x,y) =
        !> \begin{cases}
        !> \displaystyle{v_{0x} \left( \frac{x - x_\textrm{min}}{x_\textrm{max} - x_\textrm{min}} \right)^\alpha} & \mbox{(for flow in the x-direction)} \\\
        !> \displaystyle{v_{0y} \left( \frac{y - y_\textrm{min}}{y_\textrm{max} - y_\textrm{min}} \right)^\alpha} & \mbox{(for flow in the y-direction)}
        !> \end{cases} \f]
        !> where
        !> \f[ \alpha =
        !> \begin{cases}
        !> 1 & \mbox{(linear profile)} \\\
        !> 2 & \mbox{(parabolic profile)}
        !> \end{cases} \f]
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !>@param x
        !> x-coordinate
        !
        !>@param y
        !> y-coordinate
        !
        !>@return
        !> flow velocity in the y-direction at \f$\{x,y\}\f$
        !--------------------------------------------------------------
        function get_velocity_y(x,y) result(velocity_y)

          implicit none

          real(rkind), intent(in) :: x
          real(rkind), intent(in) :: y
          real(rkind)             :: velocity_y
          

          select case(flow_profile)

            case(parabolic_profile)

               select case(flow_direction)

                 case(x_direction)
                    velocity_y = v0_x_flow*((y-y_min)/(y_max-y_min))**2

                 case(y_direction)
                    velocity_y = v0_y_flow*((x-x_min)/(x_max-x_min))**2

                 case default
                    print '(''bubble_nucleation_at_wall/ic_class.f'')'
                    print '(''get_velocity_y'')'
                    print '(''flow_direction not recognized'')'
                    print '(''flow_direction: '',I2)', flow_direction
                    stop ''
               end select

            case(linear_profile)

               select case(flow_direction)

                 case(x_direction)
                    velocity_y = v0_x_flow*((y-y_min)/(y_max-y_min))

                 case(y_direction)
                    velocity_y = v0_y_flow*((x-x_min)/(x_max-x_min))

                 case default
                    print '(''bubble_nucleation_at_wall/ic_class.f'')'
                    print '(''get_velocity_y'')'
                    print '(''flow_direction not recognized'')'
                    print '(''flow_direction: '',I2)', flow_direction
                    stop ''
               end select
               
            case default

               print '(''bubble_spherical_cap/ic_class.f'')'
               print '(''get_velocity_y'')'
               print '(''flow_profile not recognized'')'
               print '(''flow_profile: '',I2)', flow_profile
               stop ''

          end select

        end function get_velocity_y


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the mass density imposed in the far-field
        !> depending on the phase in the initial conditions:
        !> bubble in saturated liquid or droplet in saturated
        !> vapor
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !>@param temperature
        !> temperature in the far-field, \f$T_\infty\f$
        !
        !>@return
        !> mass density in the far-field, \f$ \rho_\infty \f$
        !--------------------------------------------------------------
        function get_mass_far_field(temperature)
     $     result(mass)

          implicit none

          real(rkind), intent(in) :: temperature
          real(rkind)             :: mass

          select case(phase_at_center)
            case(vapor)
               mass = get_mass_density_liquid(temperature)
            case(liquid)
               mass = get_mass_density_vapor(temperature)
            case default
               print '(''bubble_next_to_wall/ic_class'')'
               print '(''get_speed_of_sound'')'
               print '(''phase at center not recognized'')'
               print '(''phase_at_center: '',I2)', phase_at_center
               stop ''
          end select

        end function get_mass_far_field


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the speed of sound in the far-field
        !> depending on the phase in the far-field
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !>@return
        !> speed of sound
        !--------------------------------------------------------------
        function get_speed_of_sound()
     $     result(c)

          implicit none

          real(rkind) :: c

          real(rkind), dimension(ne) :: nodes
          real(rkind)                :: mass

          mass = get_mass_far_field(T0)

          if(rkind.eq.8) then
             nodes(1) = mass
             nodes(2) = 0.0d0
             nodes(3) = 0.0d0
             nodes(4) = mass*(8.0d0/3.0d0*cv_r*T0-3.0d0*mass)
             
          else
             nodes(1) = mass
             nodes(2) = 0.0
             nodes(3) = 0.0
             nodes(4) = mass*(8.0/3.0*cv_r*T0-3.0*mass)

          end if

          c = speed_of_sound(nodes)

        end function get_speed_of_sound

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the smoothing function
        !> \f[ s(x) =
        !> \begin{cases}
        !> \frac{1}{2} \left[ 1.0 + \tanh \left(+\frac{2(x-x_\textrm{min})}{l_i} \right) \right] & \mbox{for} \quad \displaystyle{x < \frac{x_\textrm{min} + x_\textrm{max}}{2}} \\\
        !> \frac{1}{2} \left[ 1.0 + \tanh \left(-\frac{2(x-x_\textrm{max})}{l_i} \right) \right] & \mbox{otherwise}
        !> \end{cases}
        !> \f]
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !>@param x
        !> x-coordinate
        !
        !>@param x_borders
        !> x-coordinate borders where the smoothing function applies
        !> \f$ \{x_\textrm{min}, x_\textrm{max}\}\f$
        !
        !>@param li
        !> characteristic length for the transition between the two
        !> states of the smoothing function
        !
        !>@return
        !> smoothing factor
        !--------------------------------------------------------------
        function smoothing_fct(x,x_borders,li) result(s)

          implicit none

          real(rkind)              , intent(in) :: x
          real(rkind), dimension(2), intent(in) :: x_borders
          real(rkind)              , intent(in) :: li
          real(rkind)                           :: s

          if(x.le.(0.5d0*(x_borders(1)+x_borders(2)))) then
             s = 0.5d0*(1.0d0+Tanh(+2.0d0*(x-x_borders(1))/li))
          else
             s = 0.5d0*(1.0d0+Tanh(-2.0d0*(x-x_borders(2))/li))
          end if

        end function smoothing_fct

      end module ic_class
