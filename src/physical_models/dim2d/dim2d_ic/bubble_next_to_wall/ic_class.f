      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a bubble next to a wall
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a bubble next to a wall
      !
      !> @date
      !> 02_06_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_class

        use dim2d_dropbubble_module, only :
     $       mass_density_ellipsoid,
     $       total_energy_ellipsoid

        use dim2d_parameters, only :
     $       cv_r, we

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
     $       vapor

        use parameters_input, only :
     $       nx,ny,ne,
     $       T0,
     $       phase_at_center

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none

        private
        public :: ic


        !> @class ic
        !> class encapsulating operators to set the initial
        !> conditions and the conditions enforced at the edge of the
        !> computational domain for a domain extension test
        !---------------------------------------------------------------
        type, extends(ic_abstract) :: ic

          character(18) :: name = 'bubble_nextto_wall' !<@brief name of the initial conditions

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
          real(rkind)    :: s


          !get the mass densities corresponding to the
          !liquid and vapor phases for the initial
          !temperature field
          dliq = get_mass_density_liquid(T0)
          dvap = get_mass_density_vapor(T0)

          !get the interface length corresponding
          !to the initial temperature field
          li = get_interface_length(T0)

          !set the major and minor axes of the bubble ellipse
          a=3.0d0*li*Sqrt(2.0d0) !2.0d0*li
          b=a

          !set the center of the droplet
          xc=0.0d0
          yc=0.0d0

          !determine the flow velocities
          velocity_x = get_velocity_x()
          velocity_y = get_velocity_y()

          if(phase_at_center.eq.liquid) then
             dout = dvap
          else
             dout = dliq
          end if

          !initialize the mass, momentum and total energy fields
          do j=1, size(y_map,1)
             do i=1, size(x_map,1)

                ! coordinates
                x = x_map(i)
                y = y_map(j)

                ! bubble
                s = smoothing_fct(x,[xc-1.5d0*a,xc+1.5d0*a],li)*
     $              smoothing_fct(y,[yc-1.5d0*a,yc+1.5d0*a],li)

                nodes(i,j,1) =
     $               mass_density_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,phase_at_center)*s+
     $               dout*(1.0d0-s)

                nodes(i,j,2) =
     $               nodes(i,j,1)*velocity_x*s+
     $               dout*velocity_x*(1.0d0-s)

                nodes(i,j,3) =
     $               nodes(i,j,1)*velocity_y*s+
     $               dout*velocity_y*(1.0d0-s)
                
                nodes(i,j,4) =
     $               (total_energy_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,
     $               nodes(i,j,1),T0)
     $               +
     $               0.5d0*nodes(i,j,1)*(velocity_x**2+velocity_y**2))*s+
     $               (dout*(8.0d0/3.0d0*cv_r*T0-3.0d0*dout)+
     $                0.5d0*dout*(velocity_x**2+velocity_y**2))*(1.0d0-s)

             end do
          end do

        end subroutine apply_ic


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_ux_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical     :: side_s
          real(rkind) :: velocity_x
          real(rkind) :: c

          side_s = side

          velocity_x = get_velocity_x()
          c          = get_speed_of_sound()

          var = velocity_x/c

        end function get_mach_ux_infty


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_uy_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical     :: side_s
          real(rkind) :: velocity_y
          real(rkind) :: c

          side_s = side

          velocity_y = get_velocity_y()
          c          = get_speed_of_sound()

          var = velocity_y/c

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

          t_s    = t
          x_s    = x
          y_s    = y

          var = get_velocity_x()

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

          t_s    = t
          x_s    = x
          y_s    = y

          var = get_velocity_y()

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
          

          var = T0

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

          velocity_x  = get_velocity_x()
          velocity_y  = get_velocity_y()
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


        function get_velocity_x() result(velocity_x)

          implicit none

          real(rkind) :: velocity_x
          
          velocity_x = 0.0d0
          
        end function get_velocity_x

        
        function get_velocity_y() result(velocity_y)

          implicit none

          real(rkind) :: velocity_y
          
          velocity_y = 0.0d0

        end function get_velocity_y


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
