      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a domain extension test
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a domain extension test
      !
      !> @date
      !> 12_12_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_class
       
        use dim2d_dropbubble_module, only :
     $       mass_density_ellipsoid,
     $       total_energy_ellipsoid

        use dim2d_parameters, only :
     $       cv_r

        use dim2d_state_eq_module, only :
     $       get_mass_density_liquid,
     $       get_mass_density_vapor,
     $       get_interface_length

        use dim2d_prim_module, only :
     $       speed_of_sound

        use ic_abstract_class, only :
     $       ic_abstract

        use parameters_constant, only :
     $       liquid,
     $       vapor,
     $       x_direction,
     $       y_direction,
     $       xy_direction

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       T0,
     $       flow_direction,
     $       flow_x_side,
     $       flow_y_side,
     $       flow_velocity

        use parameters_kind, only :
     $       ikind,
     $       rkind        

        implicit none


        private
        public :: ic

        !set the phase at the center
        integer, parameter     :: phase_at_center = vapor !<@brief phase at the center (vapor: bubble in saturated liquid, liquid: droplet in saturated vapor)
        
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

          real(rkind), dimension(ne) :: far_field !<@brief governing variables imposed in the far field

          contains

          procedure, nopass :: apply_ic          !<@brief set the initial conditions                                                 
          procedure, nopass :: get_mach_ux_infty !<@brief get the Mach number along the x-direction in the far-field                 
          procedure, nopass :: get_mach_uy_infty !<@brief get the Mach number along the y-direction in the far-field                 
          procedure, nopass :: get_u_in          !<@brief get the x-component of the velocity at the edge of the computational domain
          procedure, nopass :: get_v_in          !<@brief get the y-component of the velocity at the edge of the computational domain
          procedure, nopass :: get_T_in          !<@brief get the temperature at the edge of the computational domain                
          procedure, nopass :: get_P_out         !<@brief get the pressure at the edge of the computational domain                   
                                                 
          procedure,   pass :: ini_far_field     !<@brief initialize the governing variables imposed at the edge of the computational domain
          procedure,   pass :: get_far_field     !<@brief get the governing variables imposed at the edge of the computational domain
          procedure,   pass :: set_far_field     !<@brief set the governing variables imposed at the edge of the computational domain

        end type ic


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> apply the initial conditions
        !> with a bubble located in the center and whose
        !> radius is twice the width of the interface at
        !> equilibrium
        !> \f[
        !> \begin{pmatrix} 
        !> \rho \\\ \rho u \\\ \rho v \\\ \rho E
        !> \end{pmatrix}(x,y) =
        !> \begin{pmatrix} 
        !> \rho_\textrm{bubble}(\sqrt{x^2+y^2},2L_i) \\\
        !> \rho(x,y) u_0 \\\
        !> \rho(x,y) v_0 \\\
        !> \rho(x,y) \left[ \frac{8}{3} c_v T_0 - 3 \rho(x,y) \right]
        !> + \frac{1}{2} \rho(x,y) (u_0^2 + v_0^2)
        !> + \frac{1}{2 \textrm{We}} | \nabla \rho(x,y) |^2
        !> \end{pmatrix}
        !> \f]
        !> where
        !> \f[ \rho_\textrm{bubble}(r,r_c) = \frac{\rho_\textrm{liq} + \rho_\textrm{vap}}{2}
        !> + \frac{\rho_\textrm{liq} - \rho_\textrm{vap}}{2} \tanh \left( \frac{2 (r-r_c)}{L_i} \right)\f]
        !> and \f$ L_i \f$ is the width of the interface.
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

          !local variables for the vortices
          real(rkind)    :: velocity_x, velocity_y

          !local variables for the initialization
          integer(ikind) :: i,j
          real(rkind)    :: x,y

          
          !set the center of the droplet
          xc=0.
          yc=0.


          !get the mass densities corresponding to the
          !liquid and vapor phases for the initial
          !temperature field
          dliq = get_mass_density_liquid(T0)
          dvap = get_mass_density_vapor(T0)

          !get the interface length corresponding
          !to the initial temperature field
          li = get_interface_length(T0)

          !set the major and minor axes of the bubble ellipse
          a=2.0d0*li
          b=a !a/2.0d0          


          !determine the flow velocities
          velocity_x = get_velocity_x()
          velocity_y = get_velocity_y()


          !initialize the mass, momentum and total energy fields
          do j=1, size(y_map,1)
             do i=1, size(x_map,1)

                x = x_map(i)
                y = y_map(j)

                nodes(i,j,1) = mass_density_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,phase_at_center)

                nodes(i,j,2) = nodes(i,j,1)*velocity_x

                nodes(i,j,3) = nodes(i,j,1)*velocity_y
                
                nodes(i,j,4) = total_energy_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,
     $               nodes(i,j,1),T0)
     $               +
     $               0.5d0*nodes(i,j,1)*(velocity_x**2+velocity_y**2)

             end do
          end do

        end subroutine apply_ic


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the Mach number imposed in the far-field
        !> for the velocity in the x-direction
        !> \f[ \textrm{Ma}_x = \frac{u_\infty}{c}\f]
        !> where \f$c\f$ is the speed of sound
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

          velocity_x = get_velocity_x()
          c          = get_speed_of_sound()

          var = velocity_x/c

        end function get_mach_ux_infty


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the Mach number imposed in the far-field
        !> for the velocity in the y-direction
        !> \f[ \textrm{Ma}_y = \frac{v_\infty}{c}\f]
        !> where \f$c\f$ is the speed of sound
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

          velocity_y = get_velocity_y()
          c          = get_speed_of_sound()

          var = velocity_y/c

        end function get_mach_uy_infty


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the value of the velocity
        !> in the x-direction imposed in the far-field
        !> \f[ u_\infty(t,x,y) = u_0 \f]
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

          var = get_velocity_x()

        end function get_u_in


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the value of the velocity
        !> in the y-direction imposed in the far-field
        !> \f[ v_\infty(t,x,y) = v_0 \f]
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

          var = get_velocity_y()

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
        !> \f[ P_\infty(t,x,y) = P_0 \f]
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

          mass = get_mass_far_field()

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
        !> initialize get the governing variables imposed in the far field
        !> \f[ 
        !> \begin{pmatrix}
        !> \rho_\infty \\\
        !> {\rho u}_\infty \\\
        !> {\rho v}_\infty \\\
        !> {\rho E}_\infty \\\
        !> \end{pmatrix} =
        !> \begin{pmatrix}
        !> \rho_\textrm{liq} \\\
        !> \rho_\textrm{liq} u_0 \\\
        !> \rho_\textrm{liq} v_0 \\\
        !> \rho_\textrm{liq} \left[ \frac{8}{3} c_v T_0 - 3 \rho_\textrm{liq} \right]
        !> + \frac{1}{2} \rho_\textrm{liq} \left( u_0^2 + v_0^2 \right)
        !> \end{pmatrix}
        !> \f]
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the initial conditions and
        !> the state of the conditions imposed in the far-field
        !--------------------------------------------------------------
        subroutine ini_far_field(this)

          implicit none

          class(ic), intent(inout) :: this

          real(rkind), dimension(ne) :: var
          
          real(rkind) :: mass
          real(rkind) :: velocity_x
          real(rkind) :: velocity_y
          

          mass       = get_mass_far_field()
          velocity_x = get_velocity_x()
          velocity_y = get_velocity_y()

          if(rkind.eq.8) then

             var(1) = mass
             var(2) = mass*velocity_x
             var(3) = mass*velocity_y
             var(4) = 0.5d0*mass*(velocity_x**2+velocity_y**2) +
     $                mass*(8.0d0/3.0d0*cv_r*T0-3.0d0*mass)

          else

             var(1) = mass
             var(2) = mass*velocity_x
             var(3) = mass*velocity_y
             var(4) = 0.5*mass*(velocity_x**2+velocity_y**2) +
     $                mass*(8.0/3.0*cv_r*T0-3.0*mass)

          end if

          this%far_field = var

        end subroutine ini_far_field


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
        !> \rho_\textrm{liq} \\\
        !> \rho_\textrm{liq} u_0 \\\
        !> \rho_\textrm{liq} v_0 \\\
        !> \rho_\textrm{liq} \left[ \frac{8}{3} c_v T_0 - 3 \rho_\textrm{liq} \right]
        !> + \frac{1}{2} \rho_\textrm{liq} \left( u_0^2 + v_0^2 \right)
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
        function get_far_field(this,t,x,y) result(far_field)
        
          implicit none

          class(ic)                 , intent(in) :: this
          real(rkind)               , intent(in) :: t
          real(rkind)               , intent(in) :: x
          real(rkind)               , intent(in) :: y
          real(rkind), dimension(ne)             :: far_field

          real(rkind) :: t_s,x_s,y_s

          t_s = t
          x_s = x
          y_s = y

          far_field = this%far_field

        end function get_far_field


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> set the value of the variables
        !> imposed at the edge of the computational domain
        !> depending on time and coordinates as well as the
        !> state of the object
        !> \f[ 
        !> \begin{pmatrix}
        !> \rho_\infty \\\
        !> {\rho u}_\infty \\\
        !> {\rho v}_\infty \\\
        !> {\rho E}_\infty \\\
        !> \end{pmatrix}
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !>@param this
        !> object encapsulating the initial conditions and
        !> the state of the conditions imposed in the far-field
        !
        !>@param far_field
        !> governing variables imposed in the far-field
        !--------------------------------------------------------------
        subroutine set_far_field(this,far_field)
        
          implicit none

          class(ic)                 , intent(inout) :: this
          real(rkind), dimension(ne), intent(in)    :: far_field

          this%far_field = far_field

        end subroutine set_far_field


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the initial velocity of the mean flow
        !> in the x-direction
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !>@return
        !> mean flow velocity in the x-direction
        !--------------------------------------------------------------
        function get_velocity_x() result(velocity_x)

          implicit none

          real(rkind) :: velocity_x
          
          select case(flow_direction)
            case(x_direction)
               velocity_x = u0_x_flow
            case(y_direction)
               velocity_x = u0_y_flow
            case(xy_direction)
               velocity_x = u0_xy_flow
            case default
               print '(''newgrdpt_test/ic_class.f'')'
               print '(''apply_ic'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop ''
          end select
          
        end function get_velocity_x

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the initial velocity of the mean flow
        !> in the y-direction
        !
        !> @date
        !> 03_12_2014 - initial version - J.L. Desmarais
        !
        !>@return
        !> mean flow velocity in the y-direction
        !--------------------------------------------------------------
        function get_velocity_y() result(velocity_y)

          implicit none

          real(rkind) :: velocity_y
          
          select case(flow_direction)
            case(x_direction)
               velocity_y = v0_x_flow
            case(y_direction)
               velocity_y = v0_y_flow
            case(xy_direction)
               velocity_y = v0_xy_flow
            case default
               print '(''newgrdpt_test/ic_class.f'')'
               print '(''apply_ic'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
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
        !>@return
        !> mean flow velocity in the y-direction
        !--------------------------------------------------------------
        function get_mass_far_field()
     $     result(mass)

          implicit none

          real(rkind) :: mass

          select case(phase_at_center)
            case(vapor)
               mass = get_mass_density_liquid(T0)
            case(liquid)
               mass = get_mass_density_vapor(T0)
            case default
               print '(''bubble_transported/ic_class'')'
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

          mass = get_mass_far_field()

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

      end module ic_class
