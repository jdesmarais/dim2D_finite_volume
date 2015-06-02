      !> @file
      !> ghost cells for wall boundary conditions
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute the
      !> ghost cells for the wall equilibrium boundary
      !> conditions
      !
      !> @date
      !> 02_06_2015 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module wall_xy_equilibrium_module

        use check_data_module, only :
     $       is_real_validated

        use dim2d_parameters, only :
     $       cv_r,
     $       We,
     $       Pr

        use dim2d_state_eq_module, only :
     $       get_mass_density_vapor,
     $       get_mass_density_liquid,
     $       get_surface_tension

        use parameters_kind, only :
     $       ikind,
     $       rkind

        use ridders_method_module, only :
     $       root_fct_abstract,
     $       get_root_ridder_method

        private
        public ::
     $       dmddx,
     $       temperature,
     $       md_average,
     $       temperature_average,
     $       dwallInternalEnergy_dmd,
     $       wall_x_equilibrium_root_fct,
     $       wall_x_root_fct,
     $       get_wall_x_root_brackets,
     $       get_wall_x_md_ghost_cell


        !> @class wall_x_root_fct
        !> object encapsulating the function computing 
        !> the wall equilibrium mass density ghost value
        !
        !>@param md_x
        !> mass densities at {i-1,j} and {i+1,j}
        !
        !>@param md_y
        !> mass densities at {i,j-1} and {i,j+1}
        !
        !>@param velocity_x
        !> velocity along the x-direction at {i,j}
        !
        !>@param velocity_y
        !> velocity along the y-direction at {i,j}
        !
        !>@param Ed
        !> total energy density at {i,j}
        !
        !>@param delta_x
        !> space step along the x-direction
        !
        !>@param delta_y
        !> space step along the y-direction
        !
        !>@param micro_angle
        !> micro contact angle
        !
        !>@param wall_heat_flux
        !> heat flux at the wall
        !
        !> @param ini
        !> initialization of the attributes determining
        !> the equilibrium function at the wall
        !
        !> @param f
        !> function whose root determine the equilibrium
        !> value of the mass density at the wall
        !------------------------------------------------------------
        type, extends(root_fct_abstract) :: wall_x_root_fct

          real(rkind), dimension(2) :: md_x
          real(rkind), dimension(2) :: md_y
          real(rkind)               :: velocity_x
          real(rkind)               :: velocity_y
          real(rkind)               :: Ed
          real(rkind)               :: delta_x
          real(rkind)               :: delta_y
          real(rkind)               :: micro_angle
          real(rkind)               :: wall_heat_flux

          contains

          procedure, pass :: ini => wall_x_root_ini
          procedure, pass :: f   => wall_x_root_f

        end type wall_x_root_fct


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the mass density for the ghost cell
        !> that fits to the wall equilibrium function
        !
        !> @date
        !> 02_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md_x
        !> mass densities at {i-1,j} and {i+1,j}
        !
        !>@param md_y
        !> mass densities at {i,j-1} and {i,j+1}
        !
        !>@param velocity_x
        !> velocity along the x-direction at {i,j}
        !
        !>@param velocity_y
        !> velocity along the y-direction at {i,j}
        !
        !>@param Ed
        !> total energy density at {i,j}
        !
        !>@param delta_x
        !> space step along the x-direction
        !
        !>@param delta_y
        !> space step along the y-direction
        !
        !>@param micro_angle
        !> micro contact angle
        !
        !>@param wall_heat_flux
        !> heat flux at the wall
        !--------------------------------------------------------------
        function get_wall_x_md_ghost_cell(
     $       T_guess,
     $       md_guess,
     $       md_x,
     $       md_y,
     $       velocity_x,
     $       velocity_y,
     $       Ed,
     $       delta_x,
     $       delta_y,
     $       micro_angle,
     $       wall_heat_flux)
     $       result(md)

          implicit none

          real(rkind)              , intent(in) :: T_guess
          real(rkind)              , intent(in) :: md_guess
          real(rkind), dimension(2), intent(in) :: md_x          
          real(rkind), dimension(2), intent(in) :: md_y          
          real(rkind)              , intent(in) :: velocity_x    
          real(rkind)              , intent(in) :: velocity_y    
          real(rkind)              , intent(in) :: Ed            
          real(rkind)              , intent(in) :: delta_x       
          real(rkind)              , intent(in) :: delta_y       
          real(rkind)              , intent(in) :: micro_angle   
          real(rkind)              , intent(in) :: wall_heat_flux
          real(rkind)                           :: md


          type(wall_x_root_fct)     :: wall_x_root_fct_used
          real(rkind), dimension(2) :: md_brackets


          ! initialize the parameters configurating 
          ! the wall equilibrium function whose root
          ! determines the mass density in the ghost
          ! cell
          call wall_x_root_fct_used%ini(
     $         md_x,
     $         md_y,
     $         velocity_x,
     $         velocity_y,
     $         Ed,
     $         delta_x,
     $         delta_y,
     $         micro_angle,
     $         wall_heat_flux)


          ! determine the brackets for the root
          ! solution of the wall equilibrium
          ! function
          md_brackets = get_wall_x_root_brackets(
     $         wall_x_root_fct_used,
     $         T_guess,
     $         md_guess)


          ! determine the root of the wall
          ! equilibrium function using
          ! Ridder's method
          md = get_root_ridder_method(
     $         wall_x_root_fct_used,
     $         md_brackets(1),
     $         md_brackets(2),
     $         1e-12)

        end function get_wall_x_md_ghost_cell


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the brackets around the root
        !> of the wall equilibrium function
        !
        !> @date
        !> 02_06_2015 - initial version - J.L. Desmarais
        !
        !>@param wall_x_root_fct_used
        !> wall equilibirum function
        !
        !>@param T_guess
        !> guess for the temperature at the wall
        !
        !>@param md_guess
        !> guess for the mass density at the wall
        !
        !>@param md_brackets
        !> brackets for the root
        !--------------------------------------------------------------
        function get_wall_x_root_brackets(
     $     wall_x_root_fct_used,
     $     T_guess,
     $     md_guess)
     $     result(md_brackets)

          implicit none

          class(root_fct_abstract) , intent(in) :: wall_x_root_fct_used
          real(rkind)              , intent(in) :: T_guess
          real(rkind)              , intent(in) :: md_guess
          real(rkind), dimension(2)             :: md_brackets

          
          real(rkind) :: xl
          real(rkind) :: xh
          real(rkind) :: fl
          real(rkind) :: fh
          real(rkind) :: fmd

          real(rkind) :: dmd_l
          real(rkind) :: dmd_r
          integer     :: j

          integer, parameter :: MAXIT = 99


          ! use as first potential brackets the
          ! vapor and liquid mass densities at the
          ! temperature T_guess
          xl = get_mass_density_vapor(T_guess)
          xh = get_mass_density_liquid(T_guess)
          
          fl = wall_x_root_fct_used%f(xl)
          fh = wall_x_root_fct_used%f(xh)


          ! if the brackets proposed do not lead to
          ! opposite signs of the wall equilibrium 
          ! function, new brackets should be found
          if( is_real_validated(sign(1.0d0,fl*fh),1.0d0,.false.) ) then

             dmd_l = (md_guess)/real(MAXIT+1)
             dmd_r = (3.0d0-md_guess)/real(MAXIT+1)
             fmd   = wall_x_root_fct_used%f(md_guess)

             do j=1, MAXIT

                xl = md_guess - j*dmd_l
                xh = md_guess + j*dmd_r

                fl = wall_x_root_fct_used%f(xl)
                fh = wall_x_root_fct_used%f(xh)

                if( is_real_validated(sign(1.0d0,fl*fmd),-1.0d0,.false.) ) then

                   md_brackets(1) = xl
                   md_brackets(2) = md_guess
                   return

                else if ( is_real_validated(sign(1.0d0,fh*fmd),-1.0d0,.false.) ) then

                   md_brackets(1) = md_guess
                   md_brackets(2) = xh
                   return

                end if

             end do

             print '(''wall_xy_equilibrium_module'')'
             print '(''get_wall_x_root_brackets'')'
             print '(''exceeded maximum iterations'')'
             stop ''

          else

             md_brackets(1) = xl
             md_brackets(2) = xh

          end if

        end function get_wall_x_root_brackets


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> initialize the additional parameters
        !> determining the wall equilibrium
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> additional parameters to determine the wall
        !> equilibrium
        !
        !>@param md_x
        !> mass densities at {i-1,j} and {i+1,j}
        !
        !>@param md_y
        !> mass densities at {i,j-1} and {i,j+1}
        !
        !>@param velocity_x
        !> velocity along the x-direction at {i,j}
        !
        !>@param velocity_y
        !> velocity along the y-direction at {i,j}
        !
        !>@param Ed
        !> total energy density at {i,j}
        !
        !>@param delta_x
        !> space step along the x-direction
        !
        !>@param delta_y
        !> space step along the y-direction
        !
        !>@param micro_angle
        !> micro contact angle
        !
        !>@param wall_heat_flux
        !> heat flux at the wall
        !--------------------------------------------------------------
        subroutine wall_x_root_ini(
     $       this,
     $       md_x,
     $       md_y,
     $       velocity_x,
     $       velocity_y,
     $       Ed,
     $       delta_x,
     $       delta_y,
     $       micro_angle,
     $       wall_heat_flux)

          implicit none

          class(wall_x_root_fct)   , intent(inout) :: this
          real(rkind), dimension(2), intent(in)    :: md_x          
          real(rkind), dimension(2), intent(in)    :: md_y          
          real(rkind)              , intent(in)    :: velocity_x    
          real(rkind)              , intent(in)    :: velocity_y    
          real(rkind)              , intent(in)    :: Ed            
          real(rkind)              , intent(in)    :: delta_x       
          real(rkind)              , intent(in)    :: delta_y       
          real(rkind)              , intent(in)    :: micro_angle   
          real(rkind)              , intent(in)    :: wall_heat_flux

          
          this%md_x           = md_x          
          this%md_y           = md_y          
          this%velocity_x     = velocity_x    
          this%velocity_y     = velocity_y    
          this%Ed             = Ed            
          this%delta_x        = delta_x       
          this%delta_y        = delta_y       
          this%micro_angle    = micro_angle   
          this%wall_heat_flux = wall_heat_flux

        end subroutine wall_x_root_ini


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function whose root determines the
        !> equilibrium at the wall
        !
        !> @date
        !> 02_06_2015 - initial version - J.L. Desmarais
        !
        !>@param this
        !> additional parameters needed to evaluate the
        !> wall equilibrium
        !
        !>@param x
        !> parameter for the non-linear equation
        !
        !>@param fx
        !> evaluation of the non-linear equation at x
        !--------------------------------------------------------------
        function wall_x_root_f(this,x) result(fx)

          implicit none

          class(wall_x_root_fct), intent(in) :: this
          real(rkind)           , intent(in) :: x
          real(rkind)                        :: fx

          fx = wall_x_equilibrium_root_fct(
     $         x,
     $         this%md_x,
     $         this%md_y,
     $         this%velocity_x,
     $         this%velocity_y,
     $         this%Ed,
     $         this%delta_x,
     $         this%delta_y,
     $         this%micro_angle,
     $         this%wall_heat_flux)

        end function wall_x_root_f


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the wall equilibrium function controlling
        !> the value of the mass density for the ghost cell
        !> \f$ \frac{1}{We} \frac{\partial \rho}{\partial x} 
        !> + \frac{d E^S_{\textrm{wall}}(T,\rho)}{d \rho} \f$
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density at {i,j}
        !>
        !>@param md_x
        !> mass densities at {i-1,j} and {i+1,j}
        !>
        !>@param md_y
        !> mass densities at {i,j-1} and {i,j+1}
        !>
        !>@param velocity_x
        !> velocity along the x-direction at {i,j}
        !>
        !>@param velocity_y
        !> velocity along the y-direction at {i,j}
        !>
        !>@param Ed
        !> total energy density at {i,j}
        !>
        !>@param delta_x
        !> space step along the x-direction
        !>
        !>@param delta_y
        !> space step along the y-direction
        !>
        !>@param micro_angle
        !> micro contact angle
        !>
        !>@param wall_heat_flux
        !> heat flux at the wall
        !
        !>@param equilibrium_fct
        !> function controlling the equilibrium at the wall
        !--------------------------------------------------------------
        function wall_x_equilibrium_root_fct(
     $       md,
     $       md_x,
     $       md_y,
     $       velocity_x,
     $       velocity_y,
     $       Ed,
     $       delta_x,
     $       delta_y,
     $       micro_angle,
     $       wall_heat_flux)
     $       result(equilibrium_fct)

          implicit none

          real(rkind)              , intent(in) :: md
          real(rkind), dimension(2), intent(in) :: md_x
          real(rkind), dimension(2), intent(in) :: md_y
          real(rkind)              , intent(in) :: velocity_x
          real(rkind)              , intent(in) :: velocity_y
          real(rkind)              , intent(in) :: Ed
          real(rkind)              , intent(in) :: delta_x
          real(rkind)              , intent(in) :: delta_y
          real(rkind)              , intent(in) :: micro_angle
          real(rkind)              , intent(in) :: wall_heat_flux
          real(rkind)                           :: equilibrium_fct

          real(rkind) :: md_grad_x
          real(rkind) :: md_av
          real(rkind) :: md_grad_squared
          real(rkind) :: temperature1
          real(rkind) :: T_average

          md_grad_x =
     $         dmddx(delta_x,md_x(1),md_x(2))

          md_av =
     $         md_average(md,md_x(2))

          md_grad_squared =
     $         md_grad_x**2 +
     $         dmddx(delta_y,md_y(1),md_y(2))**2

          temperature1 = temperature(
     $         md,
     $         md_grad_squared,
     $         velocity_x,
     $         velocity_y,
     $         Ed)

          T_average = temperature_average(
     $         temperature1,
     $         delta_x,
     $         wall_heat_flux)

          equilibrium_fct =
     $         1.0d0/We*md_grad_x +
     $         dwallInternalEnergy_dmd(md_av,T_average,micro_angle)

        end function wall_x_equilibrium_root_fct


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the derivative of the wall internal energy
        !> per unit surface
        !> \f$ \frac{d E_{\textrm{wall}}}{d \rho} =
        !> \frac{6 (\md - \md_{\textrm{l}})(\md - \md_{\textrm{v}})}
        !> { \left( \rho_{\textrm{l}} - \rho_{\textrm{v}} \right)^3}
        !> \sigma \Cos(\theta_m) \f$
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density
        !>
        !>@param temperature
        !> temperature
        !>
        !>@param micro_angle
        !> microscopic contact angle
        !--------------------------------------------------------------
        function dwallInternalEnergy_dmd(md,temperature,micro_angle)

          implicit none

          real(rkind), intent(in) :: md
          real(rkind), intent(in) :: temperature
          real(rkind), intent(in) :: micro_angle
          real(rkind)             :: dwallInternalEnergy_dmd

          real(rkind) :: md_vap
          real(rkind) :: md_liq
          real(rkind) :: surface_tension
          real(rkind) :: angle

          md_vap = get_mass_density_vapor(temperature)
          md_liq = get_mass_density_liquid(temperature)
          
          surface_tension = get_surface_tension(temperature)

          if(rkind.eq.8) then
             angle = micro_angle*ACOS(-1.0d0)/180.0d0
             dwallInternalEnergy_dmd = 6.0d0*(md-md_liq)*(md-md_vap)/((md_liq-md_vap)**3)*surface_tension*Cos(angle)
          else
             angle = micro_angle*ACOS(-1.0)/180.0
             dwallInternalEnergy_dmd =   6.0*(md-md_liq)*(md-md_vap)/((md_liq-md_vap)**3)*surface_tension*Cos(angle)
          end if

        end function dwallInternalEnergy_dmd


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the mass-density at midpoint $i+1/2$:
        !> \f$ \rho_{i+1/2} = \frac{\rho_i+\rho_{i+1}}{2} \f$
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param temperature0
        !> temperature at {i}
        !>
        !>@param delta_x
        !> space step
        !>
        !>@param wall heat flux
        !> heat flux at the wall
        !
        !>@param temperature_average
        !> temperature at the wall
        !--------------------------------------------------------------
        function temperature_average(temperature0,delta_x,wall_heat_flux)

          implicit none

          real(rkind), intent(in) :: temperature0
          real(rkind), intent(in) :: delta_x
          real(rkind), intent(in) :: wall_heat_flux
          real(rkind)             :: temperature_average

          
          if(rkind.eq.8) then
             temperature_average =
     $            temperature0 +
     $            0.5d0*delta_x*Pr*wall_heat_flux
          else
             temperature_average =
     $            temperature0 +
     $            0.5*delta_x*Pr*wall_heat_flux
          end if

        end function temperature_average


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the mass-density at midpoint $i+1/2$:
        !> \f$ \rho_{i+1/2} = \frac{\rho_i+\rho_{i+1}}{2} \f$
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md0
        !> mass density at {i}
        !>
        !>@param md1
        !> mass density at {i+1}
        !>
        !>@param md_average
        !> mass density at mid-point
        !--------------------------------------------------------------
        function md_average(md0,md1)

          implicit none

          real(rkind), intent(in) :: md0
          real(rkind), intent(in) :: md1
          real(rkind)             :: md_average

          if(rkind.eq.8) then
             md_average = 0.5d0*(md0+md1)
          else
             md_average = 0.5*(md0+md1)
          end if

        end function md_average


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the temperature at $i$:
        !> \f$ T_i = \frac{3}{8 cv}\left[
        !> \frac{\rho E_i}{\rho_i} - \frac{1}{2} \left(u_i^2 + v_i^2
        !> \right) - \frac{1}{2 We} \frac{(\nabla \md)^2}{\md_i}
        !> - 3 \rho_i \right] \f$
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density at {i,j}
        !>
        !>@param md_grad_squared
        !> mass density gradient at {i,j} squared
        !>
        !>@param velocity_x
        !> velocity along x-direction at {i,j}
        !>
        !>@param velocity-y
        !> velocity along y-direction at {i,j}
        !
        !>@param Ed
        !> total energy density at {i,j}
        !--------------------------------------------------------------
        function temperature(md,md_grad_squared,velocity_x,velocity_y,Ed)
        
          implicit none

          real(rkind), intent(in) :: md
          real(rkind), intent(in) :: md_grad_squared
          real(rkind), intent(in) :: velocity_x
          real(rkind), intent(in) :: velocity_y
          real(rkind), intent(in) :: Ed
          real(rkind)             :: temperature
          
          if(rkind.eq.8) then
             
             temperature = 3.0d0/(8.0d0*cv_r)*(
     $            Ed/md - 0.5d0*(velocity_x**2+velocity_y**2)
     $            - 0.5d0/We*md_grad_squared/md - 3.0d0*md)

          else
             
             temperature = 3.0/(8.0*cv_r)*(
     $            Ed/md - 0.5*(velocity_x**2+velocity_y**2)
     $            - 0.5/We*md_grad_squared/md - 3.0*md)

          end if

        end function temperature
      

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> determine the gradient of mass density at $i$:
        !> \f$ \frac{d \rho}{d x}\bigg|_i = \frac{-\rho_{i-1} + \rho_{i+1}}
        !> {2 \Delta x} \f$
        !
        !> @date
        !> 01_06_2015 - initial version - J.L. Desmarais
        !
        !>@param delta_x
        !> space step
        !>
        !>@param md0
        !> mass density evaluated at i-1
        !>
        !>@param md2
        !> mass_density evaluated at i+1
        !>
        !>@param drhodx
        !> mass density gradient at $i$
        !--------------------------------------------------------------
        function dmddx(delta_x,md0,md2)

          implicit none

          real(rkind), intent(in) :: delta_x
          real(rkind), intent(in) :: md0
          real(rkind), intent(in) :: md2
          real(rkind)             :: dmddx
          
          if(rkind.eq.8) then
             dmddx = (md2-md0)/(2.0d0*delta_x)
          else
             dmddx = (md2-md0)/(2.0*delta_x)
          end if

        end function dmddx

      end module wall_xy_equilibrium_module
