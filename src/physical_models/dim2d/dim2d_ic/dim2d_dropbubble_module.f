      !> @file
      !> module encapsulating subroutines to compute
      !> the mass density and total energy fields for
      !> initial conditions where a bubble or a droplet
      !> is present in the system
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute
      !> the mass density and total energy fields for
      !> initial conditions where a bubble or a droplet
      !> is present in the system
      !
      !> @date
      !> 27_09_2013 - initial version  - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_dropbubble_module

        use dim2d_parameters     , only : cv_r, we
        use parameters_constant  , only : liquid, vapor
        use parameters_input     , only : nx,ny,ne
        use parameters_kind      , only : rkind

        implicit none

        private
        public :: mass_density_ellipsoid,
     $            total_energy_ellipsoid


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the coordinate
        !> corresponding to an ellipsoidal 
        !> coordinate system
        !> \f[ r=\displaystyle{\sqrt{\frac{{(x-x_c)}^2}{a^2} +
        !> \frac{{(y-y_c)}^2}{b^2}}} \f]
        !> where \f$ a \f$ is the major axis and \f$ b \f$, the
        !> minor axis.
        !
        !> @date
        !> 26_09_2013 - initial version - J.L. Desmarais
        !
        !>@param x
        !> coordinate along the x-axis in a Cartesian 2D system
        !
        !>@param y
        !> coordinate along the y-axis in a Cartesian 2D system
        !
        !>@param xc
        !> coordinate along the x-axis for the center of the ellipsoid
        !
        !>@param yc
        !> coordinate along the y-axis for the center of the ellipsoid
        !
        !>@param a
        !> diameter of the ellipsoid along the x-axis (major axis)
        !
        !>@param b
        !> diameter of the ellipsoid along the y-axis (minor axis)
        !
        !>@return
        !> ellipsoidal coordinate, \f$ r \f$
        !---------------------------------------------------------------
        function get_coordinate(x,y,xc,yc,a,b) result(r)

          implicit none

          real(rkind), intent(in) :: x,y,xc,yc,a,b
          real(rkind)             :: r

          r = Sqrt((x-xc)**2/a**2+(y-yc)**2/b**2)

        end function get_coordinate


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> function giving a sign (-1 or +1) depending if there is
        !> liquid or vapor at the center of the domain
        !> \f[
        !> f(x) = 
        !> \begin{cases}
        !>   -1 & \mbox{if liquid phase in center} \\\
        !>   +1 & \mbox{if vapor phase in center}
        !> \end{cases}
        !> \f]
        !
        !> @date
        !> 27_09_2013 - initial version - J.L. Desmarais
        !
        !>@param phase_at_center
        !> choose the phase at the center of the domain: liquid or vapor
        !
        !>@return
        !> sign, \f$ f(x) \f$
        !----------------------------------------------------------------
        function choose_phase_at_center(phase_at_center) result(sign)

          implicit none

          integer, intent(in) :: phase_at_center
          integer             :: sign

          select case(phase_at_center)
            case(liquid)
               sign=-1
            case(vapor)
               sign=+1
            case default
               sign=0
               print '(''dim2d_dropbubble_module'')'
               print '(''choose_phase_at_center'')'
               stop 'phase_at_center not recognized'
          end select

        end function choose_phase_at_center


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the mass density field
        !> corresponding to an ellipsoidal droplet
        !> \f[ \rho(x,y) = \displaystyle{\frac{\rho_\textrm{liq}+\rho_\textrm{vap}}{2}
        !> + \frac{\rho_\textrm{liq}-\rho_\textrm{vap}}{2} \tanh \left(
        !> \frac{-2 \sqrt{ab} (r-r_i)}{L_i} \right)} \f]
        !> where \f$ \rho_\textrm{liq}\f$ is the mass density for the
        !> saturated liquid phase, \f$ \rho_\textrm{vap} \f$ is the mass density
        !> for the vapor phase, \f$L_i\f$ is the width of the interface,
        !> \f$r\f$ the ellipsoidal coordinate, and \f$ r_i \f$ the ellipsoidal
        !> coordinate corresponding to the location of the interface.
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param x
        !> coordinate along the x-axis in a Cartesian 2D system
        !
        !>@param y
        !> coordinate along the y-axis in a Cartesian 2D system
        !
        !>@param xc
        !> coordinate along the x-axis for the center of the ellipsoid
        !
        !>@param yc
        !> coordinate along the y-axis for the center of the ellipsoid
        !
        !>@param a
        !> diameter of the ellipsoid along the x-axis (major axis)
        !
        !>@param b
        !> diameter of the ellipsoid along the y-axis (minor axis)
        !
        !>@param li
        !> width of the interface
        !
        !>@param dliq
        !> mass density of saturated liquid water
        !
        !>@param dvap
        !> mass density of saturated vapor water
        !
        !>@param phase_at_center
        !> phase in the center        
        !
        !>@return
        !> mass density, \f$ \rho(x,y) \f$
        !---------------------------------------------------------------
        function mass_density_ellipsoid(
     $     x,y,xc,yc,a,b,li,dliq,dvap,phase_at_center)

          implicit none

          real(rkind), intent(in) :: x,y,xc,yc,a,b,li,dliq,dvap
          integer    , intent(in) :: phase_at_center
          real(rkind)             :: mass_density_ellipsoid

          real(rkind) :: scale
          real(rkind) :: r,ri
          integer     :: sign

          r     = get_coordinate(x,y,xc,yc,a,b)
          ri    = 1.0d0
          scale = Sqrt(a*b)
          sign  = choose_phase_at_center(phase_at_center)

          select case(rkind)
            case(4)
               mass_density_ellipsoid = (dliq+dvap)/2. +
     $              (dliq-dvap)/2.*tanh(sign*2*scale*(r-ri)/li)
            case(8)
               mass_density_ellipsoid = (dliq+dvap)/2.0d0 +
     $              (dliq-dvap)/2.0d0*dtanh(sign*2.0d0*scale*(r-ri)/li)
            case default
               print '(''dim2d_dropletbubble_module'')'
               print '(''mass_density_ellipsoid'')'
               stop 'rkind not supported for tanh'
           end select
        end function mass_density_ellipsoid


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the squared norm of
        !> the mass density gradient field. The mass
        !> density field corresponds to an ellipsoidal
        !> droplet
        !> \f[ {|\nabla \rho(x,y)|}^2 =
        !> {\left(\frac{\rho_\textrm{liq}-\rho_\textrm{vap}}{L_i r/\sqrt{a b}}\right)}^2
        !> {\left[ 1.0 - {\left(\textrm{tanh}\left( \frac{2 \sqrt{a b} (r - r_i)}{L_i}\right)\right)}^2\right]}^2
        !> \left[ \frac{(x-x_c)^2}{a^4} + \frac{(y-y_c)^2}{b^4} \right]
        !> \f]
        !> where \f$ \rho_\textrm{liq}\f$ is the mass density for the
        !> saturated liquid phase, \f$ \rho_\textrm{vap} \f$ is the mass density
        !> for the vapor phase, \f$L_i\f$ is the width of the interface,
        !> \f$r\f$ the ellipsoidal coordinate, \f$ r_i \f$ the ellipsoidal
        !> coordinate corresponding to the location of the interface, \f$ a \f$
        !> the major axis of the ellipse and \f$ b \f$ its minor axis.
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param x
        !> coordinate along the x-axis in a Cartesian 2D system
        !
        !>@param y
        !> coordinate along the y-axis in a Cartesian 2D system
        !
        !>@param xc
        !> coordinate along the x-axis for the center of the ellipsoid
        !
        !>@param yc
        !> coordinate along the y-axis for the center of the ellipsoid
        !
        !>@param a
        !> diameter of the ellipsoid along the x-axis (major axis)
        !
        !>@param b
        !> diameter of the ellipsoid along the y-axis (minor axis)
        !
        !>@param li
        !> width of the interface
        !
        !>@param dliq
        !> mass density of saturated liquid water
        !
        !>@param dvap
        !> mass density of saturated vapor water
        !
        !>@return
        !> squared norm of the gradient of the mass density, \f$ |\nabla \rho |^2\f$
        !---------------------------------------------------------------
        function mass_density_grad_norm_squared(
     $     x,y,xc,yc,a,b,li,dliq,dvap)

          implicit none

          real(rkind), intent(in) :: x,y,xc,yc,a,b
          real(rkind), intent(in) :: li,dliq,dvap
          real(rkind)             :: mass_density_grad_norm_squared

          real(rkind) :: r,ri
          real(rkind) :: scale
          
          scale = Sqrt(a*b)
          r     = get_coordinate(x,y,xc,yc,a,b)
          ri    = 1.0

          !< compute \f$ {| \nabla \rho |}^2 \f$
          if(r.eq.0) then
             mass_density_grad_norm_squared=0

          else
             if(rkind.eq.8) then
                mass_density_grad_norm_squared=
     $               ((dliq-dvap)/(li*r/scale))**2*
     $               (1.0d0-(dtanh(2.0d0*scale*(r-ri)/li))**2)**2*
     $               ((x-xc)**2/a**4 + (y-yc)**2/b**4)
             else
                mass_density_grad_norm_squared=
     $               ((dliq-dvap)/(li*r/scale))**2*
     $               (1.0-(tanh(2.0*scale*(r-ri)/li))**2)**2*
     $               ((x-xc)**2/a**4 + (y-yc)**2/b**4)
             end if
          end if

        end function mass_density_grad_norm_squared


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the total energy 
        !> evaluated at the cartesian point (x,y)
        !> in a uniform temperature field T0 and a
        !> mass density field corresponding to an 
        !> ellipsoidal droplet.
        !> \f[
        !> \rho E(x,y) = \rho(x,y) \, \left(\frac{8}{3}c_v
        !> T_0-3 \rho(x,y)\right) + \frac{1}{2 We} {|\nabla
        !> \rho(x,y)|}^2
        !> \f]
        !> where \f$ \rho_\textrm{liq}\f$ is the mass density for the
        !> saturated liquid phase, \f$ \rho_\textrm{vap} \f$ is the mass density
        !> for the vapor phase, \f$L_i\f$ is the width of the interface,
        !> \f$r\f$ the ellipsoidial coordinate, \f$ r_i \f$ the ellipsoidal
        !> coordinate corresponding to the location of the interface, \f$ a \f$
        !> the major axis of the ellipse and \f$ b \f$ its minor axis.
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param x
        !> coordinate along the x-axis in a Cartesian 2D system
        !
        !>@param y
        !> coordinate along the y-axis in a Cartesian 2D system
        !
        !>@param xc
        !> coordinate along the x-axis for the center of the ellipsoid
        !
        !>@param yc
        !> coordinate along the y-axis for the center of the ellipsoid
        !
        !>@param a
        !> diameter of the ellipsoid along the x-axis (major axis)
        !
        !>@param b
        !> diameter of the ellipsoid along the y-axis (minor axis)
        !
        !>@param li
        !> width of the interface
        !
        !>@param dliq
        !> mass density of saturated liquid water
        !
        !>@param dvap
        !> mass density of saturated vapor water
        !
        !>@param md
        !> mass density evaluated at the cartesian point (x,y)
        !
        !>@param T0
        !> uniform temperature in the field
        !
        !>@return
        !> total energy, \f$ \rho E(x,y)\f$
        !---------------------------------------------------------------
        function total_energy_ellipsoid(
     $     x,y,xc,yc,a,b,li,dliq,dvap,
     $     md,T0)

          implicit none

          real(rkind), intent(in) :: x,y,xc,yc,a,b,li,dliq,dvap
          real(rkind), intent(in) :: md,T0
          real(rkind)             :: total_energy_ellipsoid

          real(rkind) :: md_grad_norm_squared

          md_grad_norm_squared = mass_density_grad_norm_squared(
     $         x,y,xc,yc,a,b,li,dliq,dvap)

          if(rkind.eq.8) then
             total_energy_ellipsoid = md*(8.0d0/3.0d0*cv_r*T0-3.0d0*md)
     $            + 1.0d0/(2.0d0*we)*md_grad_norm_squared
          else
             total_energy_ellipsoid = md*(8./3.*cv_r*T0-3*md)
     $            + 1./(2*we)*md_grad_norm_squared
          end if

        end function total_energy_ellipsoid

      end module dim2d_dropbubble_module
