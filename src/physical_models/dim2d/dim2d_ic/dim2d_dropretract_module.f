      !> @file
      !> module encapsulating subroutines to compute
      !> the drop retraction initial conditions for
      !> the Diffuse Interface Model
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute
      !> the drop retraction initial conditions for
      !> the Diffuse Interface Model
      !
      !> @date
      !> 14_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_dropretract_module

        use dim2d_parameters     , only : cv_r, we
        use dim2d_state_eq_module, only : get_interface_length,
     $                                    get_mass_density_liquid,
     $                                    get_mass_density_vapor
        use field_class          , only : field
        use parameters_kind      , only : ikind, rkind

        implicit none

        private
        public :: apply_drop_retraction_ic


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> for drop retraction
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !---------------------------------------------------------------
        subroutine apply_drop_retraction_ic(field_used)

          implicit none

          class(field), intent(inout) :: field_used

          real(rkind)    :: T0
          real(rkind)    :: xc,yc,a,b,radius,rc
          real(rkind)    :: dliq,dvap,l
          integer(ikind) :: i,j


          !<set the initial temperature in the field
          T0 = 0.995

          
          !<set the center of the droplet
          xc=0.
          yc=0.


          !<get the interface length corresponding
          !>to the initial temperature field
          l = get_interface_length(T0)


          !<set the major and minor radii of the ellipse
          a=8*l**2
          b=a/2
          rc=Sqrt(a/2)


          !<get the mass densities corresponding to the
          !>liquid and vapor phases for the initial
          !>temperature field
          dliq = get_mass_density_liquid(T0)
          dvap = get_mass_density_vapor(T0)


          !<initialize the fields
          do j=1, size(field_used%y_map,1)
             do i=1, size(field_used%x_map,1)

                radius = r(
     $               field_used%x_map(i),
     $               field_used%y_map(j),
     $               xc,yc,a,b)

                field_used%nodes(i,j,1)=mass_density_ellipsoid(
     $               radius,
     $               rc,
     $               dliq, dvap,l)
                
                field_used%nodes(i,j,2)=0.0

                field_used%nodes(i,j,3)=0.0

                field_used%nodes(i,j,4)=total_energy(
     $               field_used%nodes(i,j,1),
     $               T0,
     $               field_used%x_map(i),
     $               field_used%y_map(j),
     $               xc,yc,radius,rc,a,b,dliq,dvap,l)

             end do
          end do

        end subroutine apply_drop_retraction_ic


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the coordinate
        !> corresponding to an ellipsoidial 
        !> coordinate system
        !> \f$ r=\sqrt{\frac{{(x-x_c)}^2}{a} +
        !> \frac{{(y-y_c)}^2}{b}} \f$
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
        !> diameter of the ellipsoid along the x-axis
        !
        !>@param b
        !> diameter of the ellipsoid along the y-axis
        !---------------------------------------------------------------
        function r(x,y,xc,yc,a,b)

          implicit none

          real(rkind), intent(in) :: x,y,xc,yc,a,b
          real(rkind)             :: r

          r = Sqrt((x-xc)**2/a+(y-yc)**2/b)

        end function r


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the derivative of the
        !> coordinate in an ellipsoidial 
        !> coordinate system
        !> \f$ \frac{dr}{dx}=\frac{x-x_c}{r a}\f$
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param x
        !> coordinate along the x-axis in a Cartesian 2D system
        !
        !>@param xc
        !> coordinate along the x-axis for the center of the ellipsoid
        !
        !>@param r
        !> ellipsoidal coordinate
        !
        !>@param a
        !> diameter of the ellipsoid along the x-axis
        !---------------------------------------------------------------
        function drdx(x,xc,r,a)
        
          implicit none

          real(rkind), intent(in) :: x,xc,r,a
          real(rkind)             :: drdx

          drdx = (x-xc)/(r*a)

        end function drdx

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the derivative of the
        !> coordinate in an ellipsoidial 
        !> coordinate system
        !> \f$ \frac{dr}{dy}=\frac{y-y_c}{r b}\f$
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param y
        !> coordinate along the y-axis in a Cartesian 2D system
        !
        !>@param yc
        !> coordinate along the y-axis for the center of the ellipsoid
        !
        !>@param r
        !> ellipsoidal coordinate
        !
        !>@param b
        !> diameter of the ellipsoid along the y-axis
        !---------------------------------------------------------------
        function drdy(y,yc,r,b)
        
          implicit none

          real(rkind), intent(in) :: y,yc,r,b
          real(rkind)             :: drdy

          drdy = (y-yc)/(r*b)

        end function drdy


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the mass density field
        !> corresponding to an ellipsoidal droplet
        !> \f$ \rho(x,y) = \frac{\rho_l+\rho_v}{2}
        !> + \frac{\rho_l-\rho_v}{2} \tanh \left(
        !> \frac{-2 (r-r_c)}{l} \right) \f$
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param r
        !> ellipsoidal coordinate
        !
        !>@param rc
        !> ellipsoidal coordinate of the interface center
        !
        !>@param d_liq
        !> mass density of saturated liquid water
        !
        !>@param d_vap
        !> mass density of saturated vapor water
        !
        !>@param l
        !> length of the interface
        !---------------------------------------------------------------
        function mass_density_ellipsoid(r,rc,dliq,dvap,l)

          implicit none

          real(rkind), intent(in) :: r,rc,dliq,dvap,l
          real(rkind)             :: mass_density_ellipsoid


          select case(rkind)
            case(4)
               mass_density_ellipsoid = (dliq+dvap)/2. +
     $              (dliq-dvap)/2.*tanh(-2*(r-rc)/l)
            case(8)
               mass_density_ellipsoid = (dliq+dvap)/2. +
     $              (dliq-dvap)/2.*dtanh(-2*(r-rc)/l)
            case default
               print '(''ic_drop_retractation_class'')'
               print '(''mass_density'')'
               stop 'rkind not supported for tanh'
           end select
        end function mass_density_ellipsoid


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the mass density
        !> gradient field corresponding to an ellipsoidal
        !> droplet \f$ \frac{d \rho(r)}{d z} =
        !> \frac{\rho_l-\rho_v}{l} \left( 1 - \tanh^2
        !> \left( \frac{-2 (r-r_c)}{l} \right) \right)
        !> \frac{d r}{d z}\f$ where \f$z\f$ can be either
        !> \f$x\f$ or \f$y\f$
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param r
        !> ellipsoidal coordinate
        !
        !>@param rc
        !> ellipsoidal coordinate of the interface center
        !
        !>@param r_dev
        !> derivative of the ellipsoidal coordinate along the
        !> x or the y-axis
        !
        !>@param d_liq
        !> mass density of saturated liquid water
        !
        !>@param d_vap
        !> mass density of saturated vapor water
        !
        !>@param l
        !> length of the interface
        !---------------------------------------------------------------
        function mass_density_ellipsoid_dev(r,rc,r_dev,dliq,dvap,l)

          implicit none
          
          real(rkind), intent(in) :: r,rc,r_dev,dliq,dvap,l
          real(rkind)             :: mass_density_ellipsoid_dev


          mass_density_ellipsoid_dev = (dliq-dvap)/l*
     $         (1-(tanh(-2*(r-rc)/l))**2)*r_dev

        end function mass_density_ellipsoid_dev


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the squared norm of
        !> the mass density gradient field. The mass
        !> density field corresponds to an ellipsoidal
        !> droplet \f$ {|\nabla \rho(x,y)|}^2 \f$
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param x
        !> coordinate along the x-axis in a Cartesian 2D system
        !
        !>@param xc
        !> coordinate along the x-axis for the center of the ellipsoid
        !
        !>@param y
        !> coordinate along the y-axis in a Cartesian 2D system
        !
        !>@param yc
        !> coordinate along the y-axis for the center of the ellipsoid
        !
        !>@param r
        !> ellipsoidal coordinate
        !
        !>@param rc
        !> ellipsoidal coordinate of the interface center
        !
        !>@param a
        !> diameter of the ellipsoid along the x-axis
        !
        !>@param b
        !> diameter of the ellipsoid along the y-axis
        !
        !>@param d_liq
        !> mass density of saturated liquid water
        !
        !>@param d_vap
        !> mass density of saturated vapor water
        !
        !>@param l
        !> length of the interface
        !---------------------------------------------------------------
        function mass_density_grad_norm(
     $     x,y,xc,yc,
     $     r,rc,a,b,
     $     dliq,dvap,l)

          implicit none

          real(rkind), intent(in) :: x,y,xc,yc
          real(rkind), intent(in) :: r,rc,a,b
          real(rkind), intent(in) :: dliq,dvap,l
          real(rkind)             :: mass_density_grad_norm

          real(rkind) :: rx_dev,ry_dev
          real(rkind) :: d_mde_dx, d_mde_dy


          rx_dev = drdx(x,xc,r,a)
          ry_dev = drdy(y,yc,r,b)

          d_mde_dx = mass_density_ellipsoid_dev(r,rc,rx_dev,dliq,dvap,l)
          d_mde_dy = mass_density_ellipsoid_dev(r,rc,ry_dev,dliq,dvap,l)

          
          if(r.eq.0) then
             mass_density_grad_norm=0

          else
             mass_density_grad_norm=(d_mde_dx**2+d_mde_dy**2)

          end if

        end function mass_density_grad_norm


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the total energy 
        !> evaluated at the cartesian point (x,y)
        !> in a uniform temperature field T0 and a
        !> mass density field corrsponding to an 
        !> ellisodial droplet.
        !> \f$ \rho E(x,y) = \rho(x,y)*(\frac{8}{3}*c_v
        !> *T_0-3*\rho(x,y)) + \frac{1}{2*We}*{|\nabla
        !> \rho(x,y)|}^2\f$
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param md
        !> mass density evaluated at the cartesian point (x,y)
        !
        !>@param T0
        !> uniform temperature in the field
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
        !>@param r
        !> ellipsoidal coordinate
        !
        !>@param rc
        !> ellipsoidal coordinate of the interface center
        !
        !>@param a
        !> diameter of the ellipsoid along the x-axis
        !
        !>@param b
        !> diameter of the ellipsoid along the y-axis
        !
        !>@param d_liq
        !> mass density of saturated liquid water
        !
        !>@param d_vap
        !> mass density of saturated vapor water
        !
        !>@param l
        !> length of the interface
        !---------------------------------------------------------------
        function total_energy(md,T0,x,y,xc,yc,r,rc,a,b,dliq,dvap,l)

          implicit none

          real(rkind), intent(in) :: md   !reduced mass density at (x,y)
          real(rkind), intent(in) :: T0   !reduced temperature  at (x,y)
          real(rkind), intent(in) :: x,y,xc,yc,r,rc,a,b,dliq,dvap,l
          real(rkind)             :: total_energy

          real(rkind) :: md_grad_norm


          md_grad_norm = mass_density_grad_norm(
     $         x,y,xc,yc,
     $         r,rc,a,b,
     $         dliq,dvap,l)

          total_energy = md*(8./3.*cv_r*T0-3*md) + 
     $         1./(2*we)*md_grad_norm

        end function total_energy

      end module dim2d_dropretract_module
