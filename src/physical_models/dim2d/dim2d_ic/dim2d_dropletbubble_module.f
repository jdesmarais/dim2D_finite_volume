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
      module dim2d_dropletbubble_module

        use dim2d_parameters     , only : cv_r, we
        use dim2d_state_eq_module, only : get_interface_length,
     $                                    get_mass_density_liquid,
     $                                    get_mass_density_vapor
        use field_class          , only : field
        use parameters_constant  , only : liquid, vapor
        use parameters_input     , only : nx,ny,ne
        use parameters_kind      , only : ikind, rkind

        implicit none

        private
        public :: apply_dropletbubble_ic


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
        subroutine apply_dropletbubble_ic(field_used)

          implicit none

          class(field), intent(inout) :: field_used

          integer        :: phase_at_center
          real(rkind)    :: T0,xc,yc,a,b,ri
          real(rkind)    :: dliq,dvap,li,r
          integer        :: sign
          integer(ikind) :: i,j

          
          !< choose the phase at the domain center:
          !> is it a droplet of liquid in a vapor medium ? -> vapor
          !> is it a bubble  of vapor in a liquid medium ? -> liquid
          phase_at_center = vapor


          !<set the initial temperature in the field
          T0 = 0.995

          
          !<set the center of the droplet
          xc=0.
          yc=0.


          !<get the mass densities corresponding to the
          !>liquid and vapor phases for the initial
          !>temperature field
          dliq = get_mass_density_liquid(T0)
          dvap = get_mass_density_vapor(T0)

          !<get the interface length corresponding
          !>to the initial temperature field
          li = get_interface_length(T0)

          !<set the major and minor axes of the ellipse
          a=4*li
          b=a !a/2.0d0
          ri=1.0d0


          !< depending on the phase at the center of the
          !> domain, either the sign is -1 or +1
          select case(phase_at_center)
            case(liquid)
               sign=-1
            case(vapor)
               sign=+1
          end select


          !<initialize the fields
          do j=1, ny
             do i=1, nx

                r = get_coordinate(
     $               field_used%x_map(i),
     $               field_used%y_map(j),
     $               xc,yc,a,b)

                field_used%nodes(i,j,1)=mass_density_ellipsoid(
     $               r,ri,li,dliq,dvap,sign)
                
                field_used%nodes(i,j,2)=0.0

                field_used%nodes(i,j,3)=0.0

                field_used%nodes(i,j,4)=total_energy(
     $               T0,
     $               field_used%x_map(i),field_used%y_map(j),
     $               xc,yc,r,ri,a,b,dliq,dvap,li,sign)

             end do
          end do

        end subroutine apply_dropletbubble_ic


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the coordinate
        !> corresponding to an ellipsoidial 
        !> coordinate system
        !> \f$ r=\sqrt{\frac{{(x-x_c)}^2}{a^2} +
        !> \frac{{(y-y_c)}^2}{b^2}} \f$
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
        !> diameter of the ellipsoid along the x-axis
        !
        !>@param b
        !> diameter of the ellipsoid along the y-axis
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
        function mass_density_ellipsoid(r,ri,li,dliq,dvap,sign)

          implicit none

          real(rkind), intent(in) :: r,ri,li,dliq,dvap
          integer    , intent(in) :: sign
          real(rkind)             :: mass_density_ellipsoid


          select case(rkind)
            case(4)
               mass_density_ellipsoid = (dliq+dvap)/2. +
     $              (dliq-dvap)/2.*tanh(sign*2*(r-ri)/li)
            case(8)
               mass_density_ellipsoid = (dliq+dvap)/2.0d0 +
     $              (dliq-dvap)/2.0d0*dtanh(sign*2.0d0*(r-ri)/li)
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
        !>@param ri
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
        !>@param li
        !> length of the interface
        !---------------------------------------------------------------
        function mass_density_grad_norm_squared(
     $     x,y,xc,yc,
     $     r,ri,a,b,
     $     dliq,dvap,li)

          implicit none

          real(rkind), intent(in) :: x,y,xc,yc
          real(rkind), intent(in) :: r,ri,a,b
          real(rkind), intent(in) :: dliq,dvap,li
          real(rkind)             :: mass_density_grad_norm_squared

          !< compute \f$ {| \nabla \rho |}^2 \f$
          if(r.eq.0) then
             mass_density_grad_norm_squared=0

          else
             if(rkind.eq.8) then
                mass_density_grad_norm_squared=
     $               ((dliq-dvap)/(li*r))**2*
     $               (1.0d0-(dtanh(2.0d0*(r-ri)/li))**2)**2*
     $               (x**2/a**4 + y**2/b**4)
             else
                mass_density_grad_norm_squared=
     $               ((dliq-dvap)/(li*r))**2*
     $               (1.0-(dtanh(2.0*(r-ri)/li))**2)**2*
     $               (x**2/a**4 + y**2/b**4)    
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
        !>@param ri
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
        function total_energy(T0,x,y,xc,yc,r,ri,a,b,dliq,dvap,li,sign)

          implicit none

          real(rkind), intent(in) :: T0   !reduced temperature  at (x,y)
          real(rkind), intent(in) :: x,y,xc,yc,r,ri,a,b,dliq,dvap,li
          integer    , intent(in) :: sign
          real(rkind)             :: total_energy

          real(rkind) :: md, md_grad_norm_squared


          md = mass_density_ellipsoid(r,ri,li,dliq,dvap,sign)

          md_grad_norm_squared = mass_density_grad_norm_squared(
     $         x,y,xc,yc,
     $         r,ri,a,b,
     $         dliq,dvap,li)

          if(rkind.eq.8) then
             total_energy = md*(8.0d0/3.0d0*cv_r*T0-3.0d0*md) + 
     $            1.0d0/(2.0d0*we)*md_grad_norm_squared
          else
             total_energy = md*(8./3.*cv_r*T0-3*md) + 
     $            1./(2*we)*md_grad_norm_squared
          end if

        end function total_energy

      end module dim2d_dropletbubble_module
