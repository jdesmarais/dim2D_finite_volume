      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a drop retraction
      !> (droplet with an ellipsoidal shape converging to a disc in 2-D)
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for a drop retraction
      !> (droplet with an ellipsoidal shape converging to a disc in 2-D)
      !
      !> @date
      !> 11_12_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_class
       
        use dim2d_dropbubble_module, only :
     $       mass_density_ellipsoid,
     $       total_energy_ellipsoid

        use dim2d_parameters, only :
     $       cv_r

        use dim2d_state_eq_module  , only :
     $       get_mass_density_liquid,
     $       get_mass_density_vapor,
     $       get_interface_length

        use ic_abstract_class, only :
     $       ic_abstract

        use parameters_constant, only :
     $       liquid, vapor

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       T0

        use parameters_kind, only :
     $       ikind,
     $       rkind        

        implicit none


        integer, parameter :: phase_at_center = liquid !<@brief phase at the center (vapor: bubble in saturated liquid, liquid: droplet in saturated vapor)


        private
        public :: ic


        !> @class ic
        !> class encapsulating operators to set the initial
        !> conditions and the conditions enforced at the edge of the
        !> computational domain for a domain extension test
        !---------------------------------------------------------------
        type, extends(ic_abstract) :: ic

          contains

          procedure, nopass :: apply_ic          !<@brief set the initial conditions                                                 
          procedure, nopass :: get_mach_ux_infty !<@brief get the Mach number along the x-direction in the far-field                 
          procedure, nopass :: get_mach_uy_infty !<@brief get the Mach number along the y-direction in the far-field                 
          procedure, nopass :: get_u_in          !<@brief get the x-component of the velocity at the edge of the computational domain
          procedure, nopass :: get_v_in          !<@brief get the y-component of the velocity at the edge of the computational domain
          procedure, nopass :: get_T_in          !<@brief get the temperature at the edge of the computational domain                
          procedure, nopass :: get_P_out         !<@brief get the pressure at the edge of the computational domain                   
          procedure, nopass :: get_far_field     !<@brief get the governing variables imposed at the edge of the computational domain

        end type ic


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> apply the initial conditions
        !> in the computational domain for a droplet/bubble
        !> with an ellipsoidal shape
        !> \f[
        !> \begin{pmatrix} 
        !> \rho \\\ \rho u \\\ \rho v \\\ \rho E
        !> \end{pmatrix}(x,y) =
        !> \begin{pmatrix} 
        !> \rho_\textrm{bubble}(r(x,y,a,b),r(a,0,a,b)) \\\
        !> 0 \\\
        !> 0 \\\
        !> \rho(x,y) \left[ \frac{8}{3} c_v T_0 - 3 \rho(x,y) \right]
        !> + \frac{1}{2 \textrm{We}} | \nabla \rho(x,y) |^2
        !> \end{pmatrix}
        !> \f]
        !> where
        !> \f[ \rho_\textrm{bubble}(r,r_c) = \frac{\rho_\textrm{liq} + \rho_\textrm{vap}}{2}
        !> + \frac{\rho_\textrm{liq} - \rho_\textrm{vap}}{2} \tanh \left( \frac{2 (r-r_c)}{L_i} \right)\f]
        !> \f[ r(x,y,a,b) = a b \sqrt{ \frac{x^2}{a^2} + \frac{y^2}{b^2} } \f]
        !> and \f$ L_i \f$ is the width of the interface,
        !> and \f$ a \f$ and \f$ b \f$ are the major and minor radii
        !> of the ellipsoid
        !> \f[ a = 6 L_i \f]
        !> \f[ b = 3 L_i \f]
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data    
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
          a=6.0d0*li
          b=a/2.0d0


          !initialize the fields
          do j=1, ny
             !DEC$ IVDEP
             do i=1, nx

                x = x_map(i)
                y = y_map(j)

                nodes(i,j,1)=mass_density_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,phase_at_center)

                nodes(i,j,2)=0.0d0

                nodes(i,j,3)=0.0d0
                
                nodes(i,j,4)=total_energy_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,
     $               nodes(i,j,1),T0)

             end do
          end do

        end subroutine apply_ic


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the Mach number imposed in the far-field
        !> for the velocity in the x-direction
        !> \f[ \textrm{Ma}_x = 0 \f]
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

          logical :: side_s

          side_s = side

          if(rkind.eq.8) then
             var = 0.0d0
          else
             var = 0.0
          end if

        end function get_mach_ux_infty


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the Mach number imposed in the far-field
        !> for the velocity in the y-direction
        !> \f[ \textrm{Ma}_y = 0 \f]
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

          logical :: side_s

          side_s = side

          if(rkind.eq.8) then
             var = 0.0d0
          else
             var = 0.0
          end if

        end function get_mach_uy_infty


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to obtain the value of the velocity
        !> in the x-direction imposed in the far-field
        !> \f[ u_\infty(t,x,y) = 0 \f]
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

          t_s = t
          x_s = x
          y_s = y

          if(rkind.eq.8) then
             var = 0.0d0
          else
             var = 0.0
          end if

        end function get_u_in


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to obtain the value of the velocity
        !> in the y-direction imposed in the far-field
        !> \f[ v_\infty(t,x,y) = 0 \f]
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

          t_s = t
          x_s = x
          y_s = y

          if(rkind.eq.8) then
             var = 0.0d0
          else
             var = 0.0
          end if

        end function get_v_in

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine to obtain the value of the
        !> temperature imposed in the far-field
        !> \f[ T_\infty(t,x,y) = T_0\f]
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
        !> subroutine to obtain the value of the
        !> pressure imposed in the far-field
        !> \f[ P_\infty(t,x,y) = \frac{8 \rho_\textrm{liq} T_0}{3-\rho_\textrm{liq}}
        !> - 3 \rho_\textrm{liq}^2\f]
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

          select case(phase_at_center)

            case(liquid)
               mass = get_mass_density_vapor(T0)

            case(vapor)
               mass = get_mass_density_liquid(T0)

            case default
               print '(''dim2d_ic/drop_retraction/ic_class'')'
               print '(''get_P_out'')'
               print '(''phase_at_center: '',I2)', phase_at_center
               print '(''case not recognized'')'
               stop ''

          end select

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
        !> subroutine to obtain the value of the variables
        !> imposed at the edge of the computational domain
        !> depending on time and coordinates as well as the
        !> state of the object
        !> \f[
        !> \begin{pmatrix}
        !> \rho_\infty \\\ {\rho u}_\infty \\\ {\rho v}_\infty \\\ {\rho E}_\infty
        !> \end{pmatrix} = 
        !> \begin{pmatrix}
        !> \rho_\textrm{liq}(T_0) \\\
        !> 0 \\\
        !> 0 \\\
        !> \rho_\textrm{liq}(T_0) (\frac{8}{3} c_v T_0 - 3 \rho_\textrm{liq}(T_0))
        !> \end{pmatrix}
        !> \f]
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
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

            class(ic)                 , intent(in) :: this
            real(rkind)               , intent(in) :: t
            real(rkind)               , intent(in) :: x
            real(rkind)               , intent(in) :: y
            real(rkind), dimension(ne)             :: var


            real(rkind) :: mass

            
            select case(phase_at_center)

              case(liquid)
                 mass = get_mass_density_vapor(T0)
                 
              case(vapor)
                 mass = get_mass_density_liquid(T0)
                 
              case default
                 print '(''dim2d_ic/drop_retraction/ic_class'')'
                 print '(''get_far_field'')'
                 print '(''phase_at_center: '',I2)', phase_at_center
                 print '(''case not recognized'')'
                 stop ''
                 
            end select


            if(rkind.eq.8) then

               var(1) = mass
               var(2) = var(1)*get_u_in(t,x,y)
               var(3) = var(1)*get_v_in(t,x,y)
               var(4) = 0.5d0*(var(2)**2+var(3)**2)/var(1) +
     $                  var(1)*(8.0d0/3.0d0*cv_r*T0-3.0d0*mass)

            else

               var(1) = mass
               var(2) = var(1)*get_u_in(t,x,y)
               var(3) = var(1)*get_v_in(t,x,y)
               var(4) = 0.5*(var(2)**2+var(3)**2)/var(1) +
     $                  var(1)*(8.0/3.0*cv_r*T0-3.0*mass)

            end if

          end function get_far_field

      end module ic_class
