      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for saturated liquid
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for saturated liquid
      !
      !> @date
      !> 11_12_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module ic_class
       
        use dim2d_parameters, only :
     $       cv_r

        use dim2d_state_eq_module, only :
     $       get_mass_density_liquid

        use ic_abstract_class, only :
     $       ic_abstract

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne,
     $       T0

        use parameters_kind, only :
     $       ikind,
     $       rkind

        implicit none


        private
        public :: ic

        
        !> @class ic
        !> class extending ic_abstract to encapsulate operators
        !> setting the initial conditions and the conditions
        !> enforced at the edge of the computational domain
        !> for steady state
        !---------------------------------------------------------------
        type, extends(ic_abstract) :: ic

          contains

          procedure, nopass :: apply_ic          !<@brief set the initial conditions                                                 
          procedure, nopass :: get_mach_ux_infty !<@brief get the Mach number along the x-direction in the far field                 
          procedure, nopass :: get_mach_uy_infty !<@brief get the Mach number along the y-direction in the far field                 
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
        !> subroutine to apply the initial conditions
        !> to the computational domain for steady state
        !> \f[
        !> \begin{pmatrix} 
        !> \rho \\\ \rho u \\\ \rho v \\\ \rho E
        !> \end{pmatrix}(x,y) =
        !> \begin{pmatrix} 
        !> \rho_\textrm{liq} \\\ 0 \\\ 0 \\\
        !> \rho_\textrm{liq} (\frac{8}{3} c_v T_0 - 3 \rho_\textrm{liq})
        !> \end{pmatrix}
        !> \f]
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

          
          integer(ikind) :: i,j
          real(rkind)    :: d_liq, E_liq
          real(rkind)    :: x_s, y_s

          x_s = x_map(1)
          y_s = y_map(1)

          ! compute the corresponding reduced temperature
          ! reduced liquid mas density and the total energy
          ! at the temperature asked by the user
          d_liq = get_mass_density_liquid(T0)
          E_liq = get_total_energy(d_liq,T0)

          if(rkind.eq.8) then
             do j=1, ny
                do i=1, nx
                   
                   nodes(i,j,1) =  d_liq
                   nodes(i,j,2) =  0.0d0
                   nodes(i,j,3) =  0.0d0
                   nodes(i,j,4) =  E_liq
                   
                end do
             end do
          else
             do j=1, ny
                do i=1, nx
                   
                   nodes(i,j,1) =  d_liq
                   nodes(i,j,2) =  0.0
                   nodes(i,j,3) =  0.0
                   nodes(i,j,4) =  E_liq
                   
                end do
             end do
          end if

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

          logical     :: side_s

          side_s = side

          var = 0.0d0

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

          var = 0.0d0

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

          var = 0.0d0

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

          var = 0.0d0

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

          mass = get_mass_density_liquid(T0)

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

          class(ic)     , intent(in) :: this
          real(rkind)   , intent(in) :: t
          real(rkind)   , intent(in) :: x
          real(rkind)   , intent(in) :: y
          real(rkind), dimension(ne) :: var


          real(rkind) :: t_s,x_s,y_s
          real(rkind) :: d_liq
          
          t_s = t
          x_s = x
          y_s = y

          d_liq = get_mass_density_liquid(T0)

          if(rkind.eq.8) then

             var(1) = d_liq
             var(2) = 0.0d0
             var(3) = 0.0d0
             var(4) = d_liq*(8.0d0/3.0d0*cv_r*T0-3.0d0*d_liq)

          else

             var(1) = d_liq
             var(2) = 0.0d0
             var(3) = 0.0d0
             var(4) = d_liq*(8.0/3.0*cv_r*T0-3.0*d_liq)

          end if

        end function get_far_field


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the total energy
        !> for a homogeneous liquid given its mass
        !> density and its temperature
        !> \f[ \rho E(\rho,T) = \rho \left(
        !> \frac{8}{3} c_v T - 3 \rho \right)
        !> - 3 \rho^2
        !> \f]
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !>@param mass_density
        !> mass density, \f$ \rho \f$
        !
        !>@param temperature
        !> temperature, \f$ T \f$
        !---------------------------------------------------------------
        function get_total_energy(mass_density, temperature)
     $     result(total_energy)

          implicit none

          real(rkind), intent(in) :: mass_density
          real(rkind), intent(in) :: temperature
          real(rkind)             :: total_energy


          if(rkind.eq.8) then
             total_energy = mass_density*(
     $            8.0d0/3.0d0*cv_r*temperature - 3.0d0*mass_density)
          else
             total_energy = mass_density*(
     $            8.0/3.0*cv_r*temperature - 3.0*mass_density)
          end if

        end function get_total_energy

      end module ic_class

