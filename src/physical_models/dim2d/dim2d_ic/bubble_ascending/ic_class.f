      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for bubble ascending
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for bubble ascending
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

        use dim2d_state_eq_module, only :
     $       get_mass_density_liquid,
     $       get_mass_density_vapor,
     $       get_interface_length

        use dim2d_vortex_module, only :
     $       get_vortex_velocity

        use ic_abstract_class, only :
     $       ic_abstract

        use parameters_constant, only :
     $       liquid,
     $       vapor

        use parameters_input, only :
     $       nx,
     $       ny,
     $       ne

        use parameters_kind, only :
     $       ikind,
     $       rkind        

        implicit none


        private
        public :: ic

        !set the phase at the center
        integer, parameter     :: phase_at_center = vapor
        
        !set the initial temperature in the field
        real(rkind), parameter :: T0 = 0.995d0


        !> @class ic
        !> class encapsulating operators to set the initial
        !> conditions and the conditions enforced at the edge of the
        !> computational domain for phase separation
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
        !> 11_12_2014 - initial version - J.L. Desmarais
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

          
          !< local variables for the droplet/bubble
          real(rkind)    :: xc,yc,a,b
          real(rkind)    :: dliq,dvap,li

          !< local variables for the vortices
          real(rkind)               :: x_v1, y_v1, x_v2, y_v2
          real(rkind)               :: r_v1, r_v2
          real(rkind)               :: omega_v1, omega_v2
          real(rkind), dimension(2) :: velocity_vortex1
          real(rkind), dimension(2) :: velocity_vortex2
          real(rkind)               :: velocity_x, velocity_y          

          !< local variables for the initialization
          integer(ikind) :: i,j
          real(rkind)    :: x,y

          
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

          !<set the major and minor axes of the bubble ellipse
          a=3.0d0*li
          b=a !a/2.0d0


          !< set the position of the vortex centers
          x_v1 = -1.5
          y_v1 = 0.

          x_v2 = 1.5
          y_v2 = 0.

          !< set the radius of the vortex
          r_v1 = 0.5
          r_v2 = 0.5

          !< set the strength of the vortex
          omega_v1 =  0.0d0
          omega_v2 =  0.0d0


          !<initialize the mass and total energy fields
          do j=1, size(y_map,1)
             do i=1, size(x_map,1)

                x = x_map(i)
                y = y_map(j)

                nodes(i,j,1)=mass_density_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,phase_at_center)
                
                nodes(i,j,4)=total_energy_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,
     $               nodes(i,j,1),T0)

             end do
          end do


          !< initialize the momentum and update the total energy fields
          do j=1, size(y_map,1)
             do i=1, size(x_map,1)
                
                x = x_map(i)
                y = y_map(j)

                velocity_vortex1 = get_vortex_velocity(
     $               x, y, x_v1, y_v1, r_v1, omega_v1)
                velocity_vortex2 = get_vortex_velocity(
     $               x, y, x_v2, y_v2, r_v2, omega_v2)

                velocity_x = velocity_vortex1(1)+velocity_vortex2(1)
                velocity_y = velocity_vortex1(2)+velocity_vortex2(2)

                nodes(i,j,2)=
     $               nodes(i,j,1)*velocity_x
                nodes(i,j,3)=
     $               nodes(i,j,1)*velocity_y

                if(rkind.eq.8) then
                   nodes(i,j,4)=nodes(i,j,4)+
     $                  1.0d0/2.0d0*nodes(i,j,1)*(
     $                  velocity_x**2+velocity_y**2)
                else
                   nodes(i,j,4)=nodes(i,j,4)+
     $                  1./2.*nodes(i,j,1)*(
     $                  velocity_x**2+velocity_y**2)
                end if

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

          side_s = side

          var = 0.0d0

        end function get_mach_ux_infty


        !get the variable enforced at the edge of the
        !computational domain
        function get_mach_uy_infty(side) result(var)

          implicit none

          logical, intent(in) :: side
          real(rkind)         :: var

          logical :: side_s

          side_s = side

          var = 0.0d0

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

          var = 0.0d0

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

          var = 0.0d0

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

          select case(phase_at_center)

            case(liquid)
               mass = get_mass_density_vapor(T0)

            case(vapor)
               mass = get_mass_density_liquid(T0)

            case default
               print '(''dim2d/dim2d_ic/bubble_ascending'')'
               print '(''phase at center: '',I2)', phase_at_center
               print '(''phase_at_center not recognized'')'
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
               print '(''dim2d/dim2d_ic/bubble_ascending'')'
               print '(''phase at center: '',I2)', phase_at_center
               print '(''phase_at_center not recognized'')'
               stop ''

          end select

          if(rkind.eq.8) then

             var(1) = mass
             var(2) = 0.0d0
             var(3) = 0.0d0
             var(4) = mass*(8.0d0/3.0d0*cv_r*T0-3.0d0*mass)

          else

             var(1) = mass
             var(2) = 0.0d0
             var(3) = 0.0d0
             var(4) = mass*(8.0/3.0*cv_r*T0-3.0*mass)

          end if

        end function get_far_field

      end module ic_class
