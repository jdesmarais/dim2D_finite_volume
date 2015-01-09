      !> @file
      !> class encapsulating subroutines to compute the initial
      !> conditions and the conditions enforced at the edge of
      !> the computational domain for bubble transported
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
     $       flow_y_side

        use parameters_kind, only :
     $       ikind,
     $       rkind        

        implicit none


        private
        public :: ic

        !set the phase at the center
        integer, parameter     :: phase_at_center = vapor
        
        !flow velocities for the different flow configurations
        real(rkind), parameter :: u0_x_flow  = flow_velocity*flow_x_side
        real(rkind), parameter :: u0_y_flow  = 0.0d0
        real(rkind), parameter :: u0_xy_flow = 0.5d0*SQRT(2.0d0)*flow_velocity*flow_x_side

        real(rkind), parameter :: v0_x_flow  = 0.0d0
        real(rkind), parameter :: v0_y_flow  = flow_velocity*flow_y_side
        real(rkind), parameter :: v0_xy_flow = 0.5d0*SQRT(2.0d0)*flow_velocity*flow_y_side


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
        !> 12_12_2014 - initial version - J.L. Desmarais
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
          real(rkind) :: velocity_x
          real(rkind) :: velocity_y
          
          t_s = t
          x_s = x
          y_s = y

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

        end function get_far_field


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
               print '(''bubble_transported/ic_class.f'')'
               print '(''apply_ic'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop ''
          end select
          
        end function get_velocity_x

        
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
               print '(''bubble_transported/ic_class.f'')'
               print '(''apply_ic'')'
               print '(''flow_direction not recognized'')'
               print '(''flow_direction: '',I2)', flow_direction
               stop ''
          end select

        end function get_velocity_y


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
