      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute
      !> the initial conditions of drop collision
      !
      !> @date
      !> 27_09_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_drop_collision_module

        use dim2d_dropbubble_module, only : mass_density_ellipsoid,
     $                                      total_energy_ellipsoid
        use dim2d_vortex_module,     only : get_vortex_velocity
        use dim2d_state_eq_module  , only : get_mass_density_liquid,
     $                                      get_mass_density_vapor,
     $                                      get_interface_length
        use field_class            , only : field
        use parameters_constant    , only : liquid, vapor
        use parameters_input       , only : nx,ny
        use parameters_kind        , only : ikind, rkind

        implicit none

        private
        public :: apply_drop_collision_ic


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> for drop collision
        !
        !> @date
        !> 29_10_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !---------------------------------------------------------------
        subroutine apply_drop_collision_ic(field_used)

          implicit none

          class(field), intent(inout) :: field_used


          !< local variables for the droplet/bubble
          integer        :: phase_at_center
          real(rkind)    :: T0,xc,yc,a,b
          real(rkind)    :: dliq,dvap,li

          !< local variable for the vortex
          !real(rkind)    :: vortex_xc,vortex_yc
          !real(rkind)    :: vortex_rc,vortex_omega

          !< local variables for the initialization
          integer(ikind)            :: i,j
          real(rkind)               :: x,y
          real(rkind), dimension(2) :: velocity
          real(rkind)               :: period,amp

          
          !< choose the phase at the domain center:
          !> is it a droplet of liquid in a vapor medium ? -> vapor
          !> is it a bubble  of vapor in a liquid medium ? -> liquid
          phase_at_center = liquid

          !<set the initial temperature in the field
          T0 = 0.995
          
          !<set the center of the droplet
          xc=-0.5
          yc=0.

          !<set the vortex properties
          period= 6.0
          amp   = 0.05

          !<get the mass densities corresponding to the
          !>liquid and vapor phases for the initial
          !>temperature field
          dliq = get_mass_density_liquid(T0)
          dvap = get_mass_density_vapor(T0)

          !<get the interface length corresponding
          !>to the initial temperature field
          li = get_interface_length(T0)

          !<set the major and minor axes of the bubble ellipse
          a=2.8d0*li
          b=a


          !<initialize the fields
          do j=1, ny
             !DEC$ IVDEP
             do i=1, nx

                x = field_used%x_map(i)
                y = field_used%y_map(j)

                field_used%nodes(i,j,1)=mass_density_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,phase_at_center)                
                
                velocity = get_divergence_free_sinusoidal_velocity(
     $               x,y,amp,period)

                field_used%nodes(i,j,2)=
     $               field_used%nodes(i,j,1)*velocity(1)
                field_used%nodes(i,j,3)=
     $               field_used%nodes(i,j,1)*velocity(2)
                
                field_used%nodes(i,j,4)=total_energy_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,
     $               field_used%nodes(i,j,1),T0)+
     $               0.5d0*(field_used%nodes(i,j,2)**2+
     $               field_used%nodes(i,j,3)**2)/field_used%nodes(i,j,1)

             end do
          end do

        end subroutine apply_drop_collision_ic


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial velocity
        !> field for drop collision
        !> \f{eqnarray*}{
        !> u_x     &=& - \sin \left( \frac{ 2\pi(x+y)}{T} \right)
        !>             - \sin \left( \frac{ 2\pi(x-y)}{T} \right) \\\
        !> u_y     &=&   \sin \left( \frac{ 2\pi(x+y)}{T} \right)
        !>             - \sin \left( \frac{ 2\pi(x-y)}{T} \right) \\\
        !> \f}
        !
        !> @date
        !> 30_10_2013 - initial version - J.L. Desmarais
        !
        !>@param amp:A
        !> amplitude of the sinusoidal velocity field
        !>
        !>@param period:T
        !> period of the sinusoidal velocity field
        !---------------------------------------------------------------
        function get_divergence_free_sinusoidal_velocity(
     $     x,y,amp,period)
     $     result(velocity)
        
          implicit none

          real(rkind), intent(in)   :: x
          real(rkind), intent(in)   :: y
          real(rkind), intent(in)   :: amp
          real(rkind), intent(in)   :: period
          real(rkind), dimension(2) :: velocity

          real(rkind) :: pi

          if(rkind.eq.8) then
             pi = acos(-1.0d0)

             velocity(1)=-sin((x+y)*2.0d0*pi/period)
     $                   -sin((x-y)*2.0d0*pi/period)
             velocity(2)= sin((x+y)*2.0d0*pi/period)
     $                   -sin((x-y)*2.0d0*pi/period)
             
          else             
             pi = acos(-1.0)

             velocity(1)=-sin((x+y)*2.0*pi/period)
     $                   -sin((x-y)*2.0*pi/period)
             velocity(2)= sin((x+y)*2.0*pi/period)
     $                   -sin((x-y)*2.0*pi/period)

          end if

          velocity(1)=amp*velocity(1)
          velocity(2)=amp*velocity(2)
 
        end function get_divergence_free_sinusoidal_velocity

      end module dim2d_drop_collision_module
