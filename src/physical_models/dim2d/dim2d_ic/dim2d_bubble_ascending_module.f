      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute
      !> the initial conditions of a bubble ascending
      !
      !> @date
      !> 27_09_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_bubble_ascending_module

        use dim2d_dropbubble_module, only : mass_density_ellipsoid,
     $                                      total_energy_ellipsoid
        use dim2d_vortex_module    , only : get_vortex_velocity
        use dim2d_state_eq_module  , only : get_mass_density_liquid,
     $                                      get_mass_density_vapor,
     $                                      get_interface_length
        use field_class            , only : field
        use parameters_constant    , only : liquid, vapor
        use parameters_input       , only : nx,ny,ne
        use parameters_kind        , only : ikind, rkind

        implicit none

        private
        public :: apply_bubble_ascending_ic


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> for bubble ascending
        !
        !> @date
        !> 27_09_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !---------------------------------------------------------------
        subroutine apply_bubble_ascending_ic(field_used)

          implicit none

          class(field), intent(inout) :: field_used


          !< local variables for the droplet/bubble
          integer        :: phase_at_center
          real(rkind)    :: T0,xc,yc,a,b
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
          omega_v1 =  1.0
          omega_v2 = -1.0


          !<initialize the mass and total energy fields
          do j=1, ny
             do i=1, nx

                x = field_used%x_map(i)
                y = field_used%y_map(j)

                field_used%nodes(i,j,1)=mass_density_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,phase_at_center)
                
                field_used%nodes(i,j,4)=total_energy_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,
     $               field_used%nodes(i,j,1),T0)

             end do
          end do


          !< initialize the momentum and update the total energy fields
          do j=1, ny
             do i=1, nx
                
                x = field_used%x_map(i)
                y = field_used%y_map(j)

                velocity_vortex1 = get_vortex_velocity(
     $               x, y, x_v1, y_v1, r_v1, omega_v1)
                velocity_vortex2 = get_vortex_velocity(
     $               x, y, x_v2, y_v2, r_v2, omega_v2)

                velocity_x = velocity_vortex1(1)+velocity_vortex2(1)
                velocity_y = velocity_vortex1(2)+velocity_vortex2(2)

                field_used%nodes(i,j,2)=
     $               field_used%nodes(i,j,1)*velocity_x
                field_used%nodes(i,j,3)=
     $               field_used%nodes(i,j,1)*velocity_y

                if(rkind.eq.8) then
                   field_used%nodes(i,j,4)=field_used%nodes(i,j,4)+
     $                  1.0d0/2.0d0*field_used%nodes(i,j,1)*(
     $                  velocity_x**2+velocity_y**2)
                else
                   field_used%nodes(i,j,4)=field_used%nodes(i,j,4)+
     $                  1./2.*field_used%nodes(i,j,1)*(
     $                  velocity_x**2+velocity_y**2)
                end if

             end do
          end do

        end subroutine apply_bubble_ascending_ic

      end module dim2d_bubble_ascending_module
