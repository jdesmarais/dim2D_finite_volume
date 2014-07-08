      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute
      !> the initial conditions of drop retraction
      !
      !> @date
      !> 27_09_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_nonequilibrium_module

        use dim2d_dropbubble_module, only : mass_density_ellipsoid,
     $                                      total_energy_ellipsoid
        use dim2d_state_eq_module  , only : get_mass_density_liquid,
     $                                      get_mass_density_vapor,
     $                                      get_interface_length
        use field_class            , only : field
        use parameters_constant    , only : liquid, vapor
        use parameters_input       , only : nx,ny,ne
        use parameters_kind        , only : ikind, rkind

        implicit none

        private
        public :: apply_nonequilibrium_ic


        contains

        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> for drop retraction
        !
        !> @date
        !> 27_09_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !---------------------------------------------------------------
        subroutine apply_nonequilibrium_ic(field_used)

          implicit none

          class(field), intent(inout) :: field_used


          !< local variables for the droplet/bubble
          integer        :: phase_at_center
          real(rkind)    :: T0,xc,yc,a,b
          real(rkind)    :: dliq_eq,dvap_eq
          real(rkind)    :: dliq,dvap,li


          !< local variables for the initialization
          integer(ikind) :: i,j
          real(rkind)    :: x,y

          
          !< choose the phase at the domain center:
          !> is it a droplet of liquid in a vapor medium ? -> vapor
          !> is it a bubble  of vapor in a liquid medium ? -> liquid
          phase_at_center = liquid

          !<set the initial temperature in the field
          T0 = 0.995
          
          !<set the center of the droplet
          xc=0.
          yc=0.


          !<get the mass densities corresponding to the
          !>liquid and vapor phases for the initial
          !>temperature field
          dliq_eq = get_mass_density_liquid(T0)
          dvap_eq = get_mass_density_vapor(T0)

          dliq = 2*d_liq_eq
          dvap = (d_liq+d_vap)/2.0d0

          !<get the interface length corresponding
          !>to the initial temperature field
          li = get_interface_length(T0)

          !<set the major and minor axes of the bubble ellipse
          a=6.0d0*li
          b=a


          !<initialize the fields
          do j=1, ny
             do i=1, nx

                x = field_used%x_map(i)
                y = field_used%y_map(j)

                field_used%nodes(i,j,1)=mass_density_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,phase_at_center)

                field_used%nodes(i,j,2)=0.0d0

                field_used%nodes(i,j,3)=0.0d0
                
                field_used%nodes(i,j,4)=total_energy_ellipsoid(
     $               x,y,xc,yc,a,b,li,dliq,dvap,
     $               field_used%nodes(i,j,1),T0)

             end do
          end do

        end subroutine apply_nonequilibrium_ic

      end module dim2d_nonequilibrium_module

