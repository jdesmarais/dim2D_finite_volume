      module ns1d_fluxes_module

        use field_sd_class , only : field_sd
        use parameters_kind, only : rkind, ikind
        use parameters_ns1d, only : gamma, re, pr, viscous_r,
     $                              epsilon, mach_infty

        implicit none


        private
        public :: compute_viscid2,
     $            compute_temperature,
     $            compute_viscid3,
     $            compute_pressure,
     $            compute_inviscid2,
     $            compute_inviscid3,
     $            compute_fluxes


        contains

        function compute_viscid2(velocity)

          implicit none

          real(rkind), intent(in) :: velocity
          real(rkind)             :: compute_viscid2


          if(rkind.eq.8) then
             compute_viscid2 = (2.0d0 + viscous_r)*velocity
          else
             compute_viscid2 = (2.0 + viscous_r)*velocity
          end if

        end function compute_viscid2


        function compute_temperature(mass_density, velocity, total_energy)

          implicit none
          
          real(rkind), intent(in) :: mass_density
          real(rkind), intent(in) :: velocity
          real(rkind), intent(in) :: total_energy
          real(rkind)             :: compute_temperature


          if(rkind.eq.8) then
             compute_temperature = gamma*(gamma-1.0d0)*mach_infty**2*(
     $            total_energy/mass_density - 0.5d0*velocity**2)
          else
             compute_temperature = gamma*(gamma-1.0)*mach_infty**2*(
     $            total_energy/mass_density - 0.5*velocity**2)
          end if

        end function compute_temperature


        function compute_viscid3(viscid2, velocity, temperature)

          implicit none

          real(rkind), intent(in) :: viscid2
          real(rkind), intent(in) :: velocity
          real(rkind), intent(in) :: temperature
          real(rkind)             :: compute_viscid3


          if(rkind.eq.8) then
             compute_viscid3 = viscid2*0.5d0*velocity + temperature/(Pr*(gamma-1.0d0)*mach_infty**2)
          else
             compute_viscid3 = viscid2*0.5*velocity +   temperature/(Pr*(gamma-1.0)*mach_infty**2)
          end if

        end function compute_viscid3


        function compute_pressure(mass_density,temperature) result(pressure)

          implicit none

          real(rkind), intent(in) :: mass_density
          real(rkind), intent(in) :: temperature
          real(rkind)             :: pressure


          pressure = mass_density*temperature/(gamma*mach_infty**2)

        end function compute_pressure


        function compute_inviscid2(mass_density,velocity,pressure)

          implicit none

          real(rkind), intent(in) :: mass_density
          real(rkind), intent(in) :: velocity
          real(rkind), intent(in) :: pressure
          real(rkind)             :: compute_inviscid2


          compute_inviscid2 = mass_density*velocity**2 + pressure

        end function compute_inviscid2


        function compute_inviscid3(total_energy, pressure, velocity)

          implicit none

          real(rkind), intent(in) :: total_energy
          real(rkind), intent(in) :: pressure
          real(rkind), intent(in) :: velocity
          real(rkind)             :: compute_inviscid3


          compute_inviscid3 = (total_energy+pressure)*velocity

        end function compute_inviscid3        


        subroutine compute_fluxes(
     $       field_sd_used,
     $       nodes,
     $       fluxes)


          implicit none


          class(field_sd)               , intent(in)   :: field_sd_used
          real(rkind)   , dimension(:,:), intent(in)   :: nodes
          real(rkind)   , dimension(:,:), intent(inout):: fluxes


          integer(ikind) :: bc_size
          integer(ikind) :: i
          real(rkind)    :: dx

          real(rkind), dimension(:), allocatable :: viscid2
          real(rkind), dimension(:), allocatable :: viscid3
          real(rkind), dimension(:), allocatable :: inviscid2
          real(rkind), dimension(:), allocatable :: inviscid3


          real(rkind) :: velocity, T, pressure


          !initialize the local variables
          bc_size = field_sd_used%get_bc_size()
          dx      = field_sd_used%get_dx(2)

          
          !allocate intermediate variables
          allocate(viscid2(size(nodes,1)))
          allocate(viscid3(size(nodes,1)))
          allocate(inviscid2(size(nodes,1)))
          allocate(inviscid3(size(nodes,1)))


          !compute the intermediate variable tables
          !before applying the numerical differentiation
          !-------------------------------------------------------------
          do i=1, size(nodes,1)

             velocity = nodes(i,2)/nodes(i,1)

             viscid2(i)  = compute_viscid2(velocity)
             
             T = compute_temperature(nodes(i,1),velocity,nodes(i,3))

             viscid3(i) = compute_viscid3(
     $            viscid2(i),
     $            velocity,
     $            T)

             pressure = compute_pressure(
     $            nodes(i,1),
     $            T)

             inviscid2(i) = compute_inviscid2(
     $            nodes(i,1),
     $            velocity,
     $            pressure)
             
             inviscid3(i) = compute_inviscid3(
     $            nodes(i,3),
     $            pressure,
     $            velocity)

          end do
          !-------------------------------------------------------------


          !compute the fluxes-------------------------------------------
          do i=bc_size+1, size(fluxes,1)-bc_size

             fluxes(i,1) = field_sd_used%f(nodes(:,2),i)

             fluxes(i,2) = field_sd_used%f(inviscid2,i)
     $            - epsilon*field_sd_used%dfdx(viscid2,i)

             fluxes(i,3) = field_sd_used%f(inviscid3,i)
     $            - epsilon*field_sd_used%dfdx(viscid3,i)

          end do


          !deallocate the intermediate tables
          deallocate(viscid2)
          deallocate(viscid3)
          deallocate(inviscid2)
          deallocate(inviscid3)

        end subroutine compute_fluxes

      end module ns1d_fluxes_module
