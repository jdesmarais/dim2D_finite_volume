      !> @file
      !> module encapsulating subroutines for the computation of
      !> multiphase quantities: saturated water liquid and vapor mass
      !> densities, interface length...
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines for the computation of
      !> multiphase quantities: saturated water liquid and vapor mass
      !> densities, interface length...
      !
      !> @date
      !> 14_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_state_eq_module

        use dim2d_parameters, only : we
        use parameters_kind , only : rkind

        implicit none

        private
        public :: get_interface_length,
     $            get_mass_density_liquid,
     $            get_mass_density_vapor


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get the interface length at a given temperature
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param temperature
        !> temperature at which the interface length is computed
        !---------------------------------------------------------------
        function get_interface_length(temperature)
     $       result(interface_lgh)

          implicit none

          real(rkind), intent(in) :: temperature
          real(rkind)             :: interface_lgh

          if(rkind.eq.8) then
             interface_lgh=2.0d0/SQRT(we)*(
     $            -0.19d0+1.65d0/(SQRT(1.0d0-temperature)))
          else
             interface_lgh=2./SQRT(We)*(
     $            -0.19+1.65/(SQRT(1.-temperature)))
          end if

        end function get_interface_length


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get mass density for saturated liquid water
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param T
        !> temperature at which the mass density is computed
        !---------------------------------------------------------------
        function get_mass_density_liquid(temperature) result(md_liquid)

          implicit none

          real(rkind), intent(in) :: temperature
          real(rkind)             :: md_liquid
          
          if(rkind.eq.8) then
             md_liquid = 1.0d0+2.08d0*SQRT(1.0d0-temperature)
          else
             md_liquid = 1.+2.08*SQRT(1.-temperature)
          end if

        end function get_mass_density_liquid


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> get mass density for saturated vapor water
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param T
        !> temperature at which the mass density is computed
        !---------------------------------------------------------------
        function get_mass_density_vapor(temperature) result(md_vapor)

          implicit none

          real(rkind), intent(in) :: temperature
          real(rkind)             :: md_vapor
          
          if(rkind.eq.8) then
             md_vapor = 1.0d0-1.86d0*SQRT(1.0d0-temperature)
          else
             md_vapor = 1.-1.86*SQRT(1.-temperature)
          end if

        end function get_mass_density_vapor

      end module dim2d_state_eq_module
