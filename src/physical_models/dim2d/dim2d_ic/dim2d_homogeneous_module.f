      !> @file
      !> module encapsulating subroutines to compute
      !> the homogeneous initial conditions for
      !> the Diffuse Interface Model
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute
      !> the homogeneous initial conditions for
      !> the Diffuse Interface Model
      !
      !> @date
      !> 25_09_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_homogeneous_module

        use dim2d_parameters     , only : cv_r, T_c
        use dim2d_state_eq_module, only : get_mass_density_liquid
        use field_class          , only : field
        use parameters_input     , only : nx,ny
        use parameters_kind      , only : ikind, rkind

        implicit none

        private
        public :: apply_homogeneous_ic


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> for a homogeneous liquid
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !---------------------------------------------------------------
        subroutine apply_homogeneous_ic(field_used)

          implicit none

          class(field), intent(inout) :: field_used

          integer(ikind) :: i,j
          real(rkind)    :: T0_degrees, T0, d_liq, E_liq


          !< enter the homogneous system temperature
          !> in celsius degrees
          T0_degrees   = 100

          !< compute the corresponding reduced temperature
          !> reduced liquid mas density and the total energy
          !> at the temperature asked by the user
          T0    = (T0_degrees+273.15)/T_c
          d_liq = get_mass_density_liquid(T0)
          E_liq = get_total_energy(d_liq,T0)

          if(rkind.eq.8) then
             do j=1, ny
                do i=1, nx
                   
                   field_used%nodes(i,j,1) =  d_liq
                   field_used%nodes(i,j,2) =  0.0d0
                   field_used%nodes(i,j,3) =  0.0d0
                   field_used%nodes(i,j,4) =  E_liq
                   
                end do
             end do
          else
             do j=1, ny
                do i=1, nx
                   
                   field_used%nodes(i,j,1) =  d_liq
                   field_used%nodes(i,j,2) =  0.0
                   field_used%nodes(i,j,3) =  0.0
                   field_used%nodes(i,j,4) =  E_liq
                   
                end do
             end do
          end if

        end subroutine apply_homogeneous_ic

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the total energy
        !> for a homogeneous liquid given its mass
        !> density and its temperature
        !
        !> @date
        !> 25_09_2013 - initial version - J.L. Desmarais
        !
        !>@param mass_density
        !> mass density
        !
        !>@param temperature
        !> temperature
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

      end module dim2d_homogeneous_module
