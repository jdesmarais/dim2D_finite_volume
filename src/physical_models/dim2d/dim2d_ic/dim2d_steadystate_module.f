      !> @file
      !> module encapsulating subroutines to compute
      !> the steady state initial conditions for
      !> the Diffuse Interface Model
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute
      !> the steady state initial conditions for
      !> the Diffuse Interface Model
      !
      !> @date
      !> 14_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module dim2d_steadystate_module

        use dim2d_parameters     , only : cv_r
        use dim2d_state_eq_module, only : get_mass_density_liquid
        use field_class          , only : field
        use parameters_input     , only : nx,ny
        use parameters_kind      , only : ikind, rkind

        implicit none

        private
        public :: apply_steady_state_ic


        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> subroutine computing the initial conditions
        !> for a steady state
        !
        !> @date
        !> 14_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the main variables
        !---------------------------------------------------------------
        subroutine apply_steady_state_ic(field_used)

          implicit none

          class(field), intent(inout) :: field_used

          integer(ikind) :: i,j
          real(rkind)    :: T0, d_liq


          T0   = 1.0
          d_liq = get_mass_density_liquid(T0)


          do j=1, ny
             do i=1, nx
                
                field_used%nodes(i,j,1) =  d_liq
                field_used%nodes(i,j,2) = -1.
                field_used%nodes(i,j,3) =  1.
                field_used%nodes(i,j,4) =  cv_r*T0

             end do
          end do

        end subroutine apply_steady_state_ic

      end module dim2d_steadystate_module
