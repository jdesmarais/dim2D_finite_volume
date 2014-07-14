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
        use parameters_input     , only : nx,ny,ne
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
        subroutine apply_steady_state_ic(nodes)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          integer(ikind) :: i,j
          real(rkind)    :: T0, d_liq


          T0   = 1.0
          d_liq = get_mass_density_liquid(T0)


          do j=1, size(nodes,2)
             do i=1, size(nodes,1)
                
                nodes(i,j,1) =  d_liq
                nodes(i,j,2) = -1.
                nodes(i,j,3) =  1.
                nodes(i,j,4) =  cv_r*T0

             end do
          end do

        end subroutine apply_steady_state_ic

      end module dim2d_steadystate_module
