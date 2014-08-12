      !> @file
      !> module encapsulating subroutines to compute
      !> the steady state initial conditions for
      !> the Navier-Stokes equations
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines to compute
      !> the steady state initial conditions for
      !> the Navier-Stokes equations
      !
      !> @date
      !> 12_08_2014 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module ns2d_steadystate_module

        use parameters_input    , only : nx,ny,ne
        use parameters_kind     , only : ikind, rkind

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
        !> 12_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !---------------------------------------------------------------
        subroutine apply_steady_state_ic(nodes)

          implicit none

          real(rkind), dimension(:,:,:), intent(inout) :: nodes
          integer(ikind) :: i,j

          do j=1, size(nodes,2)
             do i=1, size(nodes,1)
                
                nodes(i,j,1) =  1.0
                nodes(i,j,2) =  0.0
                nodes(i,j,3) =  0.0
                nodes(i,j,4) =  2.0

             end do
          end do

        end subroutine apply_steady_state_ic

      end module ns2d_steadystate_module
