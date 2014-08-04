      !> @file
      !> module encapsulating subroutines computing primary variables
      !> out of conservative variables for the wave 1d governing equations
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines computing primary variables
      !> out of conservative variables for the wave 1d governing equations
      !
      !> @date
      !> 04_08_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module wave1d_prim_module

        use parameters_kind , only : ikind, rkind

        implicit none

        private
        public :: position, velocity_x

        contains


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the position \f$ u \f$
        !
        !> @date
        !> 04_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ \rho \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function position(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,1)

        end function position


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the velocity along the x-axis
        !> \f$ v \f$
        !
        !> @date
        !> 04_08_2014 - initial version - J.L. Desmarais
        !
        !>@param nodes
        !> array with the grid point data
        !>
        !>@param i
        !> index along x-axis where the data is evaluated
        !>
        !>@param j
        !> index along y-axis where the data is evaluated
        !>
        !>@param var
        !> \f$ v \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function velocity_x(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,2)

        end function velocity_x

      end module wave1d_prim_module
