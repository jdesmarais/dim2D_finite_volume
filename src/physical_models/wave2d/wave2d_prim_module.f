      !> @file
      !> module encapsulating subroutines computing primary variables
      !> out of conservative variables for the wave 2d governing equations
      !
      !> @author
      !> Julien L. Desmarais
      !
      !> @brief
      !> module encapsulating subroutines computing primary variables
      !> out of conservative variables for the wave 2d governing equations
      !
      !> @date
      !> 28_07_2014 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module wave2d_prim_module

        use n_coords_module, only :
     $     get_n1_coord,
     $     get_n2_coord

        use parameters_kind , only :
     $     ikind,
     $     rkind

        implicit none

        private
        public ::
     $       position,
     $       velocity_x,
     $       velocity_y,
     $       velocity_n1,
     $       velocity_n2

        contains


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the position \f$ u \f$
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
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
        !> 28_07_2014 - initial version - J.L. Desmarais
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


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute the velocity along the y-axis
        !> \f$ w \f$
        !
        !> @date
        !> 28_07_2014 - initial version - J.L. Desmarais
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
        !> \f$ w \f$ evaluated at [i,j]
        !---------------------------------------------------------------
        function velocity_y(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var=nodes(i,j,3)

        end function velocity_y


        function velocity_n1(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var = get_n1_coord(nodes(i,j,2), nodes(i,j,3))

        end function velocity_n1


        function velocity_n2(nodes,i,j) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                               :: var

          var = get_n2_coord(nodes(i,j,2), nodes(i,j,3))

        end function velocity_n2

      end module wave2d_prim_module
