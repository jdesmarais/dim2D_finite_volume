      !> @file
      !> module implementing the abstract interface for procedures
      !> computing variable out of the conservative variables
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> abstract interface for procedures computing primary variables
      !> from conservative variables
      !
      !> @date
      !> 08_08_2013 - initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module interface_primary

        use parameters_kind, only : ikind, rkind

        implicit none

        private
        public :: get_primary_var,
     $            gradient_proc,
     $            get_secondary_var


        abstract interface


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for procedure computing any variable out of
          !> the conservative variables
          !
          !> @date
          !> 07_08_2013 - initial version - J.L. Desmarais
          !
          !>@param nodes
          !> array with the grid point data
          !
          !>@param i
          !> index along x-axis where the data is evaluated
          !
          !>@param j
          !> index along y-axis where the data is evaluated
          !
          !>@param var
          !> primary variable evaluated at [i,j]
          !--------------------------------------------------------------
          function get_primary_var(nodes,i,j) result(var)

            import ikind
            import rkind

            real(rkind), dimension(:,:,:), intent(in) :: nodes
            integer(ikind)               , intent(in) :: i
            integer(ikind)               , intent(in) :: j
            real(rkind)                               :: var

          end function get_primary_var          


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for computing the gradient
          !
          !> @date
          !> 31_07_2014 - initial version - J.L. Desmarais
          !
          !>@param nodes
          !> array with the grid point data
          !
          !>@param i
          !> index along x-axis where the gradient is evaluated
          !
          !>@param j
          !> index along y-axis where the gradient is evaluated
          !
          !>@param proc
          !> interface to compute the primitive variable from
          !> the governing variables
          !
          !>@param dn
          !> grid space step
          !
          !>@param var
          !> primary variable evaluated at [i,j]
          !--------------------------------------------------------------
          function gradient_proc(nodes,i,j,proc,dn) result(var)

            import get_primary_var
            import ikind
            import rkind

            real(rkind), dimension(:,:,:), intent(in) :: nodes
            integer(ikind)               , intent(in) :: i
            integer(ikind)               , intent(in) :: j
            procedure(get_primary_var)                :: proc
            real(rkind)                  , intent(in) :: dn
            real(rkind)                               :: var

          end function gradient_proc


          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface for procedure computing secondary variables
          !> i.e. non local variables (they require neighboring grid
          !> points to compute a gradient for example)
          !
          !> @date
          !> 07_08_2013 - initial version - J.L. Desmarais
          !
          !>@param nodes
          !> array with the grid point data
          !
          !>@param i
          !> index along x-axis where the data is evaluated
          !
          !>@param j
          !> index along y-axis where the data is evaluated
          !
          !>@param dx
          !> grid step along the x-axis
          !
          !>@param dy
          !> grid step along the y-axis
          !
          !>@param gradient_x
          !> procedure to compute the gradient along the x-direction
          !
          !>@param gradient_y
          !> procedure to compute the gradient along the y-direction
          !
          !>@param var
          !> primary variable evaluated at [i,j]
          !--------------------------------------------------------------
          function get_secondary_var(
     $     nodes,i,j,dx,dy,gradient_x,gradient_y)
     $     result(var)

            import gradient_proc
            import ikind
            import rkind            

            real(rkind), dimension(:,:,:), intent(in) :: nodes
            integer(ikind)               , intent(in) :: i
            integer(ikind)               , intent(in) :: j
            real(rkind)                  , intent(in) :: dx
            real(rkind)                  , intent(in) :: dy
            procedure(gradient_proc)                  :: gradient_x
            procedure(gradient_proc)                  :: gradient_y
            real(rkind)                               :: var

          end function get_secondary_var

        end interface

      end module interface_primary
