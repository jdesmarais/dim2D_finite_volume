      !> @file
      !> class encapsulating subroutines for the space discretization
      !> using finite difference for the Mattsson operators
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines for the space discretisation
      !
      !> @date
      !> 31_07_2014 - Initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module sd_operators_fd_module

        use interface_primary, only : get_primary_var
        use parameters_kind  , only : ikind, rkind

        implicit none

        private
        public ::
     $       gradient_x_interior,
     $       gradient_y_interior,
     $       gradient_x_x_oneside_L0,
     $       gradient_x_x_oneside_L1,
     $       gradient_x_x_oneside_R1,
     $       gradient_x_x_oneside_R0,
     $       gradient_y_y_oneside_L0,
     $       gradient_y_y_oneside_L1,
     $       gradient_y_y_oneside_R1,
     $       gradient_y_y_oneside_R0

        
        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial x}\big|_{i,j}=
        !> \frac{1}{2 \Delta x}(-u_{i-1,j} + u_{i+1,j})
        !
        !> @date
        !> 31_07_2014 - initial version  - J.L. Desmarais
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
        !>@param proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_x_interior(
     $     nodes,i,j,proc,dx)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)   , intent(in) :: dx
          real(rkind)                :: var

          if(rkind.eq.8) then

             !TAG INLINE
             var = 0.5d0/dx*(
     $            - proc(nodes,i-1,j)
     $            + proc(nodes,i+1,j))
          else
             var = 0.5d0/dx*(
     $            - proc(nodes,i-1,j)
     $            + proc(nodes,i+1,j))
          end if

        end function gradient_x_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial y}\big|_{i,j}=
        !> \frac{1}{2 \Delta y}(-u_{i,j-1} + u_{i,j+1})
        !
        !> @date
        !> 31_07_2014 - initial version  - J.L. Desmarais
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
        !>@param proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_y_interior(
     $     nodes,i,j,proc,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)   , intent(in) :: dy
          real(rkind)                :: var

          if(rkind.eq.8) then

             !TAG INLINE
             var = 0.5d0/dy*(
     $            - proc(nodes,i,j-1)
     $            + proc(nodes,i,j+1))
          else
             var = 0.5d0/dy*(
     $            - proc(nodes,i,j-1)
     $            + proc(nodes,i,j+1))
          end if

        end function gradient_y_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial x}\big|_{i,j}=
        !> \frac{1}{\Delta x}(-u_{i,j}+u_{i+1,j})\f$
        !
        !> @date
        !> 31_07_2014 - initial version  - J.L. Desmarais
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
        !>@param proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_x_x_oneside_L0(
     $     nodes,i,j,proc,dx)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !TAG INLINE
             var = 1.0d0/dx*(
     $            -proc(nodes,i,j)
     $            +proc(nodes,i+1,j))
          else
             var = 1.0/dx*(
     $            -proc(nodes,i,j)
     $            +proc(nodes,i+1,j))
          end if

        end function gradient_x_x_oneside_L0


        function gradient_x_x_oneside_L1(
     $     nodes,i,j,proc,dx)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                               :: var

          var = gradient_x_interior(nodes,i,j,proc,dx)

        end function gradient_x_x_oneside_L1


        function gradient_x_x_oneside_R1(
     $     nodes,i,j,proc,dx)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                               :: var

          var = gradient_x_interior(nodes,i,j,proc,dx)

        end function gradient_x_x_oneside_R1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial x}\big|_{i,j}=
        !> \frac{1}{\Delta x}( - u_{i-1,j} + u_{i,j} ) \f$
        !
        !> @date
        !> 31_07_2014 - initial version  - J.L. Desmarais
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
        !>@param proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dx
        !> grid step along the x-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_x_x_oneside_R0(
     $     nodes,i,j,proc,dx)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !TAG INLINE
             var = 1.0d0/dx*(
     $            - proc(nodes,i-1,j)
     $            + proc(nodes,i,j))
          else
             var = 1.0/dx*(
     $            - proc(nodes,i-1,j)
     $            + proc(nodes,i,j))
          end if

        end function gradient_x_x_oneside_R0

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial x}\big|_{i,j}=
        !> \frac{1}{\Delta y}(-u_{i,j}+u_{i,j+1})\f$
        !
        !> @date
        !> 31_07_2014 - initial version  - J.L. Desmarais
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
        !>@param proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_y_y_oneside_L0(
     $     nodes,i,j,proc,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !TAG INLINE
             var = 1.0d0/dy*(
     $            -proc(nodes,i,j)
     $            +proc(nodes,i,j+1))
          else
             var = 1.0/dy*(
     $            -proc(nodes,i,j)
     $            +proc(nodes,i,j+1))
          end if

        end function gradient_y_y_oneside_L0


        function gradient_y_y_oneside_L1(
     $     nodes,i,j,proc,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var = gradient_y_interior(nodes,i,j,proc,dy)

        end function gradient_y_y_oneside_L1


        function gradient_y_y_oneside_R1(
     $     nodes,i,j,proc,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var = gradient_y_interior(nodes,i,j,proc,dy)

        end function gradient_y_y_oneside_R1

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial y}\bigg|_{i,j}=
        !> \frac{1}{\Delta y}( - u_{i,j-1} + u_{i,j} ) \f$
        !
        !> @date
        !> 31_07_2014 - initial version  - J.L. Desmarais
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
        !>@param proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_y_y_oneside_R0(
     $     nodes,i,j,proc,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !TAG INLINE
             var = 1.0d0/dy*(
     $            - proc(nodes,i,j-1)
     $            + proc(nodes,i,j))
          else
             var = 1.0/dy*(
     $            - proc(nodes,i,j-1)
     $            + proc(nodes,i,j))
          end if

        end function gradient_y_y_oneside_R0
        
      end module sd_operators_fd_module
