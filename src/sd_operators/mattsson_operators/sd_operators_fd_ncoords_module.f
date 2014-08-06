      !> @file
      !> class encapsulating subroutines for the space discretization
      !> using finite difference for the Mattsson operators in the
      !> rotated grid ((x-y), (x+y))
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines for the space discretisation
      !
      !> @date
      !> 06_08_2014 - Initial version - J.L. Desmarais
      !-----------------------------------------------------------------
      module sd_operators_fd_ncoords_module

        use interface_primary, only : get_primary_var
        use parameters_kind  , only : ikind, rkind

        implicit none

        private
        public ::
     $       gradient_n1_interior,
     $       gradient_n2_interior,
     $       gradient_n1_oneside_L0,
     $       gradient_n1_oneside_L1,
     $       gradient_n1_oneside_R1,
     $       gradient_n1_oneside_R0,
     $       gradient_n2_oneside_L0,
     $       gradient_n2_oneside_L1,
     $       gradient_n2_oneside_R1,
     $       gradient_n2_oneside_R0

        
        contains

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\bigg|_{i,j}=
        !> \frac{1}{2 \sqrt{{\Delta x}^2 + {\Delta y}^2}}(
        !> u_{i+1,j-1} - u_{i-1,j+1})
        !
        !> @date
        !> 06_08_2014 - initial version  - J.L. Desmarais
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
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_interior(
     $     nodes,i,j,proc,dx,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !TAG INLINE
             var = 0.5d0/Sqrt(dx**2+dy**2)*(
     $            - proc(nodes,i-1,j+1)
     $            + proc(nodes,i+1,j-1))

          else
             var = 0.5/Sqrt(dx**2+dy**2)*(
     $            - proc(nodes,i-1,j+1)
     $            + proc(nodes,i+1,j-1))
          end if

        end function gradient_n1_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_2}\bigg|_{i,j}=
        !> \frac{1}{2 \sqrt{{\Delta x}^2 + {\Delta y}^2}}(
        !> u_{i+1,j+1} - u_{i-1,j-1})
        !
        !> @date
        !> 06_08_2014 - initial version  - J.L. Desmarais
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
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_interior(
     $     nodes,i,j,proc,dx,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !TAG INLINE
             var = 0.5d0/Sqrt(dx**2+dy**2)*(
     $            - proc(nodes,i-1,j-1)
     $            + proc(nodes,i+1,j+1))
          else
             var = 0.5/Sqrt(dx**2+dy**2)*(
     $            - proc(nodes,i-1,j-1)
     $            + proc(nodes,i+1,j+1))
          end if

        end function gradient_n2_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_2}\bigg|_{i,j}=
        !> \frac{1}{\sqrt{{\Delta x}^2 + {\Delta y}^2}}(
        !> u_{i+1,j-1} - u_{i,j})
        !
        !> @date
        !> 06_08_2014 - initial version  - J.L. Desmarais
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
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_oneside_L0(
     $     nodes,i,j,proc,dx,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !TAG INLINE
             var = 1.0d0/Sqrt(dx**2+dy**2)*(
     $             proc(nodes,i+1,j-1)
     $            -proc(nodes,i,j))
          else
             var = 1.0/Sqrt(dx**2+dy**2)*(
     $             proc(nodes,i+1,j-1)
     $            -proc(nodes,i,j))
          end if

        end function gradient_n1_oneside_L0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\bigg|_{i,j}=
        !> \frac{1}{2 \sqrt{{\Delta x}^2 + {\Delta y}^2}}(
        !> u_{i+1,j-1} - u_{i-1,j+1})
        !
        !> @date
        !> 06_08_2014 - initial version  - J.L. Desmarais
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
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_oneside_L1(
     $     nodes,i,j,proc,dx,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var = gradient_n1_interior(nodes,i,j,proc,dx,dy)

        end function gradient_n1_oneside_L1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\bigg|_{i,j}=
        !> \frac{1}{2 \sqrt{{\Delta x}^2 + {\Delta y}^2}}(
        !> u_{i+1,j-1} - u_{i-1,j+1})
        !
        !> @date
        !> 06_08_2014 - initial version  - J.L. Desmarais
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
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_oneside_R1(
     $     nodes,i,j,proc,dx,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var = gradient_n1_interior(nodes,i,j,proc,dx,dy)

        end function gradient_n1_oneside_R1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\bigg|_{i,j}=
        !> \frac{1}{\sqrt{{\Delta x}^2 + {\Delta y}^2}}(
        !> u_{i,j} - u_{i-1,j+1})
        !
        !> @date
        !> 06_08_2014 - initial version  - J.L. Desmarais
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
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_oneside_R0(
     $     nodes,i,j,proc,dx,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !TAG INLINE
             var = 1.0d0/Sqrt(dx**2+dy**2)*(
     $              proc(nodes,i,j)
     $            - proc(nodes,i-1,j+1))

          else
             var = 1.0/Sqrt(dx**2+dy**2)*(
     $              proc(nodes,i,j)
     $            - proc(nodes,i-1,j+1))
          end if

        end function gradient_n1_oneside_R0

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_2}\bigg|_{i,j}=
        !> \frac{1}{2 \sqrt{{\Delta x}^2 + {\Delta y}^2}}(
        !> u_{i+1,j+1} - u_{i,j})
        !
        !> @date
        !> 06_08_2014 - initial version  - J.L. Desmarais
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
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_oneside_L0(
     $     nodes,i,j,proc,dx,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          if(rkind.eq.8) then
             var = 1.0d0/sqrt(dx**2+dy**2)*(
     $            -proc(nodes,i,j)
     $            +proc(nodes,i+1,j+1))
          else
             var = 1.0/sqrt(dx**2+dy**2)*(
     $            -proc(nodes,i,j)
     $            +proc(nodes,i+1,j+1))
          end if

        end function gradient_n2_oneside_L0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_2}\bigg|_{i,j}=
        !> \frac{1}{2 \sqrt{{\Delta x}^2 + {\Delta y}^2}}(
        !> u_{i+1,j+1} - u_{i-1,j-1})
        !
        !> @date
        !> 06_08_2014 - initial version  - J.L. Desmarais
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
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_oneside_L1(
     $     nodes,i,j,proc,dx,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var = gradient_n2_interior(nodes,i,j,proc,dx,dy)

        end function gradient_n2_oneside_L1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_2}\bigg|_{i,j}=
        !> \frac{1}{2 \sqrt{{\Delta x}^2 + {\Delta y}^2}}(
        !> u_{i+1,j+1} - u_{i-1,j-1})
        !
        !> @date
        !> 06_08_2014 - initial version  - J.L. Desmarais
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
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_oneside_R1(
     $     nodes,i,j,proc,dx,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var = gradient_n2_interior(nodes,i,j,proc,dx,dy)

        end function gradient_n2_oneside_R1

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_2}\bigg|_{i,j}=
        !> \frac{1}{\sqrt{{\Delta x}^2 + {\Delta y}^2}}(
        !> u_{i,j} - u_{i-1,j-1})
        !
        !> @date
        !> 06_08_2014 - initial version  - J.L. Desmarais
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
        !>@param dy
        !> grid step along the y-axis
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_oneside_R0(
     $     nodes,i,j,proc,dx,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !TAG INLINE
             var = 1.0d0/Sqrt(dx**2+dy**2)*(
     $            - proc(nodes,i-1,j-1)
     $            + proc(nodes,i,j))
          else
             var = 1.0/Sqrt(dx**2+dy**2)*(
     $            - proc(nodes,i-1,j-1)
     $            + proc(nodes,i,j))
          end if

        end function gradient_n2_oneside_R0
        
      end module sd_operators_fd_ncoords_module
