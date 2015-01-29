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
      !> 29_01_2015 - Initial version - J.L. Desmarais
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
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{1}{2 \Delta n_1}(-u_{i-1,j+1} + u_{i+1,j-1})
        !
        !> @date
        !> 29_01_2015 - initial version  - J.L. Desmarais
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
        !> grid step along the (x-y)-axis
        !
        !>@return var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_interior(
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
     $            - proc(nodes,i-1,j+1)
     $            + proc(nodes,i+1,j-1))
          else
             var = 0.5d0/dx*(
     $            - proc(nodes,i-1,j+1)
     $            + proc(nodes,i+1,j-1))
          end if

        end function gradient_n1_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_2}\big|_{i,j}=
        !> \frac{1}{2 \Delta n_2}(-u_{i-1,j-1} + u_{i+1,j+1})
        !
        !> @date
        !> 29_01_2015 - initial version  - J.L. Desmarais
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
        !> grid step along the (x+y)-axis
        !
        !>@return var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_interior(
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
     $            - proc(nodes,i-1,j-1)
     $            + proc(nodes,i+1,j+1))
          else
             var = 0.5d0/dy*(
     $            - proc(nodes,i-1,j-1)
     $            + proc(nodes,i+1,j+1))
          end if

        end function gradient_n2_interior


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{1}{\Delta n_1}(-u_{i,j}+u_{i+1,j-1})\f$
        !
        !> @date
        !> 29_01_2015 - initial version  - J.L. Desmarais
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
        !> grid step along the (x-y)-axis
        !
        !>@return var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_oneside_L0(
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
     $            +proc(nodes,i+1,j-1))
          else
             var = 1.0/dx*(
     $            -proc(nodes,i,j)
     $            +proc(nodes,i+1,j-1))
          end if

        end function gradient_n1_oneside_L0


        function gradient_n1_oneside_L1(
     $     nodes,i,j,proc,dx)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                               :: var

          var = gradient_n1_interior(nodes,i,j,proc,dx)

        end function gradient_n1_oneside_L1


        function gradient_n1_oneside_R1(
     $     nodes,i,j,proc,dx)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                               :: var

          var = gradient_n1_interior(nodes,i,j,proc,dx)

        end function gradient_n1_oneside_R1


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_1}\big|_{i,j}=
        !> \frac{1}{\Delta n_1}( - u_{i-1,j+1} + u_{i,j} ) \f$
        !
        !> @date
        !> 29_01_2015 - initial version  - J.L. Desmarais
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
        !> grid step along the (x-y)-axis
        !
        !>@return var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n1_oneside_R0(
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
     $            - proc(nodes,i-1,j+1)
     $            + proc(nodes,i,j))
          else
             var = 1.0/dx*(
     $            - proc(nodes,i-1,j+1)
     $            + proc(nodes,i,j))
          end if

        end function gradient_n1_oneside_R0

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_2}\big|_{i,j}=
        !> \frac{1}{\Delta n_2}(-u_{i,j}+u_{i+1,j+1})\f$
        !
        !> @date
        !> 29_01_2015 - initial version  - J.L. Desmarais
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
        !> grid step along the (x+y)-axis
        !
        !>@return var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_oneside_L0(
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
     $            +proc(nodes,i+1,j+1))
          else
             var = 1.0/dy*(
     $            -proc(nodes,i,j)
     $            +proc(nodes,i+1,j+1))
          end if

        end function gradient_n2_oneside_L0


        function gradient_n2_oneside_L1(
     $     nodes,i,j,proc,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var = gradient_n2_interior(nodes,i,j,proc,dy)

        end function gradient_n2_oneside_L1


        function gradient_n2_oneside_R1(
     $     nodes,i,j,proc,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

          var = gradient_n2_interior(nodes,i,j,proc,dy)

        end function gradient_n2_oneside_R1

        
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial n_2}\bigg|_{i,j}=
        !> \frac{1}{\Delta y}( - u_{i-1,j-1} + u_{i,j} ) \f$
        !
        !> @date
        !> 29_01_2015 - initial version  - J.L. Desmarais
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
        !> grid step along the (x+y)-axis
        !
        !>@return var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function gradient_n2_oneside_R0(
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
     $            - proc(nodes,i-1,j-1)
     $            + proc(nodes,i,j))
          else
             var = 1.0/dy*(
     $            - proc(nodes,i-1,j-1)
     $            + proc(nodes,i,j))
          end if

        end function gradient_n2_oneside_R0
        
      end module sd_operators_fd_ncoords_module
