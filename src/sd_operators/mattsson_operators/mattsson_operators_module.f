      module mattsson_operators_module

        use interface_primary , only : get_primary_var,
     $                                 get_secondary_var
        use parameters_kind   , only : ikind, rkind

        implicit none

        
        private
        public :: f_x_oneside_L0,
     $            dfdx_x_oneside_L0,
     $            dfdx_x_oneside_L0_nl,
     $            dfdy_x_oneside_L0,
     $            d2fdx2_x_oneside_L0,
     $            d2fdy2_x_oneside_L0,
     $            d2fdxdy_x_oneside_L0,
     $            dgdx_x_oneside_L0,
     $            d2gdx2_x_oneside_L0,
     $            d2gdxdy_x_oneside_L0


        contains


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ u_{i+\frac{1}{2},j}=
        !>\frac{3}{2} u_{i,j} - \frac{1}{2} u_{i+1,j}\f$
        !
        !> @date
        !> 29_07_2014 - initial version  - J.L. Desmarais
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
        function f_x_oneside_L0(nodes,i,j,proc) result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                               :: var

          if(rkind.eq.8) then

             !TAG INLINE
             var = 3.0d0/2.0d0*proc(nodes,i,j)
     $            -1.0d0/2.0d0*proc(nodes,i+1,j)
          else
             var = 3.0/2.0*proc(nodes,i,j)
     $            -1.0/2.0*proc(nodes,i+1,j)
          end if

        end function f_x_oneside_L0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial x} \bigg|_{i+\frac{1}{2}
        !> ,j} = \frac{1}{\Delta x}(- 2 u_{i,j} + 3 u_{i+1,j} - u_{i+2,j})
        !> \f$
        !
        !> @date
        !> 29_07_2014 - initial version  - J.L. Desmarais
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
        function dfdx_x_oneside_L0(nodes,i,j,proc,dx) result(var)

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
     $            -2.0d0*proc(nodes,i,j)
     $            +3.0d0*proc(nodes,i+1,j)
     $            -      proc(nodes,i+2,j))
          else
             var = 1.0/dx*(
     $            -2.0*proc(nodes,i,j)
     $            +3.0*proc(nodes,i+1,j)
     $            -    proc(nodes,i+2,j))
          end if

        end function dfdx_x_oneside_L0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial x} \bigg|_{i+\frac{1}{2}
        !> ,j} = \frac{1}{\Delta x}(- 2 u_{i,j} + 3 u_{i+1,j} - u_{i+2,j})
        !> \f$
        !
        !> @date
        !> 29_07_2014 - initial version  - J.L. Desmarais
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
        function dfdx_x_oneside_L0_nl(nodes,i,j,proc,dx,dy) result(var)

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
             var = 1.0d0/dx*(
     $            -2.0d0*proc(nodes,i,j,dx,dy)
     $            +3.0d0*proc(nodes,i+1,j,dx,dy)
     $            -      proc(nodes,i+2,j,dx,dy))
          else
             var = 1.0/dx*(
     $            -2.0*proc(nodes,i,j,dx,dy)
     $            +3.0*proc(nodes,i+1,j,dx,dy)
     $            -    proc(nodes,i+2,j,dx,dy))
          end if

        end function dfdx_x_oneside_L0_nl


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial y}\bigg|
        !> _{i+\frac{1}{2},j}= \frac{1}{2 \Delta y}
        !> (-u_{i+\frac{1}{2},j-1}+u_{i+\frac{1}{2},j+1})\f$
        !
        !> @date
        !> 29_07_2014 - initial version - J.L. Desmarais
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
        function dfdy_x_oneside_L0(nodes,i,j,proc,dy)
     $     result(var)

          implicit none

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          real(rkind)                  , intent(in) :: dy
          procedure(get_primary_var)                :: proc
          real(rkind)                               :: var

          !TAG INLINE
          if(rkind.eq.8) then
             !TAG INLINE
             var = 1.0d0/(2.0d0*dy)*(
     $            - 3.0d0/2.0d0*proc(nodes,i  ,j-1)
     $            + 1.0d0/2.0d0*proc(nodes,i+1,j-1)
     $            + 3.0d0/2.0d0*proc(nodes,i  ,j+1)
     $            - 1.0d0/2.0d0*proc(nodes,i+1,j+1))
          else
             var = 1.0/(2.0*dy)*(
     $            - 3.0/2.0*proc(nodes,i  ,j-1)
     $            + 1.0/2.0*proc(nodes,i+1,j-1)
     $            + 3.0/2.0*proc(nodes,i  ,j+1)
     $            - 1.0/2.0*proc(nodes,i+1,j+1))
          end if

        end function dfdy_x_oneside_L0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial^2 u}{\partial x^2}\bigg|
        !> _{i+\frac{1}{2},j}= \frac{1}{2 {\Delta x}^2}(3 u_{i-1,j} -
        !> 5 u_{i,j} - u_{i+1,j} + 5 u_{i+2,j} -2 u_{i+3,j}) \f$
        !
        !> @date
        !> 29_07_2014 - initial version - J.L. Desmarais
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
        function d2fdx2_x_oneside_L0(
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
             var = 0.5d0/(dx**2)*(
     $            +3.0d0*proc(nodes,i-1,j)
     $            -5.0d0*proc(nodes,i  ,j)
     $            -      proc(nodes,i+1,j)
     $            +5.0d0*proc(nodes,i+2,j)
     $            -2.0d0*proc(nodes,i+3,j)
     $            )
          else
             var = 0.5/(dx**2)*(
     $            +3.0*proc(nodes,i-1,j)
     $            -5.0*proc(nodes,i  ,j)
     $            -    proc(nodes,i+1,j)
     $            +5.0*proc(nodes,i+2,j)
     $            -2.0*proc(nodes,i+3,j)
     $            )
          end if

        end function d2fdx2_x_oneside_L0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial^2 u}{\partial y^2}
        !> \bigg|_{i+\frac{1}{2},j}= \frac{1}{{\Delta y}^2}(
        !> u_{i+\frac{1}{2},j-1} - 2 u_{i+\frac{1}{2},j} +
        !> u_{i+\frac{1}{2},j+1}) \f$
        !
        !> @date
        !> 29_07_2014 - initial version - J.L. Desmarais
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
        function d2fdy2_x_oneside_L0(
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
             var = (1.0d0/dy**2)*
     $            (
     $             3.0d0/2.0d0*proc(nodes,i,j-1)
     $            -1.0d0/2.0d0*proc(nodes,i+1,j-1)
     $            -2.0d0*(
     $                 3.0d0/2.0d0*proc(nodes,i,j)
     $                -1.0d0/2.0d0*proc(nodes,i+1,j)
     $            )
     $            +3.0d0/2.0d0*proc(nodes,i,j+1)
     $            -1.0d0/2.0d0*proc(nodes,i+1,j+1)
     $            )
          else

             var = (1.0/dy**2)*
     $            (
     $             3.0/2.0*proc(nodes,i,j-1)
     $            -1.0/2.0*proc(nodes,i+1,j-1)
     $            -2.0*(
     $                 3.0/2.0*proc(nodes,i,j)
     $                -1.0/2.0*proc(nodes,i+1,j)
     $            )
     $            +3.0/2.0*proc(nodes,i,j+1)
     $            -1.0/2.0*proc(nodes,i+1,j+1)
     $            )
          end if
        end function d2fdy2_x_oneside_L0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial^2 u}{\partial x \partial y}
        !> \bigg|_{i+\frac{1}{2},j} = \frac{1}{2 \Delta y}
        !> (-\frac{\partial u}{\partial x}\big|_{i+\frac{1}{2},j-1} +
        !> \frac{\partial u}{\partial x}\big|_{i+\frac{1}{2},j+1}) \f$
        !
        !> @date
        !> 29_07_2014 - initial version - J.L. Desmarais
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
        function d2fdxdy_x_oneside_L0(
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
             var =( 2.0d0*proc(nodes,i  ,j-1)
     $            - 3.0d0*proc(nodes,i+1,j-1)
     $            +       proc(nodes,i+2,j-1)
     $            - 2.0d0*proc(nodes,i  ,j+1)
     $            + 3.0d0*proc(nodes,i+1,j+1)
     $            -       proc(nodes,i+2,j+1))*
     $            0.5d0/(dy*dx)
          else
             var =( 2.0*proc(nodes,i  ,j-1)
     $            - 3.0*proc(nodes,i+1,j-1)
     $            +     proc(nodes,i+2,j-1)
     $            - 2.0*proc(nodes,i  ,j+1)
     $            + 3.0*proc(nodes,i+1,j+1)
     $            -     proc(nodes,i+2,j+1))*
     $            0.5/(dy*dx)
          end if

        end function d2fdxdy_x_oneside_L0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial x}|_{i
        !> ,j+\frac{1}{2}}= \frac{1}{\Delta x} (-u_{i,j+\frac{1}{2}}
        !> +u_{i+1,j+\frac{1}{2}})\f$
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        function dgdx_x_oneside_L0(
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

             var = 1.0d0/(12.0d0*dx)*(
     $            +       proc(nodes,i,j-1)
     $            - 7.0d0*proc(nodes,i,j)
     $            - 7.0d0*proc(nodes,i,j+1)
     $            +       proc(nodes,i,j+2)
     $            -       proc(nodes,i+1,j-1)
     $            + 7.0d0*proc(nodes,i+1,j)
     $            + 7.0d0*proc(nodes,i+1,j+1)
     $            -       proc(nodes,i+1,j+2))
          else

             var = 1.0/(12.0*dx)*(
     $            +     proc(nodes,i,j-1)
     $            - 7.0*proc(nodes,i,j)
     $            - 7.0*proc(nodes,i,j+1)
     $            +     proc(nodes,i,j+2)
     $            -     proc(nodes,i+1,j-1)
     $            + 7.0*proc(nodes,i+1,j)
     $            + 7.0*proc(nodes,i+1,j+1)
     $            -     proc(nodes,i+1,j+2))
          end if

        end function dgdx_x_oneside_L0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial^2 u}{\partial x^2}
        !> \big|_{i,j+\frac{1}{2}}= \frac{1}{{\Delta x}^2}(
        !> u_{i,j+\frac{1}{2}}- 2 u_{i+1,j+\frac{1}{2}}
        !> + u_{i+2,j+\frac{1}{2}}) \f$
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        function d2gdx2_x_oneside_L0(
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
             var = 1.0d0/(12.0d0*dx**2) * (
     $            -      proc(nodes,i,j-1)
     $            +7.0d0*proc(nodes,i,j)
     $            +7.0d0*proc(nodes,i,j+1)
     $            -      proc(nodes,i,j+2)
     $            - 2.0d0*(
     $            -      proc(nodes,i+1,j-1)
     $            +7.0d0*proc(nodes,i+1,j)
     $            +7.0d0*proc(nodes,i+1,j+1)
     $            -      proc(nodes,i+1,j+2)
     $            )
     $            -      proc(nodes,i+2,j-1)
     $            +7.0d0*proc(nodes,i+2,j)
     $            +7.0d0*proc(nodes,i+2,j+1)
     $            -      proc(nodes,i+2,j+2)
     $            )
          else
             var = 1./(12.*dx**2) * (
     $            -  proc(nodes,i,j-1)
     $            +7*proc(nodes,i,j)
     $            +7*proc(nodes,i,j+1)
     $            -  proc(nodes,i,j+2)
     $            - 2*(
     $            -  proc(nodes,i+1,j-1)
     $            +7*proc(nodes,i+1,j)
     $            +7*proc(nodes,i+1,j+1)
     $            -  proc(nodes,i+1,j+2)
     $            )
     $            -  proc(nodes,i+2,j-1)
     $            +7*proc(nodes,i+2,j)
     $            +7*proc(nodes,i+2,j+1)
     $            -  proc(nodes,i+2,j+2)
     $            )
          end if

        end function d2gdx2_x_oneside_L0


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial^2 u}{\partial x \partial y}
        !> \bigg|_{i,j+\frac{1}{2}} = \frac{1}{\Delta x}
        !> (-\frac{\partial u}{\partial y}\big|_{i,j+\frac{1}{2}} +
        !> \frac{\partial u}{\partial y}\big|_{i+1,j+\frac{1}{2}}) \f$
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
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
        function d2gdxdy_x_oneside_L0(
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
             var =(
     $              proc(nodes,i,j)
     $            - proc(nodes,i,j+1)
     $            - proc(nodes,i+1,j)
     $            + proc(nodes,i+1,j+1))*
     $            0.5d0/(dx*dy)
          else
             var =(
     $              proc(nodes,i,j)
     $            - proc(nodes,i,j+1)
     $            - proc(nodes,i+1,j)
     $            + proc(nodes,i+1,j+1))
     $            0.5/(dx*dy)
          end if

        end function d2gdxdy_x_oneside_L0

      end module mattsson_operators_module
