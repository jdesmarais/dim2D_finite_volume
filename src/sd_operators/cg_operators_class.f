      !> @file
      !> class encapsulating subroutines for the space discretization
      !> using the operators developed by Cockburn and Gau in 
      !> “A model numerical scheme for the propagation of phase
      !> transition in solids”, SIAM J. Sci. Computing, Vol 17, No.
      !> 5, pp 1092-1121, Sep 1996
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines for the space discretisation
      !> scheme using Cockburn and Gau's operators
      !
      !> @date
      !> 30_07_2012 - initial version                   - J.L. Desmarais
      !> 04_06_2013 - functions for crossed derivatives - J.L. Desmarais
      !> 07_08_2013 - transfered in lernaeanhydra_opt   - J.L. Desmarais
      !-----------------------------------------------------------------
      module cg_operators_class

        use field_class       , only : field
        use interface_primary , only : get_primary_var
        use parameters_kind   , only : ikind, rkind
        use sd_operators_class, only : sd_operators

        implicit none

        private
        public :: cg_operators


        !> @class cg_operators
        !> class encapsulating Cockburn and Gau's spatial
        !> discretization operators
        !>
        !> @param get_bc_size
        !> get the boundary layer size
        !>
        !> @param f
        !> evaluate data at [i+1/2,j]
        !>
        !> @param dfdx
        !> evaluate \f$\frac{\partial}{\partial x}\f$ at [i+1/2,j]
        !>
        !> @param dfdy
        !> evaluate \f$\frac{\partial}{\partial y}\f$ at [i+1/2,j]
        !>
        !> @param d2fdx2
        !> evaluate \f$\frac{\partial}{\partial x^2}\f$ at [i+1/2,j]
        !>
        !> @param d2fdy2
        !> evaluate \f$\frac{\partial}{\partial y^2}\f$ at [i+1/2,j]
        !>
        !> @param d2fdxdy
        !> evaluate \f$\frac{\partial}{\partial x \partial y}\f$
        !> at [i+1/2,j]
        !>        
        !> @param g
        !> evaluate data at [i,j+1/2]
        !>
        !> @param dgdx
        !> evaluate \f$\frac{\partial}{\partial x}\f$ at [i,j+1/2]
        !>
        !> @param dgdy
        !> evaluate \f$\frac{\partial}{\partial y}\f$ at [i,j+1/2]
        !>
        !> @param d2gdx2
        !> evaluate \f$\frac{\partial}{\partial x^2}\f$ at [i,j+1/2]
        !>
        !> @param d2gdy2
        !> evaluate \f$\frac{\partial}{\partial y^2}\f$ at [i,j+1/2]
        !>
        !> @param d2gdxdy
        !> evaluate \f$\frac{\partial}{\partial x \partial y}\f$
        !> at [i,j+1/2]
        !---------------------------------------------------------------
        type, extends(sd_operators) :: cg_operators

          contains

          procedure, nopass :: get_bc_size => get_bc_size_cockburnandgau

          procedure, nopass :: f           => f_cockburnandgau
          procedure, nopass :: dfdx        => dfdx_cockburnandgau        
          procedure, nopass :: dfdy        => dfdy_cockburnandgau
          procedure, nopass :: d2fdx2      => d2fdx2_cockburnandgau
          procedure, nopass :: d2fdy2      => d2fdy2_cockburnandgau
          procedure, nopass :: d2fdxdy     => d2fdxdy_cockburnandgau

          procedure, nopass :: g           => g_cockburnandgau
          procedure, nopass :: dgdx        => dgdx_cockburnandgau
          procedure, nopass :: dgdy        => dgdy_cockburnandgau
          procedure, nopass :: d2gdx2      => d2gdx2_cockburnandgau
          procedure, nopass :: d2gdy2      => d2gdy2_cockburnandgau
          procedure, nopass :: d2gdxdy     => d2gdxdy_cockburnandgau

        end type cg_operators

        contains

        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> give the boundary layer size
        !
        !> @date
        !> 30_07_2012 - initial version - J.L. Desmarais
        !
        !>@param bc_size
        !> boundary layer size
        !---------------------------------------------------------------
        function get_bc_size_cockburnandgau() result(bc_size)
          integer :: bc_size
          bc_size = 2
        end function get_bc_size_cockburnandgau


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ u_{i+\frac{1}{2},j}=
        !>\frac{1}{12}(-u_{i-1,j}+7 u_{i,j}+ 7 u_{i+1,j} - u_{i+2,j})\f$
        !
        !> @date
        !> 07_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the data
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
        function f_cockburnandgau(field_used,i,j,proc) result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

          if(rkind.eq.8) then
             var = 1.0d0/12.0d0*(
     $            -      proc(field_used,i-1,j)
     $            +7.0d0*proc(field_used,i,j)
     $            +7.0d0*proc(field_used,i+1,j)
     $            -      proc(field_used,i+2,j)
     $            )
          else
             var = 1./12.*(
     $            -  proc(field_used,i-1,j)
     $            +7*proc(field_used,i,j)
     $            +7*proc(field_used,i+1,j)
     $            -  proc(field_used,i+2,j)
     $            )
          end if

        end function f_cockburnandgau


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial x}\big|_{i+\frac{1}{2}
        !> ,j}= \frac{1}{\Delta x}(-u_{i,j}+u_{i+1,j})\f$
        !
        !> @date
        !> 07_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the data
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
        function dfdx_cockburnandgau(
     $     field_used,i,j,proc)
     $     result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

          if(rkind.eq.8) then
             var = 1.d0/field_used%dx*(
     $            -proc(field_used,i,j)
     $            +proc(field_used,i+1,j))
          else
             var = 1./field_used%dx*(
     $            -proc(field_used,i,j)
     $            +proc(field_used,i+1,j))
          end if

        end function dfdx_cockburnandgau


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial y}\big|
        !> _{i+\frac{1}{2},j}= \frac{1}{2 \Delta y}
        !> (-u_{i+\frac{1}{2},j-1}+u_{i+\frac{1}{2},j+1})\f$
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the data
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
        function dfdy_cockburnandgau(
     $     field_used,i,j,proc)
     $     result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

          !DEC$ FORCEINLINE RECURSIVE
          if(rkind.eq.8) then
             var = 1.0d0/(24.0d0*field_used%dy)*(
     $            +       proc(field_used,i-1,j-1)
     $            - 7.0d0*proc(field_used,i  ,j-1)
     $            - 7.0d0*proc(field_used,i+1,j-1)
     $            +       proc(field_used,i+2,j-1)
     $            -       proc(field_used,i-1,j+1)
     $            + 7.0d0*proc(field_used,i  ,j+1)
     $            + 7.0d0*proc(field_used,i+1,j+1)
     $            -       proc(field_used,i+2,j+1))
          else
             var = 1./(24.*field_used%dy)*(
     $         +   proc(field_used,i-1,j-1)
     $         - 7*proc(field_used,i  ,j-1)
     $         - 7*proc(field_used,i+1,j-1)
     $         +   proc(field_used,i+2,j-1)
     $         -   proc(field_used,i-1,j+1)
     $         + 7*proc(field_used,i  ,j+1)
     $         + 7*proc(field_used,i+1,j+1)
     $         -   proc(field_used,i+2,j+1))
          end if

        end function dfdy_cockburnandgau


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial^2 u}{\partial x^2}\big|
        !> _{i+\frac{1}{2},j}= \frac{1}{2 {\Delta x}^2}(u_{i-1,j} -
        !> u_{i,j} - u_{i+1,j} + u_{i+2,j}) \f$
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the data
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
        function d2fdx2_cockburnandgau(
     $     field_used,i,j,proc)
     $     result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

          !DEC$ FORCEINLINE RECURSIVE
          if(rkind.eq.8) then
             var = 0.5d0/(field_used%dx**2)*(
     $            +proc(field_used,i-1,j)
     $            -proc(field_used,i,j)
     $            -proc(field_used,i+1,j)
     $            +proc(field_used,i+2,j)
     $            )
          else
             var = 1./(2.*(field_used%dx**2))*(
     $            +proc(field_used,i-1,j)
     $            -proc(field_used,i,j)
     $            -proc(field_used,i+1,j)
     $            +proc(field_used,i+2,j)
     $            )
          end if

        end function d2fdx2_cockburnandgau


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial^2 u}{\partial y^2}
        !> \big|_{i+\frac{1}{2},j}= \frac{1}{{\Delta y}^2}(
        !> u_{i+\frac{1}{2},j-1} - 2 u_{i+\frac{1}{2},j} +
        !> u_{i+\frac{1}{2},j+1}) \f$
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the data
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
        function d2fdy2_cockburnandgau(
     $     field_used,i,j,proc)
     $     result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

          !DEC$ FORCEINLINE RECURSIVE
          if(rkind.eq.8) then
             var = (1.0d0/(12.0d0*field_used%dy**2))*(
     $            (
     $            -       proc(field_used,i-1,j-1)
     $            + 7.0d0*proc(field_used,i  ,j-1)
     $            + 7.0d0*proc(field_used,i+1,j-1)
     $            -       proc(field_used,i+2,j-1))
     $            - 2.0d0*(
     $            -       proc(field_used,i-1,j)
     $            + 7.0d0*proc(field_used,i  ,j)
     $            + 7.0d0*proc(field_used,i+1,j)
     $            -       proc(field_used,i+2,j))
     $            + (
     $            -       proc(field_used,i-1,j+1)
     $            + 7.0d0*proc(field_used,i  ,j+1)
     $            + 7.0d0*proc(field_used,i+1,j+1)
     $            -       proc(field_used,i+2,j+1))
     $            )
          else
             var = (1./(12*field_used%dy**2))*(
     $            (
     $            -   proc(field_used,i-1,j-1)
     $            + 7*proc(field_used,i  ,j-1)
     $            + 7*proc(field_used,i+1,j-1)
     $            -   proc(field_used,i+2,j-1))
     $            - 2*(
     $            -   proc(field_used,i-1,j)
     $            + 7*proc(field_used,i  ,j)
     $            + 7*proc(field_used,i+1,j)
     $            -   proc(field_used,i+2,j))
     $            + (
     $            -   proc(field_used,i-1,j+1)
     $            + 7*proc(field_used,i  ,j+1)
     $            + 7*proc(field_used,i+1,j+1)
     $            -   proc(field_used,i+2,j+1))
     $            )
          end if
        end function d2fdy2_cockburnandgau


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial^2 u}{\partial x \partial y}
        !> \big|_{i+\frac{1}{2},j} = \frac{1}{2 \Delta y}
        !> (-\frac{\partial u}{\partial x}\big|_{i+\frac{1}{2},j-1} +
        !> \frac{\partial u}{\partial x}\big|_{i+\frac{1}{2},j+1}) \f$
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the data
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
        function d2fdxdy_cockburnandgau(
     $     field_used,i,j,proc)
     $     result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

          !DEC$ FORCEINLINE RECURSIVE
          if(rkind.eq.8) then
             var =( proc(field_used,i,j-1)
     $            - proc(field_used,i+1,j-1)
     $            - proc(field_used,i,j+1)
     $            + proc(field_used,i+1,j+1))*
     $            0.5d0/(field_used%dy*field_used%dx)
          else
             var =( proc(field_used,i,j-1)
     $            - proc(field_used,i+1,j-1)
     $            - proc(field_used,i,j+1)
     $            + proc(field_used,i+1,j+1))
     $            /(2*field_used%dy*field_used%dx)
          end if

        end function d2fdxdy_cockburnandgau


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ u_{i,j+\frac{1}{2}}=
        !>\frac{1}{12}(-u_{i,j-1}+7 u_{i,j}+ 7 u_{i,j+1} - u_{i,j+2})\f$
        !
        !> @date
        !> 07_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the data
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
        function g_cockburnandgau(
     $     field_used,i,j,proc)
     $     result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

          !DEC$ FORCEINLINE RECURSIVE
          if(rkind.eq.8) then
             var = 1.0d0/12.0d0*(
     $            -      proc(field_used,i,j-1)
     $            +7.0d0*proc(field_used,i,j)
     $            +7.0d0*proc(field_used,i,j+1)
     $            -      proc(field_used,i,j+2)
     $            )
          else
             var = 1./12.*(
     $            -  proc(field_used,i,j-1)
     $            +7*proc(field_used,i,j)
     $            +7*proc(field_used,i,j+1)
     $            -  proc(field_used,i,j+2)
     $            )
          end if

        end function g_cockburnandgau        

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial x}|_{i
        !> ,j+\frac{1}{2}}= \frac{1}{2 \Delta x} (-u_{i-1,j+\frac{1}{2}}
        !> +u_{i+1,j+\frac{1}{2}})\f$
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the data
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
        function dgdx_cockburnandgau(
     $     field_used,i,j,proc)
     $     result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

          !DEC$ FORCEINLINE RECURSIVE
          if(rkind.eq.8) then
             var = 1.0d0/(24.0d0*field_used%dx)*(
     $            +       proc(field_used,i-1,j-1)
     $            - 7.0d0*proc(field_used,i-1,j)
     $            - 7.0d0*proc(field_used,i-1,j+1)
     $            +       proc(field_used,i-1,j+2)
     $            -       proc(field_used,i+1,j-1)
     $            + 7.0d0*proc(field_used,i+1,j)
     $            + 7.0d0*proc(field_used,i+1,j+1)
     $            -       proc(field_used,i+1,j+2))
          else
             var = 1./(24.*field_used%dx)*(
     $            +   proc(field_used,i-1,j-1)
     $            - 7*proc(field_used,i-1,j)
     $            - 7*proc(field_used,i-1,j+1)
     $            +   proc(field_used,i-1,j+2)
     $            -   proc(field_used,i+1,j-1)
     $            + 7*proc(field_used,i+1,j)
     $            + 7*proc(field_used,i+1,j+1)
     $            -   proc(field_used,i+1,j+2))
          end if

        end function dgdx_cockburnandgau


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial u}{\partial y}\big|_{i
        !> ,j+\frac{1}{2}}= \frac{1}{\Delta y}(-u_{i,j}+u_{i,j+1})\f$
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the data
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
        function dgdy_cockburnandgau(
     $     field_used,i,j,proc)
     $     result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

          !DEC$ FORCEINLINE RECURSIVE
          if(rkind.eq.8) then
             var = 1.0d0/field_used%dy*(
     $            -proc(field_used,i,j)
     $            +proc(field_used,i,j+1)
     $            )
          else
             var = 1./field_used%dy*(
     $            -proc(field_used,i,j)
     $            +proc(field_used,i,j+1)
     $            )
          end if

        end function dgdy_cockburnandgau        

      
        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial^2 u}{\partial x^2}
        !> \big|_{i,j+\frac{1}{2}}= \frac{1}{{\Delta x}^2}(
        !> u_{i-1,j+\frac{1}{2}}- 2 u_{i,j+\frac{1}{2}}
        !> + u_{i+1,j+\frac{1}{2}}) \f$
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the data
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
        function d2gdx2_cockburnandgau(
     $     field_used,i,j,proc)
     $     result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

          !DEC$ FORCEINLINE RECURSIVE
          if(rkind.eq.8) then
             var = 1.0d0/(12.0d0*field_used%dx**2) * (
     $            -      proc(field_used,i-1,j-1)
     $            +7.0d0*proc(field_used,i-1,j)
     $            +7.0d0*proc(field_used,i-1,j+1)
     $            -      proc(field_used,i-1,j+2)
     $            - 2.0d0*(
     $            -      proc(field_used,i,j-1)
     $            +7.0d0*proc(field_used,i,j)
     $            +7.0d0*proc(field_used,i,j+1)
     $            -      proc(field_used,i,j+2)
     $            )
     $            -      proc(field_used,i+1,j-1)
     $            +7.0d0*proc(field_used,i+1,j)
     $            +7.0d0*proc(field_used,i+1,j+1)
     $            -      proc(field_used,i+1,j+2)
     $            )
          else
             var = 1./(12.*field_used%dx**2) * (
     $            -  proc(field_used,i-1,j-1)
     $            +7*proc(field_used,i-1,j)
     $            +7*proc(field_used,i-1,j+1)
     $            -  proc(field_used,i-1,j+2)
     $            - 2*(
     $            -  proc(field_used,i,j-1)
     $            +7*proc(field_used,i,j)
     $            +7*proc(field_used,i,j+1)
     $            -  proc(field_used,i,j+2)
     $            )
     $            -  proc(field_used,i+1,j-1)
     $            +7*proc(field_used,i+1,j)
     $            +7*proc(field_used,i+1,j+1)
     $            -  proc(field_used,i+1,j+2)
     $            )
          end if

        end function d2gdx2_cockburnandgau


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial^2 u}{\partial y^2}\big|
        !> _{i,j+\frac{1}{2}}= \frac{1}{2 {\Delta y}^2}(u_{i,j-1} -
        !> u_{i,j} - u_{i,j+1} + u_{i,j+2}) \f$
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the data
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
        function d2gdy2_cockburnandgau(
     $     field_used,i,j,proc)
     $     result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

          !DEC$ FORCEINLINE RECURSIVE
          if(rkind.eq.8) then
             var = 0.5d0/(field_used%dy**2)*(
     $            + proc(field_used,i,j-1)
     $            - proc(field_used,i,j)
     $            - proc(field_used,i,j+1)
     $            + proc(field_used,i,j+2)
     $            )
          else
             var = 1./(2.*field_used%dy**2)*(
     $            + proc(field_used,i,j-1)
     $            - proc(field_used,i,j)
     $            - proc(field_used,i,j+1)
     $            + proc(field_used,i,j+2)
     $            )
          end if

        end function d2gdy2_cockburnandgau


        !> @author
        !> Julien L. Desmarais
        !
        !> @brief
        !> compute \f$ \frac{\partial^2 u}{\partial x \partial y}
        !> \big|_{i,j+\frac{1}{2}} = \frac{1}{2 \Delta x}
        !> (-\frac{\partial u}{\partial y}\big|_{i-1,j+\frac{1}{2}} +
        !> \frac{\partial u}{\partial y}\big|_{i+1,j+\frac{1}{2}}) \f$
        !
        !> @date
        !> 08_08_2013 - initial version - J.L. Desmarais
        !
        !>@param field_used
        !> object encapsulating the data
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
        function d2gdxdy_cockburnandgau(
     $     field_used,i,j,proc)
     $     result(var)

          implicit none

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

          !DEC$ FORCEINLINE RECURSIVE
          if(rkind.eq.8) then
             var =(
     $            proc(field_used,i-1,j)
     $            - proc(field_used,i-1,j+1)
     $            - proc(field_used,i+1,j)
     $            + proc(field_used,i+1,j+1))*
     $            0.5d0/(field_used%dx*field_used%dy)
          else
             var =(
     $            proc(field_used,i-1,j)
     $            - proc(field_used,i-1,j+1)
     $            - proc(field_used,i+1,j)
     $            + proc(field_used,i+1,j+1))
     $            /(2*field_used%dx*field_used%dy)
          end if

        end function d2gdxdy_cockburnandgau
        
      end module cg_operators_class
