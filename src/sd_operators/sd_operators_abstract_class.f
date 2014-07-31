      !> @file
      !> abstract class encapsulating subroutines for the space
      !> discretisation
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines for the space discretisation
      !> scheme
      !
      !> @date
      !> 30_07_2012 - initial version                   - J.L. Desmarais
      !> 04_06_2013 - functions for crossed derivatives - J.L. Desmarais
      !> 07_08_2013 - transfered in lernaeanhydra_opt   - J.L. Desmarais
      !> 11_07_2014 - interface for erymanthianboar use - J.L. Desmarais
      !-----------------------------------------------------------------
      module sd_operators_abstract_class

        use interface_primary, only : get_primary_var, get_secondary_var
        use parameters_kind  , only : ikind, rkind

        implicit none

        private
        public :: sd_operators_abstract


        !> @class sd_operators_abstract
        !> abstract class encapsulating spatial discretization operators
        !
        !> @param get_bc_size
        !> get the boundary layer size
        !
        !> @param f
        !> evaluate data at [i-1/2,j]
        !
        !> @param dfdx
        !> evaluate \f$\frac{\partial}{\partial x}\f$ at [i-1/2,j]
        !
        !> @param dfdy
        !> evaluate \f$\frac{\partial}{\partial y}\f$ at [i-1/2,j]
        !
        !> @param d2fdx2
        !> evaluate \f$\frac{\partial}{\partial x^2}\f$ at [i-1/2,j]
        !
        !> @param d2fdy2
        !> evaluate \f$\frac{\partial}{\partial y^2}\f$ at [i-1/2,j]
        !
        !> @param d2fdxdy
        !> evaluate \f$\frac{\partial}{\partial x \partial y}\f$
        !> at [i-1/2,j]
        !
        !> @param g
        !> evaluate data at [i,j-1/2]
        !
        !> @param dgdx
        !> evaluate \f$\frac{\partial}{\partial x}\f$ at [i,j-1/2]
        !
        !> @param dgdy
        !> evaluate \f$\frac{\partial}{\partial y}\f$ at [i,j-1/2]
        !
        !> @param d2gdx2
        !> evaluate \f$\frac{\partial}{\partial x^2}\f$ at [i,j-1/2]
        !
        !> @param d2gdy2
        !> evaluate \f$\frac{\partial}{\partial y^2}\f$ at [i,j-1/2]
        !
        !> @param d2gdxdy
        !> evaluate \f$\frac{\partial}{\partial x \partial y}\f$
        !> at [i,j-1/2]
        !---------------------------------------------------------------
        type, abstract :: sd_operators_abstract

          contains
          
          procedure(get_bc_size_proc)      , nopass, deferred :: get_bc_size

          procedure(space_operator_proc)   , nopass, deferred :: f
          procedure(space_operator_proc_x) , nopass, deferred :: dfdx
          procedure(space_operator_proc_nl), nopass, deferred :: dfdx_nl
          procedure(space_operator_proc_y) , nopass, deferred :: dfdy
          procedure(space_operator_proc_x) , nopass, deferred :: d2fdx2
          procedure(space_operator_proc_y) , nopass, deferred :: d2fdy2
          procedure(space_operator_proc_xy), nopass, deferred :: d2fdxdy

          procedure(space_operator_proc)   , nopass, deferred :: g
          procedure(space_operator_proc_x) , nopass, deferred :: dgdx
          procedure(space_operator_proc_y) , nopass, deferred :: dgdy
          procedure(space_operator_proc_nl), nopass, deferred :: dgdy_nl
          procedure(space_operator_proc_x) , nopass, deferred :: d2gdx2
          procedure(space_operator_proc_y) , nopass, deferred :: d2gdy2
          procedure(space_operator_proc_xy), nopass, deferred :: d2gdxdy

        end type sd_operators_abstract


        abstract interface


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> give the boundary layer size
        !
        !> @date
        !> 30_07_2012 - initial version - J.L. Desmarais
        !
        !>@param bc_size_op
        !> boundary layer size
        !---------------------------------------------------------------
        function get_bc_size_proc() result(bc_size_op)

          integer :: bc_size_op

        end function get_bc_size_proc

        end interface


        abstract interface

        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> space operator interface
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
        !>@param proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function space_operator_proc(nodes,i,j,proc) result(var)

          import get_primary_var
          import ikind
          import rkind

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                               :: var

        end function space_operator_proc


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> space operator interface
        !
        !> @date
        !> 11_07_2014 - initial version - J.L. Desmarais
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
        function space_operator_proc_x(nodes,i,j,proc,dx) result(var)

          import get_primary_var
          import ikind
          import rkind

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                               :: var

        end function space_operator_proc_x


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> space operator interface
        !
        !> @date
        !> 11_07_2014 - initial version - J.L. Desmarais
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
        function space_operator_proc_y(nodes,i,j,proc,dy) result(var)

          import get_primary_var
          import ikind
          import rkind

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

        end function space_operator_proc_y


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> space operator interface
        !
        !> @date
        !> 11_07_2014 - initial version - J.L. Desmarais
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
        function space_operator_proc_xy(nodes,i,j,proc,dx,dy) result(var)

          import get_primary_var
          import ikind
          import rkind

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_primary_var)                :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

        end function space_operator_proc_xy


        !> @author 
        !> Julien L. Desmarais
        !
        !> @brief
        !> space operator interface
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
        !>@param proc
        !> procedure computing the special quantity evaluated at [i,j]
        !> (ex: pressure, temperature,...)
        !
        !>@param var
        !> data evaluated at [i,j]
        !---------------------------------------------------------------
        function space_operator_proc_nl(nodes,i,j,proc,dx,dy) result(var)

          import get_secondary_var
          import ikind
          import rkind

          real(rkind), dimension(:,:,:), intent(in) :: nodes
          integer(ikind)               , intent(in) :: i
          integer(ikind)               , intent(in) :: j
          procedure(get_secondary_var)              :: proc
          real(rkind)                  , intent(in) :: dx
          real(rkind)                  , intent(in) :: dy
          real(rkind)                               :: var

        end function space_operator_proc_nl

      end interface

      end module sd_operators_abstract_class
