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
      !-----------------------------------------------------------------
      module sd_operators_class

        use field_class      , only : field
        use interface_primary, only : get_primary_var
        use parameters_kind  , only : ikind, rkind

        implicit none

        private
        public :: sd_operators


        !> @class sd_operators
        !> abstract class encapsulating spatial discretization operators
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
        type, abstract :: sd_operators

          contains
          
          procedure(get_bc_size_proc)   , nopass, deferred :: get_bc_size

          procedure(space_operator_proc), nopass, deferred :: f
          procedure(space_operator_proc), nopass, deferred :: dfdx
          procedure(space_operator_proc), nopass, deferred :: dfdy
          procedure(space_operator_proc), nopass, deferred :: d2fdx2
          procedure(space_operator_proc), nopass, deferred :: d2fdy2
          procedure(space_operator_proc), nopass, deferred :: d2fdxdy

          procedure(space_operator_proc), nopass, deferred :: g
          procedure(space_operator_proc), nopass, deferred :: dgdx
          procedure(space_operator_proc), nopass, deferred :: dgdy
          procedure(space_operator_proc), nopass, deferred :: d2gdx2
          procedure(space_operator_proc), nopass, deferred :: d2gdy2
          procedure(space_operator_proc), nopass, deferred :: d2gdxdy
          

        end type sd_operators


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
        function space_operator_proc(field_used,i,j,proc) result(var)

          import field
          import get_primary_var
          import ikind
          import rkind

          class(field)  , intent(in) :: field_used
          integer(ikind), intent(in) :: i
          integer(ikind), intent(in) :: j
          procedure(get_primary_var) :: proc
          real(rkind)                :: var

        end function space_operator_proc

      end interface

      end module sd_operators_class
