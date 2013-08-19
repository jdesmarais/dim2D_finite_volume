      !> @file
      !> abstract class encapsulating subroutines to compute
      !> the time derivatives of the main variables from the
      !> governing equations using the space discretisation
      !> operators and the physical model
      !
      !> @author 
      !> Julien L. Desmarais
      !
      !> @brief
      !> class encapsulating subroutines to compute
      !> the time derivatives of the main variables
      !> in the governing equations
      !
      !> @date
      !> 13_08_2013 - initial version                   - J.L. Desmarais
      !-----------------------------------------------------------------
      module td_operators_class

        use field_class       , only : field
        use cg_operators_class, only : cg_operators
        use parameters_input  , only : nx,ny,ne
        use parameters_kind   , only : rkind
        use phy_model_eq_class, only : phy_model_eq

        implicit none


        private
        public :: td_operators


        !> @class td_operators
        !> abstract class encapsulating operators to compute
        !> the time derivatives of the main variables from 
        !> the governing equations
        !>
        !> @param compute_time_dev
        !> compute the time derivatives
        !---------------------------------------------------------------
        type, abstract :: td_operators

          contains
          procedure(time), nopass, deferred :: compute_time_dev

        end type td_operators


        abstract interface

          !> @author
          !> Julien L. Desmarais
          !
          !> @brief
          !> interface to compute the time derivatives using the
          !> space discretisation operators and the physical model
          !
          !> @date
          !> 13_08_2013 - initial version - J.L. Desmarais
          !
          !>@param field_used
          !> object encapsulating the main variables
          !
          !>@param s
          !> space discretization operators
          !
          !>@param p
          !> physical model
          !
          !>@param time_dev
          !> time derivatives
          !--------------------------------------------------------------
          function time(field_used, s, p_model) result(time_dev)

            import cg_operators
            import field
            import phy_model_eq
            import rkind
            import nx,ny,ne

            class(field)                    , intent(in)   :: field_used
            type(cg_operators)              , intent(in)   :: s
            class(phy_model_eq)             , intent(in)   :: p_model
            real(rkind), dimension(nx,ny,ne)               :: time_dev

          end function time

        end interface

      end module td_operators_class
